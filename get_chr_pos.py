import csv
import pandas as pd
from tqdm import tqdm
from io import StringIO  # 导入 io 模块中的 StringIO 类

input_sam_tsv = "mapped_reads_WT.tsv"
input_read_pos = "m6A_pos_in_reads_xgb_0.7_HMEC_WT.csv"
output_dir = "m6A_pos_on_genome_xgb_0.7_WT.tsv"

#csv.field_size_limit(1024 * 1024 * 1024)

def readintsv(filename):
    df = pd.read_csv(filename, sep='\t', engine='python', quoting = csv.QUOTE_NONE)
    df_new = df.reset_index()
    return df_new

def readincsv(filename):
    df = pd.read_csv(filename, sep=',', index_col= 0)
    df_new = df.reset_index(drop = True)
    return df_new






# 自定义分块逻辑
def custom_chunking(tsv_file, field_to_group_by):
    chunk = []
    current_chunk_key = None

    with open(tsv_file, 'r') as file:
        header = next(file).strip().split('\t')  # 读取第一行作为 header，并分割成列名
        field_index = header.index(field_to_group_by)  # 获取指定字段的索引

        for line in file:
            row = line.strip().split('\t')  # 分割每一行的值
            group = row[field_index]

            if current_chunk_key is None:
                current_chunk_key = group

            if group == current_chunk_key:
                chunk.append(row)
            else:
                yield pd.DataFrame(chunk, columns=header)
                chunk = [row]
                current_chunk_key = group

    yield pd.DataFrame(chunk, columns=header)





def getposonchrom(startpos, dict_cur):
    ref_cur = list(dict_cur.keys())[0].split("_")[0]
    chrom_pos = []
    for i in range(startpos,startpos+5):

        key_cur =  ref_cur + "_"+str(i)

        if key_cur in dict_cur.keys():

            pos_cur = dict_cur[key_cur][0]    # read pos is 0-based, ref pos is 1 based
            base_cur = dict_cur[key_cur][1]
            flag_cur = dict_cur[key_cur][2]
            if flag_cur == "16" and base_cur == "T":
                chrom_pos.append(int(pos_cur)-1) 
            elif flag_cur == "0" and base_cur == "A":
                chrom_pos.append(int(pos_cur)-1)
            else:
                continue
        
        else:
            continue
#        if pos_cur != ".":
#        if flag_cur == "16" and base_cur == "T":
#            chrom_pos.append(int(pos_cur)-1)
#        elif flag_cur == "0" and base_cur == "A":
#            chrom_pos.append(int(pos_cur)-1)
#        else:
#            continue
    # bed file is half open (left included, right exluded)
#    if len(chrom_pos) >0:
#        beginpos = min(chrom_pos) -1
#        endpos = max(chrom_pos)
#        return ref_cur, beginpos, endpos
    return ref_cur, chrom_pos
#    else:
#        return -1,-1,-1


def divideSamDataframe(df_value, df_pos_split, file_to_write):
#    df_split = df_sam.groupby('#READ_NAME')
#    df_pos_split = df_pos.groupby('read_name')
    df_value_index = df_value.index.values
#    print(df_value['#READ_NAME'])
#    with open(output_dir,'a') as file_to_write:
#        file_to_write.write('Readname' + "\t" + 'Chromosome' + "\t" + 'Start' + '\t' +"End" + "\t" + "Flag" + '\n')
    key = df_value['#READ_NAME'][df_value_index[0]]
    if key in df_pos_split.groups.keys():
        df_pos_cur = df_pos_split.get_group(key)
        cur_pos_index = df_pos_cur.index.values
        dict_cur = {}
        cur_index = df_value.index.values
        if df_value['FLAG'][cur_index[0]] == '16':
            len_read_cur = int(df_value['READ_POS'][cur_index[-1]])
            for i in cur_index:
                if df_value['OP'][i] != 'S':
                    if df_value['READ_POS'][i] != '.':
                        key_cur = df_value['CHROM'][i] + "_" + str(len_read_cur-int(df_value['READ_POS'][i]))              # flag=16, converse read_pos
#                        print(key_cur)
                        dict_cur[key_cur] = (df_value['REF_POS'][i],df_value['REF'][i],df_value['FLAG'][i])
                    else:
                        continue
                else:
                    continue
            for j in cur_pos_index:
                start = df_pos_cur['readpos_start'][j]
                ref_cur, A_pos_list = getposonchrom(start, dict_cur)
#                if ref_cur == -1 and beginpos== -1 and endpos == -1 :
                if len(A_pos_list)== 1:
                    for A_pos in A_pos_list:
                        endpos = int(A_pos)+1
                        file_to_write.write(key + '\t' + ref_cur + '\t' + str(A_pos) + '\t' + str(endpos) + '\t'+ '16' +'\n')
                else:
                    continue
        else:
            for i in cur_index:
                if df_value['OP'][i] != 'S':
                    if df_value['READ_POS'][i] != '.':
                        key_cur = df_value['CHROM'][i] + "_" + str(df_value['READ_POS'][i])                         
                        dict_cur[key_cur] = (df_value['REF_POS'][i],df_value['REF'][i],df_value['FLAG'][i])
                    else:
                        continue
                else:
                    continue
            for j in cur_pos_index:
                start = df_pos_cur['readpos_start'][j]
                ref_cur, A_pos_list = getposonchrom(start, dict_cur)
#                if ref_cur == -1 and beginpos== -1 and endpos == -1 :
                if len(A_pos_list) == 1:
                    for A_pos in A_pos_list:
                        endpos = int(A_pos)+1
                        file_to_write.write(key + '\t' + ref_cur + '\t' + str(A_pos) + '\t' + str(endpos) + '\t'+'0'+'\n')
                else:
                    continue
                            
        return


def main():
#    df_sam = readintsv(input_sam_tsv)
    df_pos = readincsv(input_read_pos)
    df_pos_split = df_pos.groupby("read_name")
#   Count the total number of chunks
    total_chunks = sum(1 for _ in custom_chunking(input_sam_tsv, '#READ_NAME'))
   # 逐块读取文件分块
    with tqdm(total=total_chunks, desc="Processing") as pbar:
        with open(output_dir,"a") as file_to_write:
            file_to_write.write('Readname' + "\t" + "Chromosome" + "\t" + "Start" + "\t" + "End" +"\t"+"Flag"+ "\n")
            for df_sam in custom_chunking(input_sam_tsv, '#READ_NAME'):
            # 在每个分块上执行操作，可以使用列名访问列
                divideSamDataframe(df_sam,df_pos_split, file_to_write)
                pbar.update(1)  # Update the progress bar
    return

if __name__=='__main__':
    main()
