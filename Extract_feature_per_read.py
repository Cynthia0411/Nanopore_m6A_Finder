import pandas as pd
import h5py
import numpy as np
import re
import os


#start_pos = 127
#end_pos = 137
kmerLength = 5
#filedir = 'D:/old/data/bioinformatics/VSG_mapped_reads_with_polyA/7656646f-042e-449a-aa0e-747bf145e6be.fast5'
#inputsamfile = 'mapped_WT2_polyA.sam'
def getStart_End(inputsam):
    with open(inputsam,'r') as file_to_read:
        with open('polyA_pos_WT2.tsv','a') as file_to_write:
            file_to_write.write('read_name'+'\t'+'reference'+'\t'+'polyA_pos_start'+'\t'+'polyA_pos_end'+'\n')
            for line in file_to_read:
                line = line.strip()
                line_list=line.split('\t')
                seq_name = line_list[0]
                reference = line_list[2]
                seq = line_list[9]
                pos_end = list(re.finditer('A{10}',seq))[-1].end()
                cigar = line_list[5]
                pos_start = len(seq) - int(list(re.finditer(r'(\d+)S',cigar))[-1].group(1))
                if pos_start <= pos_end -10:
                    file_to_write.write(seq_name+'\t'+reference+'\t'+str(pos_start)+'\t'+str(pos_end)+'\n')



def processFile(filedir,start_pos,end_pos):
    f_to_open = h5py.File(filedir)
    fastq_str = str(f_to_open['Analyses']['Basecall_1D_000']['BaseCalled_template']['Fastq'][()]).split('\\n')
    read_name = fastq_str[0][3:39]
    tra = []
    read_base = fastq_str[1]
    base = []
    read_qual = [str(int(ord(obj) - 33)) for obj in list(fastq_str[3])]
    qual = []
    tra_feature = np.array(f_to_open['Analyses']['Basecall_1D_000']['BaseCalled_template']['Trace'])
    move = f_to_open['Analyses']['Basecall_1D_000']['BaseCalled_template']['Move']
    index_move = []
    index_seq = []
    for i in range(len(move)):
        if move[i] == 1:
            index_move.append(i)
    for j in range(1,len(index_move)+1):
        index_seq.append(index_move[-j])

#    for l in range(start_pos-1,end_pos,kmerLength):
    pos_l = []
    for l in range(start_pos,end_pos-kmerLength):
        pos_l.append(l)
        tra_cur = []
        if int(l) + kmerLength <= end_pos:
            for n in range(int(l),int(l)+kmerLength):
                tra_l = (tra_feature[index_seq[n]][0:4]+tra_feature[index_seq[n]][4:])/255
                tra_l = [str(round(obj,2)) for obj in tra_l]
                tra_cur.append(','.join(tra_l))
            tra.append(','.join(tra_cur))
            base.append(read_base[l:l+kmerLength])
            qual.append(','.join(read_qual[l:l+kmerLength]))
        else:
            break
    feature_table = pd.DataFrame({'BASE':base,'QUAL':qual,'TRACE_ACGT':tra,'POS_READ':pos_l})
    feature_table.to_csv("./CDsReads_HEK293T_WT_rep1/"+read_name+'.tsv',sep='\t')

def main():
#    getStart_End(inputsamfile)
    f_ref = pd.read_csv("./CDs_pos_HEK293T_WT_rep1.tsv", sep = '\t')
    if os.path.exists("./CDsReads_HEK293T_WT_rep1")==False:
        os.mkdir("CDsReads_HEK293T_WT_rep1")
    for i in range(len(f_ref)):
        filename = f_ref['read_name'][i] + '.fast5'
        filedir = "./HEK293T-WT-rep1-fast5"+"/" + filename
        start_pos = f_ref['polyA_pos_start'][i]
        end_pos = f_ref['polyA_pos_end'][i]
        if os.path.exists(filedir):
            processFile(filedir,start_pos,end_pos)
        else:
            continue
#    getStart_End(inputsamfile)
    return

if __name__ == '__main__':
    main()


