import pandas as pd
import csv
import re

csv.field_size_limit(1024 * 1024 * 1024)

def loadfile(filename):
    df = pd.read_csv(filename, sep='\t',engine='python', quoting = csv.QUOTE_NONE)
    return df
def loadcsv(filename):
    df_ref = pd.read_csv(filename)
    return df_ref

def devide_dataframe(df,df_ref):
    df = df[df['FLAG']==0]
    df_split = list(df.groupby('#READ_NAME'))
    df_ref_cc6m_2244_t7_ecorv = df_ref.loc[df_ref['CHROM']=='cc6m_2244_t7_ecorv ']
    df_ref_cc6m_2459_t7_ecorv = df_ref.loc[df_ref['CHROM']=='cc6m_2459_t7_ecorv ']
    df_ref_cc6m_2595_t7_ecorv = df_ref.loc[df_ref['CHROM']=='cc6m_2595_t7_ecorv ']
    df_ref_cc6m_2709_t7_ecorv = df_ref.loc[df_ref['CHROM']=='cc6m_2709_t7_ecorv ']
    qual = []
    read_pos = []
    base = []
    ref = []
    read_name = []
    kmerLength = 5
    s_ref_pat = re.compile(r'^(?:[CGT]{0,4}A[CGT]{0,4})$(?<=[ACGT]{5})')
    for i in range(len(df_split)):
        print('new one!')
        curr_df = df_split[i][1]
        index_cur = df_split[i][1].index.values
        if curr_df['CHROM'][index_cur[0]] == 'cc6m_2244_t7_ecorv':
            for j in range(0,len(curr_df)-kmerLength+1):
                for k in range(len(df_ref_cc6m_2244_t7_ecorv)):
                    index_2244 = df_ref_cc6m_2244_t7_ecorv.index.values
                    if curr_df['REF_POS'][index_cur[j]] == str(df_ref_cc6m_2244_t7_ecorv['REF_POS'][index_2244[k]]+1) and curr_df['OP'][index_cur[j]]== "M":
                        s = ''.join(curr_df['BASE'][j:j+kmerLength])
                        s_ref = ''.join(curr_df['REF'][j:j+kmerLength])
                        if s_ref_pat.match(s_ref) and re.match('[ACGT]{5}', s):
                                print('========================in!')
                                qual.append(','.join([str(ord(obj)-33) for obj in curr_df['QUAL'][index_cur[j:j+kmerLength]]]))
                                read_pos.append(int(curr_df['READ_POS'][index_cur[j]]))
                                base.append(s)
                                ref.append(s_ref)
                                read_name.append(curr_df['#READ_NAME'][index_cur[j]])
        elif curr_df['CHROM'][index_cur[0]]== 'cc6m_2459_t7_ecorv':
            for j in range(0,len(curr_df)-kmerLength+1):
                for k in range(len(df_ref_cc6m_2459_t7_ecorv)):
                    index_2459 = df_ref_cc6m_2459_t7_ecorv.index.values
                    if curr_df['REF_POS'][index_cur[j]]==str(df_ref_cc6m_2459_t7_ecorv['REF_POS'][index_2459[k]]+1) and curr_df['OP'][index_cur[j]]=='M':
                        s = ''.join(curr_df['BASE'][j:j + kmerLength])
                        s_ref = ''.join(curr_df['REF'][j:j + kmerLength])
                        if s_ref_pat.match(s_ref) and re.match('[ACGT]{5}', s):
                                print("2459in!")
                                qual.append(','.join([str(ord(obj) - 33) for obj in curr_df['QUAL'][index_cur[j:j + kmerLength]]]))
                                read_pos.append(int(curr_df['READ_POS'][index_cur[j]]))
                                base.append(s)
                                ref.append(s_ref)
                                read_name.append(curr_df['#READ_NAME'][index_cur[j]])
        elif curr_df['CHROM'][index_cur[0]]== 'cc6m_2595_t7_ecorv':
            for j in range(0,len(curr_df)-kmerLength+1):
                for k in range(len(df_ref_cc6m_2595_t7_ecorv)):
                    index_2595 = df_ref_cc6m_2595_t7_ecorv.index.values
                    if curr_df['REF_POS'][index_cur[j]]== str(df_ref_cc6m_2595_t7_ecorv['REF_POS'][index_2595[k]]+1) and curr_df['OP'][index_cur[j]]=='M':
                        s = ''.join(curr_df['BASE'][j:j + kmerLength])
                        s_ref = ''.join(curr_df['REF'][j:j + kmerLength])
                        if s_ref_pat.match(s_ref) and re.match('[ACGT]{5}', s):
                                print('2595in!')
                                qual.append(','.join([str(ord(obj) - 33) for obj in curr_df['QUAL'][index_cur[j:j + kmerLength]]]))
                                read_pos.append(int(curr_df['READ_POS'][index_cur[j]]))
                                base.append(s)
                                ref.append(s_ref)
                                read_name.append(curr_df['#READ_NAME'][index_cur[j]])
        elif curr_df['CHROM'][index_cur[0]]== 'cc6m_2709_t7_ecorv':
            for j in range(0,len(curr_df)-kmerLength+1):
                for k in range(len(df_ref_cc6m_2709_t7_ecorv)):
                    index_2709 = df_ref_cc6m_2709_t7_ecorv.index.values
                    if curr_df['REF_POS'][index_cur[j]]==str(df_ref_cc6m_2709_t7_ecorv['REF_POS'][index_2709[k]]+1) and curr_df['OP'][index_cur[j]]=='M':
                        s = ''.join(curr_df['BASE'][j:j + kmerLength])
                        s_ref = ''.join(curr_df['REF'][j:j + kmerLength])
                        if s_ref_pat.match(s_ref) and re.match('[ACGT]{5}', s):
                                print('2709in!')
                                qual.append(','.join([str(ord(obj) - 33) for obj in curr_df['QUAL'][index_cur[j:j + kmerLength]]]))
                                read_pos.append(int(curr_df['READ_POS'][index_cur[j]]))
                                base.append(s)
                                ref.append(s_ref)
                                read_name.append(curr_df['#READ_NAME'][index_cur[j]])
    new = {'READ_NAME':read_name, 'READ_POS':read_pos, 'Base': base, 'REF':ref, 'QUAL':qual}
    print("Finish")
    return new

def main():
    df = loadfile('81_1_1500_top100.tsv')
    df_ref = loadcsv('curlcake_with_1_A.position.csv')
    newtable = devide_dataframe(df,df_ref)
    newdf = pd.DataFrame(newtable)
    newdf.to_csv('81_1_5mers_with_1_A.tsv', sep = '\t')
if __name__=='__main__':
    main()


