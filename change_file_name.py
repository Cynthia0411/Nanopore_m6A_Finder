import os
import h5py
import pandas as pd
import numpy as np

def change_filename(filedir):
    files = os.listdir(filedir)
    for filename in files:
        print(filename)
        portion = os.path.splitext(filename)
        if portion[1] == '.fast5':
            old_filename = filedir+'/'+filename
            f = h5py.File(old_filename)
            newname = filedir+"/"+str(f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Fastq'][()]).split(' ')[0][3:]+'.fast5'
            print(newname)
            f.close()
            os.rename(old_filename,newname)
    return
def process_files(filedir,ref_file):
    f_ref = pd.read_csv(ref_file,sep = '\t')
    ref_split = list(f_ref.groupby("READ_NAME"))
    files=os.listdir(filedir)
    kmerlength = 5
    tra = []
    for i in range(len(ref_split)):
        cur_df = ref_split[i][1]
        for filename in files:
            portion = os.path.splitext(filename)
            if portion[0] == ref_split[i][0]:
                f = h5py.File(filedir +'/'+filename)
                tra_feature = np.array(f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Trace'])
                move = f['Analyses']['Basecall_1D_001']['BaseCalled_template']['Move']
                index_m = []
                index_s = []
                for j in range(len(move)):
                    if move[j]==1:
                        index_m.append(j)
                for k in range(1,len(index_m)+1):
                    index_s.append(index_m[-k])
                for l in cur_df['READ_POS']:
                    tra_cur = []
                    for n in range(int(l),int(l)+kmerlength):
                        tra_l = (tra_feature[index_s[n]][0:4]+tra_feature[index_s[n]][4:])/255
                        tra_l = [str(round(obj,2)) for obj in tra_l]
                        print("TRACE_CUR:====================")
                        print(tra_l)
                        tra_cur.append(','.join(tra_l))
                        print('tra_cur:~~~~~~~~~~~~~~~~~~~')
                        print(tra_cur)
                    tra.append(','.join(tra_cur))

#    f_ref['trace_ACGT'] = tra
    f_ref_table = pd.DataFrame({'label':[1]*len(f_ref),'BASE':f_ref['Base'],'QUAL':f_ref['QUAL'],'TRACE_ACGT':tra})
    f_ref_table.to_csv('./select_reads/82_1_5mers_withA_trace.tsv',sep='\t')
    return 0

#process_files('./guppy_1368282_1/workspace','./select_reads/82_1_5mers_withA_order.tsv')
#change_filename(str(snakemake.input))
change_filename("./guppy_fast5/workspace")





