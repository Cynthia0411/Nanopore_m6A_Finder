import pickle
import os
from numpy import loadtxt
import pandas as pd


def main():
    xgb_model = pickle.load(open('./XGB_model_alldata.dat','rb'))
    filenames = os.listdir('./polyAFeature')
    result = []
    readName=[]
    ref = []
    tail_length = []
    df_ref = pd.read_csv('./polyA_pos_427.tsv',sep='\t')
    read_ref = {}
    for i in range(len(df_ref)):
        read_ref[df_ref['read_name'][i]] = df_ref['reference'][i]
    for files in filenames:
        filedir = './polyAFeature/'+files
        cur_pred = xgb_model.predict(loadtxt(filedir,delimiter=','))
        cur_result=sum(cur_pred)/(len(cur_pred))
        cur_readName = files.split('.')[0]
        result.append(cur_result)
        readName.append(cur_readName)
        ref.append(read_ref[cur_readName])
        tail_length.append(len(cur_pred)+4)
    polyALevelTable = pd.DataFrame({'read_name':readName,'reference':ref,'tail_length':tail_length,'polyA_level':result})
    polyALevelTable.reset_index(drop=True)
    polyALevelTable.to_csv('./polyA_level_427_tailL.csv',sep=',')

if __name__ == '__main__':
    main()


