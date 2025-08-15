import pandas as pd
import numpy as np
import os

def getReadNameList(filedir):
    files = os.listdir(filedir)
    readnameList = []
    for filename in files:
        portion = os.path.splitext(filename)
        if portion[1] == '.tsv':
            readNamedir = filedir + '/'+filename
            readnameList.append(readNamedir)
    return readnameList        
def combind2files(filename1,filename2):
#    f1 = pd.read_csv(filename1,sep='\t')
    f1 = pd.read_csv(filename1,sep=',',header = None)
    f1 = f1.sample(frac=1)
#    f2 = pd.read_csv(filename2,sep='\t')
    f2 = pd.read_csv(filename2,sep = ',',header = None)
    f2 = f2.sample(frac=1)
    resultdf = pd.concat([f1,f2],ignore_index = True)
#    resultdf.reset_index(drop=True)
    return resultdf
def getReadName(filename):
    df_read = pd.read_csv(filename,sep='\t')
    read_name = filename.split(sep='/')[-1].split(sep = '.')[0]
    df_read.reset_index(drop=True)

    return df_read,read_name


def onehotencoding(resultdf):
    seq = resultdf['BASE']
    qualScore = resultdf['QUAL']
    seqDict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "U": [0, 0, 0, 1]}
    qualDic = {"extremeLow": 0.1, 'low': 0.4, "middle": 0.6, "high": 0.7, "extremeHigh": 1}
    
    # 处理序列的one-hot编码
    seqVec = np.zeros((len(seq), 20))
    for i in range(len(seq)):
        seqList = list(seq[i])
        lineVector = []
        for j in range(5):
            lineVector.extend(seqDict[seqList[j]])
        seqVec[i] = lineVector
    
    # 处理质量值
    qualVec = np.zeros((len(qualScore), 5))
    for i in range(len(qualScore)):
        qualList = qualScore[i].split(sep=',')
        lineVector = []
        for j in range(5):
            # 只添加分类后的质量值，不再添加原始值
            if int(qualList[j]) <= 5:
                lineVector.append(qualDic["extremeLow"])
            elif 5 < int(qualList[j]) <= 8:
                lineVector.append(qualDic["low"])
            elif 8 < int(qualList[j]) <= 25:
                lineVector.append(qualDic["middle"])
            elif 25 < int(qualList[j]) <= 35:
                lineVector.append(qualDic["high"])
            elif int(qualList[j]) > 35:
                lineVector.append(qualDic["extremeHigh"])
        qualVec[i] = lineVector
    
    return seqVec, qualVec

def savefile(resultdf, qualVec, seqVec, readname):
    fr = open(readname+".csv", "w")
    
    # 获取TRACE_ACGT数据并转换为列表
    trace_data = [row.split(',') for row in resultdf['TRACE_ACGT']]
    
    for i in range(len(qualVec)):
        # 构建输出行
        # 1. 添加标签 (-1 for 81, 1 for 82)
        label = "-1" if readname == 'feature_with_1_A_81' else "1"
        
        # 2. 添加one-hot编码的序列特征
        seq_features = ",".join(str(x) for x in seqVec[i])
        
        # 3. 添加质量值特征
        qual_features = ",".join(str(x) for x in qualVec[i])
        
        # 4. 添加TRACE_ACGT特征
        trace_features = ",".join(trace_data[i])
        
        # 将所有特征组合成一行
        fr.write(f"{label},{seq_features},{qual_features},{trace_features}\n")
    
    fr.close()

def main():
    # 读取81和82的数据
    resultdf1 = pd.read_csv("81_1_5mers_with_1_A_features.tsv", sep="\t")
    resultdf2 = pd.read_csv("82_1_5mers_with_1_A_features.tsv", sep="\t")
    
    # 分别处理两个文件
    seqVec1, qualVec1 = onehotencoding(resultdf1)
    seqVec2, qualVec2 = onehotencoding(resultdf2)
    
    # 保存处理后的文件
    savefile(resultdf1, qualVec1, seqVec1, 'feature_with_1_A_81')
    savefile(resultdf2, qualVec2, seqVec2, 'feature_with_1_A_82')
    
    # 合并两个文件
    resultdf = combind2files('feature_with_1_A_81.csv', 'feature_with_1_A_82.csv')
    resultdf.to_csv('feature_with_1_A_81_82.csv', sep=',', index=False, header=False)

if __name__=='__main__':
    main()
