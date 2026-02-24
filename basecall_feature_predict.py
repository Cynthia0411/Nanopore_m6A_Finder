import pandas as pd
import numpy as np
import os

#inputnameListdir = "./m6A_reads"
inputnameListdir = "./CDsReads_HEK293T_Mettl3_KO_rep1/"
#output_dir = "./m6A_reads_feature"
output_dir = "./CDs_reads_feature_CDsReads_HEK293T_Mettl3_KO_rep1"
#inputposdir = "./m6A_sites_coverage"
#output_file = "./m6A_ref_pos.tsv"



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
    f1 = pd.read_csv(filename1,sep='\t')
    f1 = f1.sample(frac=1)
    f2 = pd.read_csv(filename2,sep='\t')
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
#    qualDic = {"low": [1, 0, 0, 0], "middle": [0, 1, 0, 0], "high": [0, 0, 1, 0], "extremHigh": [0, 0, 0, 1]}
    qualDic = {"extremeLow": 0.1,'low':0.4,"middle": 0.6, "high": 0.7, "extremeHigh": 1}
    seqVec = np.zeros((len(seq), 20))
    for i in range(len(seq)):
        seqList = list(seq[i])
        lineVector = []
        for j in range(5):
 #           print(seqList[j])
            lineVector.extend(seqDict[seqList[j]])
        seqVec[i] = lineVector
    qualVec = np.zeros((len(qualScore), 5))
    for i in range(len(qualScore)):
        qualList = qualScore[i].split(sep = ',')
        lineVector = []
        for j in range(5):
            if int(qualList[j]) <= 5:
                lineVector.append(qualDic["extremeLow"])
            elif 5 < int(qualList[j]) <= 8:
                lineVector.append(qualDic["low"])
            elif 8 < int(qualList[j]) <= 25:
                lineVector.append(qualDic["middle"])
            elif 25 < int(qualList[j]) <= 35:
                lineVector.append(qualDic["high"])
            elif int(qualList[j])>35:
                lineVector.append(qualDic["extremeHigh"])
        qualVec[i] = lineVector
    return seqVec,qualVec
def makeREF(filedir,outputfile):
    files = os.listdir(filedir)
    readnameList = {}
    for filename in files:
        portion = os.path.splitext(filename)
        if portion[1] == '.tsv':
            readNamedir = filedir + '/' + filename
            readnameList[readNamedir] = filename.split(sep=".")[0]
        else:
            continue
    with open(outputfile, "a") as file_to_write:
        file_to_write.write("readName"+"\t"+"ref" + "\t" + "pos" + "\n")
        for key in readnameList.keys():
            with open(key,"r")as file_to_read:
                for lines in file_to_read:
                    line = lines.strip()
                    readname = line.split("\t")[0]
                    ref = readnameList[key] +"_"+line.split("\t")[2]
                    pos = line.split("\t")[3]
                    file_to_write.write(readname + "\t" + ref + "\t" + pos + "\n")
            file_to_read.close()
    file_to_write.close()
    return
def savefile(resultdf,qualVec,seqVec,readname_site,outputdir):
    if os.path.exists(outputdir)==False:
        os.mkdir(outputdir)
    fr = open(outputdir+"/"+readname_site+".csv", "a")
#    fr = open('featureVec_7656646f-042e-449a-aa0e-747bf145e6be.csv','w')
#    for i in range(len(qualVec)):
#        fr.write(str(resultdf['label'][i])+','+
#            (",".join('%s' % id for id in seqVec[i])) + "," + (",".join('%s' % id for id in qualVec[i])) + "," + resultdf['TRACE_ACGT'][i] +"\n")
        
#        + ','+ str(resultdf['label'][i]) + "\n")
    for i in range(len(qualVec)):
        fr.write(
                (",".join('%s' % id for id in seqVec[i])) + "," + (",".join('%s' % id for id in qualVec[i])) + "," + resultdf['TRACE_ACGT'][i] + "\n")
#        fr.write((",".join('%s'%id for id in qualVec[i]))+","+str(classType[i])+"\n")

def main():
#    resultdf = combind2files('../select_reads/81_1_5mer_withA_trace.tsv','../select_reads/82_1_5mers_withA_trace.tsv')
#    seqVec , qualVec = onehotencoding(resultdf)
#    savefile(resultdf,qualVec,seqVec,'feature_81_82')
#    resultdf = pd.read_csv('7656646f-042e-449a-aa0e-747bf145e6be.tsv',sep = '\t',index_col=0)
#    resultdf,readname = getReadName("./polyAOfReads/000a016c-6503-46b7-b5d4-694c898161ad.tsv")
#    makeREF(inputposdir,output_file)    
    readNameList = getReadNameList(inputnameListdir)
    for filename in readNameList:
        resultdf,readname_site = getReadName(filename)
        seqVec , qualVec = onehotencoding(resultdf)
#
        savefile(resultdf,qualVec,seqVec,readname_site,output_dir)

    

if __name__=='__main__':
    main()
