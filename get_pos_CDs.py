import re
import numpy
import os
import pandas

#start pos start in number1,not0.
kmerLength = 5
#filedir = 'D:/old/data/bioinformatics/VSG_mapped_reads_with_polyA/7656646f-042e-449a-aa0e-747bf145e6be.fast5'
inputsamfile="./mapped_reads.sam"                             
#inputsamfile=str(snakemake.input)
outputtsvfile= "./CDs_pos.tsv"
#outputtsvfile=str(snakemake.output)
def getStart_End(inputsam):
    if os.path.exists(outputtsvfile):
        os.remove(outputtsvfile)
    with open(inputsam,'r') as file_to_read:
        with open(outputtsvfile,'a') as file_to_write:
            file_to_write.write('read_name'+'\t'+'reference'+'\t'+'polyA_pos_start'+'\t'+'polyA_pos_end'+'\n')
            for line in file_to_read:
                line = line.strip()
                line_list=line.split('\t')
                seq_name = line_list[0]
                flag = str(line_list[1])
                reference = line_list[2]
                seq = line_list[9]
#                pos_end = list(re.finditer('A{10}',seq))[-1].end()
                pos_start = 0
                cigar = line_list[5]
#                pos_start = len(seq) - int(list(re.finditer('(\d+)S',cigar))[-1].group(1))
                s_re = list(re.finditer('(\d+)S',cigar))
                if len(s_re) != 0:
                    if flag == "0":
                        pos_end = len(seq) - int(s_re[-1].group(1))-1
                        file_to_write.write(seq_name+'\t'+reference+'\t'+str(pos_start)+'\t'+str(pos_end)+'\n')
                    elif flag == "16":
                        pos_end = len(seq) - int(s_re[0].group(1))-1
                        file_to_write.write(seq_name + '\t' + reference + '\t' + str(pos_start) + '\t' + str(pos_end) + '\n')
    file_to_read.close()
    file_to_write.close()
    
getStart_End(inputsamfile)
