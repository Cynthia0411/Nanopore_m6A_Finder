import re
import numpy
import os
import pandas

#start pos start in number1,not0.
kmerLength = 5
#filedir = 'D:/old/data/bioinformatics/VSG_mapped_reads_with_polyA/7656646f-042e-449a-aa0e-747bf145e6be.fast5'
inputsamfile="./mapped_polyA_reads_genome.sam"
inputsam_cdsfile = "./mapped_polyA_reads_CDs.sam"
#inputsamfile=str(snakemake.input)
outputtsvfile= "./polyA_pos.tsv"
#outputtsvfile=str(snakemake.output)

def get_read_ref_dict(inputsam_cds):
    with open(inputsam_cds,"r")as file_to_read:
        readname_ref = {}
        for line in file_to_read:
            line = line.strip()
            line_list = line.split("\t")
            readname = line_list[0]
            ref = line_list[2]
            readname_ref[readname] = ref
    return readname_ref
def getStart_End(inputsam,dict_read_ref):
    if os.path.exists(outputtsvfile):
        os.remove(outputtsvfile)
    with open(inputsam,'r') as file_to_read:
        with open(outputtsvfile,'a') as file_to_write:
            file_to_write.write('read_name'+'\t'+'reference'+'\t'+'polyA_pos_start'+'\t'+'polyA_pos_end'+'\n')
            for line in file_to_read:
                line = line.strip()
                line_list=line.split('\t')
                seq_name = line_list[0]
                if seq_name in dict_read_ref.keys():
                    reference = dict_read_ref[seq_name]
                    seq_flag = line_list[1]
                    cigar = line_list[5]
                    seq = line_list[9]
                    if seq_flag == "0":
                    #    print("flag =" + seq_flag)
                        if bool(re.match('.*?\d+S$',cigar)):
                     #       print(cigar)
                            pos_end = list(re.finditer('A{10}', seq))[-1].end()  # end_pos not included
                            pos_start = len(seq) - int(list(re.finditer('(\d+)S', cigar))[-1].group(1)) - 1     # start_pos included
                            if int(pos_end) > pos_start + 5:
                                file_to_write.write(seq_name + '\t' + reference + '\t' + str(pos_start) + '\t' + str(pos_end) + '\n')
                            else:
                                continue
                        else:
                            continue
                    else:
                    #    print("flag =" + seq_flag)
                        if bool(re.match('^\d+S',cigar)):
                        #    print(cigar)
                            pos_end = len(seq)-list(re.finditer('T{10}', seq))[0].start()  # end_pos not included
                            pos_start = len(seq) - int(list(re.finditer('(\d+)S', cigar))[0].group(1))-1
                            if int(pos_end) > pos_start +5:
                                file_to_write.write(seq_name + '\t' + reference + '\t' + str(pos_start) + '\t' + str(pos_end) + '\n')
                            else:
                                continue
                        else:
                            continue
                else:
                    continue

    file_to_read.close()
    file_to_write.close()
    return

def main():
    dict_cds = get_read_ref_dict(inputsam_cdsfile)
    getStart_End(inputsamfile,dict_cds)
    return

if __name__ == '__main__':
    main()
