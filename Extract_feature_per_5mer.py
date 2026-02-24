import pandas as pd
import h5py
import numpy as np
import os

def processFile(filedir, read_positions):
    """Process a single fast5 file and extract features for specified positions"""
    f_to_open = h5py.File(filedir)
    
    # Get raw data from fast5
    fastq_str = str(f_to_open['Analyses']['Basecall_1D_000']['BaseCalled_template']['Fastq'][()]).split('\\n')
    read_name = fastq_str[0][3:39]
    read_base = fastq_str[1]
    read_qual = [str(int(ord(obj) - 33)) for obj in list(fastq_str[3])]
    
    # Get trace and move data
    tra_feature = np.array(f_to_open['Analyses']['Basecall_1D_000']['BaseCalled_template']['Trace'])
    move = f_to_open['Analyses']['Basecall_1D_000']['BaseCalled_template']['Move']
    
    # Get indices where moves occur
    index_move = []
    for i in range(len(move)):
        if move[i] == 1:
            index_move.append(i)
    
    # Reverse the indices to match sequence positions
    index_seq = list(reversed(index_move))
    
    results = []
    for pos in read_positions:
        if pos + 5 <= len(read_base):
            # Get 5-mer sequence
            kmer = read_base[pos:pos+5]
            print(kmer + "\n")
#            if kmer[2:4] != "AC":
 #               continue
            
            # Get quality scores
            qual = ','.join(read_qual[pos:pos+5])
            
            
            # Get trace features
            tra_cur = []
            for n in range(pos, pos+5):
                if n < len(index_seq):
                    tra_l = (tra_feature[index_seq[n]][0:4] + tra_feature[index_seq[n]][4:]) / 255
                    tra_l = [str(round(obj, 2)) for obj in tra_l]
                    tra_cur.extend(tra_l)
            
            trace = ','.join(tra_cur) if len(tra_cur) == 20 else None  # 4 values per position * 5 positions
            
            if trace is not None:
                results.append({
                    'label': 1,  # Default label as -1
                    'BASE': kmer,
                    'QUAL': qual,
                    'TRACE_ACGT': trace
                })
    
    f_to_open.close()
    return results

def main():
    # Read the input TSV file
    input_data = pd.read_csv('AAAAA_level_70_pos_HEK293T-WT-rep3.tsv', sep='\t')
    
    # Group by READ_NAME to process each fast5 file once
    grouped_data = input_data.groupby('ReadName')
    
    results = []
    for read_name, group in grouped_data:
        fast5_path = f"./HEK293T-WT-rep3-fast5/{read_name}.fast5"  # Adjust path as needed
        
        if os.path.exists(fast5_path):
            read_positions = group['Position'].tolist()
            file_results = processFile(fast5_path, read_positions)
            results.extend(file_results)
        else:
            print(f"File not found: {fast5_path}")
    
    # Create output DataFrame
    if results:
        output_df = pd.DataFrame(results)
        # Reorder columns to match desired format
        output_df = output_df[['label', 'BASE', 'QUAL', 'TRACE_ACGT']]
        output_df.to_csv('AAAAA_level_70_HEK293T-WT-rep3.tsv', sep='\t', index=False)
    else:
        print("No results generated")

if __name__ == '__main__':
    main()


