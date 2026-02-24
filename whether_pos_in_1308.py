import pandas as pd

def read_tsv_ref(filename):
    """读取参考tsv文件"""
    df = pd.read_csv(filename, sep='\t')
#    return df[['Chromosome', 'Start', 'End', 'Flag']]  # 只保留需要的列
    return df[['chrom', 'start', 'pos', 'strand']]

def read_bed(filename):
    """读取bed文件"""
    df = pd.read_csv(filename, sep='\t', 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand'],
                     header=None)
    return df

def find_nearby_positions(ref_df, query_bed, window_size=20):
    """找出query_bed中在ref_df位置附近window_size范围内的位置"""
    nearby_positions = []
    
    # 遍历query_bed中的每个位置
    for _, query_row in query_bed.iterrows():
        query_chrom = query_row['chrom']
        query_pos = query_row['end']  # 使用end位置作为位点位置
        
        # 在ref_df中寻找同一染色体上的位置
        ref_same_chrom = ref_df[ref_df['chrom'] == query_chrom]
        
        # 检查是否在任何参考位置的范围内
        for _, ref_row in ref_same_chrom.iterrows():
            ref_pos = ref_row['pos']
            
            # 检查是否在窗口范围内
            if abs(query_pos - ref_pos) <= window_size:
                nearby_positions.append(query_row)
                break  # 找到一个符合条件的就可以跳出内循环
    
    # 将结果转换为DataFrame
    if nearby_positions:
        result_df = pd.DataFrame(nearby_positions)
        return result_df
    else:
        return pd.DataFrame(columns=query_bed.columns)

def main():
    # 输入输出文件
    ref_tsv_file = "combined_m6A_sites_WT1_KO1_20bp_overlap.tsv"  # 第一个tsv参考文件
    query_bed_file = "predicted_m6A_sites_WT3_KO2.bed"    # 第二个bed文件
    output_file = "nearby_positions_overlap_20bp_WT3_KO2.bed"  # 输出文件
    
    # 读取输入文件
    print("Reading input files...")
    ref_df = read_tsv_ref(ref_tsv_file)
    query_bed = read_bed(query_bed_file)
    
    # 寻找附近的位置
    print("Finding nearby positions...")
    nearby_df = find_nearby_positions(ref_df, query_bed, window_size=20)
    
    # 保存结果
    print("Saving results...")
    nearby_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"Found {len(nearby_df)} positions within 20bp window")
    print("Done!")

if __name__ == "__main__":
    main()
