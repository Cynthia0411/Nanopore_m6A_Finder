def process_fastq(input_file, output_file, min_length=1500):
    """
    筛选FASTQ文件中长度大于指定值的reads
    
    参数:
    input_file: 输入的FASTQ文件路径
    output_file: 输出的FASTQ文件路径
    min_length: 最小长度阈值默认1500bp
    """
    try:
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            # FASTQ文件每4行为一组
            while True:
                # 读取4行
                header = fin.readline().strip()
                if not header:  # 文件结束
                    break
                    
                sequence = fin.readline().strip()
                plus = fin.readline().strip()
                quality = fin.readline().strip()
                
                # 检查序列长度
                if len(sequence) >= min_length:
                    # 写入符合条件的read
                    fout.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
                    
        print(f"Filtering completed. Results saved to {output_file}")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def main():
    # 设置输入输出文件路径
    input_file = "input.fastq"  # 替换为你的输入文件路径
    output_file = "filtered_1500plus.fastq"  # 替换为你想要的输出文件路径
    
    print(f"Processing {input_file}...")
    process_fastq(input_file, output_file)

if __name__ == "__main__":
    main()
