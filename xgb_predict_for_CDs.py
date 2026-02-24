import pickle
import os
from numpy import loadtxt
import pandas as pd
import numpy as np

# 输入输出路径设置
feature_dir = './CDs_reads_feature_CDsReads_HEK293T_WT_rep3'  # 请修改为您的feature文件夹路径
outputFile = './m6A_pos_in_reads_xgb_HEK293T_WT_rep3.csv'
xgb_dat = './XGB_classifier_model.pkl'
dim = 45

def get_consecutive_positions(predictions):
    """
    Convert binary predictions to list of start-end positions
    Each position should be a 5bp window: End - Start = 5
    Position intervals are half-open: [start, end)
    """
    positions = []
    for i, pred in enumerate(predictions):
        if pred == 1:
            # 每个位置都是独立的5bp窗口
            positions.append((i, i + 5))
    return positions

def process_files():
    """Process files using XGBoost model only"""
    # 加载模型
    print("Loading XGBoost model...")
    xgb_model = pickle.load(open(xgb_dat, 'rb'))
    
    # 获取所有feature文件
    print("Processing feature files...")
    result_dict = {}
    feature_files = [f for f in os.listdir(feature_dir) if f.endswith('.csv')]
    total_files = len(feature_files)
    
    for i, feature_file in enumerate(feature_files, 1):
        if i % 100 == 0:  # 打印进度
            print(f"Processing file {i}/{total_files}")
            
        try:
            # 读取feature文件
            file_path = os.path.join(feature_dir, feature_file)
            cur_feature = loadtxt(file_path, delimiter=',')
            
            # 确保feature是2D数组
            if len(cur_feature.shape) == 1:
                cur_feature = cur_feature.reshape(1, -1)
            
            # 检查feature维度
            if cur_feature.shape[1] != dim:
                print(f"Warning: Skipping {feature_file} - incorrect feature dimension")
                continue
            
            # 获取XGBoost预测
            xgb_probs = xgb_model.predict_proba(cur_feature)[:, 1]
            predictions = [1 if prob > 0.5 else 0 for prob in xgb_probs]
            
            # 获取正预测的位置
            cur_readName = feature_file[:-4]  # 移除.csv后缀
            positions = get_consecutive_positions(predictions)
            if positions:  # 只添加有正预测的reads
                result_dict[cur_readName] = positions
                
        except Exception as e:
            print(f"Error processing {feature_file}: {str(e)}")
            continue
    
    return result_dict

def main():
    # 处理文件
    print("Processing files with XGBoost model...")
    result = process_files()
    
    # 将结果转换为DataFrame并写入文件
    rows = []
    for read, positions in result.items():
        for start, end in positions:
            rows.append({
                'read_name': read,
                'readpos_start': start,
                'readpos_end': end
            })
    df = pd.DataFrame(rows)
    df.to_csv(outputFile, index=True)  # index=True 保留索引列
    
    print("Processing complete!")

if __name__ == '__main__':
    main()

