import numpy as np
import pickle
import os
from numpy import loadtxt
import pandas as pd

# 输入输出路径设置
inputDir = './polyA_reads_feature_59/'
outputFile = './m6A_polyA_per_gene_strict_with_length_TrB59_new.csv'
xgb_dat = './XGB_classifier_model_all_withA_balanced.pkl'  # 转换后的XGBoost classifier模型
rf_dat = './rf_model_balanced_all_withA.pkl'    # Random Forest模型
ref_file_dir = './polyA_pos_59.tsv'

def get_prediction_from_proba(probas, threshold=0.5):
    """
    根据预测概率确定类别
    Args:
        probas: shape (n_samples, 2) 的概率矩阵，第一列是负例概率，第二列是正例概率
        threshold: 正例概率阈值,默认0.5
    Returns:
        predictions: 预测标签数组
    """
    neg_probs = probas[:, 0]
    pos_probs = probas[:, 1]
    
    predictions = np.where(
        (pos_probs > neg_probs) & (pos_probs > threshold),
        1,
        0
    )
    return predictions

def main():
    # 加载模型
    xgb_model = pickle.load(open(xgb_dat, 'rb'))
    rf_model = pickle.load(open(rf_dat, 'rb'))
    
    filenames = os.listdir(inputDir)
    result = []
    result_count = []
    readName = []
    ref = []
    polyA_length = []
    
    # 读取参考文件
    df_ref = pd.read_csv(ref_file_dir, sep='\t')
    df_ref_index = df_ref.index.values
    read_ref = {}
    for i in df_ref_index:
        read_ref[df_ref['read_name'][i]] = df_ref['reference'][i]
    
    for files in filenames:
        filedir = inputDir + '/' + files
        cur_read = files.split(".")[0].split("_")[0]
        cur_feature = loadtxt(filedir, delimiter=',')
        
        # 确保特征是2D数组
        if np.shape(cur_feature) == (45,):
            cur_feature = np.reshape(cur_feature, (1, 45))
            
        # 获取两个模型的预测概率和预测结果
        xgb_probas = xgb_model.predict_proba(cur_feature)
        rf_probas = rf_model.predict_proba(cur_feature)
        
        xgb_predictions = get_prediction_from_proba(xgb_probas)
        rf_predictions = get_prediction_from_proba(rf_probas)
        
        # 硬投票整合结果
        counter = 0
        counter_list = []
        cur_result = 0
        
        for k in range(len(xgb_predictions)):
            # 只有两个模型都预测为正例时才认为是正例
            if xgb_predictions[k] == 1 and rf_predictions[k] == 1:
                cur_result += 1
                counter_list.append(1)
            else:
                counter_list.append(0)
            counter += 1
        
        cur_readName = files.split('.')[0]
        result_count.append(cur_result)
        result.append(cur_result/counter if counter > 0 else 0)
        readName.append(cur_read)
        ref.append(read_ref[cur_readName])
        polyA_length.append(counter_list)
    
    # 创建输出DataFrame
    polyALevelTable = pd.DataFrame({
        'read_name': readName,
        'reference': ref,
        'polyA_level': result,
        'polyA_5mer_sum': result_count,
        'polyA_length': polyA_length
    })
    
    # 保存结果
    polyALevelTable.to_csv(outputFile, sep=',')
    print(f"Results saved to: {outputFile}")

if __name__ == '__main__':
    main()


