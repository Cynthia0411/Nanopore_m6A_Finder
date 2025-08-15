import numpy as np
import matplotlib.pyplot as plt
import pickle
from sklearn.metrics import roc_curve, auc, precision_score, recall_score, f1_score, accuracy_score
from sklearn.model_selection import train_test_split
import warnings
import pandas as pd
import seaborn as sns
from sklearn.metrics import confusion_matrix

# 修改模型文件路径
xgb_model_path = './xgb_model_with_1_A_balanced.pkl'
rf_model_path = './rf_model_balanced_with_1_A.pkl'

def load_models():
    """加载XGBoost和Random Forest模型"""
    print("Loading models...")
    
    # 加载XGBoost模型
    with open(xgb_model_path, 'rb') as f:
        xgb_model = pickle.load(f)
    print("XGBoost model loaded")
    print(f"XGBoost model type: {type(xgb_model)}")  # 打印模型类型
    
    # 如果需要，将Booster转换为XGBClassifier
    if not hasattr(xgb_model, 'predict_proba'):
        import xgboost as xgb
        classifier = xgb.XGBClassifier()
        classifier._Booster = xgb_model
        classifier._le = None
        xgb_model = classifier
    
    # 加载Random Forest模型
    with open(rf_model_path, 'rb') as f:
        rf_model = pickle.load(f)
    print("Random Forest model loaded")
    
    return xgb_model, rf_model

def evaluate_model(model_name, y_pred, y_test):
    """评估模型并打印混淆矩阵"""
    # 将概率转换为二进制预测
    y_pred_binary = (y_pred > 0.5).astype(int)
    
    # 计算混淆矩阵
    cm = confusion_matrix(y_test, y_pred_binary)
    
    # 创建混淆矩阵的热图
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Negative', 'Positive'],
                yticklabels=['Negative', 'Positive'])
    plt.title(f'Confusion Matrix - {model_name}')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    
    # 保存混淆矩阵图
    plt.savefig(f'DRACH_confusion_matrix_{model_name.lower().replace(" ", "_")}.pdf', 
                bbox_inches='tight', dpi=300)
    plt.close()
    
    # 打印混淆矩阵的值
    print(f"\nConfusion Matrix for {model_name}:")
    print(f"True Negatives: {cm[0,0]}")
    print(f"False Positives: {cm[0,1]}")
    print(f"False Negatives: {cm[1,0]}")
    print(f"True Positives: {cm[1,1]}")

def get_predictions(xgb_model, rf_model, X_test, batch_size=1000):
    """获取两个模型的预测结果"""
    print("Getting predictions...")
    
    # XGBoost预测
    xgb_pred = xgb_model.predict_proba(X_test)[:, 1]
    print("XGBoost predictions completed")
    
    # Random Forest批量预测
    n_samples = len(X_test)
    rf_pred = np.zeros(n_samples)
    for i in range(0, n_samples, batch_size):
        end_idx = min(i + batch_size, n_samples)
        batch_X = X_test[i:end_idx]
        rf_pred[i:end_idx] = rf_model.predict_proba(batch_X)[:, 1]
        if (i + batch_size) % 5000 == 0:
            print(f"RF: Processed {i + batch_size}/{n_samples} samples...")
    print("Random Forest predictions completed")
    
    return xgb_pred, rf_pred

def main():
    try:
        # 加载数据
        print("Loading data...")
        df = pd.read_csv("balanced_feature_with_1_A.csv", header=None)
        x = df.iloc[:, 1:].values
        y = df.iloc[:, 0].map(lambda x: 1 if x == 1 else 0).values
        
        print("\nInitial data check:")
        print(f"Total samples: {len(y)}")
        print(f"Class distribution: {np.bincount(y)}")
        
        # 分割数据
        seed = 7
        testSize = 0.2
        X_train, X_test, y_train, y_test = train_test_split(x, y, 
                                                           test_size=testSize, 
                                                           random_state=seed,
                                                           stratify=y)
        
        # 加载模型
        xgb_model, rf_model = load_models()
        
        # 获取预测结果
        xgb_pred, rf_pred = get_predictions(xgb_model, rf_model, X_test)
        
        # 计算加权投票和硬投票
        weight_1 = 0.5
        weighted_probs = weight_1 * xgb_pred + (1-weight_1) * rf_pred
        
        hard_vote_probs = np.zeros_like(xgb_pred)
        for i in range(len(xgb_pred)):
            if xgb_pred[i] > 0.5 and rf_pred[i] > 0.5:
                hard_vote_probs[i] = 1.0
            elif xgb_pred[i] <= 0.5 and rf_pred[i] <= 0.5:
                hard_vote_probs[i] = 0.0
            else:
                hard_vote_probs[i] = 0.5
        
        # 定义模型列表
        models = [
            (xgb_pred, 'XGBoost', 'blue'),
            (rf_pred, 'Random Forest', 'green'),
            (weighted_probs, 'Weighted Vote', 'purple'),
            (hard_vote_probs, 'Hard Vote', 'orange')
        ]
        
        # 创建ROC曲线图
        plt.figure(figsize=(8, 8))
        plt.plot([0, 1], [0, 1], '--', color='red', alpha=0.8, linewidth=1)
        
        # 存储结果
        results = []
        
        # 评估每个模型
        for y_pred, model_name, color in models:
            print(f"\nEvaluating {model_name}...")
            
            if np.all(np.isfinite(y_pred)) and len(np.unique(y_test)) > 1:
                # ROC曲线
                fpr, tpr, _ = roc_curve(y_test, y_pred)
                roc_auc = auc(fpr, tpr)
                
                plt.plot(fpr, tpr,
                        color=color,
                        linewidth=2,
                        label=f'{model_name} (AUC = {roc_auc:.2f})')
                
                # 混淆矩阵
                evaluate_model(model_name, y_pred, y_test)
                
                # 计算指标
                y_pred_binary = (y_pred > 0.5).astype(int)
                metrics = {
                    'Precision': precision_score(y_test, y_pred_binary),
                    'Recall': recall_score(y_test, y_pred_binary),
                    'F1': f1_score(y_test, y_pred_binary),
                    'Accuracy': accuracy_score(y_test, y_pred_binary)
                }
            else:
                print(f"Warning: Invalid predictions for {model_name}")
                roc_auc = 0
                metrics = {'Precision': 0, 'Recall': 0, 'F1': 0, 'Accuracy': 0}
            
            results.append({
                'Model': model_name,
                'AUC': roc_auc,
                **metrics
            })
        
        # 保存结果
        df_results = pd.DataFrame(results)
        print("\nResults:")
        print(df_results)
        df_results.to_csv('4_model_result_balanced_with_1_A.tsv', sep='\t', index=False)
        
        # 自定义ROC曲线图
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.xlabel('False Positive Rate', fontsize=12)
        plt.ylabel('True Positive Rate', fontsize=12)
        plt.title('Model Comparison', fontsize=14)
        plt.legend(loc='lower right', fontsize=10)
        plt.grid(True, alpha=0.15)
        
        plt.savefig('four_roc_curves_balanced_with_1_A.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("\nDone!")
        
    except Exception as e:
        print(f"Error in main function: {str(e)}")
        import traceback
        traceback.print_exc()
        raise

if __name__ == "__main__":
    main()
