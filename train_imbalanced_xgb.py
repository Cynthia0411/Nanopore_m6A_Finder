import numpy as np
import xgboost as xgb
from xgboost import plot_importance
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import roc_auc_score, precision_recall_curve, average_precision_score, roc_curve, auc
import matplotlib.pyplot as plt
import pickle
from scipy import interp
import os
import seaborn as sns
from sklearn.metrics import confusion_matrix

def load_and_preprocess_data(file_path):
    """加载和预处理数据"""
    data = np.loadtxt(file_path, delimiter=",")
    X = data[:, 1:]  # 从第二列开始是45维特征
    y = data[:, 0]   # 第一列是标签
    
    # 将-1标签转换为0
    y = np.where(y == -1, 0, y)
    
    return X, y


def train_xgb_model(X, y):
    """训练XGBoost模型"""
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=7, stratify=y
    )
    
    # 计算类别比例
    neg_pos_ratio = np.sum(y_train == 0) / np.sum(y_train == 1)
    print(f"类别比例 (负/正): {neg_pos_ratio:.2f}")
    print(f"负样本数量: {np.sum(y_train == 0)}")
    print(f"正样本数量: {np.sum(y_train == 1)}")

    # 定义模型参数
    params = {
        'objective': 'binary:logistic',
        'eval_metric': ['auc', 'logloss'],
        'learning_rate': 0.1,
        'max_depth': 8,
        'min_child_weight': 1,
        'subsample': 0.8,
        'colsample_bytree': 0.7,
        'scale_pos_weight': neg_pos_ratio,  # 处理类别不平衡
        'tree_method': 'hist',  # 使用直方图优化算法
        'random_state': 7
    }
    
    # 创建DMatrix数据格式
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_test, label=y_test)
    
    # 设置评估列表
    evallist = [(dtrain, 'train'), (dtest, 'eval')]
    
    # 训练模型
    num_round = 1000
    model = xgb.train(
        params,
        dtrain,
        num_round,
        evallist,
        early_stopping_rounds=50,
        verbose_eval=100
    )
    
    return model, X_test, y_test, dtest

def plot_smooth_roc(y_true, y_pred_proba, n_bootstraps=1000):
    """
    绘制并保存平滑的ROC曲线包括置信区间
    """
    # 确保标签是整数类型
    y_true = y_true.astype(int)
    
    # 计算原始ROC曲线
    fpr, tpr, _ = roc_curve(y_true, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    
    # Bootstrap数组
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    
    # Bootstrap采样计算置信区间
    rng = np.random.RandomState(42)
    for i in range(n_bootstraps):
        indices = rng.randint(0, len(y_true), len(y_true))
        if len(np.unique(y_true[indices])) < 2:
            continue
        
        fpr, tpr, _ = roc_curve(y_true[indices], y_pred_proba[indices])
        tprs.append(np.interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        aucs.append(auc(fpr, tpr))
    
    # 计算平均值和置信区间
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    
    # 计算TPR的置信区间
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    
    # 绘图
    plt.figure(figsize=(10, 8))
    
    # 绘制平均ROC曲线
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=f'Mean ROC (AUC = {mean_auc:.3f} ± {std_auc:.3f})',
             lw=2, alpha=.8)
    
    # 绘制置信区间
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'± 1 std. dev.')
    
    # 绘制随机猜测的基准线
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Random Chance', alpha=.8)
    
    # 设置图形属性
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    plt.grid(True)
    
    # 保存图形
    plt.savefig('roc_curve_xgb.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return mean_auc, std_auc

def plot_confusion_matrix(y_true, y_pred):
    """绘制混淆矩阵热图"""
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Negative (0)', 'Positive (1)'],
                yticklabels=['Negative (0)', 'Positive (1)'])
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.savefig('confusion_matrix_xgb.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 计算并打印分类指标
    tn, fp, fn, tp = cm.ravel()
    print("\nDetailed Classification Metrics:")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Positives: {tp}")

def evaluate_model(model, X_test, y_test, dtest):
    """评估模型性能"""
    # 预测概率
    y_pred_proba = model.predict(dtest)
    print(y_pred_proba)
    # 将预测概率转换为类别标签（0或1）
    y_pred = (y_pred_proba > 0.5).astype(int)
    
    # 确保标签是整数类型
    y_test = y_test.astype(int)
    
    # 计算各项指标
    auc_score = roc_auc_score(y_test, y_pred_proba)
    avg_precision = average_precision_score(y_test, y_pred_proba)
    
    print(f"\nModel Performance:")
    print(f"AUC Score: {auc_score:.4f}")
    print(f"Average Precision Score: {avg_precision:.4f}")
    
    # 绘制混淆矩阵
    plot_confusion_matrix(y_test, y_pred)
    
    # 绘制并保存平滑ROC曲线
    mean_auc, std_auc = plot_smooth_roc(y_test, y_pred_proba)
    print(f"Bootstrap AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    
    # 绘制并保存PR曲线
    plt.figure(figsize=(10, 6))
    precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)
    plt.plot(recall, precision, label=f'PR curve (AP = {avg_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend()
    plt.grid(True)
    plt.savefig('pr_curve_xgb.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 绘制并保存特征重要性图
    plt.figure(figsize=(12, 6))
    plot_importance(model, max_num_features=20)
    plt.title('Feature Importance')
    plt.savefig('feature_importance_xgb.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return auc_score, avg_precision

def save_model(model, filepath):
    """使用pickle保存模型"""
    with open(filepath, 'wb') as f:
        pickle.dump(model, f)
    print(f"Model saved to {filepath}")

def load_model(filepath):
    """加载pickle格式的模型"""
    with open(filepath, 'rb') as f:
        model = pickle.load(f)
    print(f"Model loaded from {filepath}")
    return model

def main():
    # 加载数据
    X, y = load_and_preprocess_data("balanced_feature_with_1_A.csv")
    
    # 训练模型
    model, X_test, y_test, dtest = train_xgb_model(X, y)
    
    # 创建输出目录
    output_dir = 'xgb_model_evaluation_plots_with_1_A'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 保存所有图形到指定目录
    os.chdir(output_dir)
    
    # 评估模型并保存图形
    metrics = evaluate_model(model, X_test, y_test, dtest)
    
    # 返回到原始目录
    os.chdir('..')
    
    # 保存模型（使用pickle格式）
    save_model(model, "xgb_model_with_1_A_balanced.pkl")
    
    # 测试加载模型
    loaded_model = load_model("xgb_model_with_1_A_balanced.pkl")
    
    # 验证加载的模型预测结果
    test_pred_proba = loaded_model.predict(dtest)
    print("y_test:")
    print(y_test)
    test_pred = (test_pred_proba > 0.5).astype(int)
    y_test = y_test.astype(int)
    print("\nVerifying loaded model predictions...")
    print(f"AUC score with loaded model: {roc_auc_score(y_test, test_pred_proba):.4f}")

if __name__ == "__main__":
    main()
