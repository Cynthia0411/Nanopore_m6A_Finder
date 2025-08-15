import numpy as np
import hnswlib
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix, roc_curve, auc
import time
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interp
import os

class HNSWClassifier:
    def __init__(self, space='l2', dim=45):
        self.space = space
        self.dim = dim
        self.index = None
        self.scaler = StandardScaler()
        self.labels = None
        
    def fit(self, X, y, ef_construction=400, M=32):
        """
        训练HNSW模型
        
        参数:
            X: 训练特征
            y: 训练标签
            ef_construction: 构建索引时的搜索深度
            M: 每个节点的最大连接数
        """
        # 标准化特征
        X_scaled = self.scaler.fit_transform(X)
        
        # 初始化索引
        self.index = hnswlib.Index(space=self.space, dim=self.dim)
        
        # 构建索引
        self.index.init_index(
            max_elements=len(X),
            ef_construction=ef_construction,
            M=M
        )
        
        # 添加数据点
        self.index.add_items(X_scaled, np.arange(len(X)))
        
        # 保存标签
        self.labels = y
        
    def predict(self, X, k=5, ef_search=50):
        """
        预测新样本的标签
        
        参数:
            X: 待预测特征
            k: k近邻数量
            ef_search: 搜索时的候选集大小
        """
        X_scaled = self.scaler.transform(X)
        self.index.set_ef(ef_search)
        
        # 查找k近邻
        labels, distances = self.index.knn_query(X_scaled, k=k)
        
        # 对每个查询点进行多数投票
        predictions = []
        for neighbors in labels:
            neighbor_labels = self.labels[neighbors]
            # 计算正例（1）的比例
            pos_ratio = np.mean(neighbor_labels == 1)
            # 如果正例比例大于0.5，预测为1，否则为-1
            pred = 1 if pos_ratio > 0.5 else -1
            predictions.append(pred)
            
        return np.array(predictions)

def load_and_preprocess_data(file_path):
    """加载和预处理数据"""
    data = np.loadtxt(file_path, delimiter=",")
    y = data[:, 0]  # 第一列为标签
    X = data[:, 1:]  # 剩余45列为特征
    
    print(f"数据集大小: {X.shape}")
    print(f"正例数量: {np.sum(y == 1)}")
    print(f"负例数量: {np.sum(y == -1)}")
    
    return X, y

def plot_confusion_matrix(y_true, y_pred):
    """绘制混淆矩阵"""
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Negative (-1)', 'Positive (1)'],
                yticklabels=['Negative (-1)', 'Positive (1)'])
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.savefig('confusion_matrix_ANN.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 计算并打印分类指标
    tn, fp, fn, tp = cm.ravel()
    print("\nDetailed Classification Metrics:")
    print(f"True Negatives: {tn}")
    print(f"False Positives: {fp}")
    print(f"False Negatives: {fn}")
    print(f"True Positives: {tp}")

def plot_smooth_roc(y_true, y_scores, n_bootstraps=1000):
    """
    绘制平滑的ROC曲线，包括置信区间
    """
    # 将标签转换为0/1格式
    y_true_binary = (y_true + 1) / 2
    
    # 计算原始ROC曲线
    fpr, tpr, _ = roc_curve(y_true_binary, y_scores)
    roc_auc = auc(fpr, tpr)
    
    # Bootstrap数组
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)  # 用于插值的x轴点
    
    # Bootstrap采样计算置信区间
    rng = np.random.RandomState(42)
    for i in range(n_bootstraps):
        # 有放回抽样
        indices = rng.randint(0, len(y_true_binary), len(y_true_binary))
        if len(np.unique(y_true_binary[indices])) < 2:
            continue
            
        fpr, tpr, _ = roc_curve(y_true_binary[indices], y_scores[indices])
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
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=f'Mean ROC (AUC = {mean_auc:.3f} ± {std_auc:.3f})',
             lw=2, alpha=.8)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'± 1 std. dev.')
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Random Chance', alpha=.8)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.savefig('roc_curve_ANN.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return mean_auc, std_auc

def evaluate_model(y_true, y_pred, y_scores=None):
    """评估模型性能"""
    # 计算基本指标
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, pos_label=1)
    recall = recall_score(y_true, y_pred, pos_label=1)
    f1 = f1_score(y_true, y_pred, pos_label=1)
    
    print("\nModel Performance:")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"F1-score: {f1:.4f}")
    
    # 绘制混淆矩阵
    plot_confusion_matrix(y_true, y_pred)
    
    # 如果提供了预测概率/得分，绘制ROC曲线
    if y_scores is not None:
        mean_auc, std_auc = plot_smooth_roc(y_true, y_scores)
        print(f"\nROC AUC: {mean_auc:.4f} ± {std_auc:.4f}")
    
    return accuracy, precision, recall, f1

def save_model_and_labels(model, y_train, filepath_prefix):
    """保存模型和训练标签"""
    # 保存模型
    model_path = f"{filepath_prefix}_model.pkl"
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)
    
    # 保存训练标签
    labels_path = f"{filepath_prefix}_labels.npy"
    np.save(labels_path, y_train)
    
    print(f"Model saved to {model_path}")
    print(f"Training labels saved to {labels_path}")
    
def load_model_and_labels(filepath_prefix):
    """加载模型和训练标签"""
    # 加载模型
    model_path = f"{filepath_prefix}_model.pkl"
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    
    # 加载训练标签
    labels_path = f"{filepath_prefix}_labels.npy"
    y_train = np.load(labels_path)
    
    print(f"Model loaded from {model_path}")
    print(f"Training labels loaded from {labels_path}")
    return model, y_train

def main():
    # 加载数据
    X, y = load_and_preprocess_data("balanced_feature_all_within_DRACH.csv")
    
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=7, stratify=y
    )
    
    # 初始化和训练模型
    print("\nTraining HNSW model...")
    start_time = time.time()
    
    model = HNSWClassifier(space='l2', dim=45)
    model.fit(X_train, y_train, ef_construction=400, M=32)
    
    training_time = time.time() - start_time
    print(f"Training completed in {training_time:.2f} seconds")
    
    # 创建输出目录
    output_dir = 'ANN_model_evaluation_plots_all_within_DRACH_k_5'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 保存所有图形到指定目录
    os.chdir(output_dir)
    

    # 预测和评估
    print("\nMaking predictions...")
    y_pred = model.predict(X_test, k=5, ef_search=50)
    
    # 获取近邻投票的正例比例作为预测得分
    y_scores = []
    for neighbors in model.index.knn_query(model.scaler.transform(X_test), k=5)[0]:
        neighbor_labels = model.labels[neighbors]
        pos_ratio = np.mean(neighbor_labels == 1)
        y_scores.append(pos_ratio)
    y_scores = np.array(y_scores)
    
    # 评估模型（包括ROC曲线）
    metrics = evaluate_model(y_test, y_pred, y_scores)
    
    # 返回到原始目录
    os.chdir('..')
    
    # 保存模型和训练标签
    save_model_and_labels(model, y_train, "hnsw_classifier_all_DRACH")
    
    # 测试加载模型和标签
    loaded_model, loaded_y_train = load_model_and_labels("hnsw_classifier_all_DRACH")
    y_pred_loaded = loaded_model.predict(X_test, k=5, ef_search=50)
    
    # 验证加载的模型
    print("\nVerifying loaded model predictions...")
    metrics_loaded = evaluate_model(y_test, y_pred_loaded)

if __name__ == "__main__":
    main() 
