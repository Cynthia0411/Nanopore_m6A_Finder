import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path):
    # 读取数据，没有header，第一列是标签
    data = pd.read_csv(file_path, header=None)
    
    # 分离特征和标签
    X = data.iloc[:, 1:]  # 45维特征
    y = data.iloc[:, 0]   # 标签列
    
    return X, y

def train_and_evaluate(kernel='rbf', C=1.0, gamma='scale'):
    """
    参数说明：
    kernel: string, 默认='rbf'
        核函数类型：'linear', 'poly', 'rbf', 'sigmoid'
    C: float, 默认=1.0
        正则化参数，控制模型复杂度
    gamma: float or string, 默认='scale'
        'rbf', 'poly' 和 'sigmoid' 核函数的系数
    """
    
    # 加载数据
    X, y = load_data('feature_within_DRACH_81_82.csv')
    
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=7, stratify = y)
    
    # 特征标准化
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # 初始化并训练SVM模型
    svm_model = SVC(
        kernel=kernel,
        C=C,
        gamma=gamma,
        probability=True  # 启用概率估计
    )
    
    print("Training SVM model...")
    svm_model.fit(X_train_scaled, y_train)
    
    # 在测试集上进行预测
    y_pred = svm_model.predict(X_test_scaled)
    y_pred_proba = svm_model.predict_proba(X_test_scaled)[:, 1]
    
    # 打印模型参数和性能指标
    print("\nModel Parameters:")
    print(f"kernel: {kernel}")
    print(f"C: {C}")
    print(f"gamma: {gamma}")
    
    print("\nPerformance Metrics:")
    print("Accuracy:", accuracy_score(y_test, y_pred))
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))
    
    # 绘制混淆矩阵
    cm = confusion_matrix(y_test, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.savefig('svm_confusion_matrix_DRACH.png')
    plt.close()
    
    # 绘制ROC曲线
    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, 
             label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.savefig('svm_roc_curve_DRACH.png')
    plt.close()
    
    return svm_model, scaler

def main():
    try:
        print("Loading data and training SVM model...")
        
        # 可以在这里手动调整参数
        model, scaler = train_and_evaluate(
            kernel='rbf',     # 可选: 'linear', 'poly', 'rbf', 'sigmoid'
            C=1.0,           # 尝试不同的值，如 0.1, 1.0, 10.0
            gamma='scale'    # 尝试 'scale', 'auto' 或具体数值如 0.1, 0.01
        )
        
        # 保存模型和scaler
        import pickle
        with open('svm_model_DRACH.pkl', 'wb') as f:
            pickle.dump({'model': model, 'scaler': scaler}, f)
        print("Model and scaler saved as 'svm_model_DRACH.pkl'")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    main()
