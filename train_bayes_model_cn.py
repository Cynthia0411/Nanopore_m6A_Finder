import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import ComplementNB
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path):
    # 读取数据，没有header，第一列是标签
    data = pd.read_csv(file_path, header=None)
    
    # 分离特征和标签
    X = data.iloc[:, 1:]  # 45维特征
    y = data.iloc[:, 0]   # 标签列
    
    return X, y

def train_and_evaluate(alpha=1.0, norm=False, fit_prior=True):
    """
    参数说明：
    alpha: float, 默认=1.0
        加法（拉普拉斯/利德斯通）平滑参数（0表示不平滑）
    norm: boolean, 默认=False
        是否在训练后执行第二次规范化
    fit_prior: boolean, 默认=True
        是否学习类别先验概率，如果为False则使用统一先验
    """
    
    # 加载数据
    X, y = load_data('feature_within_DRACH_81_82.csv')
    
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=7)
    
    # 特征缩放到非负值（ComplementNB要求特征非负）
    scaler = MinMaxScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # 初始化并训练ComplementNB模型
    cnb_model = ComplementNB(
        alpha=alpha,
        norm=norm,
        fit_prior=fit_prior
    )
    cnb_model.fit(X_train_scaled, y_train)
    
    # 在测试集上进行预测
    y_pred = cnb_model.predict(X_test_scaled)
    y_pred_proba = cnb_model.predict_proba(X_test_scaled)[:, 1]
    
    # 计算并打印评估指标
    print("\nModel Parameters:")
    print(f"alpha: {alpha}")
    print(f"norm: {norm}")
    print(f"fit_prior: {fit_prior}")
    
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
    plt.savefig('complement_nb_confusion_matrix_DRACH.png')
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
    plt.savefig('complement_nb_roc_curve_DRACH.png')
    plt.close()
    
    return cnb_model, scaler

def main():
    try:
        print("Loading data and training ComplementNB model...")
        
        # 可以在这里手动调整参数
        model, scaler = train_and_evaluate(
            alpha=1.0,  # 尝试不同的值，如 0.1, 0.5, 1.0, 2.0
            norm=False, # 尝试 True/False
            fit_prior=True  # 尝试 True/False
        )
        
        # 保存模型和scaler
        import pickle
        with open('complement_nb_model.pkl', 'wb') as f:
            pickle.dump({'model': model, 'scaler': scaler}, f)
        print("Model and scaler saved as 'complement_nb_model.pkl'")
        
        with open('complement_nb_model.pkl', 'rb') as f:
            saved = pickle.load(f)
            model = saved['model']
            scaler = saved['scaler']
        print("Model and scaler loaded from 'complement_nb_model.pkl'")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    main()
