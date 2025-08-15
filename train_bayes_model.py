import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import GaussianNB
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

def train_and_evaluate():
    # 加载数据
    X, y = load_data('feature_within_DRACH_81_82.csv')
    
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=7)
    
    # 初始化并训练高斯朴素贝叶斯模型
    nb_model = GaussianNB()
    nb_model.fit(X_train, y_train)
    
    # 在测试集上进行预测
    y_pred = nb_model.predict(X_test)
    y_pred_proba = nb_model.predict_proba(X_test)[:, 1]
    
    # 计算并打印评估指标
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
    plt.savefig('bayes_confusion_matrix_DRACH.png')
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
    plt.savefig('bayes_roc_curve_DRACH.png')
    plt.close()
    
    # 分析每个特征的方差
    feature_variances = pd.DataFrame({
        'feature': range(X.shape[1]),
        'variance_class_-1': nb_model.var_[0],  # 标签-1的方差
        'variance_class_1': nb_model.var_[1],   # 标签1的方差
        'mean_class_-1': nb_model.theta_[0],    # 标签-1的均值
        'mean_class_1': nb_model.theta_[1]      # 标签1的均值
    })
    
    # 计算两个类别间的方差差异
    feature_variances['variance_diff'] = np.abs(feature_variances['variance_class_-1'] - 
                                              feature_variances['variance_class_1'])
    
    # 按方差差异排序
    feature_variances = feature_variances.sort_values('variance_diff', ascending=False)
    
    # 绘制前20个特征的方差对比图
    plt.figure(figsize=(15, 8))
    x = range(20)
    width = 0.35
    
    plt.bar(x, feature_variances['variance_class_-1'][:20], width, 
            label='Class -1 Variance', color='skyblue')
    plt.bar([i + width for i in x], feature_variances['variance_class_1'][:20], 
            width, label='Class 1 Variance', color='lightcoral')
    
    plt.title('Top 20 Features: Variance Comparison Between Classes')
    plt.xlabel('Feature Index')
    plt.ylabel('Variance')
    plt.legend()
    plt.xticks([i + width/2 for i in x], feature_variances['feature'][:20])
    plt.grid(True, axis='y')
    plt.savefig('bayes_feature_variance_comparison_DRACH.png')
    plt.close()
    
    # 打印前10个最具区分性的特征的详细信息
    print("\nTop 10 Most Discriminative Features:")
    print("Feature Index | Class -1 Mean | Class 1 Mean | Class -1 Variance | Class 1 Variance")
    print("-" * 75)
    for idx in feature_variances['feature'][:10]:
        print(f"{idx:12d} | {nb_model.theta_[0][idx]:12.4f} | {nb_model.theta_[1][idx]:11.4f} | "
              f"{nb_model.var_[0][idx]:15.4f} | {nb_model.var_[1][idx]:14.4f}")
    
    return nb_model

def main():
    try:
        print("Loading data and training Bayes model...")
        # 训练模型并获取评估结果
        model = train_and_evaluate()
        
        # 保存模型
        import pickle
        with open('bayes_model_DRACH.pkl', 'wb') as f:
            pickle.dump(model, f)
        print("Model training completed and saved as 'bayes_model.pkl'")
        
        # 测试加载模型
        print("Testing model loading...")
        with open('bayes_model_DRACH.pkl', 'rb') as f:
            loaded_model = pickle.load(f)
        print("Model successfully loaded!")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    main()
