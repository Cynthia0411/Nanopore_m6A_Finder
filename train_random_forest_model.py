import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.metrics import roc_curve, auc  # 添加ROC相关的导入
import matplotlib.pyplot as plt
import seaborn as sns

def load_data(file_path):
    # 读取数据，没有header，第一列是标签
    data = pd.read_csv(file_path, header=None)
    
    # 分离特征和标签
    X = data.iloc[:, 1:]  # 所有特征列（45维）
    y = data.iloc[:, 0]   # 标签列
    
    return X, y

def train_and_evaluate():
    # 加载数据
    X, y = load_data('balanced_feature_with_1_A.csv')
    
    # 划分训练集和测试集
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=7, stratify=y)
    
    # 初始化并训练随机森林模型
    rf_model = RandomForestClassifier(n_estimators=100, random_state=7)
    rf_model.fit(X_train, y_train)
    
    # 在测试集上进行预测
    y_pred = rf_model.predict(X_test)
    y_pred_proba = rf_model.predict_proba(X_test)[:, 1]  # 获取正类的概率
    
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
    plt.savefig('confusion_matrix_rf_balanced_all_with_1_A.png')
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
    plt.savefig('roc_curve_rf_balanced_all_with_1_A.png')
    plt.close()
    
    # 特征重要性分析
    feature_importance = pd.DataFrame({
        'feature': range(X.shape[1]),
        'importance': rf_model.feature_importances_
    })
    feature_importance = feature_importance.sort_values('importance', ascending=False)
    
    # 绘制特征重要性图
    plt.figure(figsize=(12, 6))
    plt.bar(range(20), feature_importance['importance'][:20])
    plt.title('Top 20 Feature Importance')
    plt.xlabel('Feature Index')
    plt.ylabel('Importance')
    plt.savefig('feature_importance_rf_balanced_all_with_1_A.png')
    plt.close()
    
    return rf_model

def main():
    try:
        print("Loading data and training model...")
        # 训练模型并获取评估结果
        model = train_and_evaluate()
        
        # 使用pickle保存模型
        import pickle
        with open('rf_model_balanced_with_1_A.pkl', 'wb') as f:
            pickle.dump(model, f)
        print("Model training completed and saved as 'rf_model_balanced_with_1_A.pkl'")
        
        # 测试加载模型（可选）
        print("Testing model loading...")
        with open('rf_model_balanced_with_1_A.pkl', 'rb') as f:
            loaded_model = pickle.load(f)
        print("Model successfully loaded!")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    main()
