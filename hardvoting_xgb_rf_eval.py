"""
Hard Voting 评估脚本：XGBoost + Random Forest 仅当两者都预测为正例时输出正例。
计算 accuracy, precision, recall, F1，并绘制 ROC 曲线与混淆矩阵热图。
"""
import argparse
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_curve,
    auc,
    confusion_matrix,
)

# 正例标签（与 CSV 中 1 对应）
POS_LABEL = 1
NEG_LABEL = -1


def load_data(csv_path: str, label_col: str = None) -> tuple:
    """加载 CSV：特征 + 标签列（1 与 -1）。"""
    df = pd.read_csv(csv_path)
    if label_col is None:
        # 默认第一列为标签
        label_col = df.columns[0]
    if label_col not in df.columns:
        raise ValueError(f"标签列 '{label_col}' 不在 CSV 中，可选: {list(df.columns)}")
    y = df[label_col].values
    X = df.drop(columns=[label_col]).values
    # 统一为 1 / -1
    y = np.where(y == POS_LABEL, POS_LABEL, NEG_LABEL)
    return X, y


def get_pos_proba(model, X: np.ndarray) -> np.ndarray:
    """获取模型对正例（1）的预测概率，形状 (n_samples,)"""
    proba = model.predict_proba(X)
    # classes_ 可能为 [0,1] 或 [-1,1] 等
    if hasattr(model, "classes_") and POS_LABEL in model.classes_:
        pos_idx = int(np.argmax(model.classes_ == POS_LABEL))
    else:
        pos_idx = 1  # 默认第二列为正例
    return proba[:, pos_idx]


def get_binary_pred(proba: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    """概率转 1/-1 预测"""
    return np.where(proba >= threshold, POS_LABEL, NEG_LABEL)


def load_xgb(path: str):
    """加载 XGBoost 模型，若为 Booster 则包装为 XGBClassifier 风格以支持 predict_proba。"""
    with open(path, "rb") as f:
        model = pickle.load(f)
    if not hasattr(model, "predict_proba"):
        try:
            import xgboost as xgb
            wrapper = xgb.XGBClassifier()
            wrapper._Booster = model
            wrapper._le = None
            # 二分类常见 classes_
            wrapper.classes_ = np.array([NEG_LABEL, POS_LABEL])
            return wrapper
        except Exception as e:
            raise RuntimeError(f"XGBoost 模型需支持 predict_proba 或为 Booster，当前无法包装: {e}")
    return model


def evaluate_and_plot(
    name_short: str,
    name_full: str,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    score: np.ndarray,
    out_dir: str,
    prefix: str,
) -> None:
    """
    计算并打印单个模型的各类指标，绘制 ROC 曲线和混淆矩阵热图。
    - name_short: 用于文件名后缀，如 'xgb'、'rf'、'hardvoting'
    - name_full:  用于图标题与打印，如 'XGBoost'、'Random Forest'
    """
    # 指标
    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred, pos_label=POS_LABEL, zero_division=0)
    rec = recall_score(y_true, y_pred, pos_label=POS_LABEL, zero_division=0)
    f1 = f1_score(y_true, y_pred, pos_label=POS_LABEL, zero_division=0)

    print(f"\n--- {name_full} 结果 ---")
    print(f"Accuracy:  {acc:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall:    {rec:.4f}")
    print(f"F1:        {f1:.4f}")

    # ROC
    fpr, tpr, _ = roc_curve(y_true, score, pos_label=POS_LABEL)
    roc_auc = auc(fpr, tpr)
    print(f"ROC AUC ({name_full}): {roc_auc:.4f}")

    # 路径
    out_dir = out_dir.rstrip("/").rstrip("\\")
    roc_path = f"{out_dir}/{prefix}_{name_short}_roc.pdf"
    cm_path = f"{out_dir}/{prefix}_{name_short}_confusion_matrix.pdf"

    # ROC 曲线
    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, color="darkorange", lw=2, label=f"{name_full} (AUC = {roc_auc:.3f})")
    plt.plot([0, 1], [0, 1], color="navy", lw=1, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate", fontsize=12)
    plt.ylabel("True Positive Rate", fontsize=12)
    plt.title(f"ROC Curve ({name_full})", fontsize=12)
    plt.legend(loc="lower right", fontsize=10)
    plt.tight_layout()
    plt.savefig(roc_path, dpi=300, bbox_inches="tight")
    plt.savefig(roc_path.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"{name_full} ROC 曲线已保存: {roc_path} 以及 PNG 版本")

    # 混淆矩阵
    cm = confusion_matrix(y_true, y_pred, labels=[NEG_LABEL, POS_LABEL])
    fig, ax = plt.subplots(figsize=(5, 4))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=[f"Pred {NEG_LABEL}", f"Pred {POS_LABEL}"],
        yticklabels=[f"True {NEG_LABEL}", f"True {POS_LABEL}"],
        ax=ax,
        cbar_kws={"label": "Count"},
    )
    ax.set_ylabel("True Label", fontsize=11)
    ax.set_xlabel("Predicted Label", fontsize=11)
    ax.set_title(f"Confusion Matrix ({name_full})", fontsize=12)
    plt.tight_layout()
    plt.savefig(cm_path, dpi=300, bbox_inches="tight")
    plt.savefig(cm_path.replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"{name_full} 混淆矩阵热图已保存: {cm_path} 以及 PNG 版本")


def main():
    parser = argparse.ArgumentParser(description="Hard Voting (XGB+RF 均预测为正才为正) 评估与绘图")
    parser.add_argument("--csv", "-i", required=True, help="带特征与标签(1/-1)的 CSV 文件")
    parser.add_argument("--xgb", required=True, help="XGBoost 模型 pkl 路径")
    parser.add_argument("--rf", required=True, help="Random Forest 模型 pkl 路径")
    parser.add_argument("--label-col", default=None, help="标签列名，默认最后一列")
    parser.add_argument("--out-dir", "-o", default=".", help="ROC 与混淆矩阵图保存目录")
    parser.add_argument("--prefix", default="hardvoting", help="输出文件名前缀")
    args = parser.parse_args()

    # 加载数据
    X, y_true = load_data(args.csv, args.label_col)
    print(f"样本数: {X.shape[0]}, 特征数: {X.shape[1]}, 正例数: {(y_true == POS_LABEL).sum()}, 负例数: {(y_true == NEG_LABEL).sum()}")

    # 加载模型
    xgb_model = load_xgb(args.xgb)
    with open(args.rf, "rb") as f:
        rf_model = pickle.load(f)

    # 两模型正例概率
    xgb_pos = get_pos_proba(xgb_model, X)
    rf_pos = get_pos_proba(rf_model, X)

    # Hard Voting：仅当两个模型都预测为正例时，最终才为正例
    xgb_pred = get_binary_pred(xgb_pos)
    rf_pred = get_binary_pred(rf_pos)
    y_pred = np.where((xgb_pred == POS_LABEL) & (rf_pred == POS_LABEL), POS_LABEL, NEG_LABEL)

    # 1）仅用 XGBoost 模型
    evaluate_and_plot(
        name_short="xgb",
        name_full="XGBoost",
        y_true=y_true,
        y_pred=xgb_pred,
        score=xgb_pos,
        out_dir=args.out_dir,
        prefix=args.prefix,
    )

    # 2）仅用 Random Forest 模型
    evaluate_and_plot(
        name_short="rf",
        name_full="Random Forest",
        y_true=y_true,
        y_pred=rf_pred,
        score=rf_pos,
        out_dir=args.out_dir,
        prefix=args.prefix,
    )

    # 3）Hard Voting（XGB 与 RF 同时为正例）
    ensemble_score = np.minimum(xgb_pos, rf_pos)
    evaluate_and_plot(
        name_short="hardvoting",
        name_full="Hard Voting (XGB ∩ RF)",
        y_true=y_true,
        y_pred=y_pred,
        score=ensemble_score,
        out_dir=args.out_dir,
        prefix=args.prefix,
    )


if __name__ == "__main__":
    main()
