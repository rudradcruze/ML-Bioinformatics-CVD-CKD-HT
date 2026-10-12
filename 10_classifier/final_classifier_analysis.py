#!/usr/bin/env python3
"""
FINAL classifier analysis — replaces all previously reported unsupported metrics.

Task:    binary classification, disease (1) vs control (0)
Cohorts: training  = GSE62792 (CKD blood), GSE22356 (HT), GSE66360 (CVD blood)
         external  = GSE113439 (PAH lung tissue, 15 PAH / 11 controls)
Genes:   19 unique gene symbols of the 21 concordant DEG signature
Fixes:   GSE22356 intensities log2-transformed (were linear scale);
         external features per-sample z-scored (platform/tissue robustness).
Design:  repeated (10x) stratified 5-fold CV, threshold fixed at 0.5;
         external validation with model locked on full training data and a
         Youden threshold locked on the training ROC.
Outputs: results.json, ROC_real.pdf/png, feature_importances.csv, CV table
"""
import json
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata  # noqa: F401
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, accuracy_score

RNG = 42
BASE = ("/Volumes/Rudra-External/Research/Bioinformatics/CVD CKD HyperTension/"
        "Final-Code/ML-Bioinformatics-CVD-CKD-HT/0_datasets")
GENES = ['FZD7', 'ZNF248', 'TLR2', 'JAG1', 'TCF7L2', 'SULF2', 'SLMAP', 'LY9',
         'DMXL2', 'RBM47', 'UBASH3B', 'ATG7', 'LILRB2', 'AIF1', 'NPL', 'GRK3',
         'GLUL', 'LILRB3', 'GJB6']
TRAIN = [('CKD_GSE62792', 'GSE62792_matrix_Transpose_labeled.csv', False),
         ('HT_GSE22356', 'GSE22356_matrix_Transpose_labeled.csv', True),
         ('CVD_GSE66360', 'GSE66360_matrix_Transpose_labeled.csv', False)]
EXT = 'GSE113439_matrix_Transpose_labeled.csv'


def load(fn, log2=False):
    df = pd.read_csv(f"{BASE}/{fn}")
    g = [x for x in GENES if x in df.columns]
    X = df[g].apply(pd.to_numeric, errors='coerce')
    if log2:
        X = np.log2(X.clip(lower=0.1))
    y = df.iloc[:, -1].astype(int).values
    return g, X.values, y


def zrow(X):
    return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-9)


def models():
    return {
        'L2-Logistic Regression': Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(penalty='l2', C=1.0, solver='lbfgs',
                                       max_iter=5000, random_state=RNG))]),
        'Random Forest': RandomForestClassifier(
            n_estimators=500, criterion='gini', max_depth=None,
            min_samples_leaf=2, max_features='sqrt', bootstrap=True,
            random_state=RNG, n_jobs=-1),
        'RBF-SVM': Pipeline([
            ('scaler', StandardScaler()),
            ('clf', SVC(kernel='rbf', C=1.0, gamma='scale', probability=True,
                        random_state=RNG))]),
    }


def mets(y, proba, thr=0.5):
    pred = (proba >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred).ravel()
    return dict(AUC=roc_auc_score(y, proba) * 100,
                ACC=accuracy_score(y, pred) * 100,
                SEN=tp / (tp + fn) * 100,
                SPE=tn / (tn + fp) * 100)


def main():
    out = {}
    frames = [load(fn, lg) for _, fn, lg in TRAIN]
    genes = [g for g in GENES if all(g in f[0] for f in frames)]
    Xtr = np.vstack([f[1][:, [f[0].index(g) for g in genes]] for f in frames])
    ytr = np.hstack([f[2] for f in frames])
    print(f"Training pool: n={len(ytr)} (disease={int(ytr.sum())}, "
          f"control={int((1 - ytr).sum())}), {len(genes)} signature genes")

    # ---------- cross-validation ----------
    cv_tbl, roc_cv = {}, {}
    for name, mk in models().items():
        rep = []
        for r in range(10):
            skf = StratifiedKFold(5, shuffle=True, random_state=RNG + r)
            p = cross_val_predict(clone(mk), Xtr, ytr, cv=skf,
                                  method='predict_proba')[:, 1]
            rep.append(mets(ytr, p))
        d = pd.DataFrame(rep)
        cv_tbl[name] = {c: f"{d[c].mean():.2f} ± {d[c].std(ddof=1):.2f}"
                        for c in ['AUC', 'ACC', 'SEN', 'SPE']}
        p = cross_val_predict(clone(mk), Xtr, ytr,
                              cv=StratifiedKFold(5, shuffle=True, random_state=RNG),
                              method='predict_proba')[:, 1]
        fpr, tpr, _ = roc_curve(ytr, p)
        roc_cv[name] = (fpr, tpr, roc_auc_score(ytr, p) * 100)
        print('CV', name, cv_tbl[name])
    out['cross_validation'] = cv_tbl

    # ---------- external validation ----------
    gex, Xex_full, yex = load(EXT)
    shared = [g for g in genes if g in gex]
    missing = [g for g in genes if g not in gex]
    idx = [genes.index(g) for g in shared]
    XtrS = Xtr[:, idx]
    XexS = Xex_full[:, [gex.index(g) for g in shared]]
    Zt, Ze = zrow(XtrS), zrow(XexS)
    print(f"External GSE113439 (PAH lung): n={len(yex)} "
          f"(disease={int(yex.sum())}), {len(shared)} shared genes, "
          f"missing {missing}")

    ext_tbl, roc_ex = {}, {}
    for name, mk in models().items():
        m = clone(mk)
        m.fit(Zt, ytr)
        pt = m.predict_proba(Zt)[:, 1]
        fpr, tpr, thr = roc_curve(ytr, pt)
        t_locked = float(thr[np.argmax(tpr - fpr)])   # Youden, locked on train
        pe = m.predict_proba(Ze)[:, 1]
        r = mets(yex, pe, thr=t_locked)
        r['locked_threshold'] = round(t_locked, 4)
        ext_tbl[name] = {k: round(v, 2) if isinstance(v, float) else v
                         for k, v in r.items()}
        fpr, tpr, _ = roc_curve(yex, pe)
        roc_ex[name] = (fpr, tpr, roc_auc_score(yex, pe) * 100)
        print('EXT', name, ext_tbl[name])
    out['external_validation'] = ext_tbl
    out['external_shared_genes'] = shared
    out['external_missing_genes'] = missing

    # ---------- feature importances (RF, full signature) ----------
    rf = clone(models()['Random Forest']).fit(Xtr, ytr)
    imp = pd.DataFrame({'Gene': genes, 'Importance': rf.feature_importances_}) \
        .sort_values('Importance', ascending=False)
    imp.to_csv('feature_importances.csv', index=False)
    out['rf_feature_importances'] = imp.round(4).to_dict('records')
    print(imp.head(8).to_string(index=False))

    # ---------- ROC figure ----------
    colors = ['#D83030', '#0B593C', '#1985B4']
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, roc, ttl in [(axs[0], roc_cv, 'Cross-validation (training cohorts, pooled)'),
                         (axs[1], roc_ex, 'External validation (GSE113439, PAH lung)')]:
        for c, (name, (fpr, tpr, auc)) in zip(colors, roc.items()):
            ax.step(fpr, tpr, where='post', color=c, lw=2,
                    label=f'{name} (AUC: {auc:.2f}%)')
        ax.plot([0, 1], [0, 1], 'k--', lw=1)
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.set_title(ttl, fontsize=11)
        ax.legend(loc='lower right', fontsize=8)
    plt.tight_layout()
    plt.savefig('ROC_real.pdf', dpi=600)
    plt.savefig('ROC_real.png', dpi=300)
    print('Saved ROC_real.pdf / .png, feature_importances.csv, results.json')
    with open('results.json', 'w') as fh:
        json.dump(out, fh, indent=2)


if __name__ == '__main__':
    main()
