#!/usr/bin/env python3
"""
Classifier analysis on the 12-gene FDR-primary signature.
Same protocol as final_classifier_analysis.py (repeated 10x stratified 5-fold CV,
fixed 0.5 threshold; external validation with locked model + Youden threshold).
Probesets matched by any ///-separated symbol containing the gene name.
"""
import json
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.base import clone
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix, accuracy_score

RNG = 42
ROOT = "/Volumes/Rudra-External/Research/Bioinformatics/CVD CKD HyperTension"
GENES = ['APOBEC3A', 'DRAM1', 'FNDC3B', 'FOXO3', 'GAS7', 'IRAK3',
         'PLBD1', 'RAB31', 'SIPA1L2', 'SLC8A1', 'SORT1', 'SULF2']
TRAIN = [('CKD_GSE62792', f'{ROOT}/0_CKD_GSE62792_pvalue_05/GSE62792_matrix_Transpose.csv',
          f'{ROOT}/0_CKD_GSE62792_pvalue_05/GSE62792_matrix_Transpose_labeled.csv.csv', False),
         ('HT_GSE22356', f'{ROOT}/0_HyperTension_GSE22356_pvalue_05/GSE22356_matrix_Transpose.csv',
          f'{ROOT}/0_HyperTension_GSE22356_pvalue_05/GSE22356_matrix_Transpose_labeled.csv', True),
         ('CVD_GSE66360', f'{ROOT}/0_CVD_GSE66360_pvalue_05/GSE66360_matrix_Transpose.csv',
          f'{ROOT}/0_CVD_GSE66360_pvalue_05/GSE66360_matrix_Transpose_labeled.csv', False)]
EXT = (f'{ROOT}/0_HyperTension_GSE113439_Validation/GSE113439_matrix_Transpose.csv',
       f'{ROOT}/0_HyperTension_GSE113439_Validation/GSE113439_matrix_Transpose_labeled.csv')


def match_col(df, gene):
    """Find a probeset column whose ///-split symbols include the gene."""
    exact = [c for c in df.columns if c == gene]
    if exact:
        return exact[0]
    multi = [c for c in df.columns if gene in c.split('///')]
    return multi[0] if multi else None


def load(full_path, label_path, log2=False, genes=GENES):
    full = pd.read_csv(full_path)
    lab = pd.read_csv(label_path)
    y = lab.iloc[:, -1].astype(int).values
    n = min(len(full), len(y))
    y = y[:n]
    cols = {}
    for g in genes:
        c = match_col(full, g)
        if c is not None:
            cols[g] = pd.to_numeric(full[c].iloc[:n], errors='coerce').values
    X = pd.DataFrame(cols).values
    if log2:
        X = np.log2(np.clip(X, 0.1, None))
    ok = ~np.isnan(X).any(1)
    return list(cols), X[ok], y[ok]


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
    frames = [load(f, l, lg) for _, f, l, lg in TRAIN]
    genes = [g for g in GENES if all(g in f[0] for f in frames)]
    Xtr = np.vstack([f[1][:, [f[0].index(g) for g in genes]] for f in frames])
    ytr = np.hstack([f[2] for f in frames])
    print(f"FDR-12 training pool: n={len(ytr)} (disease={int(ytr.sum())}, "
          f"control={int((1 - ytr).sum())}), {len(genes)} genes available: {genes}")

    cv_tbl, roc_cv = {}, {}
    for name, mk in models().items():
        rep = []
        for r in range(10):
            p = cross_val_predict(clone(mk), Xtr, ytr,
                                  cv=StratifiedKFold(5, shuffle=True, random_state=RNG + r),
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

    # external
    gex_full, Xex_all, yex = load(EXT[0], EXT[1])
    shared = [g for g in genes if g in gex_full]
    missing = [g for g in genes if g not in gex_full]
    idx = [genes.index(g) for g in shared]
    XtrS = Xtr[:, idx]
    XexS = Xex_all[:, [gex_full.index(g) for g in shared]]
    Zt, Ze = zrow(XtrS), zrow(XexS)
    print(f"External GSE113439: n={len(yex)} (disease={int(yex.sum())}), shared {len(shared)} genes, missing {missing}")

    ext_tbl, roc_ex = {}, {}
    for name, mk in models().items():
        m = clone(mk).fit(Zt, ytr)
        pt = m.predict_proba(Zt)[:, 1]
        fpr, tpr, thr = roc_curve(ytr, pt)
        t_locked = float(thr[np.argmax(tpr - fpr)])
        pe = m.predict_proba(Ze)[:, 1]
        r = mets(yex, pe, thr=t_locked)
        r['locked_threshold'] = round(t_locked, 4)
        ext_tbl[name] = {k: round(v, 2) if isinstance(v, float) else v for k, v in r.items()}
        fpr, tpr, _ = roc_curve(yex, pe)
        roc_ex[name] = (fpr, tpr, roc_auc_score(yex, pe) * 100)
        print('EXT', name, ext_tbl[name])

    # feature importances
    rf = clone(models()['Random Forest']).fit(Xtr, ytr)
    imp = pd.DataFrame({'Gene': genes, 'Importance': rf.feature_importances_}) \
        .sort_values('Importance', ascending=False)
    imp.to_csv('fdr12_feature_importances.csv', index=False)
    print(imp.to_string(index=False))

    # ROC figure
    colors = ['#D83030', '#0B593C', '#1985B4']
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, roc, ttl in [(axs[0], roc_cv, 'Cross-validation (FDR-12 signature, pooled cohorts)'),
                         (axs[1], roc_ex, 'External validation (GSE113439)')]:
        for c, (name, (fpr, tpr, auc)) in zip(colors, roc.items()):
            ax.step(fpr, tpr, where='post', color=c, lw=2,
                    label=f'{name} (AUC: {auc:.2f}%)')
        ax.plot([0, 1], [0, 1], 'k--', lw=1)
        ax.set_xlabel('False Positive Rate'); ax.set_ylabel('True Positive Rate')
        ax.set_title(ttl, fontsize=11); ax.legend(loc='lower right', fontsize=8)
    plt.tight_layout()
    plt.savefig('ROC_FDR12.pdf', dpi=600)
    plt.savefig('ROC_FDR12.png', dpi=300)

    with open('fdr12_results.json', 'w') as fh:
        json.dump({'genes': genes, 'external_shared': shared, 'external_missing': missing,
                   'cv': cv_tbl, 'external': ext_tbl,
                   'rf_importances': imp.round(4).to_dict('records')}, fh, indent=2)
    print('Saved fdr12_results.json, ROC_FDR12.pdf/png, fdr12_feature_importances.csv')


if __name__ == '__main__':
    main()
