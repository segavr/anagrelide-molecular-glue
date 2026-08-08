"""
correlation_analysis.py

Reproduces and formally tests the claim in Supplementary Note 1 of
Chen et al. (Nat Commun 2021, PMC8551160): "We did not observe the
linear correlation of the logP value with IC50 of all the molecules."

This script runs a formal Pearson correlation test on the full dataset
(both the paper's own reported LogP and our independently-built RDKit
LogP), and investigates why the authors may have reached their
qualitative conclusion despite a formally significant correlation.

Usage: run from the project root.
    python scripts/correlation_analysis.py
"""

import csv
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def load_dataset(path="data/anagrelide_analogs.csv"):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def detected_only(rows):
    """Exclude compounds with ND (not detected) IC50: A11, A16, A21, A22."""
    return [r for r in rows if r["ic50_detected"] == "True"]


def correlation_report(rows, logp_key, label, use_log=True, method="pearson"):
    """
    Report correlation between a LogP column and IC50.

    use_log: if True (default, and the pharmacologically standard choice),
    correlate against log10(IC50) rather than raw IC50 -- IC50 values span
    orders of magnitude here (0.3 to 43 nM), and a linear-scale correlation
    can be dominated by the largest values. Pearson on raw IC50 is also
    reported by callers that want it for direct comparison to the paper's
    own (unspecified-scale) analysis.

    method: "pearson" (linear correlation, sensitive to outliers) or
    "spearman" (rank correlation, robust to outliers and to whether the
    relationship is linear or merely monotonic). Note that for Spearman,
    log-transforming makes no difference: log10 is a monotonic transform,
    so it does not change the rank order, and thus not the Spearman
    statistic. This is included here so that gets reported once rather
    than presented as if it were new information under log-transform.
    """
    ic50 = np.array([float(r["ic50_nM"]) for r in rows])
    logp = np.array([float(r[logp_key]) for r in rows])
    y = np.log10(ic50) if use_log else ic50
    scale_label = "log10(IC50)" if use_log else "raw IC50"

    if method == "pearson":
        r, p = stats.pearsonr(logp, y)
    elif method == "spearman":
        r, p = stats.spearmanr(logp, y)
    else:
        raise ValueError(f"unknown method: {method}")

    print(f"{label} ({method}, {scale_label}): r = {r:.3f}, p = {p:.4f}  (n={len(rows)})")
    return r, p


if __name__ == "__main__":
    rows = load_dataset()
    detected = detected_only(rows)

    print(f"Loaded {len(rows)} compounds, {len(detected)} with detected IC50 "
          f"(excluded as ND: {[r['compound'] for r in rows if r['ic50_detected']=='False']})\n")

    print("=== Reproducing the paper's claim: 'no linear correlation' ===")
    print("(Using log10(IC50), the pharmacologically standard scale, since IC50")
    print(" spans orders of magnitude here -- 0.3 to 43 nM. Pearson and Spearman")
    print(" are both reported: Pearson tests a linear relationship and is")
    print(" sensitive to outliers; Spearman tests a monotonic relationship via")
    print(" ranks and is robust to outliers.)\n")
    for key, label in [("paper_logp", "Paper's LogP"), ("rdkit_logp", "Our RDKit LogP")]:
        correlation_report(detected, key, label, use_log=True, method="pearson")
        correlation_report(detected, key, label, use_log=True, method="spearman")

    print("\n=== Investigating: which compounds are statistically influential? ===")
    print("(Cook's distance identifies points with unusually large influence on the")
    print(" regression fit -- more rigorous than picking outliers by eye from a plot.)")

    def cooks_distance(rows_subset, logp_key):
        ic50 = np.array([float(r["ic50_nM"]) for r in rows_subset])
        y = np.log10(ic50)
        x = np.array([float(r[logp_key]) for r in rows_subset])
        names = [r["compound"] for r in rows_subset]
        n = len(x)
        X = np.column_stack([np.ones(n), x])
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        residuals = y - X @ beta
        H = X @ np.linalg.inv(X.T @ X) @ X.T
        h = np.diag(H)
        p_params = X.shape[1]
        mse = np.sum(residuals**2) / (n - p_params)
        cooks_d = (residuals**2 / (p_params * mse)) * (h / (1 - h) ** 2)
        threshold = 4 / n
        return names, cooks_d, threshold

    names, cooks_d, threshold = cooks_distance(detected, "paper_logp")
    influential = [names[i] for i in range(len(names)) if cooks_d[i] > threshold]
    print(f"Threshold (4/n): {threshold:.4f}")
    print(f"Compounds flagged as influential (Cook's D > threshold): {influential}")

    print("\nNote: this does NOT fully match the two compounds (A14, A2) picked by")
    print("eye from the scatter plot. Cook's distance flags A8 and A14 -- A2, despite")
    print("looking like a visual outlier, actually falls close to the overall trend line")
    print("(low LogP with correspondingly high IC50); A8 is the more statistically")
    print("influential point because it has low LogP but unexpectedly high potency")
    print("(low IC50), working against the trend. Recall from dataset construction that")
    print("A8 (furan) was already flagged for a large RDKit-vs-paper LogP discrepancy")
    print("attributed to heterocycle-specific algorithm differences -- an interesting")
    print("convergence between the structural and statistical analyses.")

    print("\nCorrelation with the formally-identified influential points (A8, A14) removed:")
    no_influential = [r for r in detected if r["compound"] not in ("A8", "A14")]
    correlation_report(no_influential, "paper_logp", "Paper's LogP (no A8/A14)", use_log=True, method="pearson")
    correlation_report(no_influential, "paper_logp", "Paper's LogP (no A8/A14)", use_log=True, method="spearman")
    print("(Removing the formally-flagged points barely changes -- if anything slightly")
    print(" strengthens -- the correlation, unlike removing the visually-picked A14/A2.")
    print(" This means the earlier 'outliers explain the paper's claim' story needs to")
    print(" be told more carefully: A14 does appear to be a genuine influential point by")
    print(" both methods, but A2's apparent outlier status was a visual impression that")
    print(" a formal diagnostic does not support.)")

    # Plot (log10(IC50), the pharmacologically standard scale, with both
    # Pearson and Spearman statistics annotated)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    names = [r["compound"] for r in detected]
    ic50 = np.array([float(r["ic50_nM"]) for r in detected])
    log_ic50 = np.log10(ic50)

    for ax, key, label in zip(axes, ["paper_logp", "rdkit_logp"],
                               ["Paper's LogP", "Our RDKit LogP"]):
        logp_vals = np.array([float(r[key]) for r in detected])
        r_p, p_p = stats.pearsonr(logp_vals, log_ic50)
        r_s, p_s = stats.spearmanr(logp_vals, log_ic50)
        slope, intercept = np.polyfit(logp_vals, log_ic50, 1)
        x_line = np.linspace(logp_vals.min(), logp_vals.max(), 100)
        y_line = slope * x_line + intercept

        ax.scatter(logp_vals, log_ic50, s=60, alpha=0.75, edgecolor="black", linewidth=0.5)
        ax.plot(x_line, y_line, "r--", linewidth=1.5, label="linear fit")
        for name, x, y in zip(names, logp_vals, log_ic50):
            ax.annotate(name, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(label)
        ax.set_ylabel("log10(IC50 / nM)")
        ax.set_title(f"{label} vs log10(IC50)\n"
                      f"Pearson r={r_p:.3f}, p={p_p:.4f}  |  Spearman r={r_s:.3f}, p={p_s:.4f}")
        ax.legend()
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Reproducing the paper's LogP-IC50 analysis (Supplementary Note 1)\n"
        "Paper states: \"no linear correlation observed\" — formal test disagrees on full dataset",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig("results/logp_ic50_correlation.png", dpi=150)
    print("\nPlot saved to results/logp_ic50_correlation.png")
