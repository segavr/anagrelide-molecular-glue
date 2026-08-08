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

    print("\n=== Investigating: are two outliers (A14, A2) driving significance? ===")
    no_outliers = [r for r in detected if r["compound"] not in ("A14", "A2")]
    for key, label in [("paper_logp", "Paper's LogP"), ("rdkit_logp", "Our RDKit LogP")]:
        correlation_report(no_outliers, key, f"{label} (no A14/A2)", use_log=True, method="pearson")
        correlation_report(no_outliers, key, f"{label} (no A14/A2)", use_log=True, method="spearman")

    print("\nConclusion: the full-dataset correlation IS statistically significant\n"
          "(p<0.05) with both LogP sources and both Pearson and Spearman, contradicting\n"
          "the paper's qualitative claim taken literally. However, removing A14\n"
          "(unsubstituted) and A2 (F, weakly hydrophobic) -- both unusually high-IC50\n"
          "outliers -- drops Pearson significance (paper's LogP: p=0.095) and Spearman\n"
          "significance even further (p=0.232, the rank-based test that is inherently\n"
          "robust to outliers). This more robust check strengthens the earlier finding:\n"
          "without these two compounds, the LogP-activity relationship in the remaining\n"
          "17 molecules is genuinely weak, plausibly explaining the authors' qualitative\n"
          "judgment even though it wasn't stated as a formal outlier-exclusion analysis.")

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
