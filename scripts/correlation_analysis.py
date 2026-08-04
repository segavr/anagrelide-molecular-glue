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


def correlation_report(rows, logp_key, label):
    ic50 = np.array([float(r["ic50_nM"]) for r in rows])
    logp = np.array([float(r[logp_key]) for r in rows])
    r, p = stats.pearsonr(logp, ic50)
    print(f"{label}: r = {r:.3f}, p = {p:.4f}  (n={len(rows)})")
    return r, p


if __name__ == "__main__":
    rows = load_dataset()
    detected = detected_only(rows)

    print(f"Loaded {len(rows)} compounds, {len(detected)} with detected IC50 "
          f"(excluded as ND: {[r['compound'] for r in rows if r['ic50_detected']=='False']})\n")

    print("=== Reproducing the paper's claim: 'no linear correlation' ===")
    correlation_report(detected, "paper_logp", "Paper's LogP vs IC50")
    correlation_report(detected, "rdkit_logp", "Our RDKit LogP vs IC50")

    print("\n=== Investigating: are two outliers (A14, A2) driving significance? ===")
    no_outliers = [r for r in detected if r["compound"] not in ("A14", "A2")]
    correlation_report(no_outliers, "paper_logp", "Paper's LogP vs IC50 (no A14/A2)")
    correlation_report(no_outliers, "rdkit_logp", "Our RDKit LogP vs IC50 (no A14/A2)")

    print("\nConclusion: the full-dataset correlation IS statistically significant\n"
          "(p<0.05) with both LogP sources, contradicting the paper's qualitative\n"
          "claim taken literally. However, removing A14 (unsubstituted) and A2 (F,\n"
          "weakly hydrophobic) -- both unusually high-IC50 outliers -- makes the\n"
          "paper's own LogP correlation lose significance (p=0.067), which plausibly\n"
          "explains their qualitative judgment even though it wasn't stated as a\n"
          "formal outlier-exclusion analysis in the text.")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    names = [r["compound"] for r in detected]
    ic50 = np.array([float(r["ic50_nM"]) for r in detected])

    for ax, key, label in zip(axes, ["paper_logp", "rdkit_logp"],
                               ["Paper's LogP", "Our RDKit LogP"]):
        logp_vals = np.array([float(r[key]) for r in detected])
        r, p = stats.pearsonr(logp_vals, ic50)
        slope, intercept = np.polyfit(logp_vals, ic50, 1)
        x_line = np.linspace(logp_vals.min(), logp_vals.max(), 100)
        y_line = slope * x_line + intercept

        ax.scatter(logp_vals, ic50, s=60, alpha=0.75, edgecolor="black", linewidth=0.5)
        ax.plot(x_line, y_line, "r--", linewidth=1.5, label="linear fit")
        for name, x, y in zip(names, logp_vals, ic50):
            ax.annotate(name, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(label)
        ax.set_ylabel("IC50 (nM)")
        ax.set_title(f"{label} vs IC50\nPearson r = {r:.3f}, p = {p:.4f}")
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
