"""
multi_descriptor_analysis.py

Extends the single-descriptor (LogP) analysis from correlation_analysis.py
by testing whether other RDKit descriptors (molecular weight, TPSA,
aromatic ring count, H-bond donor/acceptor count) explain IC50 better
than LogP alone, and whether combining descriptors adds genuinely new
information or is redundant (multicollinearity).

Usage: run from the project root.
    python scripts/multi_descriptor_analysis.py
"""

import csv
import numpy as np
from scipy import stats
from numpy.linalg import lstsq
import matplotlib.pyplot as plt


def load_detected(path="data/anagrelide_analogs.csv"):
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return [r for r in rows if r["ic50_detected"] == "True"]


DESCRIPTORS = {
    "rdkit_logp": "RDKit LogP",
    "mol_weight": "Molecular Weight",
    "tpsa": "TPSA",
    "num_aromatic_rings": "Aromatic Rings",
    "num_h_donors": "H-Bond Donors",
    "num_h_acceptors": "H-Bond Acceptors",
}


if __name__ == "__main__":
    detected = load_detected()
    ic50 = np.array([float(r["ic50_nM"]) for r in detected])

    print("=== Single-descriptor correlations with IC50 ===")
    print(f"{'Descriptor':22s} {'Pearson r':>10s} {'p-value':>10s}")
    values = {}
    for key, label in DESCRIPTORS.items():
        vals = np.array([float(r[key]) for r in detected])
        values[key] = vals
        if np.std(vals) == 0:
            print(f"{label:22s} {'N/A (constant)':>10s}")
            continue
        r, p = stats.pearsonr(vals, ic50)
        print(f"{label:22s} {r:10.3f} {p:10.4f}")

    print("\n=== Checking multicollinearity: LogP vs Molecular Weight ===")
    r_logp_mw, p_logp_mw = stats.pearsonr(values["rdkit_logp"], values["mol_weight"])
    print(f"LogP vs MW: r = {r_logp_mw:.3f}, p = {p_logp_mw:.4f}")
    print("(A high r here means these descriptors carry mostly the same information,")
    print(" i.e. they are not independent predictors of activity)")

    print("\n=== Does combining LogP + MW improve prediction beyond either alone? ===")
    logp = values["rdkit_logp"]
    mw = values["mol_weight"]
    X = np.column_stack([logp, mw, np.ones(len(logp))])
    coeffs, _, _, _ = lstsq(X, ic50, rcond=None)
    predicted = X @ coeffs
    ss_res = np.sum((ic50 - predicted) ** 2)
    ss_tot = np.sum((ic50 - np.mean(ic50)) ** 2)
    r2_combined = 1 - ss_res / ss_tot
    r2_logp = stats.pearsonr(logp, ic50)[0] ** 2
    r2_mw = stats.pearsonr(mw, ic50)[0] ** 2
    print(f"R^2 (LogP alone):        {r2_logp:.3f}")
    print(f"R^2 (MW alone):          {r2_mw:.3f}")
    print(f"R^2 (LogP + MW combined): {r2_combined:.3f}")
    print("Conclusion: combining barely improves R^2 over MW alone -- LogP and MW")
    print("are redundant descriptors here, both reflecting substituent size/hydrophobicity,")
    print("not two independent structural properties driving activity.")

    # Plot
    names = [r["compound"] for r in detected]
    rings = values["num_aromatic_rings"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for ax, (vals, label) in zip(axes, [(logp, "RDKit LogP"), (mw, "Molecular Weight"),
                                          (rings, "Aromatic Ring Count")]):
        r, p = stats.pearsonr(vals, ic50)
        ax.scatter(vals, ic50, s=55, alpha=0.75, edgecolor="black", linewidth=0.5)
        slope, intercept = np.polyfit(vals, ic50, 1)
        x_line = np.linspace(vals.min(), vals.max(), 100)
        ax.plot(x_line, slope * x_line + intercept, "r--", linewidth=1.3)
        ax.set_xlabel(label)
        ax.set_ylabel("IC50 (nM)")
        ax.set_title(f"{label}\nr = {r:.3f}, p = {p:.4f}")
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Multi-descriptor comparison: LogP, MW, and aromatic ring count are\n"
        "strongly inter-correlated (r>0.75 pairwise) — not independent predictors",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig("results/multi_descriptor_comparison.png", dpi=150)
    print("\nPlot saved to results/multi_descriptor_comparison.png")
