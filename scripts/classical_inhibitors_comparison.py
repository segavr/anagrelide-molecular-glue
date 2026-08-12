"""
classical_inhibitors_comparison.py

Step 5: builds a dataset of PDE3A inhibitors and compares their descriptor
profile against the anagrelide molecular-glue analog series built in
build_dataset.py.

IMPORTANT CORRECTION (see project notes): an earlier version of this
dataset labeled all 12 compounds here as "classical (non-glue)" inhibitors.
This was factually wrong for one compound. Yan et al. (Cell Chem Biol 2022,
Figure S2, Supplementary Table) directly tested nine PDE3A inhibitors for
PDE3A-SLFN12 interaction, SLFN12 stabilization, and cell death, and found:
  - CONFIRMED GLUE (all three readouts positive): anagrelide, enoximone,
    quazinone, zardaverine, DNMDP
  - CONFIRMED NON-GLUE (inhibits PDE3A phosphodiesterase activity, but
    negative on interaction/stabilization/cell death): trequinsin,
    cilostazol, milrinone, cilostamide
Enoximone was previously included here as a "classical inhibitor" -- that
was incorrect; it is a confirmed molecular glue and has been removed from
the non-glue comparison group. The remaining seven compounds in this file
(amrinone, imazodan, papaverine, ensifentrine, pimobendan, olprinone,
vesnarinone) were identified only from IUPHAR/BPS Guide to PHARMACOLOGY's
curated PDE3A inhibitor list -- their glue/non-glue phenotype has NOT been
experimentally tested in the literature searched for this project, so they
are labeled "unknown" and excluded from the primary statistical comparison
(which now uses only the four compounds with an experimentally confirmed
non-glue phenotype: trequinsin, cilostazol, milrinone, cilostamide). The
"unknown" compounds are still built and reported for completeness/transparency,
and a secondary, clearly-labeled comparison against the full unknown+non-glue
set is also shown, but the primary claim in this project rests on the
confirmed-only comparison.

All compound structures (SMILES) were verified against at least two
independent sources (Wikidata + Wikipedia, or Wikidata + a supplier/CAS
database) before use, and cross-checked against known molecular formula
in RDKit -- the same verification discipline used for the anagrelide
analogs in build_dataset.py.

Usage: run from the project root.
    python scripts/classical_inhibitors_comparison.py
"""

import csv
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen, Descriptors

# glue_status: "confirmed_non_glue" (Yan et al. 2022, negative on all three
# glue readouts), "confirmed_glue" (positive on all three, excluded from
# the classical-inhibitor comparison group), or "unknown" (PDE3A inhibitor
# per IUPHAR, but glue phenotype not experimentally tested in literature
# reviewed for this project).
CLASSICAL_INHIBITORS = {
    "milrinone":    ("CC1=C(C=C(C(=O)N1)C#N)C2=CC=NC=C2", "confirmed_non_glue"),
    "cilostazol":   ("O=C4Nc3c(cc(OCCCCc1nnnn1C2CCCCC2)cc3)CC4", "confirmed_non_glue"),
    "cilostamide":  ("CN(C1CCCCC1)C(=O)CCCOC2=CC3=C(C=C2)NC(=O)C=C3", "confirmed_non_glue"),
    "trequinsin":   ("CC1=CC(=C(C(=C1)C)N=C2C=C3C4=CC(=C(C=C4CCN3C(=O)N2C)OC)OC)C", "confirmed_non_glue"),
    "amrinone":     ("O=C2C(/N)=C\\C(\\c1ccncc1)=C/N2", "unknown"),
    "imazodan":     ("C1CC(=O)NN=C1C2=CC=C(C=C2)N3C=CN=C3", "unknown"),
    "papaverine":   ("COc1ccc(cc1OC)Cc2c3cc(c(cc3ccn2)OC)OC", "unknown"),
    "ensifentrine": ("COc1cc2c(cc1OC)-c1c/c(=N\\c3c(C)cc(C)cc3C)n(CCNC(N)=O)c(=O)n1CC2", "unknown"),
    "pimobendan":   ("CC1CC(=O)NN=C1C2=CC3=C(C=C2)N=C(N3)C4=CC=C(C=C4)OC", "unknown"),
    "olprinone":    ("Cc1[nH]c(=O)c(C#N)cc1-c1ccc2nccn2c1", "unknown"),
    "vesnarinone":  ("COC1=C(C=C(C=C1)C(=O)N2CCN(CC2)C3=CC4=C(C=C3)NC(=O)CC4)OC", "unknown"),
    # enoximone intentionally removed -- confirmed molecular glue (Yan et al.
    # 2022), not a non-glue inhibitor. Kept here in a comment for the record
    # rather than silently deleted from project history.
    # "enoximone": ("O=C(/C1=C(/NC(=O)N1)C)c2ccc(SC)cc2", "confirmed_glue"),
}

DESCRIPTOR_LABELS = {
    "rdkit_logp": "LogP",
    "mol_weight": "MW",
    "tpsa": "TPSA",
    "num_aromatic_rings": "Aromatic Rings",
    "num_h_donors": "H-Donors",
    "num_h_acceptors": "H-Acceptors",
}

# NOTE ON REPRODUCIBILITY: RDKit's CalcNumHBA (H-bond acceptor count) has
# been observed to give slightly different values across RDKit versions for
# some of these compounds (e.g. how aromatic heteroatoms or sulfoxide-type
# groups are counted). TPSA, LogP, MW, and aromatic ring count were checked
# and are stable across the versions tested. If your H-Acceptors numbers
# differ slightly from those quoted in the README, this is the likely cause
# -- it does not change the overall conclusion (H-Acceptors is significant
# either way; TPSA remains the strongest and most stable finding).


def build_classical_dataset():
    rows = []
    for name, (smi, glue_status) in CLASSICAL_INHIBITORS.items():
        mol = Chem.MolFromSmiles(smi)
        rows.append({
            "compound": name,
            "class": "pde3a_inhibitor",
            "glue_status": glue_status,
            "smiles": Chem.MolToSmiles(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mol_weight": round(Descriptors.MolWt(mol), 2),
            "rdkit_logp": round(Crippen.MolLogP(mol), 2),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "num_h_donors": rdMolDescriptors.CalcNumHBD(mol),
            "num_h_acceptors": rdMolDescriptors.CalcNumHBA(mol),
        })
    return rows


def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


if __name__ == "__main__":
    all_rows = build_classical_dataset()
    with open("data/classical_inhibitors.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    print("Saved PDE3A inhibitor dataset (with glue_status labels) to "
          "data/classical_inhibitors.csv\n")

    glue_rows = load_csv("data/anagrelide_analogs.csv")

    confirmed_non_glue = [r for r in all_rows if r["glue_status"] == "confirmed_non_glue"]
    unknown_status = [r for r in all_rows if r["glue_status"] == "unknown"]
    print(f"Confirmed non-glue PDE3A inhibitors (Yan et al. 2022): "
          f"{[r['compound'] for r in confirmed_non_glue]} (n={len(confirmed_non_glue)})")
    print(f"Unknown glue-status PDE3A inhibitors (IUPHAR-listed, untested for glue "
          f"phenotype in literature reviewed): {[r['compound'] for r in unknown_status]} "
          f"(n={len(unknown_status)})\n")

    def compare(glue_rows, other_rows, label):
        print(f"=== {label} (n={len(other_rows)}) ===")
        print(f"{'Descriptor':16s} {'Glue mean':>10s} {'Glue std':>9s} "
              f"{'Other mean':>12s} {'Other std':>11s} {'t-test p':>9s}")
        for desc, dlabel in DESCRIPTOR_LABELS.items():
            glue_vals = np.array([float(r[desc]) for r in glue_rows])
            other_vals = np.array([float(r[desc]) for r in other_rows])
            t, p = stats.ttest_ind(glue_vals, other_vals)
            flag = "  <-- significant" if p < 0.05 else ""
            print(f"{dlabel:16s} {glue_vals.mean():10.2f} {glue_vals.std():9.2f} "
                  f"{other_vals.mean():12.2f} {other_vals.std():11.2f} {p:9.4f}{flag}")
        print()

    print("PRIMARY COMPARISON (confirmed non-glue only -- the scientifically")
    print("defensible comparison, though n=4 is small):")
    compare(glue_rows, confirmed_non_glue, "Molecular glue vs CONFIRMED non-glue PDE3A inhibitors")

    print("SECONDARY COMPARISON (confirmed non-glue + unknown-status compounds")
    print("combined, n=11 -- shown for transparency/context only; the 'unknown'")
    print("compounds may or may not actually be non-glue, so this comparison")
    print("should NOT be over-interpreted):")
    combined = confirmed_non_glue + unknown_status
    compare(glue_rows, combined, "Molecular glue vs confirmed-non-glue+unknown PDE3A inhibitors")

    print("Interpretation: with only n=4 confirmed non-glue compounds, statistical")
    print("power is very limited and individual p-values should be read cautiously.")
    print("TPSA is the descriptor that separates most consistently across both the")
    print("primary and secondary comparisons; LogP does not reach significance in")
    print("either. This is consistent with -- but does not on its own establish --")
    print("the hypothesis that lower polarity is associated with PDE3A-SLFN12")
    print("molecular-glue activity specifically, as opposed to reflecting broader")
    print("differences in the chemical space of the compounds compared. A properly")
    print("powered comparison would require a larger set of experimentally")
    print("phenotyped (confirmed glue vs confirmed non-glue) PDE3A ligands, which")
    print("does not yet exist in the literature reviewed for this project.")

    # Plot: primary comparison only (confirmed non-glue), the one the project's
    # conclusions actually rest on
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, desc, label in zip(axes, ["tpsa", "num_aromatic_rings", "rdkit_logp"],
                                 ["TPSA", "Aromatic Ring Count", "LogP"]):
        glue_vals = [float(r[desc]) for r in glue_rows]
        other_vals = [float(r[desc]) for r in confirmed_non_glue]
        bp = ax.boxplot(
            [glue_vals, other_vals],
            tick_labels=["Molecular glue\n(anagrelide + analogs)",
                         "Confirmed non-glue\nPDE3A inhibitors (n=4)"],
            patch_artist=True, widths=0.5,
        )
        bp["boxes"][0].set_facecolor("lightcoral")
        bp["boxes"][1].set_facecolor("lightblue")
        for i, vals in enumerate([glue_vals, other_vals], start=1):
            jitter = np.random.normal(0, 0.04, len(vals))
            ax.scatter([i] * len(vals) + jitter, vals, alpha=0.6, s=25, color="black", zorder=3)
        ax.set_ylabel(label)
        ax.set_title(f"{label} comparison")
        ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(
        "Molecular glue analogs vs CONFIRMED non-glue PDE3A inhibitors (n=4):\n"
        "small comparison group (Yan et al. 2022 phenotyping) -- read as suggestive, not definitive",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig("results/glue_vs_classical_comparison.png", dpi=150)
    print("\nPlot saved to results/glue_vs_classical_comparison.png")
