"""
build_dataset.py

Builds the full anagrelide analog dataset (anagrelide + 22 analogs, A1-A22)
from Supplementary Table 2 of Chen et al., Nature Communications 2021,
12:6204 (PMC8551160), using RDKit to programmatically construct each
structure and compute descriptors.

METHOD NOTES (read before modifying):

1. Structures are built via RDKit RWMol.AddBond on a verified scaffold,
   NEVER by hand-editing SMILES text. A hand-edited SMILES string was
   found to silently produce a wrong structure (missing one ring
   chlorine) that still parsed without any RDKit error -- this is why
   direct atom-level construction is used throughout.

2. The anagrelide core SMILES was verified against PubChem (CID 2182,
   which is the same compound as CID 135409400 -- same InChIKey
   OTBXOEAOVRKTNQ-UHFFFAOYSA-N; PubChem simply has two CIDs for one
   structure). Formula, molecular weight, and RDKit LogP all matched
   independently-sourced values before this scaffold was trusted.

3. R1 attaches at ring atom index 6, R2 at ring atom index 9, in the
   scaffold SMILES "O=C1CN2Cc3cc(Cl)ccc3NC2=N1" (freshly parsed each
   time -- RDKit atom indexing is stable for a given SMILES string but
   should not be assumed stable across edits).

4. Two structural corrections were made during construction after
   comparing built structures' RDKit LogP against the paper's reported
   LogP, and re-inspecting the table images at higher resolution when
   discrepancies appeared:
     - R2 position was initially confused with a fixed ring chlorine;
       resolved by checking that anagrelide itself has R2=H (per the
       table's own text data) and testing candidate atoms against A5
       (R1=Cl, R2=Cl, paper LogP=2.64).
     - A18 was first misread as a benzyl group (Ph-CH2-) from a
       low-resolution image; a higher-resolution re-crop showed it is
       actually a directly-attached ortho-tolyl (2-methylphenyl) group.

5. Some discrepancies between our RDKit LogP and the paper's reported
   LogP remain even after structural correction, specifically for
   N- and O-containing heterocyclic substituents: A8 (furan, diff
   +0.96), A10 (pyridine, diff +0.67), A18 (diff -0.25, smaller but
   still notable). These are consistent in direction and are most
   likely due to RDKit's Crippen LogP algorithm handling heteroatom
   contributions differently than whatever tool the paper's authors
   used (likely ChemDraw or similar) -- not structural errors, since
   ring-attachment isomer choice (e.g. thiophen-2-yl vs 3-yl) was
   shown to make zero difference to LogP for symmetric-enough cases
   (A7, A11). This is documented here as a known limitation, not
   silently corrected by forcing agreement.

6. Three compounds have "ND" (not detected) IC50 in the original assay:
   A11, A16, A21, A22 -- handled here as Python None, not silently
   dropped or treated as zero.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen, Descriptors
import csv

SCAFFOLD_SMILES = "O=C1CN2Cc3cc(Cl)ccc3NC2=N1"
R1_ATOM_IDX = 6
R2_ATOM_IDX = 9


def build_analog(r1_smiles, r2_smiles):
    """
    Build an anagrelide analog by attaching r1_smiles at R1 (ring
    position 7) and r2_smiles at R2 (ring position 8), via direct
    RDKit atom bonding. Either may be None (meaning H, i.e. no
    substituent added at that position).

    Returns (mol, error_message).
    """
    scaffold = Chem.MolFromSmiles(SCAFFOLD_SMILES)
    combined = Chem.RWMol(scaffold)

    if r1_smiles is not None:
        frag1 = Chem.MolFromSmiles(r1_smiles)
        if frag1 is None:
            return None, f"bad R1 fragment SMILES: {r1_smiles}"
        offset1 = combined.GetNumAtoms()
        combined.InsertMol(frag1)
        combined.AddBond(R1_ATOM_IDX, offset1, Chem.BondType.SINGLE)

    if r2_smiles is not None:
        frag2 = Chem.MolFromSmiles(r2_smiles)
        if frag2 is None:
            return None, f"bad R2 fragment SMILES: {r2_smiles}"
        offset2 = combined.GetNumAtoms()
        combined.InsertMol(frag2)
        combined.AddBond(R2_ATOM_IDX, offset2, Chem.BondType.SINGLE)

    mol = combined.GetMol()
    try:
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as e:
        return None, str(e)


# Full dataset from Supplementary Table 2.
# Format: name -> (R1 SMILES or None, R2 SMILES or None,
#                   IC50 in nM or None if ND, paper LogP, paper CLogP)
DATASET = {
    "Anagrelide": ("Cl", None, 6.67, 2.08, 1.02),
    "A1":  ("Br", None, 2.5, 2.36, 1.152),
    "A2":  ("F", None, 34.2, 1.68, 0.572),
    "A3":  ("C=C", None, 10.40, 2.16, 1.15),
    "A4":  ("c1ccccc1", None, 0.56, 3.20, 2.067),
    "A5":  ("Cl", "Cl", 0.70, 2.64, 1.62),
    "A6":  ("c1ccc(C)cc1", None, 0.30, 3.69, 2.57),
    "A7":  ("c1ccsc1", None, 1.77, 3.13, 1.7214),
    "A8":  ("c1ccoc1", None, 0.60, 1.76, 1.24),
    "A9":  ("c1ccc(Cl)cc1", None, 1.81, 3.76, 2.78),
    "A10": ("c1ccncc1", None, 6.22, 1.86, 0.587),
    "A11": ("c1ccc2ccccc2c1", None, None, 4.2, 3.24),  # ND
    "A12": ("c1ccc(C#N)cc1", None, 2.38, 3.23, 1.50),
    "A13": ("c1ccc(O)cc1", None, 12.30, 2.81, 1.46),
    "A14": (None, None, 42.88, 1.53, 0.43),
    "A15": ("c1ccccc1", "Cl", 1.14, 3.76, 2.53),
    "A16": ("c1ccc(C(=O)O)cc1", None, None, 2.76, 1.83),  # ND
    "A17": ("Br", "Cl", 1.30, 2.91, 1.725),
    "A18": ("c1ccccc1C", None, 1.09, 3.69, 2.27),  # ortho-tolyl, corrected
    "A19": ("c1cccc(C)c1", None, 1.37, 3.69, 2.57),
    "A20": ("c1ccc(CC)cc1", None, 1.80, 4.11, 3.10),
    "A21": ("c1ccc(CCC)cc1", None, None, 4.52, 3.62),  # ND
    "A22": ("c1ccc(C(C)C)cc1", None, None, 4.44, 3.49),  # ND
}


def build_full_dataset():
    """Build all compounds and compute descriptors. Returns a list of dicts."""
    rows = []
    for name, (r1, r2, ic50, paper_logp, paper_clogp) in DATASET.items():
        mol, err = build_analog(r1, r2)
        if mol is None:
            print(f"WARNING: failed to build {name}: {err}")
            continue
        rows.append({
            "compound": name,
            "smiles": Chem.MolToSmiles(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "mol_weight": round(Descriptors.MolWt(mol), 2),
            "rdkit_logp": round(Crippen.MolLogP(mol), 2),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "num_h_donors": rdMolDescriptors.CalcNumHBD(mol),
            "num_h_acceptors": rdMolDescriptors.CalcNumHBA(mol),
            "ic50_nM": ic50,
            "ic50_detected": ic50 is not None,
            "paper_logp": paper_logp,
            "paper_clogp": paper_clogp,
        })
    return rows


if __name__ == "__main__":
    rows = build_full_dataset()

    print(f"Built {len(rows)} / {len(DATASET)} compounds\n")
    print(f"{'Compound':12s} {'RDKit LogP':>11s} {'Paper LogP':>11s} {'Diff':>7s} {'IC50 (nM)':>10s}")
    for row in rows:
        diff = row["rdkit_logp"] - row["paper_logp"]
        ic50_str = f"{row['ic50_nM']:.2f}" if row["ic50_nM"] is not None else "ND"
        flag = "  <-- known heterocycle gap" if abs(diff) > 0.5 else ""
        print(f"{row['compound']:12s} {row['rdkit_logp']:11.2f} {row['paper_logp']:11.2f} {diff:7.2f} {ic50_str:>10s}{flag}")

    out_path = "data/anagrelide_analogs.csv"
    fieldnames = list(rows[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nSaved to {out_path}")
