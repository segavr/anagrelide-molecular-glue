"""
build_analogs.py

Programmatically builds anagrelide analog structures by substituting the R1
group (position 7 on the quinazoline ring) using RDKit molecule editing
(RWMol + AddBond), NOT text-based SMILES editing.

IMPORTANT LESSON LEARNED (see project notes): manually editing a SMILES
string by hand (e.g. deleting/replacing characters) silently produced a
WRONG structure (lost one of the two chlorines) that still parsed
successfully in RDKit without any error. Text-editing SMILES is NOT safe
even when the result parses. Always verify against an independent source
(here: PubChem CID 2182/135409400 for anagrelide itself) before trusting
a hand-built or hand-edited SMILES string.

Verified against Supplementary Table 2 (Chen et al., Nat Commun 2021,
PMC8551160):
  - Anagrelide (R1=Cl): formula C10H7Cl2N3O matches PubChem exactly;
    RDKit LogP = 2.12 vs paper's reported LogP = 2.08 (close; different
    LogP algorithms are expected to differ slightly)
  - A1 (R1=Br): RDKit LogP = 2.23 vs paper 2.36
  - A2 (R1=F): RDKit LogP = 1.60 vs paper 1.68
  - A3 (R1=vinyl): RDKit LogP = 2.11 vs paper 2.16
  - A14 (R1=H): RDKit LogP = 1.46 vs paper 1.53
All within ~0.1-0.15 LogP units, consistent with known algorithm variance.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Crippen

# Anagrelide scaffold with R1 position (ring position 7) set to H.
# This SMILES was verified programmatically: removing atom index 7 (the
# variable chlorine) from the canonical anagrelide SMILES
# "C1C2=C(C=CC(=C2Cl)Cl)NC3=NC(=O)CN31" and re-sanitizing gives this
# exact string, with formula C10H8ClN3O (one Cl remaining at the fixed
# ring position, R1 position now bearing an aromatic H).
SCAFFOLD_R1_IS_H = "O=C1CN2Cc3cc(Cl)ccc3NC2=N1"

# In this scaffold (as parsed fresh each time), atom index 6 is the
# aromatic carbon bearing the free H at the R1 attachment point.
R1_ATTACHMENT_ATOM_IDX = 6


def build_r1_analog(r1_group_smiles):
    """
    Build an anagrelide analog by attaching r1_group_smiles at the R1
    position via direct atom bonding (RWMol.AddBond), not string editing.

    r1_group_smiles: SMILES of the substituent as a standalone fragment,
    written starting from its attachment atom (e.g. 'Cl', 'Br', 'F',
    'C=C' for vinyl, 'c1ccccc1' for phenyl -- RDKit will bond the first
    atom of the parsed fragment to the scaffold).

    Returns (mol, error_message). mol is None if construction failed.
    """
    scaffold = Chem.MolFromSmiles(SCAFFOLD_R1_IS_H)
    frag = Chem.MolFromSmiles(r1_group_smiles)
    if frag is None:
        return None, f"Could not parse R1 fragment SMILES: {r1_group_smiles}"

    combined = Chem.RWMol(scaffold)
    frag_offset = combined.GetNumAtoms()
    combined.InsertMol(frag)
    combined.AddBond(R1_ATTACHMENT_ATOM_IDX, frag_offset, Chem.BondType.SINGLE)

    mol = combined.GetMol()
    try:
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as e:
        return None, str(e)


def descriptors_for(mol):
    """Return (formula, RDKit LogP) for a molecule."""
    formula = rdMolDescriptors.CalcMolFormula(mol)
    logp = Crippen.MolLogP(mol)
    return formula, logp


if __name__ == "__main__":
    # Sanity-check against known compounds from the paper before trusting
    # this on the full 22-analog dataset.
    test_cases = {
        "Anagrelide": ("Cl", 2.08),
        "A1": ("Br", 2.36),
        "A2": ("F", 1.68),
        "A3": ("C=C", 2.16),
        "A14": ("[H]", 1.53),
    }

    print(f"{'Compound':12s} {'R1':8s} {'Formula':18s} {'RDKit LogP':>11s} {'Paper LogP':>11s} {'Diff':>7s}")
    for name, (r1, paper_logp) in test_cases.items():
        mol, err = build_r1_analog(r1)
        if mol is None:
            print(f"{name:12s} ERROR: {err}")
            continue
        formula, logp = descriptors_for(mol)
        print(f"{name:12s} {r1:8s} {formula:18s} {logp:11.2f} {paper_logp:11.2f} {logp - paper_logp:7.2f}")
