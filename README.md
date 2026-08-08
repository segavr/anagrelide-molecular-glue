# Anagrelide as a PDE3A-SLFN12 Molecular Glue: A Structure-Activity and Descriptor Analysis

## Research Question

Can the structure-activity relationships of the anagrelide molecular-glue series be quantitatively verified and extended using RDKit descriptors, and do they structurally distinguish this class from classical (non-glue) PDE3A inhibitors?

## Motivation

Anagrelide, alongside DNMDP and nauclefine, induces apoptosis by acting as a molecular glue between PDE3A and SLFN12 — a rare, structurally characterized example of this drug mechanism (Chen et al., Nature Communications 2021, PMC8551160), with a real cryo-EM structure available (PDB 7EG0). This makes it possible to cross-check 2D descriptor-based analysis directly against 3D experimental structure.

## Data Source

Chen et al., Nat Commun 2021, 12:6204, Supplementary Table 2 (22 anagrelide analogs + anagrelide itself, with measured IC50 against PDE3A-SLFN12-dependent apoptosis).

## Methods

- All 23 compound structures built programmatically in RDKit (direct atom bonding, not text-edited SMILES), independently validated against the paper's reported LogP
- Statistical analysis: Pearson correlation, multiple linear regression, t-tests (SciPy)
- Classical PDE3A inhibitor structures (milrinone, amrinone, enoximone, cilostazol, cilostamide, imazodan) each verified against 2+ independent sources
- Real 3D structure analysis using PDB 7EG0 (cryo-EM structure of the anagrelide-PDE3A-SLFN12 complex)

## Key Results

1. **Structure verification caught real errors**: three structural ambiguities were found and corrected during dataset construction by cross-checking against the paper's reported LogP — see notebooks/01_analysis_and_results.ipynb for details.

2. **The paper's own "no correlation" claim doesn't fully hold up to a formal test**: both Pearson (linear) and Spearman (rank-based, outlier-robust) correlation between LogP and log10(IC50) are statistically significant on the full dataset (Pearson r = -0.623 to -0.750, p < 0.005; Spearman r = -0.504 to -0.606, p < 0.03, depending on LogP source), contradicting the qualitative statement in Supplementary Note 1. Removing two outlier compounds (A14, A2) drops both statistics well below significance for the paper's own LogP (Pearson p = 0.095, Spearman p = 0.232) — the fact that even the outlier-robust Spearman test loses significance strengthens the case that the authors' qualitative judgment reflects a genuine weak relationship in the remaining 17 compounds, not just a visual artifact of two extreme points.

3. **Within this series, lipophilicity shows a stronger association with activity than the evaluated polarity-related descriptors**: LogP and molecular weight are redundant predictors (r = 0.879 with each other); TPSA and H-bonding descriptors show no significant correlation with IC50. This finding is specific to the ~19-compound anagrelide analog series analyzed here and should not be generalized beyond it.

4. **Molecular glues are structurally distinguishable from classical PDE3A inhibitors**: significantly higher LogP (p = 0.035) and much lower TPSA (p = 0.0002, near-complete separation).

5. **Independent 3D structural confirmation**: real atomic distances in PDB 7EG0 confirm the paper's stated binding mechanism (H-bonds to H961/Q1001, pi-stacking with F1004, hydrophobic contact with SLFN12's I557/I558) — corroborating the 2D descriptor findings from a completely different type of evidence.

See [notebooks/01_analysis_and_results.ipynb](notebooks/01_analysis_and_results.ipynb) for the full analysis with tables and plots.

## Repository Structure

```
├── data/
│   ├── anagrelide_analogs.csv       # 23 molecular glue structures + descriptors
│   ├── classical_inhibitors.csv     # 6 classical PDE3A inhibitors + descriptors
│   └── 7EG0.pdb                     # Cryo-EM structure (downloaded from RCSB)
├── scripts/
│   ├── build_dataset.py             # Builds and validates all 23 structures
│   ├── correlation_analysis.py      # Reproduces/tests paper's LogP-IC50 claim
│   ├── multi_descriptor_analysis.py # Multi-descriptor comparison
│   ├── classical_inhibitors_comparison.py  # Glue vs classical inhibitor comparison
│   ├── pdb_structural_analysis.py   # Real 3D structure analysis
│   └── run_all.py                   # Runs the full pipeline in sequence
├── notebooks/
│   └── 01_analysis_and_results.ipynb
├── results/
│   ├── logp_ic50_correlation.png
│   ├── multi_descriptor_comparison.png
│   └── glue_vs_classical_comparison.png
└── README.md
```

## Reproducing the Results

```bash
# Requires RDKit, scipy, matplotlib, pandas
python scripts/run_all.py
```

Note: pdb_structural_analysis.py requires data/7EG0.pdb, which must be downloaded manually from https://files.rcsb.org/download/7EG0.pdb (this step is optional and will be skipped automatically if the file isn't present).

## Limitations

- The classical-inhibitor comparison group is small (n=6); p-values should be read as suggestive, not conclusive
- Four compounds (A11, A16, A21, A22) have no detected IC50 and were excluded from correlation analyses
- RDKit's Crippen LogP shows systematic discrepancies (~0.6-1.0 units) from the paper's reported LogP specifically for N/O-heterocyclic substituents, most likely reflecting different LogP algorithms rather than structural errors
- This is a retrospective analysis of published data, not new experimental or computational (docking/DFT) work

## References

Chen, J., Liu, N., Huang, Y. et al. Structure of PDE3A-SLFN12 complex and structure-based design for a potent apoptosis inducer of tumor cells. Nat Commun 12, 6204 (2021). https://doi.org/10.1038/s41467-021-26546-8
