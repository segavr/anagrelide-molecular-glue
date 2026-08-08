"""
run_all.py

Master script that runs the full analysis pipeline in order, from raw
structure-building through the final comparison with classical PDE3A
inhibitors and the real-structure verification.

This does not reimplement the individual analyses -- it calls each
already-validated script in sequence, so each step remains independently
readable, testable, and re-runnable on its own. See each script's own
docstring for the detailed methodology and validation notes.

Pipeline order:
  1. build_dataset.py               -- builds and validates all 23 molecular
                                        glue structures (anagrelide + analogs)
  2. correlation_analysis.py        -- reproduces and critically tests the
                                        paper's LogP-IC50 claim
  3. multi_descriptor_analysis.py   -- tests whether other descriptors add
                                        information beyond LogP
  4. classical_inhibitors_comparison.py -- builds classical inhibitor
                                        dataset, compares descriptor profiles
  5. pdb_structural_analysis.py     -- analyzes real 3D structure (7EG0)
                                        NOTE: requires data/7EG0.pdb to be
                                        present already (downloaded manually
                                        from files.rcsb.org -- see that
                                        script's docstring for why this
                                        can't be automated here)

Usage: run from the project root.
    python scripts/run_all.py
"""

import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent

PIPELINE = [
    "build_dataset.py",
    "correlation_analysis.py",
    "multi_descriptor_analysis.py",
    "classical_inhibitors_comparison.py",
    "pdb_structural_analysis.py",
]


def run_step(script_name):
    script_path = SCRIPTS_DIR / script_name
    print(f"\n{'=' * 70}")
    print(f"Running: {script_name}")
    print("=" * 70)

    if script_name == "pdb_structural_analysis.py" and not Path("data/7EG0.pdb").exists():
        print("SKIPPED: data/7EG0.pdb not found. This optional step requires")
        print("manually downloading the structure from:")
        print("  https://files.rcsb.org/download/7EG0.pdb")
        print("and saving it to data/7EG0.pdb before running this step.")
        return True

    result = subprocess.run([sys.executable, str(script_path)], capture_output=False)
    return result.returncode == 0


if __name__ == "__main__":
    print("Anagrelide molecular glue analysis -- full pipeline")

    failed = []
    for script_name in PIPELINE:
        ok = run_step(script_name)
        if not ok:
            failed.append(script_name)

    print(f"\n{'=' * 70}")
    if failed:
        print(f"Pipeline finished with errors in: {', '.join(failed)}")
        sys.exit(1)
    else:
        print("Pipeline completed successfully. Outputs:")
        print("  data/anagrelide_analogs.csv")
        print("  data/classical_inhibitors.csv")
        print("  results/logp_ic50_correlation.png")
        print("  results/multi_descriptor_comparison.png")
        print("  results/glue_vs_classical_comparison.png")
