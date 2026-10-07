# AffiTox GitHub example (16 targets × 5 ligands × 5 methods)

Same six information blocks as the Zenodo records, truncated to `ligand_0001` … `ligand_0000005`
for every target and every method. Seed 42. If a method failed one of these five ligands, the
row stays with a missing score; do not substitute `ligand_0006`.

```bash
source config/project.env.sh
# EXAMPLE_N=5  # first five ligand_id rows; omit the variable for the full panel
bash pipeline/prepare/run_curate.sh
PYTHONPATH=src python src/analysis/ligand_identity.py
bash pipeline/postprocess/run_merge.sh   # full panel → analysis/excluding_2z5x_3mjg/tables
EXAMPLE_N=5 bash pipeline/postprocess/run_merge.sh --output-dir example/04_merged
```

Full-panel PoseBusters conversion and docking dumps are on Zenodo records 03 and 05
(not the historical DOIs 10.5281/zenodo.20825057–067). Boltz-2 / DynamicBind folders
for this 16×5 slice stay on disk locally; they are gitignored (too large for GitHub).
