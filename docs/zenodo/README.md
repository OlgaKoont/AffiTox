# AffiTox data release — six related Zenodo records

Upload **six sibling records** (one information block each), then link them from
one Zenodo community/project. Do **not** split by protein (the old DOIs
10.5281/zenodo.20825057–067 were per-target zips).

Canonical GitHub: https://github.com/OlgaKoont/AffiTox

| Record | Content | Workspace tree |
|--------|---------|----------------|
| 01_raw | ChEMBL API snapshots + 16 PDB files | `release/zenodo/01_raw/` |
| 02_prepared | curated `ligand_id` tables, maps, boxes, receptors, weights | `release/zenodo/02_prepared/` |
| 03_docking | frozen docking dumps; sharded 03a–03e (Zenodo 50 GiB limit) | `release/zenodo/03_docking/` + `payloads/03/` |
| 04_merged | `merged_ligands_docking_<pdb>.csv` (join on ligand_id) | `release/zenodo/04_merged/` |
| 05_posebusters | full-panel PoseBusters CSVs + converted SDF | `release/zenodo/05_posebusters/` |
| 06_metrics_figures | correlation / nEF / inferential tables + PNG | `release/zenodo/06_metrics_figures/` |

Each record must include `NOTICE`, `README.md`, and `SHA256SUMS.txt`.

Panel PDB order: 1g5m, 2v5z, 3eyg, 3jy9, 3lxk, 11ue, 4ase, 4f65, 4tz4, 4zau, 5jkv, 5mo4, 6gqj, 6jok, 5lf3, 7kk3.

Canonical analysis root: `analysis/excluding_2z5x_3mjg/`. Publication heatmaps:
`figures/by_protein/correlations/heatmaps/pearson_heatmap_ultramarine_ivory_vermilion_square_cells.png`
(and Spearman / Kendall / combined in the same ultramarine square-cell style). PNG only.

GPU docking (Boltz-2, DynamicBind) is not bit-identical on rerun. Record 03 is the
scientific snapshot. Rebuild 04–06 from those files for matching tables.

Fill `N_panel.txt` after the ligand_id merge (do not freeze 10497).
