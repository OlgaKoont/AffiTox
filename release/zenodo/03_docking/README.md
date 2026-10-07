# AffiTox Zenodo record 03 — docking outputs

Frozen poses, logs, and native engine JSON/SDF for the 16-target panel.

Re-running Boltz-2 or DynamicBind on a GPU is **not bit-identical**. This record is the
scientific snapshot. Layers 4–6 must be rebuilt from these files, not from a new GPU run,
if you need the same numbers.

Run `pack_from_workspace.sh` on the cluster to zip `results/<pdb>/`.
