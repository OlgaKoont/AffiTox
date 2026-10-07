# AffiTox pipeline wrappers

Executable stages from `run_pipeline.sh`. Scientific interpretation belongs in `analysis/` outputs, not in stage logs.

```text
pipeline/
├── install/setup_environment.sh
├── prepare/run_curate.sh      # layer 1 → 2, assigns ligand_id
├── prepare/run_prepare.sh
├── dock/run_docking.sh
├── hpc/                       # optional SLURM examples, not the public interface
└── postprocess/
    ├── run_merge.sh
    ├── run_posebusters.sh
    ├── run_article_analysis.sh
    └── run_si_full.sh
```

Tracked extras: `rename_zenodo_affitox.py`, `upload_zenodo.py` (Zenodo packaging; not part of analysis reproduction).

```bash
source config/project.env.sh          # TOXAFFINITY_ROOT, TARGETS, METHODS, paths
bash run_pipeline.sh --help           # stage list; default stage is analysis
bash run_pipeline.sh analysis         # manuscript analysis from deposited tables
```

`bash run_pipeline.sh all` is prepare, dock, merge, posebusters, analysis. It does not install packages and does not run Boltz-2.

Stage contracts, gaps, and configs: [docs/PIPELINE.md](../docs/PIPELINE.md). Per-stage notes: [install](install/README.md), [prepare](prepare/README.md), [dock](dock/README.md), [postprocess](postprocess/README.md).
