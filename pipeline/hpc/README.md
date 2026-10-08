# Optional HPC wrappers

The public interface is portable bash:

```bash
source config/project.env.sh
bash run_pipeline.sh curate
bash run_pipeline.sh merge
```

Scripts in this folder are **examples** for the aichem/aihub cluster (`sbatch`,
partition `aichem`). On another machine, keep the Python commands and replace
only the scheduler header (partition, account, modules, GPU request).

Never use these SLURM files as the only documented entrypoint.
