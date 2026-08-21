# Install stage

```bash
source config/project.env.sh
bash pipeline/install/setup_environment.sh
# or
bash run_pipeline.sh install
```

Runs `pip install -e ".[analysis]"` with fallbacks to `.[rdkit]` then `-e .` (`pipeline/install/setup_environment.sh`). This does **not** install GNINA, QVina2, Boltz-2, DynamicBind, PLAPT, or PoseBusters.

Set `PYTHON` to the interpreter that should receive the package. Details: [docs/INSTALL.md](../../docs/INSTALL.md).
