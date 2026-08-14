# Pipeline Install Stage

This directory contains bootstrap wrappers for setting up a runnable project environment.  
It is intended for users who start from a fresh machine/session and need deterministic package resolution before running benchmark stages.

## Repository layout

```text
pipeline/install/
└── setup_environment.sh      editable install wrapper
```

## Minimal run

```bash
source config/project.env.sh
bash pipeline/install/setup_environment.sh
```

## Output contract

| Script | Input | Output | Use case |
|---|---|---|---|
| `setup_environment.sh` | `PYTHON`, project root, optional `environment.yml` | editable package install (`pip install -e .`) | enable CLI/modules for all pipeline stages |

## Interpretation scope

- This stage validates runtime availability, not scientific quality.
- Successful install does not guarantee external tools (GNINA/QVina/DynamicBind/etc.) are installed/configured.

## Reproducibility

Key parameters and effects:
- `PYTHON=/path/to/python`  
  - Controls interpreter and dependency versions; changing it can change numeric outputs downstream.
- Active conda/venv before running script  
  - Determines where packages are installed; wrong env often causes missing binary/module failures later.

Examples:

```bash
# Install into the Python that is already on PATH (conda/venv)
PYTHON="$(command -v python)" \
bash pipeline/install/setup_environment.sh
```

```bash
# Standard default from project.env.sh
bash run_pipeline.sh install
```

## Evaluation status

- Implemented now:
  - editable install bootstrap with optional conda hint.
- Planned / evolving:
  - optional automated validation of external binaries per method.

## Citation, license, contribution

### Citation

No method-specific citation required for install stage itself; cite methods used in later stages.

### License

Follow repository license and package/tool licenses installed in the environment.

### Contribution

- Keep install logic minimal and non-destructive.
- If adding checks, make failures actionable (exact package/binary and fix).
