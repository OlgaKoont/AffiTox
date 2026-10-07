# Data availability (manuscript template)

Fill the two placeholders after (1) the GitHub `v1.0.0` software archive exists on
Zenodo and (2) the umbrella data-manifest record exists. Do **not** create the
GitHub release until CI is green, including the analysis smoke test.

The source code, curated benchmark inputs, and figure-generation workflows are
publicly available at https://github.com/OlgaKoont/AffiTox and archived as
version 1.0.0 at Zenodo (software DOI: [DOI]). Data are six **related** Zenodo
records (raw ChEMBL+PDB, prepared inputs and weights, docking dumps, merged
tables, PoseBusters/conversion, metrics and PNG figures), indexed from one
community/project (data DOI: [DOI]). Layout: [`docs/zenodo/README.md`](zenodo/README.md)
and `release/zenodo/`.

The analysis workflow can be reproduced from record 04–06 without rerunning
docking. GPU re-runs of Boltz-2 / DynamicBind are not bit-identical; record 03
is the frozen scientific snapshot. Full docking requires third-party methods and
checkpoints (NOTICE; [`SOFTWARE_REGISTRY.md`](SOFTWARE_REGISTRY.md)). ChEMBL and
PDB redistributions remain under their original terms.

The previous per-protein zips `10.5281/zenodo.20825057`–`10.5281/zenodo.20825067`
are **historical only** (old panel including 7awe/2z5x/3mjg). Do not announce
them as the current AffiTox datasets.

Do not add `.zenodo.json` alongside `CITATION.cff`: if both exist, Zenodo ignores
`CITATION.cff`.
