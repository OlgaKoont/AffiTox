# Data availability (manuscript template)

Fill the two placeholders after (1) the GitHub `v1.0.0` software archive exists on
Zenodo and (2) the umbrella data-manifest record exists. Do **not** create the
GitHub release until CI is green, including the analysis smoke test.

The source code, curated benchmark inputs, manuscript-level analysis tables, and
figure-generation workflows are publicly available in the AffiTox GitHub repository
and archived as version 1.0.0 at Zenodo (software DOI: [DOI]). Raw docking and
prediction outputs are available as six open Zenodo records indexed by the AffiTox
data-release manifest (data DOI: [DOI]). The analysis-level workflow can be
reproduced from the deposited tables without rerunning the computationally intensive
docking calculations. Full pipeline execution requires separately installed
third-party methods and model checkpoints, whose versions and configurations are
documented in the repository. Third-party source data remain subject to their
original licenses.

Until the software and umbrella DOIs exist, cite the six open records
`10.5281/zenodo.20825057` to `10.5281/zenodo.20825067` and
https://github.com/OlgaKoont/AffiTox.

Do not add `.zenodo.json` alongside `CITATION.cff`: if both exist, Zenodo ignores
`CITATION.cff`.
