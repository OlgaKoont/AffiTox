# AffiTox v1.0.0 data-release manifest

Umbrella index for the six open Zenodo raw-data records. Upload this directory
(`README.md` + `affitox_data_manifest.tsv`) as a small Zenodo **Dataset** record
and link the six parts with **Has part**. Do not re-upload the ~202 GB archives.

## Canonical code

https://github.com/OlgaKoont/AffiTox

## Open data parts

| Part | DOI | Status |
|------|-----|--------|
| 1/6 | [10.5281/zenodo.20825057](https://doi.org/10.5281/zenodo.20825057) | open |
| 2/6 | [10.5281/zenodo.20825059](https://doi.org/10.5281/zenodo.20825059) | open |
| 3/6 | [10.5281/zenodo.20825061](https://doi.org/10.5281/zenodo.20825061) | open |
| 4/6 | [10.5281/zenodo.20825063](https://doi.org/10.5281/zenodo.20825063) | open |
| 5/6 | [10.5281/zenodo.20825065](https://doi.org/10.5281/zenodo.20825065) | open |
| 6/6 | [10.5281/zenodo.20825067](https://doi.org/10.5281/zenodo.20825067) | open |

Filenames, sizes, and MD5 checksums: `affitox_data_manifest.tsv` (queried from the
Zenodo API; about 202 GB in total). All 16 AffiTox targets are covered.

## License caveat

The CC BY 4.0 license applies to the original AffiTox curation, organization, and
derived results. Redistributed third-party source records remain subject to their
original terms and licenses.

## Metadata edits on the six existing records

Do not re-upload binaries. In the Zenodo web UI, set:

- related GitHub URL: `https://github.com/OlgaKoont/AffiTox` (not `docking-benchmark`);
- version: `1.0.0`;
- creator ORCID for Olga Konovalova: `0000-0003-0391-3211`;
- community: `affitox` (remove unused `toxdock-bench` if still attached);
- creators/contributors matching actual data-generation credit.
