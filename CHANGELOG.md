# Changelog

## v2 (19 October 2026) - label correction

- The v1 dataset filtered lipidation notes with the substring `"palmitoyl"`, admitting
  2,381 N-/O-palmitoyl sites (35% of positives). v2 keeps only S-palmitoyl cysteine on a
  cysteine residue: **4,424 positives** (449 experimental, 3,975 inferred), 67,266 sites,
  21,568 proteins, 786 organisms, 557 human positives.
- Splits: every protein (including background negatives) is clustered; groups sharing an
  exact window are merged. Cross-fold identical windows: **0** (v1: 751 windows).
- Results corrected throughout: cluster CV AUC 0.944 (0.940-0.948), MCC 0.611; human
  0.914 (0.900-0.927); species holdout 0.885; hard negatives 0.956.
- The "UniProt annotation error" claim for the 66 non-cysteine sites is withdrawn: those
  were legitimate N6-/N-/O-palmitoyl annotations admitted by our filter.
- Cross-species matrix reduced to six species; bacterial N-palmitoyl sites removed.
- Proteome scan rescored with the forest view (ESM view not tractable at proteome scale
  on the available CPU host); tiers reported as candidate rankings.
- ESM-2 view omitted from the scan; the released prediction service still uses the full
  sequence ensemble.
