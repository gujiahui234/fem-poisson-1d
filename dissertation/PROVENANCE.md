# Reuse and packaging changes

The solver was developed from earlier interval and square-domain implementations, then used for the final dissertation. It is a separate experiment suite from the historical v0.1.0 root examples; both remain available.

| Resource | Origin and release treatment |
|---|---|
| `fem_experiments.py` | Existing final-dissertation solver. Only the optional private historical-output comparison was removed. |
| `audit_results.py` | Existing independent result audit. Accepts either the original 151 checks or the standalone 150; manifest paths resolve relative to each isolated run. |
| `audit_rate_bounds.py`, `check_integral_limits.py` | Reused without source changes. |
| `make_template_tables.py`, `make_readability_figures.py` | Reused without source changes from the latest readability revision. |
| `make_profile.py` | Profile plotting block extracted from the earlier `make_assets.py`; it does not assemble or solve a system. |
| `reproduce.py` | New public wrapper: isolated runs, stages, reference comparison and summary-rate verification. |
| `reference_results/` | All 14 numerical result files copied byte-for-byte; `reference_manifest.json` records their SHA-256 values. |

`source_provenance.json` records the original and published hashes of the six principal scripts. No private absolute paths are needed to interpret that record.

Local thesis builders, declaration/template checks and revision scripts are intentionally not execution dependencies of the public numerical package. They refer to the private dissertation and protected-source manifests. Historical generation scripts that could overwrite the document are also excluded. The eight exported tables retain symbolic manuscript cross-references for insertion into the paper; no complete paper is published.

The new wrapper compares recomputed numeric CSV fields with the bundled reference using declared tolerances and reruns independent identities. Figure byte hashes are not cross-platform numerical tests. Release verification reports distinguish saved-data audits from actual FEM recomputation.
