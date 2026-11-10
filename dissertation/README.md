# Dissertation experiments: regularity and maximum-norm errors

This directory reproduces the numerical resources used in the final dissertation. The original smooth 1D examples remain available at the repository root. This extension studies P1 Poisson elements on an interval, with the unit square as a numerical extension.

- Interval solutions: `sin(pi*x)` and `x**(5/3)-x`; meshes 8, 16, 32, 64, 128.
- Square solutions: each interval solution multiplied by `sin(pi*y)`; grids 8, 16, 32, 64, split southwest to northeast.
- Triangle quadrature degrees 2 through 5 mean polynomial exactness, with 3, 4, 6, 7 points respectively. Degree 5 supplies the main comparison.
- Four interval errors, square integral errors, sampled maxima, nodal discrepancies, independent reference quadrature and sampling sensitivity are retained.

[中文说明](README.zh-CN.md) · [Reuse and provenance](PROVENANCE.md)

## Run (Python 3.12)

From the repository root, install `dissertation/requirements-lock.txt` in a virtual environment. No course files, local credentials, LaTeX installation or Codex installation are required.

Windows Command Prompt:

```cmd
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r dissertation\requirements-lock.txt
.venv\Scripts\python.exe dissertation\reproduce.py
.venv\Scripts\python.exe dissertation\reproduce.py --recompute
```

Linux/macOS:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r dissertation/requirements-lock.txt
.venv/bin/python dissertation/reproduce.py
.venv/bin/python dissertation/reproduce.py --recompute
```

The default command audits the included saved results and rebuilds eight table fragments and three figures. It does **not** solve the FEM systems. `--recompute` assembles and solves all systems, then runs the same audits and compares with the saved reference. Both modes end with `REPRODUCTION PASSED`. Each creates a new directory in `dissertation/runs/`; `--output NEW_DIRECTORY` selects another new directory and refuses an existing one.

## Inspect the output

| Path inside a run | Contents |
|---|---|
| `results/` | CSV errors and rates, quadrature definitions, eight NPZ solution vectors, experiment validation |
| `qa/results_audit/audit.json` | 144 checks, including independent high-precision norm identities and stored-vector diagnostics |
| `qa/rate_bound_audit.json` | 8 checks of reference sampled-rate intervals |
| `qa/integral_limits_audit.json` | 15 checks of integral constants and normalized values |
| `qa/readability_figures.json` | Reference slopes and observed final rates |
| `tables/` | Eight LaTeX table fragments; cross-references target the dissertation and are not a standalone document |
| `figures/` | Three figures in PNG and PDF; the profile sampling is illustrative, not the reported maximum-search algorithm |
| `reproduction_summary.json` | Independent reference comparison, table-rate checks and pass/fail status |
| `logs/` | Logs for each stage |

The distributed reference preserves 151 historical experiment assertions. A standalone recomputation runs 150: one comparison with a private historical smooth-2D output has been removed. The numerical algorithms and tolerances are otherwise retained. The 144/8/15 audits are rerun in both modes; these numbers count assertions, not independent experiments.

Reference comparisons use relative tolerance `2e-6` and absolute tolerance `1e-12`; raw NPZ/PDF bytes need not agree between platforms. Immutable bundled reference hashes are checked separately. Dependencies are pinned to the tested versions; CI runs the recomputation on Windows and Ubuntu. Consult the current Actions result rather than assuming a historical pass applies to a new commit.

## Interpretation limits

Exact interval Galerkin integration gives nodal exactness for this constant-coefficient, zero-Dirichlet P1 problem. The power solution has maximum-error powers 5/3 and 2/3; the integral norms approach powers 2 and 1. The numerical scripts validate examples; the mathematical proof is not replaced by a table.

The square maximum is sampled on a barycentric lattice. The observed power-case rates near 1.83 are finite-grid observations, not a proved asymptotic exponent. Degree-five quadrature is not exact integration of a fractional-power source. Two-dimensional nodal discrepancies may persist under accurate integration. Dashed 1D plot guides have arbitrary vertical scales and are neither fits nor error bounds; 2D legend rates are measured.

## References and attribution

- S. C. Brenner and L. R. Scott, *The Mathematical Theory of Finite Element Methods*, third edition, Springer, 2008. DOI: [10.1007/978-0-387-75934-0](https://doi.org/10.1007/978-0-387-75934-0).
- D. A. Dunavant, “High degree efficient symmetrical Gaussian quadrature rules for the triangle,” *International Journal for Numerical Methods in Engineering* 21(6), 1129–1148, 1985. DOI: [10.1002/nme.1620210612](https://doi.org/10.1002/nme.1620210612).
- MATH 4890 course notes and assigned Poisson problems informed the study; they are not redistributed.

Development used AI assistance in code, derivations, documentation and verification. Code and numerical reference data use the repository's [MIT license](../LICENSE). Third-party packages retain their own licenses. The dissertation, personal correspondence, recordings, course PDFs, university artwork and local fonts are excluded.
