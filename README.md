# P1 finite elements for the 1D Poisson equation

This code solves `-u'' = f` on `[0,1]` with `u(0) = u(1) = 0`, using continuous piecewise-linear finite elements. It reproduces two examples:

* `u = x(1-x)/2`, with `f = 1`.
* `u = sin(pi*x)`, with `f = pi^2 sin(pi*x)`.

Both use 10, 20, 40, 80 and 160 elements. The program assembles the stiffness matrix and load vector, solves the interior system, and then evaluates the errors. It does not substitute exact nodal values for the computed solution.

[中文说明](README.zh-CN.md)

## Run the release

Install **Python 3.12**. Download and extract `fem-poisson-1d-v0.1.0.zip` from [Releases](https://github.com/gujiahui234/fem-poisson-1d/releases), then open a terminal in the extracted folder.

Windows Command Prompt:

```cmd
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.venv\Scripts\python.exe reproduce.py
```

macOS or Linux:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
.venv/bin/python reproduce.py
```

The final line should include `REPRODUCTION PASSED`. Installation needs internet access; the calculations run offline. No Anaconda, Word, Codex, course files or API keys are required.

`requirements-lock.txt` pins the tested dependencies, including transitive dependencies. `requirements.txt` lists only the four direct dependencies. Python 3.12 on Windows is tested locally; GitHub Actions runs the same commands on Windows and Linux. Check the Actions result for those remote runs.

## Results

The default run writes to `results/`:

| File | Contents |
| --- | --- |
| `errors.csv` | Maximum, L2 and energy-seminorm errors, with observed rates |
| `diagnostics.csv` | Nodal errors, residuals, quadrature and maximum-search comparisons |
| `validation.json` | Numerical checks with measured values and tolerances |
| `nodal_solutions.npz` | Meshes, computed nodal values and interior loads |
| `figures/` | Two simple figures, each in PNG and PDF |
| `metadata.json` | Runtime versions, numerical settings and file hashes |
| `reproduction_summary.json` | Comparison with the included reference data |

For `N=160`, the maximum errors are approximately `4.8828125e-6` (quadratic) and `4.8188718e-5` (sine). The observed rates tend to 2 for the maximum and L2 errors, and 1 for the energy seminorm. The full reference table is [examples/expected_errors.csv](examples/expected_errors.csv).

`reproduce.py` runs 118 numerical checks and 23 tests, then compares the ten result rows with the reference table. Small floating-point differences are allowed (`rtol=2e-6`, `atol=1e-12`). Figure bytes may differ between platforms and Matplotlib versions; numerical values are the comparison target.

## Other runs

With the virtual environment activated:

```sh
python code/run_experiments.py --case sine --meshes 10 20 40 80 160 --output-dir results_sine
python code/run_tests.py
```

Available cases are `quadratic`, `sine` and `both`. Mesh sizes must be increasing integers of at least 2. The acceptance thresholds and order checks were selected for the default smooth examples; a different mesh sequence may fail these checks even if the solver is behaving correctly.

## Numerical details

Loads and squared integral errors use five-point Gauss quadrature, checked against ten points. Maximum errors are found from element endpoints and roots of `u' - u_h'`. This search is complete for the two supplied concave exact solutions, up to root-finding and rounding error. Adding an arbitrary oscillatory solution requires a more general maximum search.

Checks include the two-element result `U=1/8`, maximum error `1/32`, exact quadratic error formulas and independent sine integral-error formulas evaluated using Decimal. A deliberately underintegrated one-point load demonstrates why nodal exactness depends on the exact Galerkin conditions.

The examples are limited to 1D, constant coefficients, homogeneous Dirichlet data and uniform meshes. These calculations do not establish a general theorem or a 2D result.

## Files and provenance

`code/run_experiments.py` contains the numerical calculation and plotting. `code/reference_errors.py` contains the independent sine formulas. `code/test_experiments.py` checks the solver; `code/test_evidence.py` checks stale or inconsistent evidence. Development used AI assistance; the numerical checks and raw outputs are included so the calculations can be inspected and rerun.

This repository contains code and reference numerical data. The student's report and course notes are not distributed here. Code is released under the [MIT license](LICENSE).
