"""Reproducible first-week P1 Poisson experiments (no changes to source files).

Run with the existing Anaconda interpreter. All results are computed from the
assembled Galerkin system; exact nodal values are used only for diagnostics.
The maximum search is complete for these two concave manufactured solutions:
on each element u' is decreasing, so u' - u_h' has at most one interior zero.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np
import scipy
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from reference_errors import sine_integral_reference
from evidence_utils import RUN_SOURCES, RUN_OUTPUTS, hashes


@dataclass(frozen=True)
class Case:
    name: str
    u: Callable
    du: Callable
    f: Callable


CASES = {
    "quadratic": Case("quadratic", lambda x: x * (1 - x) / 2,
                      lambda x: 0.5 - x, lambda x: np.ones_like(x)),
    "sine": Case("sine", lambda x: np.sin(np.pi * x),
                 lambda x: np.pi * np.cos(np.pi * x),
                 lambda x: np.pi**2 * np.sin(np.pi * x)),
}


@dataclass
class Solution:
    case: Case
    x: np.ndarray
    U: np.ndarray
    A: csr_matrix
    b: np.ndarray

    @property
    def n(self):
        return len(self.x) - 1


def assemble_and_solve(case: Case, n: int, quadrature_order: int = 5) -> Solution:
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError("n must be an integer >= 2 (number of elements)")
    if quadrature_order < 1:
        raise ValueError("quadrature_order must be >= 1")
    x = np.linspace(0.0, 1.0, n + 1)
    xi, weights = leggauss(quadrature_order)
    phi = np.array([(1 - xi) / 2, (1 + xi) / 2])
    ii, jj, vv = [], [], []
    full_b = np.zeros(n + 1)
    for k in range(n):
        h = x[k + 1] - x[k]
        local_A = np.array([[1.0, -1.0], [-1.0, 1.0]]) / h
        qx = (x[k] + x[k + 1]) / 2 + h * xi / 2
        local_b = h / 2 * (phi @ (weights * case.f(qx)))
        for a in range(2):
            full_b[k + a] += local_b[a]
            for b in range(2):
                ii.append(k + a)
                jj.append(k + b)
                vv.append(local_A[a, b])
    full_A = coo_matrix((vv, (ii, jj)), shape=(n + 1, n + 1)).tocsr()
    # Homogeneous Dirichlet values are eliminated, preserving symmetry.
    A = full_A[1:-1, 1:-1].tocsr()
    b = full_b[1:-1].copy()
    U = np.zeros(n + 1)
    U[1:-1] = spsolve(A, b)
    return Solution(case, x, U, A, b)


def maximum_error(sol: Solution) -> tuple[float, float]:
    """Endpoints plus every stationary point, to numerical root precision."""
    best, location = -1.0, 0.0
    for k in range(sol.n):
        left, right = sol.x[k:k + 2]
        slope = (sol.U[k + 1] - sol.U[k]) / (right - left)
        g = lambda z: float(sol.case.du(z) - slope)
        candidates = [left, right]
        if g(left) * g(right) < 0:
            candidates.append(brentq(g, left, right, xtol=5e-15,
                                      rtol=4 * np.finfo(float).eps))
        for z in candidates:
            uh = sol.U[k] + slope * (z - left)
            error = float(abs(sol.case.u(z) - uh))
            if error > best:
                best, location = error, float(z)
    return best, location


def integral_errors(sol: Solution, quadrature_order: int = 5) -> tuple[float, float]:
    xi, w = leggauss(quadrature_order)
    e0_sq, ee_sq = 0.0, 0.0
    for k in range(sol.n):
        left, right = sol.x[k:k + 2]
        h = right - left
        qx = (left + right) / 2 + h * xi / 2
        slope = (sol.U[k + 1] - sol.U[k]) / h
        uh = sol.U[k] + slope * (qx - left)
        e0_sq += h / 2 * np.dot(w, (sol.case.u(qx) - uh)**2)
        ee_sq += h / 2 * np.dot(w, (sol.case.du(qx) - slope)**2)
    return math.sqrt(e0_sq), math.sqrt(ee_sq)


def dense_maximum(sol: Solution, samples: int = 1001) -> float:
    if samples < 2:
        raise ValueError("samples must include both endpoints")
    t = np.linspace(0.0, 1.0, samples)
    best = 0.0
    for k in range(sol.n):
        xx = sol.x[k] + (sol.x[k + 1] - sol.x[k]) * t
        uh = sol.U[k] * (1 - t) + sol.U[k + 1] * t
        best = max(best, float(np.max(np.abs(sol.case.u(xx) - uh))))
    return best


def error_values(sol: Solution, order: int) -> dict:
    ei, xm = maximum_error(sol)
    e0, ee = integral_errors(sol, order)
    return {"E_inf": ei, "E_L2": e0, "E_H1sem": ee, "x_at_E_inf": xm}


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)  # float repr preserves full double precision


def make_figures(rows: list[dict], output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "savefig.dpi": 240})
    output.mkdir(parents=True, exist_ok=True)
    names = list(dict.fromkeys(r["case"] for r in rows))
    fig, ax = plt.subplots(figsize=(6.4, 2.65), layout="constrained")
    for name in names:
        rr = [r for r in rows if r["case"] == name]
        color, marker = ("#2463A0", "o") if name == "quadratic" else ("#333333", "s")
        ax.loglog([r["N"] for r in rr], [r["E_inf"] for r in rr], marker + "-",
                  color=color, markersize=5, linewidth=1.5, label=name.capitalize())
    ticks = sorted({r["N"] for r in rows})
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: str(int(x))))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set(xlabel="Number of elements N", ylabel="Maximum error")
    ax.grid(axis="y", which="major", alpha=0.15)
    ax.legend(frameon=False, fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"convergence.{suffix}", metadata={"Creator": "P1 Poisson experiment"})
    plt.close(fig)

    sol = assemble_and_solve(CASES["quadratic"], 2)
    x = np.linspace(0, 1, 1001)
    uh = np.interp(x, sol.x, sol.U)
    fig, ax = plt.subplots(figsize=(6.4, 2.4), layout="constrained")
    ax.plot(x, sol.case.u(x), label="Exact solution", color="#2463A0", linewidth=1.7)
    ax.plot(x, uh, "--", label="Two-element solution", color="#333333", linewidth=1.5)
    ax.scatter(sol.x, sol.U, color="#333333", s=22, zorder=3)
    ax.set(xlabel="x", ylabel="u", xticks=[0, 0.5, 1], yticks=[0, 0.0625, 0.125])
    ax.legend(fontsize=9, frameon=False)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"two_element_check.{suffix}", metadata={"Creator": "P1 Poisson experiment"})
    plt.close(fig)


def run(cases: list[str], meshes: list[int], output: Path) -> dict:
    if not meshes or any(n < 2 for n in meshes) or meshes != sorted(set(meshes)):
        raise ValueError("meshes must be unique, increasing integers >= 2")
    output.mkdir(parents=True, exist_ok=True)
    rows, diagnostics, checks, arrays = [], [], [], {}

    def check(name, value, limit, comparison="<="):
        passed = value <= limit if comparison == "<=" else value >= limit
        checks.append({"name": name, "value": float(value), "comparison": comparison,
                       "limit": float(limit), "passed": bool(passed)})

    hand = assemble_and_solve(CASES["quadratic"], 2)
    hand_errors = error_values(hand, 5)
    for name, value, exact in [("A", hand.A[0, 0], 4), ("b", hand.b[0], 0.5),
                               ("U", hand.U[1], 0.125), ("E_inf", hand_errors["E_inf"], 1 / 32)]:
        check(f"hand_{name}_absolute_difference", abs(value - exact), 1e-14)
    # A deliberately underintegrated load distinguishes exact Galerkin theory
    # from a numerical method whose load functional has been changed.
    underintegrated = assemble_and_solve(CASES["sine"], 10, 1)
    underintegrated_nodal = float(np.max(np.abs(underintegrated.U - underintegrated.case.u(underintegrated.x))))
    check("one_point_load_exposes_nodal_error", underintegrated_nodal, 1e-4, ">=")
    for name in cases:
        previous = None
        for n in meshes:
            sol = assemble_and_solve(CASES[name], n, 5)
            errors = error_values(sol, 5)
            row = {"case": name, "N": n, "h": 1 / n,
                   **{k: errors[k] for k in ("E_inf", "E_L2", "E_H1sem")}}
            for error_key, rate_key in [("E_inf", "p_inf"), ("E_L2", "p_L2"),
                                         ("E_H1sem", "p_H1sem")]:
                row[rate_key] = (math.log(previous[error_key] / row[error_key]) /
                                 math.log(previous["h"] / row["h"])) if previous else None
            rows.append(row)
            previous = row
            residual = sol.A @ sol.U[1:-1] - sol.b
            residual_abs = float(np.linalg.norm(residual, ord=np.inf))
            norm_A = float(np.max(np.asarray(abs(sol.A).sum(axis=1))))
            norm_b = float(np.linalg.norm(sol.b, ord=np.inf))
            backward = residual_abs / (norm_A * np.max(np.abs(sol.U)) + norm_b)
            nodal = float(np.max(np.abs(sol.U - sol.case.u(sol.x))))
            dense = dense_maximum(sol)
            sol10 = assemble_and_solve(CASES[name], n, 10)
            values10 = error_values(sol10, 10)
            same_solution10 = integral_errors(sol, 10)
            q_change = max(abs(errors[k] - values10[k]) / errors[k]
                           for k in ("E_inf", "E_L2", "E_H1sem"))
            diag = {"case": name, "N": n, "nodal_error": nodal,
                    "boundary_error": float(max(abs(sol.U[0]), abs(sol.U[-1]))),
                    "residual_inf": residual_abs, "relative_residual_b": residual_abs / norm_b,
                    "scaled_backward_error": float(backward),
                    "max_location": errors["x_at_E_inf"], "dense_max_1001": dense,
                    "dense_relative_gap": (errors["E_inf"] - dense) / errors["E_inf"],
                    "load_5_to_10_max_change": float(np.max(np.abs(sol.b - sol10.b))),
                    "nodal_5_to_10_max_change": float(np.max(np.abs(sol.U - sol10.U))),
                    "all_errors_5_to_10_max_relative_change": q_change,
                    "L2_integration_only_relative_change": abs(errors["E_L2"] - same_solution10[0]) / errors["E_L2"],
                    "energy_integration_only_relative_change": abs(errors["E_H1sem"] - same_solution10[1]) / errors["E_H1sem"]}
            diagnostics.append(diag)
            prefix = f"{name}_N{n}"
            arrays[prefix + "_nodes"] = sol.x
            arrays[prefix + "_U"] = sol.U
            arrays[prefix + "_load_interior"] = sol.b
            check(prefix + "_nodal", nodal, 1e-10)
            check(prefix + "_boundary", diag["boundary_error"], 0)
            check(prefix + "_relative_residual", diag["relative_residual_b"], 1e-10)
            check(prefix + "_backward_error", backward, 5e-15)
            check(prefix + "_quadrature_change", q_change, 1e-7)
            check(prefix + "_dense_gap", abs(diag["dense_relative_gap"]), 2e-6)
            check(prefix + "_dense_cannot_exceed_max", dense - errors["E_inf"], 1e-13)
            if name == "quadratic":
                h = 1 / n
                for key, exact in {"E_inf": h * h / 8, "E_L2": h * h / math.sqrt(120),
                                   "E_H1sem": h / math.sqrt(12)}.items():
                    check(prefix + "_formula_" + key, abs(errors[key] / exact - 1), 2e-8)
            else:
                e0_reference, ee_reference = sine_integral_reference(n)
                for key, exact in (("E_L2", e0_reference), ("E_H1sem", ee_reference)):
                    difference = abs(errors[key] / exact - 1)
                    diag["sine_formula_" + key + "_relative_difference"] = difference
                    check(prefix + "_formula_" + key, difference, 2e-8)
            # Keep one rectangular diagnostics table for both cases.
            for key in ("E_L2", "E_H1sem"):
                diag.setdefault("sine_formula_" + key + "_relative_difference", None)
        case_rows = [r for r in rows if r["case"] == name]
        for row in case_rows[-3:]:
            for key, expected in [("p_inf", 2), ("p_L2", 2), ("p_H1sem", 1)]:
                if row[key] is not None:
                    check(f"{name}_N{row['N']}_{key}_difference", abs(row[key] - expected), 0.1)
    write_csv(output / "errors.csv", rows)
    write_csv(output / "diagnostics.csv", diagnostics)
    np.savez_compressed(output / "nodal_solutions.npz", **arrays)
    validation = {"passed": all(c["passed"] for c in checks), "number_of_checks": len(checks),
                  "negative_control": {"case": "sine", "N": 10, "load_quadrature": 1,
                                       "nodal_error": underintegrated_nodal,
                                       "purpose": "Show that underintegration can destroy nodal exactness"},
                  "hand_calculation": {"A": hand.A.toarray().tolist(), "b": hand.b.tolist(),
                                       "U_including_boundary": hand.U.tolist(), **hand_errors},
                  "checks": checks}
    write_json(output / "validation.json", validation)
    import matplotlib
    metadata = {"recorded_at_utc": datetime.now(timezone.utc).isoformat(),
                "python_executable": sys.executable, "python_version": sys.version,
                "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__,
                "matplotlib": matplotlib.__version__, "command_arguments": sys.argv,
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "cases": cases, "meshes": meshes, "arithmetic": "IEEE 754 float64",
                "load_quadrature": 5, "error_quadrature": 5, "comparison_quadrature": 10,
                "dense_points_per_element": 1001,
                "max_method": "endpoints and unique stationary point per element (brentq)",
                "max_root_xtol": 5e-15,
                "rates": "log(E_previous/E_current)/log(h_previous/h_current), before rounding",
                "validation_passed": validation["passed"]}
    make_figures(rows, output / "figures")
    metadata["sine_reference"] = "Closed-form exact interpolant integrals, Decimal precision 65"
    metadata["code_sha256"] = hashes(Path(__file__).resolve().parents[1], RUN_SOURCES)
    metadata["artifact_sha256"] = hashes(output, RUN_OUTPUTS)
    write_json(output / "metadata.json", metadata)
    print(f"case       N      E_inf         p_inf      E_L2          p_L2       E_H1sem      p_H1sem")
    for row in rows:
        rates = ["-" if row[k] is None else f"{row[k]:.6f}" for k in ("p_inf", "p_L2", "p_H1sem")]
        print(f"{row['case']:10} {row['N']:4d} {row['E_inf']:.8e} {rates[0]:>9} "
              f"{row['E_L2']:.8e} {rates[1]:>9} {row['E_H1sem']:.8e} {rates[2]:>9}")
    print(f"Checks: {sum(c['passed'] for c in checks)}/{len(checks)} passed. Output: {output}")
    return validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=["both", *CASES], default="both")
    parser.add_argument("--meshes", type=int, nargs="+", default=[10, 20, 40, 80, 160])
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "runs" / "primary")
    args = parser.parse_args()
    if any(n < 2 for n in args.meshes) or args.meshes != sorted(set(args.meshes)):
        parser.error("--meshes must contain unique, increasing integers >= 2")
    names = list(CASES) if args.case == "both" else [args.case]
    result = run(names, args.meshes, args.output_dir.resolve())
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
