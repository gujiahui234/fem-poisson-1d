"""Independent small checks for the research quantities and assembly."""
import math
import numpy as np
import pytest
from run_experiments import CASES, assemble_and_solve, maximum_error, integral_errors, dense_maximum
from reference_errors import sine_integral_reference


def test_two_element_hand_calculation():
    sol = assemble_and_solve(CASES["quadratic"], 2)
    np.testing.assert_allclose(sol.A.toarray(), [[4]], atol=1e-14, rtol=0)
    np.testing.assert_allclose(sol.b, [0.5], atol=1e-14, rtol=0)
    np.testing.assert_allclose(sol.U, [0, 0.125, 0], atol=1e-14, rtol=0)
    assert abs(maximum_error(sol)[0] - 1 / 32) < 1e-14
    assert max(abs(sol.U - sol.case.u(sol.x))) < 1e-14


@pytest.mark.parametrize("n", [3, 10, 40, 160])
def test_quadratic_closed_form_norms(n):
    sol = assemble_and_solve(CASES["quadratic"], n)
    measured = [maximum_error(sol)[0], *integral_errors(sol)]
    h = 1 / n
    expected = [h**2 / 8, h**2 / math.sqrt(120), h / math.sqrt(12)]
    np.testing.assert_allclose(measured, expected, rtol=2e-8, atol=1e-14)


def test_three_interior_nodes_match_notes_example_2_5():
    sol = assemble_and_solve(CASES["quadratic"], 4)
    np.testing.assert_allclose(sol.A.toarray(), [[8, -4, 0], [-4, 8, -4], [0, -4, 8]])
    np.testing.assert_allclose(sol.b, [0.25] * 3, atol=1e-14)
    np.testing.assert_allclose(sol.U[1:-1], [3 / 32, 1 / 8, 3 / 32], atol=1e-14)
    assert np.linalg.eigvalsh(sol.A.toarray()).min() > 0


def test_sine_extremum_search_against_independent_analytic_stationary_points():
    sol = assemble_and_solve(CASES["sine"], 10)
    values = []
    for k in range(sol.n):
        a, b = sol.x[k:k + 2]
        slope = (sol.U[k + 1] - sol.U[k]) / (b - a)
        z = math.acos(slope / math.pi) / math.pi
        assert a < z < b
        values.append(abs(math.sin(math.pi * z) - (sol.U[k] + slope * (z - a))))
    np.testing.assert_allclose(maximum_error(sol)[0], max(values), rtol=1e-12)
    assert abs(dense_maximum(sol) / maximum_error(sol)[0] - 1) < 2e-6


def test_quadrature_refinement_and_residual():
    sol5 = assemble_and_solve(CASES["sine"], 40, 5)
    sol10 = assemble_and_solve(CASES["sine"], 40, 10)
    np.testing.assert_allclose(sol5.U, sol10.U, rtol=0, atol=1e-12)
    np.testing.assert_allclose(integral_errors(sol5, 5), integral_errors(sol5, 10), rtol=1e-8)
    assert np.max(np.abs(sol5.A @ sol5.U[1:-1] - sol5.b)) < 1e-12


@pytest.mark.parametrize("n", [10, 20, 40, 80, 160])
def test_sine_integrals_against_closed_forms(n):
    sol = assemble_and_solve(CASES["sine"], n, 5)
    np.testing.assert_allclose(integral_errors(sol, 5), sine_integral_reference(n),
                               rtol=2e-8, atol=1e-14)


def test_underintegrated_load_does_not_have_exact_nodal_values():
    coarse = assemble_and_solve(CASES["sine"], 10, 1)
    accurate = assemble_and_solve(CASES["sine"], 10, 5)
    assert np.max(np.abs(coarse.U - coarse.case.u(coarse.x))) > 1e-4
    assert np.max(np.abs(accurate.U - accurate.case.u(accurate.x))) < 1e-12


@pytest.mark.parametrize("n", [0, 1, -4, 2.5])
def test_invalid_element_count(n):
    with pytest.raises(ValueError):
        assemble_and_solve(CASES["quadratic"], n)
