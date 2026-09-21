"""Independent integral-error references for the uniform-mesh sine interpolant.

These formulas use exact element integrals and trigonometric sums, not Gauss
quadrature or the assembled solver. They describe the exact Galerkin solution
only under the nodal-exactness hypotheses stated in the report.
"""
from decimal import Decimal, localcontext


def sine_integral_reference(n: int) -> tuple[float, float]:
    """Return (L2 error, energy seminorm) with stable high-precision evaluation.

    theta = pi/N; S = 2*(1-cos(theta))/theta**2.
    E0**2 = (5+cos(theta))/6 - S; EE**2 = pi**2*(1-S)/2.
    Direct float64 subtraction is inaccurate for the small L2 squared error.
    Decimal and the convergent cosine Taylor series avoid that cancellation.
    """
    if not isinstance(n, int) or isinstance(n, bool) or n < 2:
        raise ValueError("n must be an integer >= 2")
    with localcontext() as ctx:
        ctx.prec = 65
        pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406")
        theta = pi / n
        term = cosine = Decimal(1)
        for k in range(1, 100):
            term *= -theta * theta / ((2*k - 1) * (2*k))
            cosine += term
            if abs(term) < Decimal("1e-67"):
                break
        else:
            raise ArithmeticError("cosine reference series did not converge")
        s = 2 * (1 - cosine) / theta**2
        e0_sq = (5 + cosine) / 6 - s
        ee_sq = pi**2 * (1 - s) / 2
        return float(e0_sq.sqrt()), float(ee_sq.sqrt())
