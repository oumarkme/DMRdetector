"""R-compatible statistics used by the v2.2 algorithm.

v2.2 was written in R. These helpers reproduce the R semantics that differ from
numpy/scipy defaults, so Python results match v2.2.
"""

from __future__ import annotations

from collections.abc import Iterable
import decimal
from fractions import Fraction
import math

import numpy as np
from scipy import stats


_EPS = np.finfo(float).eps  # R's .Machine$double.eps
# R's M_LOG10_2; math.log10 is correctly rounded for this argument, so the
# constant does not depend on the libm of the host
_LOG10_2 = math.log10(2.0)


def exp_rounded(x: float) -> float:
    """Returns exp(x) correctly rounded to float64, on every CPU and libm.

    np.exp and the libm exp are accurate to about 1 ulp, not exact, and the
    value for the same argument depends on the SIMD loop that numpy selects at
    run time and on the libm variant (FMA or not) that glibc selects. This
    function evaluates exp with `decimal` (50 digits, correctly rounded
    results, integer arithmetic only) and converts the decimal to the nearest
    float64. The result is recomputed with 80 digits; a double-rounding hazard
    (a 50-digit result that lies within 10^-50 relative of a midpoint between
    two doubles) would show as a mismatch and raises, which has probability
    below 1e-34 per call. It costs about 20 microseconds per call, so use it
    for tables that are computed once.

    Args:
      x: A finite or infinite float; nan gives nan.

    Returns:
      exp(x) rounded to nearest, ties to even; 0.0 on underflow and inf on
      overflow.
    """
    if math.isnan(x):
        return math.nan
    if math.isinf(x):
        return 0.0 if x < 0 else math.inf
    if x == 0.0:
        return 1.0
    val = []
    for prec in (50, 80):
        ctx = decimal.Context(prec=prec, Emin=-999999, Emax=999999)
        val.append(float(ctx.exp(decimal.Decimal(x))))
    if val[0] != val[1]:
        raise RuntimeError(f"exp({x!r}) is not stable between 50 and 80 "
                           f"decimal digits: {val[0]!r} and {val[1]!r}")
    return val[0]


def t_test_one_sample_p(x: Iterable[float] | np.ndarray,
                        mu: float = 0.0) -> float:
    """Two-sided one-sample t-test p-value, as R's t.test(x, mu = mu)$p.value.

    v2.2 wraps t.test in tryCatch and maps errors to NA. R stops when fewer than
    two observations are given ("not enough 'x' observations") or when stderr <
    10 * .Machine$double.eps * abs(mean) ("data are essentially constant"). Both
    return NaN here.
    """
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    n = x.size
    if n < 2:
        return float("nan")
    mx = float(x.mean())
    stderr = float(np.sqrt(np.var(x, ddof=1) / n))
    if stderr < 10 * _EPS * abs(mx):
        return float("nan")
    t = (mx - mu) / stderr
    if np.isnan(t):
        return float("nan")
    return float(2.0 * stats.t.sf(abs(t), n - 1))


def p_adjust_bh(p: Iterable[float] | np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjustment, as R's p.adjust(p, method = "BH").

    NaN values stay NaN and are not counted in n (R drops NA before adjusting
    and uses n = number of non-NA values).
    """
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    # decreasing, as R order(decreasing=TRUE)
    order = np.argsort(-pv, kind="stable")
    ranks = np.arange(n, 0, -1)
    adj = np.minimum.accumulate((n / ranks) * pv[order])
    adj = np.minimum(adj, 1.0)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def nextn(n: int, factors: tuple[int, ...] = (2, 3, 5)) -> int:
    """Returns the smallest integer >= n whose prime factors are in `factors`.

    This is R's nextn.
    """
    m = max(int(n), 1)
    while True:
        k = m
        for f in factors:
            while k % f == 0:
                k //= f
        if k == 1:
            return m
        m += 1


def r_round(x: float | np.ndarray, digits: int = 0) -> np.ndarray:
    """Rounds as R >= 4.0.0 does.

    round() picks the closer of the two decimal candidates, ties to even.
    """
    x = np.asarray(x, dtype=float)
    out = np.array(x, copy=True)
    if digits == 0:
        return np.rint(x)
    ok = np.isfinite(x) & (x != 0)
    # a Fraction is rounded once, so no libm pow is involved
    p10 = float(Fraction(10)**digits)
    ax = np.abs(x[ok])
    x10 = p10 * ax
    i10 = np.floor(x10)
    xd = i10 / p10
    xu = np.ceil(x10) / p10
    dd, du = ax - xd, xu - ax
    pick_d = (dd < du) | ((np.fmod(i10, 2.0) == 0) & (dd == du))
    res = np.where(pick_d, xd, xu)
    # R returns x unchanged when digits exceed double precision
    # R: digits + M_LOG10_2 * (0.5 + logb(x)) > DBL_DIG. frexp gives logb
    # exactly
    # (np.log2 can round up to the next integer just below a power of 2, and
    # numpy may evaluate it with an AVX512 loop)
    logb = np.frexp(ax)[1] - 1
    big = (digits + (logb + 0.5) * _LOG10_2) > 15
    res = np.where(big, ax, res)
    out[ok] = np.sign(x[ok]) * res
    return out
