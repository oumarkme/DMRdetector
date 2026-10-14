import numpy as np

from dmrdetector.util import rstats


def test_t_test_constant_and_short():
    assert np.isnan(rstats.t_test_one_sample_p([0.5, 0.5, 0.5]))
    assert np.isnan(rstats.t_test_one_sample_p([0.5]))
    assert rstats.t_test_one_sample_p([0.1, 0.2, 0.3, 0.4]) < 0.05


def test_p_adjust_bh_matches_r():
    # R: p.adjust(c(0.01, 0.02, 0.03, NA, 0.5), "BH") -> 0.04 0.04 0.04 NA 0.5
    got = rstats.p_adjust_bh([0.01, 0.02, 0.03, np.nan, 0.5])
    np.testing.assert_allclose(got[[0, 1, 2, 4]], [0.04, 0.04, 0.04, 0.5])
    assert np.isnan(got[3])


def test_nextn():
    assert [rstats.nextn(n) for n in (1, 7, 13, 1021)] == [1, 8, 15, 1024]


def test_exp_rounded_special_values_and_known_roundings():
    assert rstats.exp_rounded(0.0) == 1.0 and rstats.exp_rounded(-0.0) == 1.0
    assert np.isnan(rstats.exp_rounded(float("nan")))
    assert rstats.exp_rounded(float("-inf")) == 0.0
    assert rstats.exp_rounded(float("inf")) == float("inf")
    assert rstats.exp_rounded(710.0) == float("inf")
    assert rstats.exp_rounded(-800.0) == 0.0
    # R 4.4 and glibc: exp(1) and exp(-1), as printed by sprintf("%a")
    assert rstats.exp_rounded(1.0).hex() == float.fromhex(
        "0x1.5bf0a8b145769p+1").hex()
    assert rstats.exp_rounded(-1.0).hex() == float.fromhex(
        "0x1.78b56362cef38p-2").hex()
    # smallest subnormal: exp(-744.4400719213812) = 4.94e-324 rounds to 5e-324
    assert rstats.exp_rounded(-744.4400719213812) == 5e-324


def test_gaussian_kernel_equals_r_bit_for_bit():
    from dmrdetector.steps import atlas
    # R 4.4: k <- exp(-(x^2) / (2 * 150^2)); sprintf("%a", k[x == i])
    want = {
        0: "0x1p+0",
        1: "0x1.fffd165aaf821p-1",
        7: "0x1.ff715ad6456d8p-1",
        100: "0x1.99fa40bc6c5f7p-1",
        333: "0x1.5c79922e1c972p-4",
        450: "0x1.6c0504695c417p-7",
    }
    k = atlas.gaussian_kernel(150)
    assert k.size == 901 and k.dtype == np.float64
    for i, h in want.items():
        assert k[450 + i].hex() == float.fromhex(h).hex(), i
        assert k[450 - i].hex() == float.fromhex(h).hex(), -i


def test_gaussian_kernel_does_not_call_numpy_exp_or_libm(monkeypatch):
    import math
    from dmrdetector.steps import atlas

    def refuse(*args, **kwargs):
        raise AssertionError("the kernel must not use np.exp or math.exp")

    monkeypatch.setattr(np, "exp", refuse)
    monkeypatch.setattr(math, "exp", refuse)
    k = atlas.gaussian_kernel(23)  # a sigma that no other test caches
    assert k.size == 139 and k[69] == 1.0 and (k == k[::-1]).all()
    assert 0 < k[0] < k[1] < 1


def test_gaussian_kernel_returns_a_private_array():
    from dmrdetector.steps import atlas
    k = atlas.gaussian_kernel(31)
    k[:] = 0.0
    assert atlas.gaussian_kernel(31).max() == 1.0


def test_r_round_large_values_use_logb_as_r_does():
    # digits + log10(2) * (0.5 + logb(x)) > 15: round() returns x unchanged
    assert rstats.r_round(np.array([3.0e13 + 0.04]), 1)[0] == 3.0e13
    big = 1.5 * 2.0**47 + 0.25
    assert rstats.r_round(np.array([big]), 1)[0] == big
    # just below 2^47: np.log2 rounds up to exactly 47.0 and would say "too
    # large to round"; logb is 46 (R: 14.998 <= 15), so the value is rounded
    below = np.nextafter(2.0**47, 0.0)
    assert rstats.r_round(np.array([below]), 1)[0] == 2.0**47
