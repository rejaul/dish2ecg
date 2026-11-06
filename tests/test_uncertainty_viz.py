"""Confidence bands and Bland-Altman helpers in dish2ecg.viz."""
import numpy as np
import pandas as pd
import pytest

from dish2ecg import viz

DATA = __import__("pathlib").Path(__file__).resolve().parents[1] / "results"


def test_band_is_narrowest_at_the_centring_point():
    x = np.linspace(-2, 2, 101)
    _, lo, hi = viz.confidence_band(x, slope=10, intercept=50, slope_se=2,
                                    intercept_se=3, centre=0.0, df=20)
    w = hi - lo
    assert np.argmin(w) == 50                       # the centre of the grid
    assert w[0] > w[50] and w[-1] > w[50]


def test_band_at_the_centre_reproduces_the_reported_interval():
    """The cascade reports effect +- t(0.95, n_groups-1) * intercept_se; the band must match."""
    from scipy import stats
    df = 14
    _, lo, hi = viz.confidence_band([0.0], slope=100.0, intercept=128.7, slope_se=10.2,
                                    intercept_se=21.99, centre=0.0, df=df)
    t = stats.t.ppf(0.95, df)
    assert lo[0] == pytest.approx(128.7 - t * 21.99)
    assert hi[0] == pytest.approx(128.7 + t * 21.99)


def test_band_matches_the_published_invitro_interval():
    iv = pd.read_csv(DATA / "invitro_concentration_response.csv")
    r = iv[(iv.drug == "Dofetilide") & (iv.analysis == "all_concentrations") &
           (iv.endpoint == "ddFPDc")].iloc[0]
    _, lo, hi = viz.confidence_band([r.at_concentration], slope=r.slope, intercept=r.intercept,
                                    slope_se=r.slope_se, intercept_se=r.intercept_se,
                                    centre=r.at_concentration, df=int(r.n_groups) - 1)
    assert lo[0] == pytest.approx(r.effect_ci90_lower, abs=0.05)
    assert hi[0] == pytest.approx(r.effect_ci90_upper, abs=0.05)


def test_clustering_widens_the_limits_of_agreement():
    """Ignoring repeated measures understates the spread; that is the whole point."""
    rng = np.random.default_rng(0)
    g = np.repeat(np.arange(12), 8)
    offs = rng.normal(0, 14, 12)[g]
    a = rng.normal(100, 3, len(g))
    b = a - (offs + rng.normal(0, 3, len(g)))
    naive = viz.bland_altman(a, b)
    clustered = viz.bland_altman(a, b, group=g)
    assert clustered["sd"] > naive["sd"] * 0.9
    assert clustered["var_between"] > 0


def test_variance_components_recover_known_variances():
    rng = np.random.default_rng(1)
    g = np.repeat(np.arange(60), 20)
    d = rng.normal(0, 5, 60)[g] + rng.normal(0, 2, len(g))
    vb, vw = viz.variance_components(d, g)
    assert vb == pytest.approx(25, rel=0.45)
    assert vw == pytest.approx(4, rel=0.2)


def test_bland_altman_reproduces_the_shipped_optical_comparison():
    arm = pd.read_csv(DATA / "arm_bycell.csv").dropna(subset=["MEA", "VSO"])
    stats_row = pd.read_csv(DATA / "bland_altman_stats.csv").iloc[0]
    s = viz.bland_altman(arm.VSO, arm.MEA, group=arm.drug)
    assert s["n"] == int(stats_row.n)
    assert s["bias"] == pytest.approx(stats_row.bias_ms, abs=0.01)
    assert s["sd"] == pytest.approx(stats_row.sd_diff_ms, abs=0.05)
    assert s["prop_bias_slope"] == pytest.approx(stats_row.prop_bias_slope, abs=0.01)


def test_optical_and_electrode_are_not_interchangeable():
    """The source study calls them equivalent measures; the limits of agreement do not."""
    s = pd.read_csv(DATA / "bland_altman_stats.csv").iloc[0]
    assert s.loa_upper_ms - s.loa_lower_ms > 400        # against drug effects of 10-130 ms
    assert abs(s.prop_bias_slope) > 0.3                 # and the gap grows with the effect


def test_rate_correction_choice_can_flip_the_ich_e14_decision():
    """Verapamil: Fridericia says below 10 ms, Bazett says above. Documented, not hidden."""
    m = pd.read_csv(DATA / "methodology_sensitivity.csv")
    v = m[(m.drug == "Verapamil") & (m.question == "rate correction")].set_index("variant")
    assert not bool(v.loc["Fridericia (used)", "exceeds_10ms"])
    assert bool(v.loc["Bazett", "exceeds_10ms"])


def test_preformed_double_delta_is_close_to_the_white_paper_model():
    """Our simplification should cost a few ms, not a conclusion."""
    m = pd.read_csv(DATA / "methodology_sensitivity.csv")
    p = m[m.question == "placebo adjustment"]
    wide = p.pivot_table(index="drug", columns="variant", values="effect_ms")
    gap = (wide["pre-formed double delta (used)"] - wide["white-paper model (estimated)"]).abs()
    assert gap.max() < 5.0
