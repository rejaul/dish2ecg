"""Rate corrections against closed-form ground truth."""
import numpy as np
import pandas as pd
import pytest

from dish2ecg import (bazett, fridericia, framingham, hodges, jtpeakc_hnatkova2019,
                      power_correction, rate_correct, fit_power_exponent,
                      fit_individual_exponents, individual_correction)
from dish2ecg.provenance import ProvenanceError


def test_bazett_matches_formula():
    qt, rr = np.array([400.0, 380.0]), np.array([1000.0, 800.0])
    got = bazett(qt, rr).values.to_numpy()
    want = qt / np.sqrt(rr / 1000.0)
    np.testing.assert_allclose(got, want)
    # at RR = 1000 ms exactly, every power correction is the identity
    assert got[0] == pytest.approx(400.0)


def test_fridericia_matches_formula():
    qt, rr = np.array([400.0, 420.0]), np.array([1000.0, 1200.0])
    got = fridericia(qt, rr).values.to_numpy()
    np.testing.assert_allclose(got, qt / np.cbrt(rr / 1000.0))


def test_bazett_overcorrects_relative_to_fridericia_at_high_rate():
    # RR = 600 ms is 100 bpm; Bazett must give the larger corrected value
    qt, rr = 350.0, 600.0
    assert bazett([qt], [rr]).values[0] > fridericia([qt], [rr]).values[0]


def test_framingham_and_hodges_are_identity_at_60_bpm():
    assert framingham([400.0], [1000.0]).values[0] == pytest.approx(400.0)
    assert hodges([400.0], [1000.0]).values[0] == pytest.approx(400.0)


def test_jtpeakc_hnatkova2019_exact_constant():
    # JTpc = JTp + 0.150 * (1 - RR) in seconds; RR = 0.8 s -> +30 ms
    got = jtpeakc_hnatkova2019([200.0], [800.0]).values[0]
    assert got == pytest.approx(230.0)
    assert jtpeakc_hnatkova2019([200.0], [1000.0]).values[0] == pytest.approx(200.0)


def test_power_correction_recovers_known_exponent():
    rng = np.random.default_rng(0)
    rr = rng.uniform(700, 1300, 4000)
    alpha_true = 0.37
    qt = 400.0 * (rr / 1000.0) ** alpha_true
    alpha, diag = fit_power_exponent(qt, rr)
    assert alpha == pytest.approx(alpha_true, abs=1e-6)
    assert diag["r_squared"] == pytest.approx(1.0, abs=1e-9)
    # applying it removes the rate dependence entirely
    c = power_correction(qt, rr, alpha).values.to_numpy()
    assert c.std() == pytest.approx(0.0, abs=1e-8)


def test_fit_power_exponent_with_noise():
    rng = np.random.default_rng(1)
    rr = rng.uniform(700, 1300, 20000)
    qt = 400.0 * (rr / 1000.0) ** 0.33 * np.exp(rng.normal(0, 0.01, rr.size))
    alpha, diag = fit_power_exponent(qt, rr)
    assert alpha == pytest.approx(0.33, abs=0.01)
    assert diag["n"] == 20000


def test_individual_correction_recovers_per_subject_exponents():
    rng = np.random.default_rng(2)
    rows, truth = [], {}
    for s in range(12):
        a = 0.25 + 0.03 * s
        truth[s] = a
        rr = rng.uniform(700, 1300, 200)
        rows.append(pd.DataFrame({"subj": s, "rr": rr, "qt": 400.0 * (rr / 1000.0) ** a}))
    df = pd.concat(rows, ignore_index=True)
    ex = fit_individual_exponents(df, interval_col="qt", cycle_length_col="rr",
                                  subject_col="subj")
    for _, r in ex.iterrows():
        assert r["alpha"] == pytest.approx(truth[r["subj"]], abs=1e-6)
    ep = individual_correction(df, interval_col="qt", cycle_length_col="rr",
                               subject_col="subj", exponents=ex)
    np.testing.assert_allclose(ep.values.to_numpy(), 400.0, atol=1e-8)


def test_individual_correction_requires_fallback_for_unfitted_subject():
    df = pd.DataFrame({"subj": [0, 1], "qt": [400.0, 400.0], "rr": [900.0, 900.0]})
    ex = pd.DataFrame({"subj": [0, 1], "alpha": [0.33, np.nan]})
    with pytest.raises(ValueError, match="fallback_alpha"):
        individual_correction(df, interval_col="qt", cycle_length_col="rr",
                              subject_col="subj", exponents=ex)
    ep = individual_correction(df, interval_col="qt", cycle_length_col="rr",
                               subject_col="subj", exponents=ex, fallback_alpha=0.33)
    assert ep.provenance[-1].parameters["n_records_on_fallback"] == 1


def test_invitro_path_uses_the_same_algebra():
    """FPD against the interbeat interval is the same operation as QT against RR."""
    fpd, ibi = np.array([350.0]), np.array([1500.0])
    a = power_correction(fpd, ibi, 1 / 3, name="FPDcF", domain="invitro")
    b = fridericia(fpd, ibi, name="QTcF")
    assert a.values[0] == pytest.approx(b.values[0])
    assert a.domain == "invitro" and b.domain == "clinical"


def test_rate_correct_dispatch_and_unknown_method():
    assert rate_correct([400.0], [900.0], "fridericia").name == "QTcF"
    with pytest.raises(ValueError, match="unknown correction"):
        rate_correct([400.0], [900.0], "nonsense")
    with pytest.raises(ValueError, match="requires alpha"):
        rate_correct([400.0], [900.0], "power")


def test_zero_cycle_length_is_nan_not_inf():
    v = power_correction([400.0, 400.0], [900.0, 0.0], 0.5).values
    assert np.isfinite(v[0]) and np.isnan(v[1])


def test_length_mismatch_raises():
    with pytest.raises(ValueError, match="length mismatch"):
        power_correction([1.0, 2.0], [1000.0], 0.5)


def test_provenance_is_recorded():
    ep = fridericia([400.0], [900.0])
    assert ep.steps == ["rate_correction"]
    assert ep.provenance[0].method == "fridericia"
    assert "ICH E14" in ep.provenance[0].citation
    assert "QTcF" in ep.describe()
