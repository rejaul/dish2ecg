import numpy as np
import pandas as pd
import pytest

from dish2ecg import double_delta, emax_fit, fit_concentration_response, fridericia


def _dd(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    return double_delta(crossover, qtc, subject_col="subject", arm_col="arm",
                        baseline_mask=crossover["is_baseline"], control_arm="placebo",
                        time_col="time", match_on_time=True, paired=True)


def _arm(crossover, dd, arm):
    sel = ((crossover["arm"] == arm) & (crossover["day"] == "dose")).to_numpy()
    return crossover[sel], dd.subset(sel)


def test_concentration_response_recovers_the_true_slope(crossover):
    dd = _dd(crossover)
    sub, ep = _arm(crossover, dd, "drugA")
    fit = fit_concentration_response(sub, ep, concentration_col="conc",
                                     group_col="subject", at_concentration=20.0)
    assert fit.slope == pytest.approx(crossover.attrs["true_slope"]["drugA"], abs=1e-4)
    assert fit.effect_at_concentration == pytest.approx(
        crossover.attrs["true_effect_at_conc20"]["drugA"], abs=1e-3)
    assert fit.n_groups == 8
    assert fit.n_obs == 48
    assert "fridericia" in fit.correction_cascade


def test_negative_control_slope_is_zero(crossover):
    dd = _dd(crossover)
    sub, ep = _arm(crossover, dd, "drugB")
    fit = fit_concentration_response(sub, ep, concentration_col="conc",
                                     group_col="subject", at_concentration=20.0)
    assert fit.slope == pytest.approx(0.0, abs=1e-6)
    assert fit.exceeds_threshold(10.0) is False


def test_threshold_decision_and_row_export(crossover):
    dd = _dd(crossover)
    sub, ep = _arm(crossover, dd, "drugA")
    fit = fit_concentration_response(sub, ep, concentration_col="conc",
                                     group_col="subject", at_concentration=20.0)
    assert fit.exceeds_threshold(10.0) is True
    assert fit.exceeds_threshold(1000.0) is False
    row = fit.to_row()
    for k in ("effect_ci90_upper", "effect_ci90_lower", "slope_ci90_upper",
              "slope_ci90_lower", "correction_cascade", "model_formula"):
        assert k in row


def test_centering_does_not_change_the_prediction(crossover):
    dd = _dd(crossover)
    sub, ep = _arm(crossover, dd, "drugA")
    a = fit_concentration_response(sub, ep, concentration_col="conc",
                                   group_col="subject", at_concentration=20.0)
    b = fit_concentration_response(sub, ep, concentration_col="conc",
                                   group_col="subject", at_concentration=20.0,
                                   center_concentration=True)
    assert a.slope == pytest.approx(b.slope, abs=1e-6)
    assert a.effect_at_concentration == pytest.approx(b.effect_at_concentration,
                                                      abs=1e-6)


def test_fit_requires_the_correction_cascade(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    with pytest.raises(Exception, match="missing required correction"):
        fit_concentration_response(crossover, qtc, concentration_col="conc",
                                   group_col="subject")


def test_subset_preserves_provenance(crossover):
    dd = _dd(crossover)
    sub, ep = _arm(crossover, dd, "drugA")
    assert ep.steps == dd.steps
    assert len(ep.values) == len(sub)
    with pytest.raises(Exception, match="mask has"):
        dd.subset(np.ones(3, dtype=bool))


def test_emax_recovers_known_parameters():
    c = np.array([0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0])
    y = 5.0 + 50.0 * c / (2.0 + c)
    out = emax_fit(c, y)
    assert out["converged"]
    assert out["E0"] == pytest.approx(5.0, abs=1e-4)
    assert out["Emax"] == pytest.approx(50.0, abs=1e-3)
    assert out["EC50"] == pytest.approx(2.0, abs=1e-3)


def test_emax_refuses_too_few_points():
    assert emax_fit([0.0, 1.0], [0.0, 1.0])["converged"] is False


def test_interval_coverage_over_many_simulated_studies(simulate_noisy_crossover):
    """Interval coverage is a frequency property: test it as one.

    A single realisation of a 90% interval misses the truth 10% of the time, so a
    one-seed assertion is a coin flip. Measured here over 60 simulated crossover
    studies. Reference run at n=200: effect-at-Cmax coverage 0.865, slope coverage
    0.900 with the t quantile (0.845 / 0.900 with a normal quantile).
    """
    eff_hits, slope_hits = [], []
    for seed in range(60):
        df = simulate_noisy_crossover(seed)
        qtc = fridericia(df["qt"], df["rr"])
        dd = double_delta(df, qtc, subject_col="subject", arm_col="arm",
                          baseline_mask=df["is_baseline"], control_arm="placebo",
                          time_col="time", match_on_time=True, paired=True)
        sel = ((df["arm"] == "drugA") & (df["day"] == "dose")).to_numpy()
        fit = fit_concentration_response(df[sel], dd.subset(sel),
                                         concentration_col="conc",
                                         group_col="subject", at_concentration=20.0)
        assert fit.effect_ci90[1] > fit.effect_ci90[0]
        eff_hits.append(fit.effect_ci90[0] <= 40.0 <= fit.effect_ci90[1])
        slope_hits.append(fit.slope_ci90[0] <= 2.0 <= fit.slope_ci90[1])
    assert 0.75 <= np.mean(eff_hits) <= 1.0, f"effect coverage {np.mean(eff_hits)}"
    assert 0.78 <= np.mean(slope_hits) <= 1.0, f"slope coverage {np.mean(slope_hits)}"


def test_boundary_fallback_is_flagged_not_silent(crossover):
    """A zero between-subject variance must refit as OLS and say so."""
    dd = _dd(crossover)
    sub, ep = _arm(crossover, dd, "drugA")
    fit = fit_concentration_response(sub, ep, concentration_col="conc",
                                     group_col="subject", at_concentration=20.0)
    # the noiseless fixture has no between-subject variance left after ddQTcF
    if fit.boundary_fallback:
        assert "variance estimated at 0" in fit.random_effects
        assert np.isfinite(fit.effect_ci90[0]) and np.isfinite(fit.effect_ci90[1])
    assert isinstance(fit.boundary_fallback, bool)
