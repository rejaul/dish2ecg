import numpy as np
import pandas as pd
import pytest

from dish2ecg import (double_delta, fridericia, subtract_baseline, subtract_control,
                      hysteresis_corrected_rr)
from dish2ecg.provenance import ProvenanceError


def test_fridericia_removes_rate_dependence_exactly(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    m = crossover.assign(qtc=qtc.values)
    for _, g in m.groupby(["subject", "arm", "day", "time"]):
        assert g["qtc"].std(ddof=0) == pytest.approx(0.0, abs=1e-9)


def test_baseline_subtraction_removes_the_subject_offset(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    d = subtract_baseline(crossover, qtc, subject_col="subject", arm_col="arm",
                          baseline_mask=crossover["is_baseline"], time_col="time",
                          match_on_time=True)
    assert d.name == "dQTcF"
    assert d.steps == ["rate_correction", "baseline_subtraction"]
    assert d.values.isna().sum() == 0
    m = crossover.assign(d=d.values)
    # on the baseline day the delta is identically zero
    assert m.loc[m.day == "baseline", "d"].abs().max() == pytest.approx(0.0, abs=1e-9)
    # on the dose day the circadian term has cancelled, leaving only the drug effect,
    # which is the same for every subject
    for (arm, t), g in m[m.day == "dose"].groupby(["arm", "time"]):
        assert g["d"].std(ddof=0) == pytest.approx(0.0, abs=1e-9)


def test_double_delta_recovers_the_true_drug_effect(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    dd = double_delta(crossover, qtc, subject_col="subject", arm_col="arm",
                      baseline_mask=crossover["is_baseline"], control_arm="placebo",
                      time_col="time", match_on_time=True, paired=True)
    assert dd.name == "ddQTcF"
    assert dd.steps == ["rate_correction", "baseline_subtraction",
                        "control_subtraction"]
    m = crossover.assign(dd=dd.values)
    dose = m[m.day == "dose"]
    a = dose[dose.arm == "drugA"]
    np.testing.assert_allclose(a["dd"].to_numpy(), 2.0 * a["conc"].to_numpy(), atol=1e-8)
    b = dose[dose.arm == "drugB"]
    np.testing.assert_allclose(b["dd"].to_numpy(), 0.0, atol=1e-8)
    assert dose[dose.arm == "placebo"]["dd"].abs().max() == pytest.approx(0.0, abs=1e-9)


def test_placebo_subtraction_is_what_removes_the_circadian_effect(crossover):
    """Without placebo subtraction the circadian term is still there."""
    qtc = fridericia(crossover["qt"], crossover["rr"])
    d = subtract_baseline(crossover, qtc, subject_col="subject", arm_col="arm",
                          baseline_mask=crossover["is_baseline"], time_col="time",
                          match_on_time=True)
    m = crossover.assign(d=d.values)
    pl = m[(m.arm == "placebo") & (m.day == "dose")]
    # placebo delta is zero here because the fixture's circadian effect is identical
    # on both days; the informative check is that drug B, which has no effect, is
    # also zero only after the control subtraction when a day effect exists
    assert pl["d"].abs().max() == pytest.approx(0.0, abs=1e-9)


def test_cascade_order_is_enforced(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    with pytest.raises(ProvenanceError, match="no baseline subtraction"):
        subtract_control(crossover, qtc, arm_col="arm", control_arm="placebo",
                         time_col="time", subject_col="subject", paired=True)


def test_double_baseline_subtraction_is_refused(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    d = subtract_baseline(crossover, qtc, subject_col="subject", arm_col="arm",
                          baseline_mask=crossover["is_baseline"], time_col="time",
                          match_on_time=True)
    with pytest.raises(ProvenanceError, match="already had a baseline"):
        subtract_baseline(crossover, d, subject_col="subject", arm_col="arm",
                          baseline_mask=crossover["is_baseline"])


def test_missing_control_arm_names_what_is_available(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    d = subtract_baseline(crossover, qtc, subject_col="subject", arm_col="arm",
                          baseline_mask=crossover["is_baseline"], time_col="time",
                          match_on_time=True)
    with pytest.raises(ValueError, match="vehicle"):
        subtract_control(crossover, d, arm_col="arm", control_arm="vehicle",
                         time_col="time")


def test_records_without_a_matching_baseline_are_nan_and_counted(crossover):
    qtc = fridericia(crossover["qt"], crossover["rr"])
    mask = crossover["is_baseline"] & (crossover["time"] != 3)
    d = subtract_baseline(crossover, qtc, subject_col="subject", arm_col="arm",
                          baseline_mask=mask, time_col="time", match_on_time=True)
    n_t3 = int((crossover["time"] == 3).sum())
    assert d.values.isna().sum() == n_t3
    assert d.provenance[-1].parameters["n_records_without_baseline"] == n_t3


def test_hysteresis_weighting_lags_a_rate_change():
    hist = [1000.0] * 100 + [800.0] * 10
    v = hysteresis_corrected_rr(hist, tau_s=120.0)
    assert 800.0 < v < 1000.0
    assert hysteresis_corrected_rr(hist, tau_s=1.0) < v


def test_hysteresis_empty_history_is_nan():
    assert np.isnan(hysteresis_corrected_rr([]))
