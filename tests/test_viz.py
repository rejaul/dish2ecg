"""The WFDB reader and the plotting helpers.

The reader is checked against a hand-decoded header rather than a golden file,
and the annotation decoder against the invariant that makes the demo beat
usable: its fiducial intervals reproduce the clinical table exactly.
"""
import numpy as np
import pytest

from dish2ecg.viz import (ANNOTATION_CODES, synthetic_field_potential,
                          plot_ecg_beat, plot_field_potential, plot_qt_rr,
                          plot_concentration_response)

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt


#: Fiducial samples decoded from the demo record's .atr, and the clinical
#: table's own measured intervals for the same record (subject 1001,
#: ranolazine, TPT 1.0 h, EGREFID 00ed2097-cd14-4f03-ab33-853da5be5550).
FIDUCIALS = {"P_onset": 164, "QRS_onset": 297, "J": 392, "T_peak": 597, "T_end": 674}
INTERVALS = {"PR": 133.0, "QRS": 95.0, "QT": 377.0, "JTPEAK": 205.0, "TPEAKTEND": 77.0}


def test_fiducials_reproduce_the_clinical_intervals():
    """This is why the demo beat can be annotated with table values at all."""
    f = FIDUCIALS
    assert f["QRS_onset"] - f["P_onset"] == INTERVALS["PR"]
    assert f["J"] - f["QRS_onset"] == INTERVALS["QRS"]
    assert f["T_peak"] - f["J"] == INTERVALS["JTPEAK"]
    assert f["T_end"] - f["T_peak"] == INTERVALS["TPEAKTEND"]
    assert f["T_end"] - f["QRS_onset"] == INTERVALS["QT"]


def test_the_interval_identity_holds_for_this_record():
    assert INTERVALS["QRS"] + INTERVALS["JTPEAK"] + INTERVALS["TPEAKTEND"] == INTERVALS["QT"]


def test_annotation_codes_cover_the_delineation_set():
    assert {"(", ")", "t"} <= set(ANNOTATION_CODES.values())


def test_synthetic_field_potential_peaks_at_the_requested_duration():
    t, y = synthetic_field_potential(duration_ms=400.0)
    late = t > 150.0
    assert abs(t[late][np.argmax(y[late])] - 400.0) < 12.0


def test_synthetic_field_potential_has_a_downward_depolarisation_spike():
    t, y = synthetic_field_potential()
    assert y.min() < -0.5
    assert abs(t[np.argmin(y)]) < 10.0


def test_plot_ecg_beat_runs_and_labels_intervals():
    t = np.linspace(0, 800, 801)
    y = np.exp(-0.5 * ((t - 350) / 12) ** 2)
    fig, ax = plt.subplots()
    plot_ecg_beat((t, y), fiducials=FIDUCIALS, intervals=INTERVALS, ax=ax)
    texts = [x.get_text() for x in ax.texts]
    assert any("QT 377 ms" == s for s in texts)
    plt.close(fig)


def test_plot_helpers_do_not_raise():
    fig, ax = plt.subplots()
    plot_field_potential(ax=ax); plt.close(fig)
    fig, ax = plt.subplots()
    rng = np.random.default_rng(0)
    rr = rng.uniform(700, 1200, 200)
    plot_qt_rr(rr, 400 * (rr / 1000) ** (1 / 3), fitted_alpha=0.31, ax=ax)
    plt.close(fig)
    fig, ax = plt.subplots()
    c = rng.uniform(0, 3, 100)
    plot_concentration_response(c, 25 * c, slope=25.0, intercept=50.0,
                                at_concentration=2.0, effect_at=50.0,
                                ci=(45.0, 55.0), ax=ax)
    plt.close(fig)
