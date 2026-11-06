"""The two waveform results, pinned against the committed tables."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

RES = Path(__file__).resolve().parents[1] / "results"
pytest.importorskip("scipy")
from scipy import stats  # noqa: E402

WINDOW_S, QT_ADAPTATION_S = 10.0, (60.0, 120.0)


def _rr():
    return pd.read_csv(RES / "waveform_rr_window.csv")


def _dl():
    d = pd.read_csv(RES / "waveform_delineation.csv")
    return d[d.qt_tangent.between(250, 650)]


def test_r_detector_agrees_with_the_published_rr():
    """If this fails the beat detection is wrong and nothing downstream is safe."""
    a = _rr()
    err = a.rr_mean - a.table_RR
    assert abs(err.median()) < 2.0
    assert (err.abs() > 20).mean() < 0.05


def test_the_window_is_far_shorter_than_qt_adaptation():
    """Why hysteresis correction is unidentifiable here, not merely unimplemented."""
    assert WINDOW_S / QT_ADAPTATION_S[0] < 0.2


def test_within_window_rr_variation_is_not_negligible():
    a = _rr()
    assert a.rr_sd.median() > 10.0
    assert a.rr_range.median() > 40.0


def test_implied_qtcf_spread_exceeds_the_threshold_for_most_records():
    a = _rr()
    assert (a.qtcf_spread > 10).mean() > 0.5


def test_tangent_tend_differs_from_the_study_rule_by_a_few_ms():
    d = _dl()
    assert -12.0 < d.dQT.median() < 0.0


def test_the_delineation_difference_is_rate_dependent():
    """The reason an extraction-method difference does not cancel in the cascade."""
    d = _dl()
    slope, _, _, p, _ = stats.linregress(d.table_RR, d.dQT)
    assert slope < 0
    assert p < 0.05
    assert abs(slope * 300) > 3.0      # ms across a 300 ms RR change
