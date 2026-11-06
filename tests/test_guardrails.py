import numpy as np
import pandas as pd
import pytest

from dish2ecg import (Endpoint, fridericia, power_correction, provenance_columns,
                      to_free_cmax_multiples, to_percent_change, to_standardized_effect)
from dish2ecg.provenance import DomainError, ProvenanceError


def test_invitro_endpoint_refuses_to_be_a_qt_interval():
    fpdc = power_correction([350.0], [1500.0], 1 / 3, name="FPDcF", domain="invitro")
    with pytest.raises(DomainError, match="cannot be expressed as a QT interval"):
        fpdc.as_qt_milliseconds()


def test_clinical_endpoint_allows_milliseconds():
    qtc = fridericia([400.0], [900.0])
    assert qtc.as_qt_milliseconds().iloc[0] == pytest.approx(
        400.0 / (0.9 ** (1 / 3)))


def test_percent_endpoint_is_not_milliseconds():
    qtc = fridericia([410.0], [900.0])
    pct = to_percent_change(qtc, [400.0])
    assert pct.units == "percent"
    with pytest.raises(ProvenanceError, match="not milliseconds"):
        pct.as_qt_milliseconds()


def test_require_steps_names_what_is_missing():
    qtc = fridericia([400.0], [900.0])
    with pytest.raises(ProvenanceError) as e:
        qtc.require_steps("rate_correction", "control_subtraction")
    assert "control_subtraction" in str(e.value)
    assert qtc.require_steps("rate_correction") is qtc


def test_free_cmax_multiples_rejects_bad_cmax():
    with pytest.raises(ValueError, match="positive and finite"):
        to_free_cmax_multiples([1e-9], 0.0)
    v = to_free_cmax_multiples([2e-9, 1e-8], 2e-9)
    assert v.tolist() == pytest.approx([1.0, 5.0])


def test_standardized_effect_requires_positive_sd():
    qtc = fridericia([400.0], [900.0])
    with pytest.raises(ValueError, match="positive and finite"):
        to_standardized_effect(qtc, 0.0)


def test_provenance_travels_into_the_results_table():
    qtc = fridericia([400.0], [900.0])
    cols = provenance_columns(qtc)
    assert cols["endpoint"] == "QTcF"
    assert "fridericia" in cols["correction_cascade"]
    assert "ICH E14" in cols["citations"]
    frame = qtc.to_frame()
    assert frame.attrs["domain"] == "clinical"
    assert frame.attrs["provenance"]


def test_length_mismatch_between_df_and_endpoint_is_caught():
    from dish2ecg import subtract_baseline
    df = pd.DataFrame({"s": [1, 2], "a": ["x", "x"]})
    ep = fridericia([400.0], [900.0])
    with pytest.raises(ProvenanceError, match="different lengths"):
        subtract_baseline(df, ep, subject_col="s", arm_col="a",
                          baseline_mask=pd.Series([True, False]))
