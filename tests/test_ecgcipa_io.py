"""The ADEG pivot must not silently drop rows on sparse index keys."""
import numpy as np
import pandas as pd
import pytest

from dish2ecg.io import ecgcipa_wide


def _adeg():
    """Minimal ADEG with the two sparse keys that broke the naive pivot:
    ABLFL is 'Y' only on the derived baseline row, EGREPNUM is missing there."""
    rows = []
    for rep in (1, 2, 3):
        for param, val in (("QT", 400.0 + rep), ("RR", 900.0), ("JTP", 200.0)):
            rows.append(dict(USUBJID="S1", TRTA="Drug", APERIOD=1, ATPTN=5,
                             ATPT="2 hrs", EGREPNUM=float(rep), ABLFL=np.nan,
                             PARAMCD=param, AVAL=val))
    for param, val in (("QT", 395.0), ("RR", 905.0), ("JTP", 198.0)):
        rows.append(dict(USUBJID="S1", TRTA="Drug", APERIOD=1, ATPTN=1,
                         ATPT="Baseline", EGREPNUM=np.nan, ABLFL="Y",
                         PARAMCD=param, AVAL=val))
    return pd.DataFrame(rows)


def test_pivot_keeps_every_record():
    w = ecgcipa_wide(_adeg())
    assert len(w) == 4, "3 replicates + 1 derived baseline row"
    assert {"QT", "RR", "JTP"} <= set(w.columns)
    assert w["QT"].notna().all()


def test_baseline_flag_is_usable_as_a_boolean():
    w = ecgcipa_wide(_adeg())
    assert w["ABLFL"].isin(["Y", "N"]).all()
    assert int(w["ABLFL"].eq("Y").sum()) == 1


def test_replicate_number_is_filled_not_dropped():
    w = ecgcipa_wide(_adeg())
    assert set(w["EGREPNUM"]) == {0, 1, 2, 3}
    assert w.loc[w["ABLFL"].eq("Y"), "EGREPNUM"].iloc[0] == 0


def test_empty_pivot_raises_rather_than_returning_nothing():
    bad = _adeg().assign(PARAMCD=np.nan)
    with pytest.raises(Exception):
        ecgcipa_wide(bad)
