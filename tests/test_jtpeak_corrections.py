"""The two J-Tpeak corrections are different, and the difference matters.

Reference numbers measured on the 5740 complete records of PhysioNet ECGCIPA
v1.0.0, whose ADEG table ships the study's own derived JTPC:

* their JTPC is a pure power law, exponent recovered as 0.5799 with intercept
  0.000 and residual SD 0.15 ms;
* our Fridericia QTcF matches their QTCF with median difference 0.000 ms,
  SD 0.190 ms, and only 3 of 5740 records differing by more than 0.05 ms;
* the Hnatkova 2019 linear correction differs from their power-law JTPC by a
  median of 0.44 ms, SD 3.56 ms, maximum 21.1 ms.
"""
import numpy as np
import pytest

from dish2ecg import jtpeakc_power, jtpeakc_hnatkova2019, rate_correct
from dish2ecg.corrections import JTPEAK_POWER_ALPHA_CIPA


def test_power_correction_matches_its_formula():
    # RR = 0.8 s, alpha = 0.58 -> JTp / 0.8**0.58
    got = jtpeakc_power([200.0], [800.0]).values.iloc[0]
    assert got == pytest.approx(200.0 / 0.8 ** 0.58)
    assert jtpeakc_power([200.0], [1000.0]).values.iloc[0] == pytest.approx(200.0)


def test_the_cipa_exponent_is_the_documented_one():
    assert JTPEAK_POWER_ALPHA_CIPA == pytest.approx(0.58)


def test_the_two_corrections_disagree_materially():
    """They are not interchangeable: at 100 bpm they differ by several ms."""
    rr = np.array([600.0, 800.0, 1000.0, 1200.0])
    jtp = np.array([180.0, 200.0, 220.0, 235.0])
    a = jtpeakc_power(jtp, rr).values.to_numpy()
    b = jtpeakc_hnatkova2019(jtp, rr).values.to_numpy()
    assert np.max(np.abs(a - b)) > 3.0
    # both are the identity at RR = 1 s
    assert a[2] == pytest.approx(jtp[2])
    assert b[2] == pytest.approx(jtp[2])


def test_power_correction_is_steeper_at_fast_rates():
    rr, jtp = 600.0, 180.0
    assert jtpeakc_power([jtp], [rr]).values.iloc[0] > \
           jtpeakc_hnatkova2019([jtp], [rr]).values.iloc[0]


def test_dispatch_knows_both():
    assert rate_correct([200.0], [800.0], "jtpeak_power").name == "JTpeakc"
    assert rate_correct([200.0], [800.0], "hnatkova2019_jtpeak").name == "JTpeakc"


def test_provenance_records_which_one_was_used():
    p = jtpeakc_power([200.0], [800.0]).provenance[0]
    assert p.method == "jtpeak_power"
    assert "ECGCIPA" in p.citation
    v = jtpeakc_hnatkova2019([200.0], [800.0]).provenance[0]
    assert v.method == "hnatkova2019_jtpeak"
    assert "s41598-019-51491-4" in v.citation
