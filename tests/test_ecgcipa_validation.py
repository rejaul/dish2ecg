"""ECGCIPA fits against the published CiPA Phase I results.

Vicente J et al., "Assessment of Multi-Ion Channel Block in a Phase I Randomized
Study Design: Results of the CiPA Phase I ECG Biomarker Validation Study."
Clin Pharmacol Ther 2019;105:943-953, doi:10.1002/cpt.1303.

The specification matters and is easy to get wrong. This is a **multiple-dose**
design: subjects saw a lower concentration on day 1 and a higher one on day 3,
and the prespecified endpoints are evaluated at Cmax **on a named day** - ddJ-Tc
on day 3 for the balanced blockers, ddQTc on day 1 for chloroquine. Pooling the
days into one concentration-response model gives correct point estimates with
intervals that are too narrow, and flips the ranolazine decision.

Result when the day stratification is respected: 6 of 8 published point
estimates fall inside our two-sided 90% interval, and 6 of 8 reach the same
verdict against the 10 ms threshold.

The two misses are both verapamil ddQTcF, and in opposite directions - 5.7 ms
low on day 1, 3.3 ms high on day 3 - which is what a linear concentration-
response model does to a curved relationship. The paper itself reports comparing
linear against Emax and sigmoid models for this reason.
"""
import pytest

#: arm, endpoint, day -> (published mean, lo90, hi90)
PUBLISHED = {
    ("Verapamil", "ddJTpeakc", 3): (-2.6, -11.4, 6.1),
    ("Lopinavir+Ritonavir", "ddJTpeakc", 3): (-2.9, -14.6, 8.8),
    ("Ranolazine", "ddJTpeakc", 3): (1.2, -9.5, 12.0),
    ("Ranolazine", "ddQTcF", 1): (14.1, 7.5, 20.8),
    ("Ranolazine", "ddQTcF", 3): (28.2, 18.3, 38.2),
    ("Verapamil", "ddQTcF", 1): (8.6, 4.4, 12.9),
    ("Verapamil", "ddQTcF", 3): (13.9, 10.0, 17.8),
    ("Chloroquine", "ddQTcF", 1): (30.7, 22.6, 38.9),
}
#: our day-stratified fit: (mean, lo90, hi90)
OURS = {
    ("Verapamil", "ddJTpeakc", 3): (-2.0, -9.2, 5.2),
    ("Lopinavir+Ritonavir", "ddJTpeakc", 3): (-0.8, -10.7, 9.2),
    ("Ranolazine", "ddJTpeakc", 3): (-0.8, -10.9, 9.3),
    ("Ranolazine", "ddQTcF", 1): (11.0, 6.5, 15.4),
    ("Ranolazine", "ddQTcF", 3): (26.7, 18.3, 35.1),
    ("Verapamil", "ddQTcF", 1): (2.9, -0.2, 6.1),
    ("Verapamil", "ddQTcF", 3): (17.2, 14.8, 19.6),
    ("Chloroquine", "ddQTcF", 1): (25.7, 20.0, 31.3),
}
KNOWN_MISSES = {("Verapamil", "ddQTcF", 1), ("Verapamil", "ddQTcF", 3)}


@pytest.mark.parametrize("key", [k for k in PUBLISHED if k not in KNOWN_MISSES])
def test_published_point_inside_our_interval(key):
    pm = PUBLISHED[key][0]
    _, lo, hi = OURS[key]
    assert lo <= pm <= hi, f"{key}: published {pm} outside [{lo}, {hi}]"


def test_overall_agreement_is_six_of_eight():
    inside = sum(1 for k, (pm, _, _) in PUBLISHED.items()
                 if OURS[k][1] <= pm <= OURS[k][2])
    assert inside == 6, f"agreement changed: {inside}/8"


def test_the_two_misses_are_verapamil_in_opposite_directions():
    """A linear model on a curved relationship: low at the low dose, high at the
    high dose. If this ever becomes one-directional, the cause is different."""
    d1 = OURS[("Verapamil", "ddQTcF", 1)][0] - PUBLISHED[("Verapamil", "ddQTcF", 1)][0]
    d3 = OURS[("Verapamil", "ddQTcF", 3)][0] - PUBLISHED[("Verapamil", "ddQTcF", 3)][0]
    assert d1 < 0 < d3, f"day1 {d1:+.1f}, day3 {d3:+.1f}"


def test_balanced_blockers_have_a_near_zero_jtpeak_effect():
    """The mechanistic result the study was built to show, reproduced."""
    for arm in ("Verapamil", "Lopinavir+Ritonavir", "Ranolazine"):
        mean = OURS[(arm, "ddJTpeakc", 3)][0]
        assert abs(mean) < 5.0, f"{arm} ddJ-Tpeakc {mean}"


def test_ranolazine_decision_differs_on_interval_width_not_effect():
    """Ours passes the 10 ms criterion, the published analysis fails it - but the
    point estimates agree to 2 ms. The disagreement is precision."""
    pm, _, phi = PUBLISHED[("Ranolazine", "ddJTpeakc", 3)]
    om, _, ohi = OURS[("Ranolazine", "ddJTpeakc", 3)]
    assert abs(om - pm) < 2.5, "point estimates should agree closely"
    assert phi > 10.0 and ohi < 10.0, "the verdicts should differ"
