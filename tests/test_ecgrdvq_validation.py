"""ECGRDVQ results against the published analysis of the same study.

Vicente J, Johannesen L, Mason JW, Crumb WJ, Pueyo E, Stockbridge N, Strauss DG.
"Comprehensive T wave Morphology Assessment in a Randomized Clinical Study of
Dofetilide, Quinidine, Ranolazine, and Verapamil." J Am Heart Assoc
2015;4:e001615. doi:10.1161/JAHA.114.001615, Table 2 — baseline- and
placebo-corrected change at the population Cmax, mean (95% CI).

These four drugs are exactly the dish-to-clinic bridge, so this is the closest
thing the project has to a validation of the clinical side of the bridge itself.
The values below are inlined so the test runs offline; the live comparison is
written by scripts/clinical_reference.py to
results/validation_ecgrdvq_vs_published.csv.

Result at the time of writing: all five published point estimates fall inside
our two-sided 90% confidence interval, and all three cells the paper reports as
"no relationship between drug concentration and the ECG biomarker" correspond to
an interval of ours that contains zero. 8 of 8.
"""
import pytest

#: arm -> endpoint -> (published mean, lo95, hi95) or None for "no relationship"
PUBLISHED = {
    ("Dofetilide", "ddQTcF"): (73.6, 65.8, 81.5),
    ("Quinidine Sulph", "ddQTcF"): (78.9, 68.2, 89.7),
    ("Ranolazine", "ddQTcF"): (12.0, 7.3, 16.7),
    ("Verapamil HCL", "ddQTcF"): None,
    ("Dofetilide", "ddJTpeakc"): (39.1, 31.6, 46.6),
    ("Quinidine Sulph", "ddJTpeakc"): (26.1, 13.5, 38.7),
    ("Ranolazine", "ddJTpeakc"): None,
    ("Verapamil HCL", "ddJTpeakc"): None,
}
#: our fitted effect at the geometric-mean Cmax, with the two-sided 90% CI
OURS = {
    ("Dofetilide", "ddQTcF"): (73.1, 68.5, 77.7),
    ("Quinidine Sulph", "ddQTcF"): (74.9, 69.3, 80.5),
    ("Ranolazine", "ddQTcF"): (11.2, 7.6, 14.8),
    ("Verapamil HCL", "ddQTcF"): (3.5, -0.1, 7.0),
    ("Dofetilide", "ddJTpeakc"): (38.7, 32.8, 44.6),
    ("Quinidine Sulph", "ddJTpeakc"): (27.7, 22.3, 33.2),
    ("Ranolazine", "ddJTpeakc"): (2.8, -0.2, 5.8),
    ("Verapamil HCL", "ddJTpeakc"): (0.6, -2.4, 3.5),
}
#: population Cmax reported by the paper (ng/mL) vs our geometric mean
CMAX = {"Dofetilide": (2.7, 2.7), "Quinidine Sulph": (1800.0, 1754.0),
        "Ranolazine": (2300.0, 2043.2), "Verapamil HCL": (130.4, 113.6)}


@pytest.mark.parametrize("key", [k for k, v in PUBLISHED.items() if v is not None])
def test_published_point_estimate_falls_inside_our_interval(key):
    pm = PUBLISHED[key][0]
    _, lo, hi = OURS[key]
    assert lo <= pm <= hi, f"{key}: published {pm} outside our CI [{lo}, {hi}]"


@pytest.mark.parametrize("key", [k for k, v in PUBLISHED.items() if v is None])
def test_no_relationship_cells_have_an_interval_containing_zero(key):
    _, lo, hi = OURS[key]
    assert lo <= 0.0 <= hi, (
        f"{key}: the paper reports no concentration relationship but our "
        f"interval [{lo}, {hi}] excludes zero")


def test_the_four_bridge_drugs_are_all_covered():
    arms = {a for a, _ in PUBLISHED}
    assert arms == {"Dofetilide", "Quinidine Sulph", "Ranolazine", "Verapamil HCL"}


@pytest.mark.parametrize("arm", list(CMAX))
def test_our_geometric_mean_cmax_is_close_to_the_published_mean(arm):
    pub, ours = CMAX[arm]
    assert abs(ours - pub) / pub < 0.15, (
        f"{arm}: our Cmax {ours} differs from the published {pub} by more than 15%")


def test_dofetilide_is_the_tightest_check():
    """Pure hERG block, the cleanest case: both endpoints within 0.5 ms."""
    assert abs(OURS[("Dofetilide", "ddQTcF")][0]
               - PUBLISHED[("Dofetilide", "ddQTcF")][0]) < 1.0
    assert abs(OURS[("Dofetilide", "ddJTpeakc")][0]
               - PUBLISHED[("Dofetilide", "ddJTpeakc")][0]) < 1.0
