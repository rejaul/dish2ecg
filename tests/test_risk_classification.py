"""qNet vs the dish endpoints as TdP-risk predictors on the 12 CiPA training drugs.

NOTE: these 12 drugs are the set the ORd-CiPA model and qNet were optimised on
(Dutta 2017), so qNet's separation here is IN-SAMPLE. These tests assert that the
implementation reproduces the published training result, not that qNet generalises.
"""
import numpy as np, pandas as pd, pytest, pathlib

RES = pathlib.Path(__file__).resolve().parents[1] / "results"


@pytest.fixture(scope="module")
def X():
    return pd.read_csv(RES / "risk_classification_cipa12.csv")


def test_all_twelve_drugs_ran_at_four_doses(X):
    assert len(X) == 12 and X.n_ok.eq(4).all()
    assert sorted(X.CiPA.value_counts().tolist()) == [4, 4, 4]


def test_control_qnet_and_apd_match_published_ord_cipa():
    raw = pd.read_csv(RES / "qnet_raw.csv")
    assert raw.qnet_control.nunique() == 1 and raw.apd90_control.nunique() == 1
    assert raw.qnet_control.iloc[0] == pytest.approx(0.0785, abs=0.002)
    assert raw.apd90_control.iloc[0] == pytest.approx(286.7, abs=1.0)
    assert int((~raw.ok).sum()) == 0


def test_qnet_bands_are_completely_non_overlapping(X):
    g = X.groupby("CiPA").qnet_mean_1to4.agg(["min", "max"])
    assert g.loc[2, "max"] < g.loc[1, "min"]      # high  below intermediate
    assert g.loc[1, "max"] < g.loc[0, "min"]      # intermediate below low


def test_qnet_is_monotone_in_risk(X):
    from scipy import stats
    assert stats.spearmanr(-X.qnet_mean_1to4, X.CiPA).statistic > 0.9


def test_the_dish_finds_high_risk_but_not_low_risk(X):
    """Blinova's best predictor flags every high-risk drug and grades the rest poorly."""
    def auc(score, pos):
        s1, s0 = score[pos], score[~pos]
        return (np.greater.outer(s1, s0).sum() + 0.5 * np.equal.outer(s1, s0).sum()) / (len(s1) * len(s0))
    ead, risk = X.ead_fraction.values, X.CiPA.values
    assert auc(ead, risk == 2) == 1.0            # all four high-risk drugs on top
    assert auc(-ead, risk == 0) < 0.8            # cannot pick out the low-risk ones
    # and the dish's prolongation endpoints are weaker at the high end than its events
    assert auc(X.dish_max_ddFPDc.values, risk == 2) < 1.0


def test_cmax_outside_tested_range_is_flagged_not_silently_extrapolated(X):
    flagged = set(X.loc[~X.cmax_in_tested_range, "drug"])
    assert flagged == {"chlorpromazine", "cisapride"}
