"""Out-of-sample TdP-risk test: all 28 CiPA drugs, 16 held out from qNet."""
import numpy as np, pandas as pd, pytest, pathlib
from scipy import stats

RES = pathlib.Path(__file__).resolve().parents[1] / "results"


def auc(score, pos):
    s1, s0 = score[pos], score[~pos]
    return (np.greater.outer(s1, s0).sum() + 0.5 * np.equal.outer(s1, s0).sum()) / (len(s1) * len(s0))


@pytest.fixture(scope="module")
def X():
    return pd.read_csv(RES / "risk28_full.csv")


def test_panel_is_the_full_cipa_set_split_12_16(X):
    assert len(X) == 28
    assert X.split.value_counts().to_dict() == {"validation": 16, "training": 12}
    assert sorted(X.CiPA.value_counts().tolist()) == [8, 9, 11]


def test_every_drug_dose_simulation_succeeded():
    raw = pd.read_csv(RES / "qnet28_raw.csv")
    assert len(raw) == 2240 and int((~raw.ok).sum()) == 0
    assert raw.qnet_control.iloc[0] == pytest.approx(0.0785, abs=0.002)


def test_bootstrap_median_reproduces_the_published_point_fit():
    """The 12 training drugs have a published optimal hERG fit; the ensemble agrees."""
    a = pd.read_csv(RES / "qnet28_by_drug.csv")[["drug", "qnet_median"]]
    b = pd.read_csv(RES / "qnet_by_drug.csv")[["drug", "qnet_mean_1to4"]]
    m = a.merge(b, on="drug")
    assert len(m) == 12
    d = (m.qnet_median - m.qnet_mean_1to4).abs()
    assert d.median() < 0.001 and d.max() < 0.005
    assert stats.spearmanr(m.qnet_median, m.qnet_mean_1to4).statistic > 0.95


def test_qnet_separates_high_risk_perfectly_including_held_out(X):
    assert X.loc[X.CiPA == 2, "qnet_median"].max() < X.loc[X.CiPA != 2, "qnet_median"].min()
    assert auc(-X.qnet_median.values, X.CiPA.values == 2) == 1.0
    V = X[X.split == "validation"]
    assert auc(-V.qnet_median.values, V.CiPA.values == 2) == 1.0


def test_the_dish_matches_qnet_on_clearing_low_risk_out_of_sample(X):
    """Withdraws the 12-drug claim that the dish cannot grade safety."""
    V = X[X.split == "validation"]
    low = V.CiPA.values == 0
    # higher qNet = safer, so qNet itself is the low-risk score;
    # higher ddFPDc = riskier, so the dish score is negated.
    qnet_low = auc(V.qnet_median.values, low)
    dish_low = auc(-V.dish_at_cmax.values, low)
    assert dish_low >= qnet_low          # dish is at least as good at low-risk, held out
    assert dish_low > 0.9


def test_qnet_false_alarms_are_named(X):
    """Two low-risk drugs sit inside the intermediate band - both false alarms, not misses."""
    lo, hi = X.loc[X.CiPA == 1, "qnet_median"].min(), X.loc[X.CiPA == 1, "qnet_median"].max()
    inside = set(X[(X.CiPA == 0) & X.qnet_median.between(lo, hi)].drug_name)
    assert inside == {"metoprolol", "tamoxifen"}
    assert X.loc[X.CiPA == 0, "qnet_median"].min() > X.loc[X.CiPA == 2, "qnet_median"].max()
