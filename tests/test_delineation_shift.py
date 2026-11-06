"""The T-end rule does not explain the 6-7 ms between-study shift."""
import numpy as np, pandas as pd, pytest, pathlib

RES = pathlib.Path(__file__).resolve().parents[1] / "results"


@pytest.fixture(scope="module")
def dd():
    return pd.read_csv(RES / "delineation_cascade_effect.csv")


def test_both_studies_were_delineated_with_the_same_code():
    D = pd.read_csv(RES / "delineation_both.csv")
    assert set(D.study) == {"ecgrdvq", "ecgcipa"}
    assert (D.groupby("study").size() > 2000).all()
    # the two rules differ per record, and differ BETWEEN studies
    med = D[D.plausible].groupby("study").dQT.median()
    assert med["ecgrdvq"] == pytest.approx(-4.4, abs=0.6)
    assert med["ecgcipa"] == pytest.approx(-7.0, abs=0.6)


def test_the_cascade_absorbs_almost_all_of_the_delineation_difference(dd):
    """Per record the rules differ by 4-7 ms; after baseline+placebo, under 2 ms."""
    assert dd.ddelta.abs().max() < 2.5
    for r in dd.itertuples():
        assert r.lo < r.ddelta < r.hi


def test_delineation_cannot_account_for_the_between_study_shift(dd):
    w = dd.pivot(index="drug", columns="study", values="ddelta")
    between = (w["ecgrdvq"] - w["ecgcipa"]).abs()
    assert between.max() < 2.5          # against a 6-7 ms discrepancy
    assert between["Ranolazine"] < 1.0
