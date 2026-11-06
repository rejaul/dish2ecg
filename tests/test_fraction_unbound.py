"""Fraction unbound for bepridil, chlorpromazine, cisapride and terfenadine.

These four are NOT in the dish-clinic bridge, so none of the repository's published results
depends on them. The tests pin the values to their stated sources and to independent
cross-checks, so a later edit cannot silently change them.
"""
import json
import numpy as np
import pandas as pd
import pytest
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
FOUR = ["bepridil", "chlorpromazine", "cisapride", "terfenadine"]


@pytest.fixture(scope="module")
def ref():
    return pd.read_csv(DATA / "drug_reference.csv").set_index("drug_name")


def test_all_four_drugs_now_have_a_value_and_a_source(ref):
    for d in FOUR:
        assert 0 < ref.loc[d, "fraction_unbound"] < 0.1, d      # all four are >90% bound
        j = json.loads(ref.loc[d, "source_for_each_value"])
        assert j["fraction_unbound"] != "not found" and len(j["fraction_unbound"]) > 80


def test_point_estimates_follow_the_stated_convention(ref):
    # label range midpoint for chlorpromazine (92-97 %) and cisapride (97.5-98 %)
    assert ref.loc["chlorpromazine", "fraction_unbound"] == pytest.approx(1 - (0.92 + 0.97) / 2, abs=1e-9)
    assert ref.loc["cisapride", "fraction_unbound"] == pytest.approx(1 - (0.975 + 0.98) / 2, abs=1e-9)
    assert ref.loc["terfenadine", "fraction_unbound"] == pytest.approx(0.03)          # 97 %
    # bepridil is a BOUND, not a measurement: geometric mean of 0.003 and 0.010
    assert ref.loc["bepridil", "fraction_unbound"] == pytest.approx(np.sqrt(0.003 * 0.010), abs=0.0006)


def test_point_lies_inside_its_tabulated_range(ref):
    for d in ("bepridil", "chlorpromazine", "cisapride"):
        r = ref.loc[d]
        assert r.fu_low <= r.fraction_unbound <= r.fu_high, d


def test_bepridil_fu_matches_the_value_implied_by_label_total_cmax(ref):
    """Crumb free Cmax / Vascor label steady-state total Cmax (2332 ng/mL, MW 366.5)."""
    crumb = pd.read_csv(DATA / "crumb2016_free_cmax.csv").assign(k=lambda x: x.drug.str.lower()).set_index("k")
    total_nM = 2332 / 366.5 * 1000
    implied = crumb.loc["bepridil", "free_cmax_nM"] / total_nM
    assert implied == pytest.approx(0.005, abs=0.0005)
    r = ref.loc["bepridil"]
    assert r.fu_low <= implied <= r.fu_high


def test_derived_total_cmax_is_free_over_fu(ref):
    for d in FOUR:
        r = ref.loc[d]
        assert r.total_cmax_molar == pytest.approx(r.free_cmax_molar / r.fraction_unbound, rel=1e-6)


def test_terfenadine_free_cmax_divergence_is_flagged_not_hidden(ref):
    """CiPA's 4 nM is 13x Crumb/Blinova's 0.3 nM; the note must say why it is not comparable."""
    crumb = pd.read_csv(DATA / "crumb2016_free_cmax.csv").assign(k=lambda x: x.drug.str.lower()).set_index("k")
    ratio = ref.loc["terfenadine", "free_cmax_molar"] * 1e9 / crumb.loc["terfenadine", "free_cmax_nM"]
    assert ratio == pytest.approx(13.3, abs=0.2)
    note = json.loads(ref.loc["terfenadine", "source_for_each_value"])["total_cmax_molar"]
    assert "CAUTION" in note and "inhibitor" in note


def test_values_used_by_the_bridge_are_unchanged(ref):
    """The four bridge drugs' fu feed published results; this table must agree with them."""
    expected = {"dofetilide": 0.35, "quinidine": 0.16, "ranolazine": 0.38, "verapamil": 0.10}
    for d, v in expected.items():
        assert ref.loc[d, "fraction_unbound"] == pytest.approx(v)
