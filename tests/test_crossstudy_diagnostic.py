"""Pins the cross-study diagnostic so the documented explanation cannot drift.

These read the committed results tables; they do not need the PhysioNet data.
"""
from pathlib import Path

import pandas as pd

RES = Path(__file__).resolve().parents[1] / "results"


def _t():
    return pd.read_csv(RES / "placebo_commonmode_test.csv").set_index("drug")


def test_ecgcipa_reads_higher_for_all_three_shared_drugs():
    assert (_t()["shift_ddQTcF"] > 0).all()


def test_removing_placebo_subtraction_shrinks_every_shift():
    t = _t()
    assert (t["shift_dQTcF_only"] < t["shift_ddQTcF"]).all()


def test_placebo_group_explains_verapamil_but_not_ranolazine_or_dofetilide():
    t = _t()
    assert abs(t.loc["Verapamil", "shift_dQTcF_only"]) < 2.0
    assert t.loc["Ranolazine", "shift_dQTcF_only"] > 4.0
    assert t.loc["Dofetilide", "shift_dQTcF_only"] > 4.0


def test_replication_table_labels_dofetilide_as_a_crossover():
    r = pd.read_csv(RES / "cross_study_replication.csv").set_index("drug")
    assert "crossover" in r.loc["Dofetilide", "ecgcipa_design"]
    assert "parallel" in r.loc["Ranolazine", "ecgcipa_design"]
    assert "parallel" in r.loc["Verapamil", "ecgcipa_design"]


def test_power_law_swap_moves_ecgrdvq_jtpeakc_by_at_most_1p5_ms():
    j = pd.read_csv(RES / "ecgrdvq_jtpeakc_power_sensitivity.csv").set_index("arm")
    linear = {"Dofetilide": 38.7, "Quinidine Sulph": 27.7,
              "Ranolazine": 2.8, "Verapamil HCL": 0.6}
    for arm, v in linear.items():
        assert abs(j.loc[arm, "ddJTpeakc_power"] - v) <= 1.5, arm
