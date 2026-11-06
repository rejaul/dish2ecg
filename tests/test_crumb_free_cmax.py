"""Crumb 2016 Table 1 free Cmax, and the invariance it does not threaten."""
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parents[1] / "data" / "crumb2016_free_cmax.csv"


def test_table_has_the_full_panel():
    d = pd.read_csv(DATA)
    assert len(d) == 30
    assert d["free_cmax_nM"].gt(0).all()


def test_matches_the_blinova_reference_for_shared_drugs_except_quinidine():
    d = pd.read_csv(DATA).set_index(lambda i: None) if False else pd.read_csv(DATA)
    cm = dict(zip(d.drug.str.lower(), d.free_cmax_nM))
    blinova = {"dofetilide": 2.0, "ranolazine": 1900.0, "verapamil": 50.0,
               "terfenadine": 0.30, "sotalol": 15000.0, "cisapride": 2.6,
               "diltiazem": 130.0, "mexiletine": 2500.0, "bepridil": 30.0,
               "chlorpromazine": 34.5, "ondansetron": 370.0}
    for k, v in blinova.items():
        assert 0.85 <= cm[k] / v <= 1.15, f"{k}: Crumb {cm[k]} vs Blinova {v}"
    # the one documented exception
    assert cm["quinidine"] / 3000.0 < 0.35
