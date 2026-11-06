"""The APD90/FPD arm comparison in Blinova 2018, and why it is confounded."""
import numpy as np, pandas as pd, pytest, pathlib

RES = pathlib.Path(__file__).resolve().parents[1] / "results"


@pytest.fixture(scope="module")
def loso():
    return pd.read_csv(RES / "loso_site_profiles.csv")


def test_the_optical_arm_is_a_single_site():
    L = pd.read_csv(RES / "loso_site_profiles.csv")
    vso = L[L.arm.str.startswith("VSO")]
    assert len(vso) == 1 and int(vso.site.iloc[0]) == 2
    assert vso.platform.iloc[0] == "CLY"


def test_every_site_profile_uses_the_full_drug_by_concentration_grid(loso):
    # 28 drugs x 4 concentrations = 112 cells; sites cover most of them
    assert loso.n_cells.min() > 80 and loso.n_cells.max() <= 112


def test_optical_gain_exceeds_one_and_every_mea_site(loso):
    v = loso[loso.site == 2].iloc[0]
    mea = loso[loso.arm == "MEA (FPD)"]
    assert v.gain == pytest.approx(1.436, abs=0.01)
    assert v.gain - 1.96 * v.gain_se > 1.0          # CI excludes unity
    assert v.gain > mea.gain.max()                   # highest of all ten sites
    assert (v.gain - mea.gain.mean()) / mea.gain.std() > 2.0


def test_paired_arm_slope_is_about_one_and_a_half():
    from scipy import stats
    t = pd.read_csv(RES / "arm_bycell.csv")
    lr = stats.linregress(t.MEA, t.VSO)
    assert lr.slope == pytest.approx(1.507, abs=0.01)
    assert lr.slope - 1.96 * lr.stderr > 1.0
    assert len(t) == 97


def test_insilico_gets_all_four_signs_right_where_the_dish_gets_three():
    M = pd.read_csv(RES / "insilico_vs_dish_vs_clinic.csv")
    assert len(M) == 4
    sil = (np.sign(M.dQT_sur) == np.sign(M.clin_ddQTcF_ms)).sum()
    dish = (np.sign(M.dish_ddFPDc_ms) == np.sign(M.clin_ddQTcF_ms)).sum()
    assert sil == 4 and dish == 3
    # verapamil is the one the dish inverts
    v = M[M.Drug == "Verapamil"].iloc[0]
    assert v.dish_ddFPDc_ms < -50 and v.dQT_sur > 0 and v.clin_ddQTcF_ms > 0


def test_uncoupled_transmural_summation_degenerates_to_the_m_cell():
    """Without electrotonic coupling the QT surrogate is just the M-cell APD90.

    True wherever the M cell repolarises last, which is 121 of 136 conditions.
    The 15 exceptions are all supratherapeutic (>= 3.78x free Cmax), where the
    M cell loses a stable single repolarisation (two fail outright) - the M cell
    destabilises first, so the surrogate falls back to the endocardial cell.
    """
    A = pd.read_csv(RES / "ord_apd90.csv")
    T = pd.read_csv(RES / "ord_transmural.csv")
    assert int((~A.ok).sum()) == 0 and len(A) == 408

    same = np.isclose(T.qt_surrogate, T["mid"])
    assert same.sum() == 121 and (~same).sum() == 15
    assert T.loc[~same, "multiple"].min() >= 3.7          # all supratherapeutic
    assert T.loc[same, "multiple"].max() <= 30            # sweep ceiling

    # at every therapeutic match point used by the bridge, it IS the M cell
    for mu in (1.074, 0.288, 0.956, 0.500):
        row = T[np.isclose(T.multiple, mu)]
        assert np.allclose(row.qt_surrogate, row["mid"])


def test_ord_control_apd90_matches_the_published_model():
    A = pd.read_csv(RES / "ord_apd90.csv")
    ctl = A.groupby("cell").apd90_ctl.first()
    assert ctl["endo"] == pytest.approx(265.5, abs=1.0)
    assert ctl["epi"] == pytest.approx(228.0, abs=1.0)
    assert ctl["mid"] == pytest.approx(329.9, abs=1.0)
    assert ctl["mid"] > ctl["endo"] > ctl["epi"]     # canonical ORd ordering
