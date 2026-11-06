"""The in-vitro rate correction against an independently published one.

Supplemental Table 1 of the CiPA MEA pilot (Millard D et al., Toxicol Sci
2018;164:550-562, doi:10.1093/toxsci/kfy110) reports mean baseline beat period,
field potential duration and rate-corrected field potential duration for all 18
participating studies. The values below were read verbatim out of the EMF text
records of that table by scripts/validate_invitro_correction.py; they are
inlined here so the test runs offline.

Interpretation: the published FPDc is the mean of per-well corrected values,
while the value computable from the table is the correction applied to the
means. They differ by Jensen's inequality, so the test asserts agreement to a
few milliseconds and that the residual grows with beat-period spread - not exact
equality.
"""
import numpy as np
import pytest
from scipy import stats

from dish2ecg import power_correction
from dish2ecg.io import emf_table, emf_text_records

# study, cell type, platform, BP mean, BP sd, FPD mean, FPDc published
MILLARD_TABLE_S1 = [
    (1, "AXG", "AXN", 861.0, 42.0, 258.0, 271.0),
    (2, "AXG", "AXN", 939.0, 27.0, 276.0, 282.0),
    (3, "AXG", "AXN", 1091.0, 95.0, 327.0, 317.0),
    (4, "AXG", "MCS", 1042.0, 118.0, 346.0, 341.0),
    (5, "CDI", "AXN", 2103.0, 164.0, 658.0, 513.0),
    (6, "CDI", "AXN", 1664.0, 133.0, 572.0, 482.0),
    (7, "CDI", "MCS", 1733.0, 215.0, 658.0, 546.0),
    (8, "CDI", "MCS", 1364.0, 307.0, 580.0, 524.0),
    (9, "CDI", "ACA", 2111.0, 128.0, 712.0, 577.0),
    (10, "SCI", "AXN", 1671.0, 156.0, 586.0, 494.0),
    (11, "AXG", "MCS", 764.0, 112.0, 268.0, 293.0),
    (12, "AXG", "ACA", 905.0, 118.0, 253.0, 263.0),
    (13, "AXG", "MCS", 954.0, 303.0, 326.0, 327.0),
    (14, "CDI", "ACA", 1493.0, 145.0, 718.0, 627.0),
    (15, "CDI", "MCS", 1848.0, 767.0, 588.0, 482.0),
    (16, "CDI", "ACA", 1646.0, 72.0, 569.0, 482.0),
    (17, "CDI", "AXN", 1918.0, 122.0, 711.0, 572.0),
    (18, "GEH", "MCS", 3372.0, 1135.0, 652.0, 437.0),
]
#: Study 9 is an unexplained outlier: 22 ms off with a LOW beat-period CV of
#: 0.061, so averaging order does not account for it. Excluded from the
#: agreement assertion and asserted separately so the exclusion stays visible.
OUTLIER_STUDY = 9


def _arrays(exclude=()):
    rows = [r for r in MILLARD_TABLE_S1 if r[0] not in exclude]
    bp = np.array([r[3] for r in rows])
    bp_sd = np.array([r[4] for r in rows])
    fpd = np.array([r[5] for r in rows])
    pub = np.array([r[6] for r in rows])
    return bp, bp_sd, fpd, pub


def test_fridericia_reproduces_published_fpdc():
    bp, _, fpd, pub = _arrays(exclude=(OUTLIER_STUDY,))
    ours = power_correction(fpd, bp, 1 / 3, name="FPDcF",
                            domain="invitro").values.to_numpy()
    err = np.abs(ours - pub)
    assert err.max() < 5.0, f"worst study off by {err.max():.1f} ms"
    assert np.median(err) < 1.5
    assert (100 * err / pub).max() < 2.0


def test_the_outlier_study_is_still_an_outlier():
    """Study 9 does not reproduce, and that is recorded rather than hidden."""
    row = next(r for r in MILLARD_TABLE_S1 if r[0] == OUTLIER_STUDY)
    ours = power_correction([row[5]], [row[3]], 1 / 3,
                            domain="invitro").values.iloc[0]
    assert abs(ours - row[6]) > 15.0
    assert row[4] / row[3] < 0.10, "its beat-period spread is small, so Jensen's " \
                                   "inequality does not explain the gap"


def test_residual_grows_with_beat_period_spread():
    bp, bp_sd, fpd, pub = _arrays(exclude=(OUTLIER_STUDY,))
    ours = power_correction(fpd, bp, 1 / 3, domain="invitro").values.to_numpy()
    rho, p = stats.spearmanr(bp_sd / bp, np.abs(ours - pub))
    assert rho > 0.4, f"rho={rho:.3f}"
    assert p < 0.10


def test_implied_exponent_is_fridericia():
    bp, _, fpd, pub = _arrays(exclude=(OUTLIER_STUDY,))
    alpha = float(np.polyfit(np.log(bp / 1000.0), np.log(fpd / pub), 1)[0])
    assert alpha == pytest.approx(1 / 3, abs=0.02), f"implied alpha {alpha:.4f}"


def test_emf_reader_on_a_synthetic_record():
    """The metafile reader is parsing the format, not guessing at bytes."""
    import struct
    text = "FPDc"
    payload = text.encode("utf-16-le")
    head = 8 + 16 + 4 + 4 + 4 + 8 + 4 + 4
    off_string = head
    size = head + len(payload)
    rec = (struct.pack("<II", 84, size) + b"\x00" * 16
           + struct.pack("<I", 2) + struct.pack("<ff", 1.0, 1.0)
           + struct.pack("<ii", 11, 22)
           + struct.pack("<I", len(text)) + struct.pack("<I", off_string)
           + payload)
    recs = emf_text_records(rec)
    assert recs == [(22, 11, "FPDc")]
    assert emf_table(rec) == [["FPDc"]]
