import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# allow running the suite from a checkout without installing
_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))


@pytest.fixture
def crossover():
    """Synthetic crossover study with a known truth.

    Design mirrors ECGRDVQ: 8 subjects x 3 arms (placebo, drugA, drugB), each arm
    with a pre-dose BASELINE day and a DOSE day, both sampled at the same 6 clock
    times, so a time-matched baseline exists for every post-dose record.

    Truth built in
    --------------
    * QT depends on RR with exponent exactly 1/3, so Fridericia removes the rate
      dependence completely.
    * A subject offset (380 + 10*s ms) is present on both days.
    * A circadian effect varies with clock time and is present on both days and in
      every arm, so only placebo subtraction removes it.
    * Drug A adds 2.0 ms per ng/mL; drug B adds nothing.
    """
    rng = np.random.default_rng(42)
    arms = {"placebo": 0.0, "drugA": 2.0, "drugB": 0.0}
    times = [0, 1, 2, 3, 4, 5]
    circadian = {0: 0.0, 1: 4.0, 2: 7.0, 3: 5.0, 4: 2.0, 5: -1.0}
    conc_profile = {0: 0.0, 1: 10.0, 2: 20.0, 3: 15.0, 4: 8.0, 5: 3.0}
    rows = []
    for s in range(8):
        offset = 380.0 + 10.0 * s
        for arm, slope in arms.items():
            for day in ("baseline", "dose"):
                for t in times:
                    on_drug = (day == "dose") and arm != "placebo"
                    conc = conc_profile[t] if on_drug else 0.0
                    rr = rng.uniform(800.0, 1200.0)
                    qtc_true = offset + circadian[t] + slope * conc
                    rows.append(dict(
                        subject=s, arm=arm, day=day, time=t, rr=rr,
                        qt=qtc_true * (rr / 1000.0) ** (1 / 3),
                        conc=conc, is_baseline=(day == "baseline"),
                    ))
    df = pd.DataFrame(rows)
    df.attrs["true_slope"] = {"drugA": 2.0, "drugB": 0.0}
    df.attrs["true_effect_at_conc20"] = {"drugA": 40.0, "drugB": 0.0}
    return df


@pytest.fixture
def crossover_noisy():
    """Same design as ``crossover`` but with realistic measurement noise.

    The noiseless fixture gives a zero-width confidence interval, which cannot test
    interval coverage. Here each ECG carries 6 ms of residual noise and each subject
    a 12 ms random intercept, so the mixed model has something to estimate.
    """
    rng = np.random.default_rng(7)
    arms = {"placebo": 0.0, "drugA": 2.0}
    times = [1, 2, 3, 4, 5]
    circadian = {1: 4.0, 2: 7.0, 3: 5.0, 4: 2.0, 5: -1.0}
    conc_profile = {1: 10.0, 2: 20.0, 3: 15.0, 4: 8.0, 5: 3.0}
    rows = []
    for s in range(24):
        offset = 380.0 + rng.normal(0, 12.0)
        for arm, slope in arms.items():
            for day in ("baseline", "dose"):
                for t in [0] + times:
                    on_drug = (day == "dose") and arm != "placebo"
                    conc = conc_profile.get(t, 0.0) if on_drug else 0.0
                    rr = rng.uniform(800.0, 1200.0)
                    qtc = (offset + circadian.get(t, 0.0) + slope * conc
                           + rng.normal(0, 6.0))
                    rows.append(dict(subject=s, arm=arm, day=day, time=t, rr=rr,
                                     qt=qtc * (rr / 1000.0) ** (1 / 3), conc=conc,
                                     is_baseline=(day == "baseline")))
    df = pd.DataFrame(rows)
    df.attrs["true_slope"] = {"drugA": 2.0}
    df.attrs["true_effect_at_conc20"] = {"drugA": 40.0}
    return df


@pytest.fixture
def simulate_noisy_crossover():
    """Factory: build a noisy crossover study for a given seed.

    24 subjects x (placebo, drugA), baseline and dose days, 6 ms residual noise and
    a 12 ms subject random intercept, so the mixed model has variance to estimate
    and the confidence interval has non-zero width.
    """
    def _make(seed: int) -> pd.DataFrame:
        rng = np.random.default_rng(seed)
        arms = {"placebo": 0.0, "drugA": 2.0}
        circ = {1: 4.0, 2: 7.0, 3: 5.0, 4: 2.0, 5: -1.0}
        prof = {1: 10.0, 2: 20.0, 3: 15.0, 4: 8.0, 5: 3.0}
        rows = []
        for s in range(24):
            off = 380.0 + rng.normal(0, 12.0)
            for arm, slope in arms.items():
                for day in ("baseline", "dose"):
                    for t in [0, 1, 2, 3, 4, 5]:
                        on = (day == "dose") and arm != "placebo"
                        conc = prof.get(t, 0.0) if on else 0.0
                        rr = rng.uniform(800.0, 1200.0)
                        qtc = (off + circ.get(t, 0.0) + slope * conc
                               + rng.normal(0, 6.0))
                        rows.append(dict(subject=s, arm=arm, day=day, time=t, rr=rr,
                                         qt=qtc * (rr / 1000.0) ** (1 / 3),
                                         conc=conc, is_baseline=(day == "baseline")))
        return pd.DataFrame(rows)
    return _make
