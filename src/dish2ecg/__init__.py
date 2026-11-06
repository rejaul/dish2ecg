"""dish2ecg - corrected repolarization endpoints from dish and clinic, on one scale.

The package implements one correction cascade and applies it symmetrically to
hiPSC-CM field-potential / action-potential recordings and to clinical 12-lead ECG
data, so that a concentration-response relationship measured in a dish can be
compared with a clinical concentration-QTc relationship.

It does not convert a field potential duration into a QT interval, and it refuses to
try - see :meth:`dish2ecg.provenance.Endpoint.as_qt_milliseconds`.
"""
from .provenance import (CorrectionRecord, DomainError, Endpoint, ProvenanceError,
                         provenance_columns)
from .corrections import (bazett, fit_individual_exponents, fit_power_exponent,
                          framingham, fridericia, hodges, individual_correction,
                          jtpeakc_power, jtpeakc_hnatkova2019, jtpeakc_vicente2019,
                          power_correction, rate_correct)
from .cascade import (double_delta, hysteresis_corrected_rr, subtract_baseline,
                      subtract_control)
from .models import ConcentrationResponseFit, emax_fit, fit_concentration_response
from .bridge import (BridgeComparison, compare_slopes, decision_agreement,
                     to_free_cmax_multiples, to_percent_change,
                     to_standardized_effect)

__version__ = "0.1.0"

__all__ = [
    "CorrectionRecord", "Endpoint", "ProvenanceError", "DomainError",
    "provenance_columns",
    "power_correction", "bazett", "fridericia", "framingham", "hodges",
    "jtpeakc_hnatkova2019",
    "jtpeakc_vicente2019",  # deprecated alias, see corrections.py
    "jtpeakc_vicente2019",  # deprecated alias "jtpeakc_power", "fit_power_exponent",
    "fit_individual_exponents",
    "individual_correction", "rate_correct",
    "subtract_baseline", "subtract_control", "double_delta",
    "hysteresis_corrected_rr",
    "fit_concentration_response", "ConcentrationResponseFit", "emax_fit",
    "to_free_cmax_multiples", "to_percent_change", "to_standardized_effect",
    "compare_slopes", "decision_agreement", "BridgeComparison",
    "__version__",
]

# `viz` is imported on demand: it needs matplotlib, which the core does not.
from . import viz as viz  # noqa: E402,F401
__all__ = list(__all__) + ["viz"]
