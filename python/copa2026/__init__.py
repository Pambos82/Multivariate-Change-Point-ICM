"""Feature-based conformal martingale benchmark and figure workflows."""

__version__ = "0.2.0"

from .betting import combination, prepare_knn, evaluate_knn
from .experiment import (
    conformal_pvalues, pvalues, ellipsoid_density, ncms_from_arrays,
    run_experiment, run_replay, run_knn_cm_experiment,
)

__all__ = [
    "combination", "prepare_knn", "evaluate_knn", "pvalues",
    "conformal_pvalues", "ellipsoid_density", "ncms_from_arrays",
    "run_experiment", "run_replay", "run_knn_cm_experiment",
]
