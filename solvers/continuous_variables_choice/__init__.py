"""
Implements global function to return values vector
"""

from .cnmcts import cnmcts
from .crbnmcts import crbnmcts
from .cnrpa import cnrpa
from .cgnrpa import cgnrpa
from .cabgnrpa import cabgnrpa
from .baselines import *

SAMPLING_FUNCTIONS = {
    "cmaes": pygmo_baseline,
    "sade": pygmo_baseline,
    "sga": pygmo_baseline,
    "simulated_annealing": pygmo_baseline,
    "pso": pygmo_baseline,
    "gaco": pygmo_baseline,
    "uniform": uniform_variables_values_vector,
    "gaussian_cma_es": gaussian_variables_values_vector,
    "cnmcts": cnmcts,
    "cnrpa": cnrpa,
    "cgnrpa": cgnrpa,
    "crbnmcts": crbnmcts,
    "cabgnrpa": cabgnrpa,
    "genetic": genetic_algorithm,
}


def get_variables_values(inputs_values: dict):
    sampling_function = inputs_values.get("sampling_function", None)
    assert (
        not sampling_function is None
    ), "Specified sampling_function is unrecognizable"

    return SAMPLING_FUNCTIONS[sampling_function](**inputs_values)
