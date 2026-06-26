"""
Implements global function to return values vector
"""

from solvers.continuous_variables_choice.cnmcts import cnmcts
from solvers.continuous_variables_choice.crbnmcts import crbnmcts
from solvers.continuous_variables_choice.cnrpa import cnrpa
from solvers.continuous_variables_choice.cgnrpa import cgnrpa
from solvers.continuous_variables_choice.one_shot_vector_choice import (
    uniform_variables_values_vector,
    gaussian_variables_values_vector,
    genetic_algorithm,
)

SAMPLING_FUNCTIONS = {
    "uniform": uniform_variables_values_vector,
    "gaussian_cma_es": gaussian_variables_values_vector,
    "cnmcts": cnmcts,
    "cnrpa": cnrpa,
    "cgnrpa": cgnrpa,
    "crbnmcts": crbnmcts,
    "genetic": genetic_algorithm,
}


def get_variables_values(inputs_values: dict):
    sampling_function = inputs_values.get("sampling_function", None)
    assert (
        not sampling_function is None
    ), "Specified sampling_function is unrecognizable"

    return SAMPLING_FUNCTIONS[sampling_function](
        **inputs_values
    )
