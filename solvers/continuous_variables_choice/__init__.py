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
)
from solvers.continuous_variables_choice.values_separators import separate_values

SAMPLING_FUNCTIONS = {
    "uniform": uniform_variables_values_vector,
    "gaussian_cma_es": gaussian_variables_values_vector,
    "cnmcts": cnmcts,
    "cnrpa": cnrpa,
    "cgnrpa": cgnrpa,
    "crbnmcts": crbnmcts,
}


def get_variables_values(inputs_values: dict):
    sampling_function = inputs_values.get("sampling_function", None)
    assert (
        not sampling_function is None
    ), "Specified sampling_function is unrecognizable"
    input_vector, total_delta_velocity = SAMPLING_FUNCTIONS[sampling_function](
        **inputs_values
    )
    departure_epoch, times_of_flight = separate_values(input_vector)

    return (
        departure_epoch,
        times_of_flight,
        total_delta_velocity,
    )
