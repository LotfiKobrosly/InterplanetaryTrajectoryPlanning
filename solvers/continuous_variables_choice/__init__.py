"""
Implements global function to return values vector
"""
from solvers.continuous_variables_choice.cnmcts import cnmcts
from solvers.continuous_variables_choice.cnrpa import cnrpa
from solvers.continuous_variables_choice.one_shot_vector_choice import uniform_variables_values_vector, gaussian_variables_values_vector
from solvers.continuous_variables_choice.values_separators import *

SAMPLING_FUNCTIONS = {
    "uniform": uniform_variables_values_vector,
    "gaussian_cma_es": gaussian_variables_values_vector,
    "cnmcts": cnmcts,
}
SEPERATION_FUNCTIONS = {
    "uniform": separate_vector_values,
    "gaussian_cma_es": separate_vector_values,
    "cnmcts": separate_sequence_values,
    "cnrpa": separate_sequence_values,
    "cgnrpa": separate_sequence_values,
}

def get_variables_values(inputs_values: dict):
    sampling_function = inputs_values.get("sampling_function", None)
    assert not sampling_function is None, "Specified sampling_function is unrecognizable"
    input_vector, total_delta_velocity = SAMPLING_FUNCTIONS[sampling_function](**inputs_values)
    departure_epoch, times_of_flight, flyby_parameters = SEPERATION_FUNCTIONS[sampling_function](
        input_vector, len(inputs_values["planets_sequence"])
    )

    return (
        departure_epoch,
        times_of_flight,
        flyby_parameters,
        total_delta_velocity,
    )