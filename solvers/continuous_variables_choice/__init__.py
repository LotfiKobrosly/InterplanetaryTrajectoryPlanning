"""
Implements global function to return values vector
"""
from solvers.continuous_variables_choice.cnmcts import cnmcts
from solvers.continuous_variables_choice.one_shot_vector_choice import uniform_variables_values_vector, gaussian_variables_values_vector
from solvers.continuous_variables_choice.values_vector_utils import separate_values

def get_variables_values(inputs_values: dict):
    sampling_function = inputs_values["sampling_function"]
    sampling_functions_dict = {
        "uniform": uniform_variables_values_vector,
        "gaussian_cma_es": gaussian_variables_values_vector,
        "cnmcts": cnmcts,
    }
    input_vector, total_delta_velocity = sampling_functions_dict[sampling_function](**inputs_values)
    departure_epoch, times_of_flight, flyby_parameters = separate_values(
        input_vector, len(inputs_values["planets_sequence"])
    )
    return (
        departure_epoch,
        times_of_flight,
        flyby_parameters,
        total_delta_velocity,
    )