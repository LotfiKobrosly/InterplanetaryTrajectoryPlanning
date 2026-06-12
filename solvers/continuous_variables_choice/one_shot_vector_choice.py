"""
Implements a one-shot choice of the control variables in a vector.

For the gaussian-based implementation, we update the means and covariance with the CMA-ES (Covariance Matrix Adaptation Evolution Strategy)
"""

import warnings
import numpy as np
import cma
from utils.trajectory_evaluation import evaluate_mga_trajectory
from utils.constants import VARIABLES_BOUNDS

warnings.filterwarnings("ignore")


def separate_values(input_vector: np.ndarray, planets_sequence_length: int):
    """
    Simply seperates the values vector into the desired variables compatible with Trajectory class from classes.trajectory
    """
    departure_epoch = input_vector[0]
    times_of_flight = input_vector[1:planets_sequence_length]
    flyby_parameters = input_vector[planets_sequence_length:]
    flyby_parameters = [
        (flyby_parameters[i], flyby_parameters[i + 1])
        for i in range(0, len(flyby_parameters), 2)
    ]
    return departure_epoch, times_of_flight, flyby_parameters


def uniform_variables_values_vector(
    bounds: list, planets_sequence: list, n_iterations: int = 500, *args, **kwargs
):
    """
    Returns a uniformly sampled vector from the given bounds
    """
    best_vector = None
    best_value = 1e30
    for iteration in range(n_iterations):
        # Imposing a strict while loop here without the max iterations condition might
        # result in an endless loop, given the pure random aspect of this algorithm
        vector = np.random.uniform(
            np.array([bound[0] for bound in bounds]),
            np.array([bound[1] for bound in bounds]),
        )
        result = evaluate_mga_trajectory(
            planets_sequence, *separate_values(vector, len(planets_sequence))
        )

        if not (result is None) and (result[0] < best_value):
            best_value = result[0]
            best_vector = vector
    if best_vector is None:
        best_vector = vector
    return best_vector


def gaussian_variables_values_vector(
    bounds: list, planets_sequence: list, n_iterations: int = 500, *args, **kwargs
):
    """
    Returns a sample of a fitted normal distribution through a CMA-ES.
    """
    lower_bounds = np.array([bound[0] for bound in bounds])
    upper_bounds = np.array([bound[1] for bound in bounds])

    def normalize(x, low, high):
        return [(xi - l) / (h - l) for xi, l, h in zip(x, low, high)]

    def denormalize(x, low, high):
        return [xi * (h - l) + l for xi, l, h in zip(x, low, high)]

    def objective_function(normalized_vector: np.ndarray):
        result = evaluate_mga_trajectory(
            planets_sequence,
            *separate_values(
                denormalize(normalized_vector, lower_bounds, upper_bounds),
                len(planets_sequence),
            )
        )
        if result is None:
            return 1e10  # penalty for invalid trajectory

        return result[0] / 1000

    # --- Initial mean: center of search space (normalized) ---
    initial_input = [0.95] * len(bounds)
    sigma0 = 10  # initial step size (in normalized space)

    # --- CMA-ES options ---
    options = cma.CMAOptions()
    options["bounds"] = [[0.0] * len(bounds), [1.0] * len(bounds)]  # normalized bounds
    options["maxiter"] = n_iterations
    options["popsize"] = 20  # λ : samples per iteration
    # options['CMA_diagonal'] = True
    options["tolx"] = 1e-6  # convergence on x
    options["tolfun"] = 1e-6  # convergence on f
    options["verbose"] = -9

    # --- Run ---
    estimator = cma.CMAEvolutionStrategy(initial_input, sigma0, options)

    while not estimator.stop():
        candidates = estimator.ask()  # sample λ candidates
        fitnesses = [objective_function(x) for x in candidates]
        estimator.tell(candidates, fitnesses)  # update means, covariance, sigma
        estimator.logger.add()
        estimator.disp()

    result_normalized = estimator.result.xbest
    final_vector = denormalize(result_normalized, lower_bounds, upper_bounds)
    # print("Chosen vector:", final_vector)
    return final_vector


def get_variables_values(bounds: list, planets_sequence: list, sampling_function: str):
    sampling_functions_dict = {
        "uniform": uniform_variables_values_vector,
        "gaussian_cma_es": gaussian_variables_values_vector,
    }
    input_vector = sampling_functions_dict[sampling_function](bounds, planets_sequence)
    return separate_values(input_vector, len(planets_sequence))
