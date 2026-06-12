"""
Implements a Continous Nested Monte Carlo Search for the continuous decision variables
"""
import numpy as np
from utils.trajectory_evaluation import evaluate_mga_trajectory
from solvers.continuous_variables_choice.values_vector_utils import separate_values

def cnmcts(values_sequence: list = None, bounds: list = None, planets_sequence: list = None, level: int = 0, bandwidth: int = 10, *args, **kwargs):
    while len(values_sequence) < len(bounds):
        if level == 0:
            advancement = len(values_sequence)
            values_sequence.append(
                np.random.uniform(*bounds[advancement])
            )

        else:
            best_value = None
            best_delta_v = np.inf
            advancement = len(values_sequence)

            # Ensure we have the same number of possibilities (prevent redundant values)
            list_of_values = list()
            for _ in range(bandwidth):
                new_value = np.random.uniform(*bounds[advancement])
                if not (new_value in list_of_values):
                    list_of_values.append(new_value)

            # Iterating
            for new_value in list_of_values:
                temporary_values_sequence, delta_v = cnmcts(
                    values_sequence=values_sequence + [new_value],
                    bounds=bounds,
                    planets_sequence=planets_sequence,
                    level=level - 1,
                    bandwidth=bandwidth,
                )
                if delta_v < best_delta_v:
                    best_delta_v = delta_v
                    best_value = new_value

            values_sequence.append(best_value)

    return (
        values_sequence,
        evaluate_mga_trajectory(
            planets_sequence, *separate_values(values_sequence, len(planets_sequence))
        )[0],
    )