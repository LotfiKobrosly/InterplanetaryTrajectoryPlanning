"""
Implements a Continous Nested Monte Carlo Search for the continuous decision variables
Here, the values sequence is represented as follows:
values_sequence = [departure_epoch, time_of_flight_0, ..., time_of_flight_n]
"""

import numpy as np
from utils.trajectory_evaluation import evaluate_mga_trajectory
from utils.constants import RANDOM_GENERATOR
from solvers.continuous_variables_choice.values_separators import (
    separate_values,
)


def cnmcts(
    values_sequence: list = list(),
    bounds: list = None,
    planets_sequence: list = None,
    level: int = 0,
    bandwidth: int = 10,
    *args,
    **kwargs
):
    if level == 0:
        while len(values_sequence) < len(planets_sequence):
            values_sequence.append(
                RANDOM_GENERATOR.uniform(*bounds[len(values_sequence)])
            )

    else:
        best_result = np.inf
        best_sequence = None
        while len(values_sequence) < len(planets_sequence):
            temporary_values_sequences = [values_sequence[:] for _ in range(bandwidth)]
            for new_values_sequence in temporary_values_sequences:

                new_values_sequence.append(
                    RANDOM_GENERATOR.uniform(*bounds[len(values_sequence)])
                )
                new_values_sequence, result = cnmcts(
                    values_sequence=new_values_sequence,
                    bounds=bounds,
                    planets_sequence=planets_sequence,
                    level=level - 1,
                    bandwidth=bandwidth,
                )
                if result < best_result:
                    best_sequence = new_values_sequence
                    best_result = result
            values_sequence.append(best_sequence[len(values_sequence)])

    return (
        values_sequence,
        evaluate_mga_trajectory(
            planets_sequence,
            *separate_values(values_sequence),
        )[0],
    )
