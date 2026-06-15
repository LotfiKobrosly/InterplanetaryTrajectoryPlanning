"""
Implements a Continous Nested Monte Carlo Search for the continuous decision variables
Here, the values sequence is represented as follows:
values_sequence = [departure_epoch, (time_of_flight_0, radius_0, angle_0), ..., (time_of_flight_(n-1), radius_(n-1), angle_(n-1), time_of_flight_n]
"""
import numpy as np
from utils.trajectory_evaluation import evaluate_mga_trajectory
from solvers.continuous_variables_choice.values_separators import separate_sequence_values

def cnmcts(values_sequence: list = None, bounds: list = None, planets_sequence: list = None, level: int = 0, bandwidth: int = 10, *args, **kwargs):
    while len(values_sequence) < len(planets_sequence):
        advancement = len(values_sequence)
        if advancement == 0 or (advancement == (len(planets_sequence) - 1)):
            # Departure epoch or final time of flight
            low_bound, high_bound = bounds[advancement]
        else:
            # Time of flight and flyby parameters
            low_bound = np.array([bound[0] for bound in bounds[advancement]])
            high_bound = np.array([bound[1] for bound in bounds[advancement]])
        
        if level == 0:
            chosen_move = np.random.uniform(low_bound, high_bound)

        else:
            best_result = np.inf
            chosen_move = None
            sampled_values_list = list()
            for _ in range(bandwidth):
                sampled_value = np.random.uniform(low=low_bound, high=high_bound)
                sampled_values_list.append(sampled_value)
            for value in sampled_values_list:
                new_values_sequence = values_sequence[:]
                new_values_sequence.append(value)
                new_values_sequence, result = cnmcts(
                    new_values_sequence,
                    bounds,
                    planets_sequence,
                    level - 1,
                    bandwidth,
                )
                if result < best_result:
                    chosen_move = value
                    best_result = result
            if chosen_move is None:
                chosen_move = np.random.uniform(low_bound, high_bound)
        values_sequence.append(chosen_move)
    return (
        values_sequence,
        evaluate_mga_trajectory(
            planets_sequence, *separate_sequence_values(values_sequence, len(planets_sequence))
        )[0],
    )
