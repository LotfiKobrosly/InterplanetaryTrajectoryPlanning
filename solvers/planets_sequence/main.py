import time
from copy import deepcopy
import numpy as np
import pykep as pk
from utils.constants import (
    PLANETS,
    RANDOM_GENERATOR,
    DEPARTURE_EPOCH_BOUNDS,
    UNFEASIBILITY_VALUE,
    NEXT_PLANETS,
)
from utils.basic_functions import get_tof_bounds, insert_sorted
from utils.trajectory_validator import trajectory_is_valid
from solvers.continuous_variables_choice import (
    get_variables_values,
)
from solvers.planets_sequence.nrpa import run_mga_nrpa

if __name__ == "__main__":
    # Hyperparameters
    level = 2
    n_policies = 40
    sequence_size = 6
    learning_rate = 0.1
    policy_types_list = ["simple", "last_planet", "sequence_dependent"]

    # Planets
    start_body = "Earth"
    target_body = "Saturn"

    # Continuous solver
    continuous_values_sampling = "sade"
    timeout = 2

    # Best sequences archive
    sequences_archive = list()
    values_vectors_archive = list()
    scores_archive = list()
    max_archive_size = 100

    # Running
    start_time = time.time()

    results = dict()

    for policy_type in policy_types_list:
        print("\n\nRunning for policy type:", policy_type)
        sequence, values_vector, best_dv = run_mga_nrpa(
            level=level,
            n_policies=n_policies,
            policy_type=policy_type,
            start_body=start_body,
            target_body=target_body,
            sequence_size=sequence_size,
            learning_rate=learning_rate,
            continuous_values_sampling=continuous_values_sampling,
            sequences_archive=sequences_archive,
            values_vectors_archive=values_vectors_archive,
            scores_archive=scores_archive,
            max_archive_size=max_archive_size,
            timeout=timeout,
        )

        # Solving time
        total_time = time.time() - start_time
        hours = int(total_time) // 3600
        minutes = (int(total_time) % 3600) // 60
        seconds = total_time - (minutes * 60 + hours * 3600)

        # Travel time
        travel_time = int(np.sum(values_vector[1:]))
        years = travel_time // 365
        days = travel_time - years * 365

        print("Best found sequence:", sequence)
        print(f"Total DV: {best_dv / 1000:.3f} km/s")
        print(f"Total travel time: {years} years and {days} days")
        print("\nNow verifying constraint satisfaction:")

        # Constraint violation verification
        valid_sequence_not_found = True
        for sequence_index, sequence in enumerate(sequences_archive):
            if (
                not (
                    (sequence is None)
                    or (values_vectors_archive[sequence_index] is None)
                )
            ) and trajectory_is_valid(
                planets_sequence=sequence,
                values_vector=values_vectors_archive[sequence_index],
            ):
                print("Best valid sequence:", sequence)
                print(f"Total DV: {scores_archive[sequence_index] / 1000:.3f} km/s")
                travel_time = int(np.sum(values_vectors_archive[sequence_index][1:]))
                years = travel_time // 365
                days = travel_time - years * 365
                print(f"Total travel time: {years} years and {days} days")
                print(f"Solving time: {hours}h, {minutes}min and {seconds:.3f}s")
                results[policy_type] = {
                    "Sequence": sequence,
                    "Travel time": str(years) + " years and " + str(days) + " days",
                    "Total DV": f"{scores_archive[sequence_index] / 1000:.3f} km/s",
                    "Solve time": f"Solving time: {hours}h, {minutes}min and {seconds:.3f}s",
                }
                valid_sequence_not_found = False
                break

    print("\n\nRecap:")
    for policy_type, result in results.items():
        print("\nFor policy type: " + policy_type)
        for key, value in result.items():
            print(key + ": " + str(value))
