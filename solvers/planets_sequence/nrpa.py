"""
Implement a Nested Rollout Policy Adaptation for the planets' sequence setting.
In this specific setting, we fix the size of the planets sequence.
"""

import time
from copy import deepcopy
import numpy as np
import pykep as pk
from utils.constants import (
    PLANETS,
    RANDOM_GENERATOR,
    DEPARTURE_EPOCH_BOUNDS,
    UNFEASIBILITY_VALUE,
)
from utils.basic_functions import get_tof_bounds, insert_sorted
from utils.trajectory_validator import trajectory_is_valid
from solvers.continuous_variables_choice import (
    get_variables_values,
)

def policy_rollout(
    start_body: str = "Earth",
    target_body: str = "Saturn",
    sequence_size: int = 6,
    policy: dict = dict(),
    continuous_values_sampling: str = "sade",
    *args,
    **kwargs,
):
    # Sequence initialization
    planets_sequence = [start_body]

    # Construct sequence
    for planet_index in range(sequence_size - 2):
        choices = list()
        probabilities = list()
        for planet, weight in policy[planet_index].items():
            choices.append(planet)
            probabilities.append(np.exp(weight))
        probabilities = np.array(probabilities)
        probabilities /= np.sum(probabilities)
        next_planet = RANDOM_GENERATOR.choice(choices, p=probabilities)
        planets_sequence.append(next_planet)

    planets_sequence.append(target_body)

    # Set TOF bounds
    tof_bounds = list()
    for planet_index, planet in enumerate(planets_sequence[:-1]):
        tof_bounds.append(get_tof_bounds(planet, planets_sequence[planet_index + 1]))

    # Sequence evaluation
    trajectory = pk.trajopt.mga(
        seq=[pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence],
        t0=DEPARTURE_EPOCH_BOUNDS,
        tof=tof_bounds,
        vinf=3000,
    )
    try:
        values_vector, total_dv, _, _ = get_variables_values(
            {
                "sampling_function": continuous_values_sampling,
                "timeout": kwargs.get("timeout", 5),
                "evaluator": trajectory,
                "population_size": kwargs.get("population_size", 100),
            }
        )
    except ValueError:
        values_vector, total_dv = None, UNFEASIBILITY_VALUE

    # Results
    return planets_sequence, values_vector, total_dv


def adapt_policy(policy: list, planets_sequence: list, learning_rate: float):
    for id in range(len(planets_sequence[:-2])):
        normalizer = 0
        policy[id][planets_sequence[id + 1]] += learning_rate
        for weight in policy[id].values():
            normalizer += np.exp(weight)
        for planet, weight in policy[id].items():
            policy[id][planet] -= learning_rate * np.exp(weight) / normalizer

def mga_nrpa(
    level: int = 1,
    n_policies: int = 100,
    policy: list = None,
    start_body: str = "Earth",
    target_body: str = "Saturn",
    sequence_size: int = 6,
    planets_pool: list = PLANETS,
    learning_rate: float = 0.1,
    global_best_sequence: list = None,
    global_best_values_vector: list = None,
    global_best_dv: float = np.inf,
    continuous_values_sampling: str = "sade",
    sequences_archive: list = None,
    values_vectors_archive: list = None,
    scores_archive: list = None,
    max_archive_size: int = 100,
    *args,
    **kwargs,
):
    assert not (policy is None), "policy needs to be initialized at 0"
    if level == 0:
        return policy_rollout(
            start_body=start_body,
            target_body=target_body,
            sequence_size=sequence_size,
            policy=policy,
            continuous_values_sampling=continuous_values_sampling,
            timeout=5,
        )
    else:
        local_best_dv = np.inf
        local_best_sequence = None
        local_best_values_vector = None
        for _ in range(n_policies):
            new_policy = deepcopy(policy)
            planets_sequence, values_vector, total_dv = mga_nrpa(
                policy=new_policy,
                level=level - 1,
                n_policies=n_policies,
                start_body=start_body,
                target_body=target_body,
                sequence_size=sequence_size,
                planets_pool=planets_pool,
                learning_rate=learning_rate,
                global_best_sequence=local_best_sequence,
                global_best_dv=local_best_dv,
                continuous_values_sampling=continuous_values_sampling,
                sequences_archive=sequences_archive,
                values_vectors_archive=values_vectors_archive,
                scores_archive=scores_archive,
                max_archive_size=max_archive_size,
                args=args,
                kwargs=kwargs,
            )

            # Update archives
            insert_sorted(sequences_archive, scores_archive, planets_sequence, total_dv)
            insert_sorted(values_vectors_archive, scores_archive, values_vector, total_dv)
            while len(sequences_archive) > max_archive_size:
                sequences_archive.pop(-1)
                values_vectors_archive.pop(-1)

            # Update lower level's best sequence
            if total_dv < local_best_dv:
                local_best_dv = total_dv
                local_best_sequence = planets_sequence
                local_best_values_vector = values_vector
            if local_best_dv < UNFEASIBILITY_VALUE:
                adapt_policy(new_policy, local_best_sequence, learning_rate)

        # Update higher level's best sequence 
        if local_best_dv < global_best_dv:
            global_best_dv = local_best_dv
            global_best_sequence = local_best_sequence
            global_best_values_vector = local_best_values_vector

        # Only adapt if no constraints are violated
        if global_best_dv < UNFEASIBILITY_VALUE:
            adapt_policy(policy, global_best_sequence, learning_rate)

        return global_best_sequence, global_best_values_vector, global_best_dv

def run_mga_nrpa(
    level: int = 1,
    n_policies: int = 100,
    start_body: str = "Earth",
    target_body: str = "Saturn",
    sequence_size: int = 6,
    learning_rate: float = 0.1,
    continuous_values_sampling: str = "sade",
    sequences_archive: list = None,
    values_vectors_archive: list = None,
    scores_archive: list = None,
    max_archive_size: int = 100,
    *args,
    **kwargs,
):
    # Verifications
    assert start_body != target_body, "Target and start bodies are identical"
    assert sequence_size > 2, "Sequence size <= 2: there are no intermediary flybys, which is a single Lambert Problem, not MGA"

    # Set planets pool
    planets_pool = PLANETS[:max(PLANETS.index(start_body), PLANETS.index(target_body))]

    # Initialize policy
    policy = [
        {planet: 0 for planet in planets_pool}
        for _ in range(sequence_size - 2)
    ]

    # Run NRPA
    return mga_nrpa(
        policy=policy,
        level=level,
        n_policies=n_policies,
        start_body=start_body,
        target_body=target_body,
        sequence_size=sequence_size,
        planets_pool=planets_pool,
        learning_rate=learning_rate,
        global_best_sequence=None,
        global_best_dv=np.inf,
        continuous_values_sampling=continuous_values_sampling,
        sequences_archive=sequences_archive,
        values_vectors_archive=values_vectors_archive,
        scores_archive=scores_archive,
        max_archive_size=max_archive_size,
        args=args,
        kwargs=kwargs,
    )

if __name__ == "__main__":

    # Hyperparameters
    level = 1
    n_policies = 5
    sequence_size = 6
    learning_rate = 0.1

    # Planets
    start_body = "Earth"
    target_body = "Saturn"

    # Continuous solver
    continuous_values_sampling = "sade"
    timeout = 5

    # Best sequences archive
    sequences_archive = list()
    values_vectors_archive = list()
    scores_archive = list()
    max_archive_size = 100

    # Running
    start_time = time.time()
    sequence, values_vector, best_dv = run_mga_nrpa(
        level=level,
        n_policies=n_policies,
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
    total_time = time.time() - start_time
    minutes = int(total_time) // 60
    seconds = total_time - minutes * 60

    print("Best found sequence:", sequence)
    print(f"Total DV: {best_dv / 1000:.3f} km/s")
    print(f"Total travel time: {int(np.sum(values_vector[1:]))} days")
    print(f"Solving time: {minutes}:{seconds:.3f}s")

    # Constraint violation verification
    valid_sequence_not_found = True
    for sequence_index, sequence in enumerate(sequences_archive):
        if trajectory_is_valid(planets_sequence=sequence, values_vector=values_vectors_archive[sequence_index]):
            print("Best valid sequence:", sequence)
            print(f"Total DV: {scores_archive[sequence_index] / 1000:.3f} km/s")
            print(f"Total travel time: {int(np.sum(values_vectors_archive[1:]))} days")
            valid_sequence_not_found = False
            break

    if valid_sequence_not_found:
        print("No valid values vector found")
