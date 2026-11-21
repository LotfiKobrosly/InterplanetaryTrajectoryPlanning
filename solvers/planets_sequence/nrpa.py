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
    NEXT_PLANETS,
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
    policy_type: str = "simple",
    legal_moves: dict = None,
    continuous_values_sampling: str = "sade",
    *args,
    **kwargs,
):
    # Sequence initialization
    planets_sequence = [start_body]

    # Construct sequence
    for planet_index in range(sequence_size - 2):

        ## Initialization
        choices = list()
        probabilities = list()
        last_planet = planets_sequence[-1]
        choices = legal_moves[last_planet]
        ## Possible actions and their probabilities

        if policy_type == "simple":
            current_policy = policy[planet_index]
        else:

            step_policy = policy.get(planet_index, dict())
            if policy_type == "last_planet":
                current_policy = step_policy.get(last_planet, None)
            else:
                current_policy = step_policy.get(tuple(planets_sequence), None)
        if current_policy is None:
            probabilities = np.ones(len(choices))
        else:
            probabilities = [np.exp(current_policy[planet]) for planet in choices]

        ## Normalize and sample
        probabilities = np.array(probabilities)
        probabilities /= np.sum(probabilities)
        try:
            next_planet = RANDOM_GENERATOR.choice(choices, p=probabilities)
        except ValueError:
            print("Last planet: ", last_planet)
            print("Advancement: ", len(planets_sequence))
            print("Choices: ", choices)
            raise ValueError("Empty sequence")

        ## Extend sequence
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


def adapt_policy(
    policy: dict = None,
    policy_type: str = "simple",
    legal_moves: dict = None,
    planets_sequence: list = None,
    learning_rate: float = 0.1,
):
    # Code correcteness
    assert not (policy is None), "Undefined policy"
    assert not (planets_sequence is None), "Empty planets_sequence"
    assert not (learning_rate is None), "Error in argument passing of legal_moves"

    # Adaptation per step
    for id in range(len(planets_sequence[:-2])):
        normalizer = 0
        last_planet = planets_sequence[id]

        # Verify policy type
        if policy_type == "simple":
            current_policy = policy[id]
        else:
            step_policy = policy.get(id, dict())
            policy[id] = step_policy
            if policy_type == "last_planet":
                key = planets_sequence[id]
            else:
                key = tuple(planets_sequence[: id + 1])
            current_policy = step_policy.get(key, None)

            if current_policy is None:
                current_policy = {planet: 0 for planet in legal_moves[last_planet]}
            step_policy[key] = current_policy

        # Adapt policy
        current_policy[planets_sequence[id + 1]] += learning_rate
        for weight in current_policy.values():
            normalizer += np.exp(weight)
        for planet, weight in current_policy.items():
            current_policy[planet] -= learning_rate * np.exp(weight) / normalizer


def mga_nrpa(
    level: int = 1,
    n_policies: int = 100,
    policy: list = None,
    policy_type: str = "simple",
    start_body: str = "Earth",
    target_body: str = "Saturn",
    sequence_size: int = 6,
    legal_moves: dict = None,
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
            legal_moves=legal_moves,
            policy=policy,
            policy_type=policy_type,
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
                policy_type=policy_type,
                start_body=start_body,
                target_body=target_body,
                sequence_size=sequence_size,
                legal_moves=legal_moves,
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
            insert_sorted(
                values_vectors_archive, scores_archive, values_vector, total_dv
            )
            while len(sequences_archive) > max_archive_size:
                sequences_archive.pop(-1)
                values_vectors_archive.pop(-1)

            # Update lower level's best sequence
            if total_dv < local_best_dv:
                local_best_dv = total_dv
                local_best_sequence = planets_sequence
                local_best_values_vector = values_vector
            if local_best_dv < UNFEASIBILITY_VALUE:
                adapt_policy(
                    policy=new_policy,
                    policy_type=policy_type,
                    legal_moves=legal_moves,
                    planets_sequence=local_best_sequence,
                    learning_rate=learning_rate,
                )

        # Update higher level's best sequence
        if local_best_dv < global_best_dv:
            global_best_dv = local_best_dv
            global_best_sequence = local_best_sequence
            global_best_values_vector = local_best_values_vector

        # Only adapt if no constraints are violated
        if global_best_dv < UNFEASIBILITY_VALUE:
            adapt_policy(
                policy=policy,
                policy_type=policy_type,
                planets_sequence=global_best_sequence,
                legal_moves=legal_moves,
                learning_rate=learning_rate,
            )

        return global_best_sequence, global_best_values_vector, global_best_dv


def run_mga_nrpa(
    level: int = 1,
    n_policies: int = 100,
    start_body: str = "Earth",
    target_body: str = "Saturn",
    sequence_size: int = 6,
    policy_type: str = "simple",
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
    # assert start_body != target_body, "Target and start bodies are identical"
    assert (
        sequence_size > 2
    ), "Sequence size <= 2: there are no intermediary flybys, which is a single Lambert Problem, not MGA"

    # Set planets pool
    planets_pool = PLANETS[: max(PLANETS.index(start_body), PLANETS.index(target_body))]
    legal_moves = {
        planet: list(
            (set(planets_pool).intersection(set(NEXT_PLANETS[planet]))).difference(
                {target_body}
            )
        )
        for planet in planets_pool
    }

    # Initialize policy
    if policy_type == "simple":
        policy = {
            step: {planet: 0 for planet in planets_pool}
            for step in range(sequence_size - 2)
        }
    elif policy_type == "last_planet":
        policy = {0: {start_body: {planet: 0 for planet in legal_moves[start_body]}}}
    elif policy_type == "sequence_dependent":
        policy = {0: {(start_body,): {planet: 0 for planet in legal_moves[start_body]}}}
    else:
        raise ValueError("Unrecognized policy type: " + policy_type)

    # Run NRPA
    return mga_nrpa(
        policy=policy,
        level=level,
        n_policies=n_policies,
        policy_type=policy_type,
        legal_moves=legal_moves,
        start_body=start_body,
        target_body=target_body,
        sequence_size=sequence_size,
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
    n_policies = 10
    sequence_size = 5
    learning_rate = 0.1
    policy_type = "last_planet"  # "simple", "last_planet", or "sequence_dependent"

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
    print(f"Solving time: {hours}h, {minutes}min and {seconds:.3f}s")

    # Constraint violation verification
    valid_sequence_not_found = True
    for sequence_index, sequence in enumerate(sequences_archive):
        if (
            not ((sequence is None) or (values_vectors_archive[sequence_index] is None))
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
            valid_sequence_not_found = False
            break

    if valid_sequence_not_found:
        print("No valid values vector found")
