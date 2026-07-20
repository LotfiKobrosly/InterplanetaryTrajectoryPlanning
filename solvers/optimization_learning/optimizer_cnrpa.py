"""
Implements a GcNRPA aiming to learn how to optimize the search for better values
"""

import time
from copy import deepcopy
import numpy as np
import pykep as pk
import matplotlib.pyplot as plt
from utils.gaussian_kernel import GaussianKernel
from utils.constants import RANDOM_GENERATOR, GAUSSIAN_KERNEL_THRESHOLD
from utils.basic_functions import *
from utils.udp_wrapper import CountingEvaluator

SCORE_FUNCTION_LIST = [
    "best_delta_v",
    "weighted_differences_sum",
    "differences_sum",
    "largest_diff",
]

INITIAL_STATE_STRATEGIES = [
    "random",
    "mixed",
    "best_current_value",
]


def score_function(
    sequence: list = None,
    score_type: str = "best_delta_v",
    evaluator: pk.trajopt.mga = None,
    cumsum_weights: np.ndarray = None,
    *args,
    **kwargs,
):
    delta_v_list = [
        evaluator.fitness(denormalize(values_vector, *evaluator.get_bounds()))[0]
        for values_vector in sequence
    ]
    best_vector = sequence[np.argmin(delta_v_list)]
    if score_type == "best_delta_v":
        score = min(delta_v_list)

    elif score_type == "weighted_differences_sum":
        assert (
            cumsum_weights is not None
        ), "cumsum_weights must be speciefied for score_type = weighted_differences_sum"
        differences = [
            delta_v - delta_v_list[index - 1]
            for index, delta_v in enumerate(delta_v_list[1:])
        ]
        score = np.sum(cumsum_weights * np.array(differences))

    elif score_type == "differences_sum":
        score = delta_v_list[-1] - delta_v_list[0]

    elif score_type == "largest_diff":
        score = min(delta_v_list) - max(delta_v_list)

    else:
        raise ValueError("Score type unknown")

    return score, best_vector


def get_initial_state(
    initial_state_strategy: str = "random",
    mixture_probability: float = None,
    best_values_archive: list = None,
    best_values_scores: list = None,
    vector_size: int = None,
    *args,
    **kwargs,
):
    if initial_state_strategy == "mixed":
        assert mixture_probability is not None, "Undefined mixture_probability"
        if RANDOM_GENERATOR.uniform() < mixture_probability:
            initial_state_strategy = "random"
        else:
            initial_state_strategy = "best_current_value"

    if initial_state_strategy == "random":
        assert (
            vector_size is not None
        ), "vector_size unspecified for random initial state_strategy"
        return RANDOM_GENERATOR.uniform(0, 1, size=vector_size)
    elif initial_state_strategy == "best_current_value":
        assert (best_values_archive is not None) and (best_values_scores is not None), (
            "best_values_archive is None for initial state strategy = "
            + initial_state_strategy
        )
        if best_values_archive:
            probabilities = 1 / np.array(
                [delta_v / sum(best_values_scores) for delta_v in best_values_scores]
            )
            probabilities /= np.sum(probabilities)
            chosen_index = RANDOM_GENERATOR.choice(
                list(range(len(best_values_archive))), p=probabilities
            )
            return best_values_archive[chosen_index]
        else:
            return RANDOM_GENERATOR.uniform(0, 1, size=vector_size)
    else:
        raise ValueError("Wrong initial_state_strategy = " + initial_state_strategy)


def policy_playout(
    initial_state: np.ndarray = None,
    policy: dict = None,
    max_steps: int = 100,
    std_factor: float = 0.1,
    movement_range: tuple = (-0.2, 0.2),
    *args,
    **kwargs,
):

    sequence = [initial_state]
    vector_size = len(initial_state)
    action_sequence = list()

    for _ in range(max_steps):
        if policy:
            gaussian_kernel = GaussianKernel(sequence[-1], sigma=std_factor)
            values, weights = list(), list()
            for key in policy.keys():
                value = np.array(policy[key])
                weight = gaussian_kernel.pdf(key)
                if weight >= GAUSSIAN_KERNEL_THRESHOLD:
                    weights.append(weight)
                    values.append(value)
            if values:
                weights = np.array(weights)
                weights /= np.sum(weights)
                try:
                    action = RANDOM_GENERATOR.normal(
                        weights @ np.array(values), std_factor
                    )
                except:
                    print("Weights:", np.array(weights).shape)
                    print("Values:", np.array(values).shape)
                    raise ValueError("Shape mismatch")
            else:
                action = RANDOM_GENERATOR.uniform(*movement_range, size=vector_size)
        else:
            action = RANDOM_GENERATOR.uniform(*movement_range, size=vector_size)
        action_sequence.append(action)
        sequence.append(
            truncate(sequence[-1] + action, [0] * vector_size, [1] * vector_size)
        )

    return sequence, action_sequence


def adapt_policy(
    policy: dict = dict(),
    learning_rate: float = 0.01,
    best_sequence: list = None,
    best_actions: list = None,
    *args,
    **kwargs,
):
    if policy:
        # Adapt visited states
        best_sequence = [code(list(vector)) for vector in best_sequence]
        for vector_index, vector in enumerate(best_sequence[:-1]):
            if vector in policy.keys():
                policy[vector] += learning_rate * (
                    best_actions[vector_index] - policy[vector]
                )

        # Adapt nearby states
        for state, action in policy.items():
            if state not in best_sequence:
                kernel = GaussianKernel(state, 0.1)
                weights, values = list(), list()
                for vector_index, vector in enumerate(best_sequence[:-1]):
                    weight = kernel.pdf(vector)
                    if weight >= GAUSSIAN_KERNEL_THRESHOLD:
                        weights.append(weight)
                        values.append(best_actions[vector_index])
                if weights:
                    try:
                        new_value = np.array(weights) @ np.array(values)
                    except:
                        print("Weights:", np.array(weights).shape)
                        print("Values:", np.array(values).shape)
                        raise ValueError("Shape mismatch")
                    policy[state] = action + learning_rate * (new_value - action)

    else:
        for vector_index, vector in enumerate(best_sequence[:-1]):
            policy[code(list(vector))] = best_actions[vector_index]

    return policy


def run_optimizer_cnrpa(
    level: int = 0,
    n_policies: int = 100,
    learning_rate: float = 0.01,
    evaluator: pk.trajopt.mga = None,
    policy: dict = None,
    initial_state_strategy: str = "random",
    score_type: str = "best_delta_v",
    current_iteration: int = 0,
    best_score: float = None,
    best_sequence: list = None,
    best_actions: list = None,
    best_values_archive: list = None,
    best_values_scores: list = None,
    score_evolution: list = None,
    time_list: list = None,
    archive_size: int = None,
    max_steps: int = 100,
    start_time: float = 0,
    timeout: float = 10,
    cumsum_weights: np.ndarray = None,
    movement_range: tuple = (-0.2, 0.2),
    *args,
    **kwargs,
):
    lower_bounds, upper_bounds = evaluator.get_bounds()
    current_time = time.time() - start_time
    if level == 0:
        mixture_probability = 1 - current_iteration / n_policies
        initial_state = get_initial_state(
            initial_state_strategy=initial_state_strategy,
            mixture_probability=mixture_probability,
            best_values_archive=best_values_archive,
            best_values_scores=best_values_scores,
            vector_size=len(lower_bounds),
        )
        sequence, actions_sequence = policy_playout(
            initial_state=initial_state,
            policy=policy,
            max_steps=max_steps,
            movement_range=movement_range,
            std_factor=0.01 + 1 / np.sqrt(current_iteration + 1),
        )
        score, best_vector = score_function(
            sequence=sequence,
            score_type=score_type,
            evaluator=evaluator,
            cumsum_weights=cumsum_weights,
        )
        return (
            sequence,
            actions_sequence,
            score,
            best_vector,
        )
    else:
        current_policy = deepcopy(policy)
        for current_iteration in range(n_policies):
            sequence, actions_sequence, score, best_vector = run_optimizer_cnrpa(
                level=level - 1,
                n_policies=n_policies,
                policy=current_policy,
                learning_rate=learning_rate,
                evaluator=evaluator,
                initial_state_strategy=initial_state_strategy,
                score_type=score_type,
                best_score=best_score,
                best_sequence=best_sequence,
                best_actions=best_actions,
                current_iteration=current_iteration,
                best_values_archive=best_values_archive,
                best_values_scores=best_values_scores,
                score_evolution=score_evolution,
                time_list=time_list,
                archive_size=archive_size,
                max_steps=max_steps,
                start_time=start_time,
                timeout=timeout,
                cumsum_weights=cumsum_weights,
                movement_range=movement_range,
            )
            current_delta_v = evaluator.fitness(
                denormalize(best_vector, lower_bounds, upper_bounds)
            )[0]

            # Potentially add new found value to the archive
            if best_values_archive:
                insert_sorted(
                    best_values_archive,
                    best_values_scores,
                    best_vector,
                    current_delta_v,
                )
                while len(best_values_archive) > archive_size:
                    best_values_archive.pop(-1)
                    best_values_scores.pop(-1)
            else:
                best_values_archive.append(best_vector)
                best_values_scores.append(current_delta_v)

            # Adapt policy
            if score < best_score:
                best_score = score
                best_sequence = sequence[:]
                best_actions = actions_sequence[:]
            current_policy = adapt_policy(
                policy=current_policy,
                learning_rate=learning_rate,
                best_actions=best_actions,
                best_sequence=best_sequence,
            )
            current_time = time.time() - start_time
            score_evolution.append(best_score)
            time_list.append(current_time)
            if current_time > timeout:
                break

        return (
            best_sequence,
            best_actions,
            best_score,
            best_values_archive[0],
        )


def optimizer_cnrpa(
    level: int = 1,
    n_policies: int = 100,
    learning_rate: float = 0.01,
    evaluator: pk.trajopt.mga = None,
    initial_state_strategy: str = "random",
    score_type: str = "best_delta_v",
    archive_size: int = None,
    max_steps: int = 100,
    movement_range: tuple = (-0.2, 0.2),
    timeout: float = 10,
    *args,
    **kwargs,
):
    assert evaluator is not None, "Undefined evaluator"
    assert level >= 0, "level is negative"
    assert n_policies > 0, "n_policies is not positive"
    cumsum_weights = None
    if score_type == "weighted_differences_sum":
        cumsum_weights = np.array(
            [1 / (max_steps - counter) for counter in range(max_steps)]
        )
        cumsum_weights /= np.sum(cumsum_weights)
    score_evolution, time_list = list(), list()

    start_time = time.time()
    best_sequence, best_actions, best_score, best_vector = run_optimizer_cnrpa(
        level=level,
        n_policies=n_policies,
        policy=dict(),
        learning_rate=learning_rate,
        evaluator=evaluator,
        initial_state_strategy=initial_state_strategy,
        score_type=score_type,
        best_score=np.inf,
        best_sequence=None,
        best_actions=None,
        current_iteration=0,
        best_values_archive=list(),
        best_values_scores=list(),
        score_evolution=score_evolution,
        time_list=time_list,
        archive_size=archive_size,
        max_steps=max_steps,
        start_time=start_time,
        timeout=timeout,
        cumsum_weights=cumsum_weights,
        movement_range=(-movement_range, movement_range),
    )
    total_time = time.time() - start_time
    best_vector = denormalize(best_vector, *evaluator.get_bounds())
    return best_vector, evaluator.fitness(best_vector)[0], score_evolution, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = CountingEvaluator(pk.trajopt.gym.cassini1)

    # Variables bounds
    bounds = [
        (low_bound, high_bound)
        for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    ]

    # General input values
    inputs_values = {
        "evaluator": udp,
        "bounds": bounds,
        "timeout": 180,
        "level": 3,
        "learning_rate": 0.5744087681488131,
        "n_policies": 297,
        "initial_state_strategy": "mixed",
        "score_type": "differences_sum",
        "archive_size": 81,
        "max_steps": 50,
        "movement_range": 0.00064558763425047,
    }
    values_sequence, best_value, scores_list, time_list = optimizer_cnrpa(
        **inputs_values
    )
    print(f"Best Delta V: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
    # print(f"Total number of evaluations: {udp.count}")

    figure = plt.figure(figsize=(10, 10))
    plt.plot(time_list, scores_list)
    plt.title(
        f"Best value {min(scores_list) / 1000:.3f} km/s found first after {time_list[scores_list.index(min(scores_list))]:.3f} s"
    )
    plt.show()

    axe = udp.plot(values_sequence, figsize=(20, 20))
    # figure = axe.figure
    axe.view_init(90, 0)
    axe.axis("off")
    axe.set_title(
        "Optimizing with GcNRPA" + r": $\Delta$V = " + f"{best_value / 1000:.3f} km/s"
    )
    plt.show()
