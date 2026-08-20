"""
Implements a cNMCTS aiming to learn how to optimize the search for better values
"""

import time
from copy import deepcopy
import numpy as np
import pykep as pk
import matplotlib.pyplot as plt
from utils.constants import RANDOM_GENERATOR
from utils.basic_functions import *
from solvers.continuous_variables_choice.baselines import pygmo_baseline
from solvers.optimization_learning.general import score_function, get_initial_state


def run_optimizer_cnmcts(
    evaluator: pk.trajopt.mga = None,
    best_results: dict = None,
    level: int = 1,
    bandwidth: int = 20,
    sequence: list = None,
    score_type: str = "best_delta_v",
    max_steps: int = 100,
    movement_range: float = 0.2,
    timeout: float = 10,
    start_time: float = 0,
    cumsum_weights: np.ndarray = None,
    score_evolution: list = None,
    time_list: list = None,
    *args,
    **kwargs,
):
    current_time = time.time() - start_time
    vector_size = len(sequence[0])
    while len(sequence) <= max_steps:
        if level == 0:
            sequence.append(
                truncate(
                    sequence[-1]
                    + RANDOM_GENERATOR.uniform(*movement_range, size=vector_size),
                    [0] * vector_size,
                    [1] * vector_size,
                )
            )

        else:
            actions_pool = [
                RANDOM_GENERATOR.uniform(*movement_range, size=vector_size)
                for _ in range(bandwidth)
            ]
            for action in actions_pool:
                temporary_sequence = deepcopy(sequence)
                temporary_sequence.append(
                    truncate(
                        sequence[-1]
                        + RANDOM_GENERATOR.uniform(*movement_range, size=vector_size),
                        [0] * vector_size,
                        [1] * vector_size,
                    )
                )
                temporary_sequence = run_optimizer_cnmcts(
                    sequence=temporary_sequence,
                    evaluator=evaluator,
                    level=level - 1,
                    bandwidth=bandwidth,
                    best_results=best_results,
                    score_type=score_type,
                    max_steps=max_steps,
                    movement_range=movement_range,
                    timeout=timeout,
                    start_time=start_time,
                    cumsum_weights=cumsum_weights,
                    score_evolution=score_evolution,
                    time_list=time_list,
                )

                current_time = time.time() - start_time
                if current_time > timeout:
                    break
                score, best_current_vector = score_function(
                    sequence=temporary_sequence,
                    score_type=score_type,
                    evaluator=evaluator,
                    cumsum_weights=cumsum_weights,
                )
                best_current_vector_score = evaluator.fitness(
                    denormalize(best_current_vector, *evaluator.get_bounds())
                )[0]
                if score < best_results["best_sequence_score"]:
                    best_results["best_sequence_score"] = score
                    best_results["best_sequence"] = deepcopy(temporary_sequence)
                if best_current_vector_score < best_results["best_vector_score"]:
                    best_results["best_vector_score"] = best_current_vector_score
                    best_results["best_vector"] = deepcopy(best_current_vector)
            sequence.append(best_results["best_sequence"][len(sequence)])
            current_time = time.time() - start_time
            score_evolution.append(best_results["best_sequence_score"])
            time_list.append(current_time)
            if current_time > timeout:
                break
    return sequence


def optimizer_cnmcts(
    level: int = 1,
    bandwidth: int = 20,
    evaluator: pk.trajopt.mga = None,
    initial_state_strategy: str = "random",
    score_type: str = "best_delta_v",
    max_steps: int = 100,
    movement_magnitude: float = 0.2,
    warm_starter: str = None,
    timeout: float = 10,
    *args,
    **kwargs,
):
    assert evaluator is not None, "Undefined evaluator"
    assert (isinstance(bandwidth, int)) and (
        bandwidth > 0
    ), "bandwidth must be a positive integer"
    cumsum_weights = None
    if score_type == "weighted_differences_sum":
        cumsum_weights = np.array(
            [1 / (max_steps - counter) for counter in range(max_steps)]
        )
        cumsum_weights /= np.sum(cumsum_weights)
    score_evolution, time_list = list(), list()
    lower_bounds, upper_bounds = evaluator.get_bounds()
    if initial_state_strategy == "random":
        sequence = [RANDOM_GENERATOR.uniform(0, 1, size=len(lower_bounds))]
    elif initial_state_strategy == "middle":
        sequence = [[0.5] * len(lower_bounds)]
    elif initial_state_strategy == "warm_start":
        assert (warm_starter is not None), "Unspecified warm starting algorithm"
        warm_start_parameters = {
            "evaluator": deepcopy(evaluator),
            "solver": warm_starter,
            "timeout": 10,
            "population_size": 50,
        }
        start_vector, _, _, _ = pygmo_baseline(**warm_start_parameters)
        sequence = [normalize(start_vector, *evaluator.get_bounds())]
    else:
        raise ValueError("Unknown initial state")
    start_time = time.time()
    best_results = {
        "best_sequence": None,
        "best_sequence_score": np.inf,
        "best_vector": None,
        "best_vector_score": np.inf,
    }
    sequence = run_optimizer_cnmcts(
        evaluator=evaluator,
        sequence=sequence,
        best_results=best_results,
        level=level,
        bandwidth=bandwidth,
        score_type=score_type,
        max_steps=max_steps,
        movement_range=(-movement_magnitude, movement_magnitude),
        timeout=timeout,
        start_time=start_time,
        cumsum_weights=cumsum_weights,
        score_evolution=score_evolution,
        time_list=time_list,
    )
    best_vector = denormalize(best_results["best_vector"], *evaluator.get_bounds())
    return best_vector, evaluator.fitness(best_vector)[0], score_evolution, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = pk.trajopt.gym.cassini1

    # General input values
    inputs_values = {
        "evaluator": udp,
        "timeout": 300,
        "level": 2,
        "bandwidth": 400,
        "initial_state_strategy": "warm_start",
        "score_type": "weighted_differences_sum",
        "max_steps": 50,
        "movement_range": 0.001,
        "warm_starter": "pso",
    }
    values_sequence, best_value, scores_list, time_list = optimizer_cnmcts(
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
        "Optimizing with cNMCTS" + r": $\Delta$V = " + f"{best_value / 1000:.3f} km/s"
    )
    plt.show()
