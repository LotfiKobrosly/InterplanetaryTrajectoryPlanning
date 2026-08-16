"""
Useful functions for Optimization Learning algorithms
"""

import numpy as np
import pykep as pk
from utils.constants import RANDOM_GENERATOR
from utils.basic_functions import denormalize

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
