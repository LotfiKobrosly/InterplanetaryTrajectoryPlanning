"""
Implements a GcNRPA aiming to learn how to optimize the search for better values
"""

import time
from copy import deepcopy
import numpy as np
import pykep as pk
import matplotlib.pyplot as plt
from utils.gaussian_kernel import GaussianKernel
from utils.constants import RANDOM_GENERATOR, RANDOM_SEED, GAUSSIAN_KERNEL_THRESHOLD
from utils.basic_functions import *
from utils.udp_wrapper import CountingEvaluator
from solvers.continuous_variables_choice.baselines import *
from solvers.optimization_learning.optimizer_cnrpa import adapt_policy
from solvers.optimization_learning.general import (
    score_function,
    get_initial_state,
)


def biased_policy_playout(
    initial_state: np.ndarray = None,
    policy: dict = None,
    max_steps: int = 100,
    std_factor: float = 0.1,
    bias_value: np.ndarray = None,
    bias_std: float = 1.0,
    movement_range: tuple = (-0.2, 0.2),
    tau: float = 1.5,
    *args,
    **kwargs,
):

    sequence = [initial_state]
    vector_size = len(initial_state)
    action_sequence = list()

    for _ in range(max_steps):

        # Policy component
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

        # Bias component
        bias_direction = (bias_value - sequence[-1]) / movement_range[1]
        # print("Shape:", np.shape(action))

        # Mixing policy and bias
        action_sequence.append(
            sample_mixture_1d(
                n_samples=vector_size,
                mu1=action,
                sigma1=std_factor,
                mu2=bias_direction,
                sigma2=bias_std,
                weight1=1 / tau,
            )[0]
        )

        sequence.append(
            truncate(
                sequence[-1] + action_sequence[-1], [0] * vector_size, [1] * vector_size
            )
        )

    return sequence, action_sequence


def run_optimizer_cabgnrpa(
    level: int = 0,
    n_policies: int = 100,
    learning_rate: float = 0.01,
    evaluator: pk.trajopt.mga = None,
    policy: dict = None,
    biases_values: list = None,
    bias_handler: pg.algorithm = None,
    zeta: float = 0.2,
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
        # GACO iteration
        biases_values = bias_handler.evolve(biases_values)
        bias_values = biases_values.get_x().copy()
        bias_fitness = biases_values.get_f().copy()

        # Get bias value
        density_values = 1 / np.array(bias_fitness).flatten()
        density_values /= density_values.sum()
        bias = list()
        bias_std = list()
        for dimension in range(len(lower_bounds)):
            values = bias_values[:, dimension].reshape(-1)
            values = [
                normalize(value, lower_bounds[dimension], upper_bounds[dimension])
                for value in values
            ]

            bias_center, bias_sigma = fit_gaussian_from_density(
                values,
                density_values,
                zeta,
            )
            bias.append(bias_center)
            bias_std.append(bias_sigma)
            if bias_sigma < 0:
                print("Here:", bias_sigma)

        # Get initial state
        ## For mixed initial state
        mixture_probability = 1 - current_iteration / n_policies
        initial_state = get_initial_state(
            initial_state_strategy=initial_state_strategy,
            mixture_probability=mixture_probability,
            best_values_archive=best_values_archive,
            best_values_scores=best_values_scores,
            vector_size=len(lower_bounds),
        )

        # Playout
        sequence, actions_sequence = biased_policy_playout(
            initial_state=initial_state,
            policy=policy,
            max_steps=max_steps,
            movement_range=movement_range,
            bias_value=bias,
            bias_std=bias_std,
            std_factor=0.01 + 1 / np.sqrt(current_iteration + 1),
        )

        # Compute score
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
        current_biases_values = deepcopy(biases_values)
        current_bias_handler = deepcopy(bias_handler)
        for current_iteration in range(n_policies):
            sequence, actions_sequence, score, best_vector = run_optimizer_cabgnrpa(
                level=level - 1,
                n_policies=n_policies,
                policy=current_policy,
                learning_rate=learning_rate,
                evaluator=evaluator,
                initial_state_strategy=initial_state_strategy,
                score_type=score_type,
                biases_values=current_biases_values,
                bias_handler=current_bias_handler,
                zeta=zeta,
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


def optimizer_cabgnrpa(
    level: int = 1,
    n_policies: int = 100,
    learning_rate: float = 0.01,
    evaluator: pk.trajopt.mga = None,
    initial_state_strategy: str = "random",
    score_type: str = "best_delta_v",
    archive_size: int = None,
    max_steps: int = 100,
    movement_range: float = 0.2,
    timeout: float = 10,
    solver_parameters: dict = None,
    solver: str = "gaco",
    tau: float = 1.5,
    zeta: float = 0.2,
    *args,
    **kwargs,
):
    assert evaluator is not None, "Undefined evaluator"
    assert level >= 0, "level is negative"
    assert n_policies > 0, "n_policies is not positive"

    # If weighted cumulative sum
    cumsum_weights = None
    if score_type == "weighted_differences_sum":
        cumsum_weights = np.array(
            [1 / (max_steps - counter) for counter in range(max_steps)]
        )
        cumsum_weights /= np.sum(cumsum_weights)
    score_evolution, time_list = list(), list()

    # Define GACO bias generator
    solver = PYGMO_SOLVERS[solver]["function"]

    # kernel_size = solver_parameters["kernel_size"]
    # n_generations = solver_parameters["n_generations"]
    # elitism_factor = solver_parameters["elitism_factor"]
    problem = pg.problem(evaluator)
    bias_handler = pg.algorithm(solver(**solver_parameters))
    biases_values = pg.population(problem, size=archive_size, seed=RANDOM_SEED)

    start_time = time.time()
    best_sequence, best_actions, best_score, best_vector = run_optimizer_cabgnrpa(
        level=level,
        n_policies=n_policies,
        policy=dict(),
        biases_values=biases_values,
        bias_handler=bias_handler,
        zeta=zeta,
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
        tau=tau,
        cumsum_weights=cumsum_weights,
        movement_range=(-movement_range, movement_range),
    )
    total_time = time.time() - start_time
    best_vector = denormalize(best_vector, *evaluator.get_bounds())
    return best_vector, evaluator.fitness(best_vector)[0], score_evolution, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = CountingEvaluator(pk.trajopt.gym.cassini1)

    # # Variables bounds
    # bounds = [
    #     (low_bound, high_bound)
    #     for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    # ]

    # General input values
    inputs_values = {
        "evaluator": udp,
        "timeout": 300,
        "level": 2,
        "tau": 1.5,
        "learning_rate": 0.05,
        "n_policies": 10000,
        "initial_state_strategy": "mixed",
        "score_type": "best_delta_v",
        "archive_size": 100,
        "max_steps": 1,
        "zeta": 0.1,
        "movement_range": 0.001,
        "solver": "cmaes",
        "solver_parameters": {
            "gen": 10,
            "force_bounds": True,
            "sigma0": 0.5,
            "ftol": 1e-4,
            "memory": True,
        },
    }
    values_sequence, best_value, scores_list, time_list = optimizer_cabgnrpa(
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
        "Optimizing with GcABGNRPA"
        + r": $\Delta$V = "
        + f"{best_value / 1000:.3f} km/s"
    )
    plt.show()
