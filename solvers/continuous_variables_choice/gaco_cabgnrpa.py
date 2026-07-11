"""
Implements a modified version of a GcABGNRPA, aided with GACO (see pygmo documentation)
"""

import time
from copy import deepcopy
import numpy as np
import pykep as pk
import pygmo as pg
from utils.constants import (
    GAUSSIAN_KERNEL_THRESHOLD,
    RANDOM_SEED,
    RANDOM_GENERATOR,
    VELOCITY_NORMALIZING_FACTOR,
    DV_LAUNCHER,
    UNFEASIBILITY_VALUE,
    N_CANDIDATES,
)
from utils.basic_functions import *
from utils.gaussian_kernel import GaussianKernel
from solvers.continuous_variables_choice.cnrpa import adapt_policy


def gaco_cabgnrpa_playout(
    policy: dict = dict(),
    bias_value: np.ndarray = None,
    bias_std: float = 1.0,
    bounds: list = None,
    states_sequence: list = list(),
    std_factor: float = 1,
    tau: float = 10,
):
    values_sequence, states_sequence = list(), list()
    for advancement, bound in enumerate(bounds):
        # Get bounds
        low_bound, high_bound = bound

        # Choose departure epoch
        if advancement == 0:
            # No state defined
            if policy:
                chosen_value = truncate(
                    sample_mixture_1d(
                        n_samples=1,
                        mu1=denormalize(
                            RANDOM_GENERATOR.normal(policy[0], std_factor),
                            low_bound,
                            high_bound,
                        ),
                        sigma1=std_factor * (high_bound - low_bound),
                        mu2=bias_value[0],
                        sigma2=bias_std[0],
                        weight1=1 / tau,
                    )[0],
                    low_bound,
                    high_bound,
                )
            else:
                chosen_value = RANDOM_GENERATOR.uniform(low_bound, high_bound)
            states_sequence.append(normalize(chosen_value, low_bound, high_bound))

        else:
            # Policy candidate
            if policy:
                current_policy = policy[advancement]
                gaussian_kernel = GaussianKernel(states_sequence, sigma=std_factor)
                values, weights = list(), list()
                for key in current_policy.keys():
                    if isinstance(current_policy[key], (float, int)):
                        value = current_policy[key]
                    else:
                        value = np.array(current_policy[key])
                    weight = gaussian_kernel.pdf(value)
                    if weight >= GAUSSIAN_KERNEL_THRESHOLD:
                        weights.append(weight)
                        values.append(value)
                if values:
                    weights = np.array(weights)
                    weights /= np.sum(weights)
                    policy_weights_candidate = denormalize(
                        RANDOM_GENERATOR.normal(
                            weights @ np.array(values).T, std_factor
                        ),
                        low_bound,
                        high_bound,
                    )
                else:
                    policy_weights_candidate = RANDOM_GENERATOR.uniform(
                        low_bound, high_bound
                    )
            else:
                policy_weights_candidate = RANDOM_GENERATOR.uniform(
                    low_bound, high_bound
                )

            # Choose value
            chosen_value = truncate(
                sample_mixture_1d(
                    n_samples=1,
                    mu1=policy_weights_candidate,
                    sigma1=std_factor * (high_bound - low_bound) / 2,
                    mu2=bias_value[advancement],
                    sigma2=bias_std[advancement],
                    weight1=1 / tau,
                )[0],
                low_bound,
                high_bound,
            )
            states_sequence.append(
                normalize(
                    chosen_value,
                    low_bound,
                    high_bound,
                )
            )

        values_sequence.append(chosen_value)
    return values_sequence, states_sequence


def run_gaco_cabgnrpa(
    evaluator: pk.trajopt.mga,
    policy: dict = dict(),
    biases_values: list = None,
    bias_handler: pg.algorithm = None,
    level: int = 0,
    n_policies: int = 10,
    zeta: float = 0.2,
    bounds: list = None,
    current_iteration: int = 0,
    learning_rate: float = 0.01,
    timeout: float = 10,
    start_time: float = 0,
    best_values_sequence: list = None,
    best_states_sequence: list = None,
    best_value: float = UNFEASIBILITY_VALUE,
    best_values_list: list = None,
    time_list: list = None,
    tau: float = 10,
    *args,
    **kwargs,
):
    assert not (bounds is None), "bounds is None"
    assert not (bias_handler is None), "No algorithm specified or given as argument"
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
        for dimension in range(len(bounds)):

            bias_center, bias_sigma = fit_gaussian_from_density(
                bias_values[:, dimension],
                density_values,
                zeta,
            )
            bias.append(bias_center)
            bias_std.append(bias_sigma)
        values_sequence, states_sequence = gaco_cabgnrpa_playout(
            policy=policy,
            bias_value=bias,
            bias_std=bias_std,
            bounds=bounds,
            std_factor=0.01 + 1 / np.sqrt(current_iteration + 1),
            tau=tau,
        )
        return (
            values_sequence,
            states_sequence,
            evaluator.fitness(values_sequence)[0],
        )
    else:
        # Save current policy
        current_policy = deepcopy(policy)
        current_biases_values = deepcopy(biases_values)

        for current_iteration in range(n_policies):
            values_sequence, states_sequence, total_delta_v = run_gaco_cabgnrpa(
                evaluator=evaluator,
                policy=current_policy,
                biases_values=current_biases_values,
                bias_handler=bias_handler,
                level=level - 1,
                n_policies=n_policies,
                zeta=zeta,
                bounds=bounds,
                current_iteration=current_iteration,
                timeout=timeout,
                start_time=start_time,
                best_values_sequence=best_values_sequence,
                best_states_sequence=best_states_sequence,
                best_value=best_value,
                best_values_list=best_values_list,
                time_list=time_list,
                tau=tau,
            )
            if total_delta_v < best_value:
                best_value = total_delta_v
                best_values_sequence = values_sequence[:]
                best_states_sequence = states_sequence[:]
            current_time = time.time() - start_time
            if best_value < UNFEASIBILITY_VALUE:
                current_policy = adapt_policy(
                    best_values_sequence=best_values_sequence,
                    best_states_sequence=best_states_sequence,
                    policy=current_policy,
                    learning_rate=learning_rate,
                    bounds=bounds,
                )
                best_values_list.append(best_value)
                time_list.append(current_time)
            if current_time > timeout:
                break
        return (
            best_values_sequence,
            best_states_sequence,
            evaluator.fitness(best_values_sequence)[0],
        )


def gaco_cabgnrpa(
    evaluator: pk.trajopt.mga,
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    learning_rate: float = 0.01,
    timeout: float = 10,
    zeta: float = 0.2,
    kernel_size: int = 63,
    n_generations: int = 10,
    elitism_factor: float = 1.0,
    tau: float = 10,
    *args,
    **kwargs,
):
    # Define GACO bias generator
    problem = pg.problem(evaluator)
    bias_handler = pg.algorithm(
        pg.gaco(
            gen=n_generations,
            ker=kernel_size,
            q=elitism_factor,
            seed=RANDOM_SEED,
            memory=True,
        )
    )
    biases_values = pg.population(problem, size=kernel_size, seed=RANDOM_SEED)

    # Launch stopwatch
    start_time = time.time()
    best_values_list, time_list = list(), list()

    # Launch solver
    best_values_sequence, best_states_sequence, best_value = run_gaco_cabgnrpa(
        evaluator=evaluator,
        policy=dict(),
        biases_values=biases_values,
        bias_handler=bias_handler,
        level=level,
        zeta=zeta,
        n_policies=n_policies,
        bounds=bounds,
        current_iteration=0,
        learning_rate=learning_rate,
        timeout=timeout,
        start_time=start_time,
        best_values_sequence=None,
        best_states_sequence=None,
        best_value=UNFEASIBILITY_VALUE,
        best_values_list=best_values_list,
        time_list=time_list,
        tau=tau,
    )
    return best_values_sequence, best_value, best_values_list, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = pk.trajopt.gym.cassini1

    # Variables bounds
    bounds = [
        (low_bound, high_bound)
        for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    ]

    # General input values
    inputs_values = {
        "evaluator": udp,
        "bounds": bounds,
        "timeout": 60,
        "level": 2,
        "learning_rate": 0.25,
        "n_policies": 100,
        "tau": 1.3,
        "zeta": 0.2,
        "kernel_size": 63,
        "n_generations": 10,
        "elitism_factor": 0.2,
    }
    values__sequence, best_value, values_list, time_list = gaco_cabgnrpa(
        **inputs_values
    )
    print(f"Best Delta V: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
