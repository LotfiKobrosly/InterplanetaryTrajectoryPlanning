"""
Implements Continuous Generalized Nested Rollout Policy for continuous variables' values choice
Here, the values sequence is represented as follows:
values_sequence = [departure_epoch, time_of_flight_0, ..., time_of_flight_n]

We would do the lambert leg values computation as well as the gravity assist impulses as a way to get heuristic values for GNRPA:
- At each step for action selection (except first step), generate candidates from the policy and store their weights
- Compute for each candidate TOF the Lambert leg and the flyby result
- Use these results to generate the biases values
- Add the biases to the candidates' weights and do the selection

Policy adaptation will then occur in the same manner as in cNRPA
"""

import time
from copy import deepcopy
import numpy as np
import pykep as pk
from utils.constants import (
    GAUSSIAN_KERNEL_THRESHOLD,
    RANDOM_GENERATOR,
    VELOCITY_NORMALIZING_FACTOR,
    DV_LAUNCHER,
    UNFEASIBILITY_VALUE,
    N_CANDIDATES,
)
from utils.basic_functions import *
from utils.gaussian_kernel import (
    GaussianKernel,
)
from utils.heuristic_functions import porkchop_scan
from utils.udp_wrapper import CountingEvaluator
from solvers.continuous_variables_choice.cnrpa import adapt_policy


def biased_policy_playout(
    policy: dict,
    bounds: list,
    planets_sequence: list = None,
    std_factor: float = 1,
    tau: float = 10,
):
    values_sequence, states_sequence = list(), list()
    epoch_list, planets_radii_list, planets_velocities_list = list(), list(), list()
    last_arrival_velocity = None
    for advancement, planet in enumerate(planets_sequence):
        # Get bounds
        low_bound, high_bound = bounds[advancement]

        # Choose departure epoch
        if advancement == 0:
            # No biases at this stage
            if policy:
                chosen_value = round(
                    denormalize(
                        RANDOM_GENERATOR.normal(policy[0], std_factor),
                        low_bound,
                        high_bound,
                    )
                )
            else:
                chosen_value = RANDOM_GENERATOR.uniform(low_bound, high_bound)

            epoch_list.append(chosen_value)
            radius, velocity = planets_sequence[0].eph(chosen_value)
            planets_radii_list.append(radius)
            planets_velocities_list.append(velocity)
            states_sequence.append(normalize(chosen_value, low_bound, high_bound))
            last_arrival_velocity = velocity

        else:
            # The choice of the next value will depend on a candidate generated
            # from the probability weights, and a candidate from the heuristic function
            # Then, we sample according to a mixture of two normal distributions
            # around these two candidates, with 1/tau being the weight of the policy candidate

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
                    weight = gaussian_kernel.pdf(key)
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
            epoch_list.append(epoch_list[-1] + chosen_value)
            planet_radius, planet_velocity = planet.eph(epoch_list[-1])
            planets_velocities_list.append(planet_velocity)
            if (advancement > 1) and (advancement < (len(planets_sequence) - 1)):
                # Heuristic function candidate
                ## Get candidates
                _, bias_candidates, local_delta_v_list = porkchop_scan(
                    planets_sequence[advancement - 1],
                    planets_sequence[advancement],
                    [epoch_list[-1], epoch_list[-1] + 1],
                    [low_bound, high_bound],
                    last_arrival_velocity,
                    n_t0=1,
                    n_tof=N_CANDIDATES,
                )
                ## Keep valid values
                local_delta_v_list = local_delta_v_list.flatten()
                bias_candidates = bias_candidates[
                    np.where(local_delta_v_list < UNFEASIBILITY_VALUE)
                ]
                local_delta_v_list = local_delta_v_list[
                    np.where(local_delta_v_list < UNFEASIBILITY_VALUE)
                ]
                if list(bias_candidates):
                    bias_candidate, bias_sigma = fit_gaussian_from_density(
                        bias_candidates, 1000 / local_delta_v_list
                    )

                    # Choosing the best candidate
                    chosen_value = truncate(
                        sample_mixture_1d(
                            n_samples=1,
                            mu1=policy_weights_candidate,
                            sigma1=denormalize(std_factor, low_bound, high_bound),
                            mu2=bias_candidate,
                            sigma2=bias_sigma,
                            weight1=1 / tau,
                        )[0],
                        low_bound,
                        high_bound,
                    )
                else:
                    chosen_value = truncate(
                        RANDOM_GENERATOR.normal(
                            policy_weights_candidate,
                            (high_bound - low_bound) * std_factor,
                        ),
                        low_bound,
                        high_bound,
                    )
            else:
                chosen_value = truncate(
                    RANDOM_GENERATOR.normal(
                        policy_weights_candidate,
                        (high_bound - low_bound) * std_factor,
                    ),
                    low_bound,
                    high_bound,
                )

            # Computing next state
            lambert_leg = pk.lambert_problem(
                tof=chosen_value * pk.DAY2SEC,
                r0=planets_radii_list[-1],
                r1=planet_radius,
                mu=pk.MU_SUN,
                cw=False,
                multi_revs=0,
            )

            arrival_velocity = np.array(lambert_leg.v1[0])
            planets_radii_list.append(planet_radius)

            last_arrival_velocity = arrival_velocity
            states_sequence.append(
                normalize(
                    chosen_value,
                    low_bound,
                    high_bound,
                )
            )

        values_sequence.append(chosen_value)
    return values_sequence


def run_cgnrpa(
    evaluator: pk.trajopt.mga,
    policy: dict = dict(),
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    planets_sequence: list = None,
    current_iteration: int = 0,
    learning_rate: float = 0.01,
    timeout: float = 10,
    start_time: float = 0,
    best_values_sequence: list = None,
    best_value: float = UNFEASIBILITY_VALUE,
    best_values_list: list = None,
    time_list: list = None,
    tau: float = 10,
    *args,
    **kwargs,
):
    assert not (planets_sequence is None), "planets_sequence is None"
    assert not (bounds is None), "bounds is None"
    current_time = time.time() - start_time
    if level == 0:

        values_sequence = biased_policy_playout(
            policy=policy,
            bounds=bounds,
            planets_sequence=planets_sequence,
            std_factor=0.01 + 1 / np.sqrt(current_iteration + 1),
            tau=tau,
        )
        return (
            values_sequence,
            evaluator.fitness(values_sequence)[0],
        )
    else:
        current_policy = deepcopy(policy)
        for current_iteration in range(n_policies):
            values_sequence, total_delta_v = run_cgnrpa(
                evaluator=evaluator,
                policy=current_policy,
                level=level - 1,
                n_policies=n_policies,
                learning_rate=learning_rate,
                planets_sequence=planets_sequence,
                bounds=bounds,
                current_iteration=current_iteration,
                timeout=timeout,
                start_time=start_time,
                best_values_sequence=best_values_sequence,
                best_value=best_value,
                best_values_list=best_values_list,
                time_list=time_list,
                tau=tau,
            )
            if total_delta_v < best_value:
                best_value = total_delta_v
                best_values_sequence = values_sequence[:]
            current_time = time.time() - start_time
            if best_value < UNFEASIBILITY_VALUE:
                current_policy = adapt_policy(
                    best_values_sequence=best_values_sequence,
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
            best_value,
        )


def cgnrpa(
    evaluator: pk.trajopt.mga,
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    planets_sequence: list = None,
    learning_rate: float = 0.01,
    timeout: float = 10,
    tau: float = 10,
    *args,
    **kwargs,
):
    start_time = time.time()
    best_values_list, time_list = list(), list()
    best_values_sequence, best_value = run_cgnrpa(
        evaluator=evaluator,
        policy=dict(),
        level=level,
        n_policies=n_policies,
        bounds=bounds,
        planets_sequence=planets_sequence,
        current_iteration=0,
        learning_rate=learning_rate,
        timeout=timeout,
        start_time=start_time,
        best_values_sequence=None,
        best_value=UNFEASIBILITY_VALUE,
        best_values_list=best_values_list,
        time_list=time_list,
        tau=tau,
    )
    return best_values_sequence, best_value, best_values_list, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = CountingEvaluator(pk.trajopt.gym.cassini1)
    planets_sequence = [
        pk.planet(pk.udpla.jpl_lp("Earth")),
        pk.planet(pk.udpla.jpl_lp("Venus")),
        pk.planet(pk.udpla.jpl_lp("Venus")),
        pk.planet(pk.udpla.jpl_lp("Earth")),
        pk.planet(pk.udpla.jpl_lp("Jupiter")),
        pk.planet(pk.udpla.jpl_lp("Saturn")),
    ]

    # Variables bounds
    bounds = [
        (low_bound, high_bound)
        for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    ]

    # General input values
    inputs_values = {
        "evaluator": udp,
        "planets_sequence": planets_sequence,
        "bounds": bounds,
        "timeout": 120,
        "level": 3,
        "n_policies": 17908,
        "learning_rate": 0.154937318595788,
        "tau": 2.545153120408655,
    }
    values_sequence, best_value, values_list, time_list = cgnrpa(**inputs_values)
    print(f"Delta V: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
    print(f"Total number of evaluations: {udp.count}")
