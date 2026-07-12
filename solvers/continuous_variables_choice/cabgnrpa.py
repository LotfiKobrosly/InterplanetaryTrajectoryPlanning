"""
Implements Continuous Adaptive Bias Generalized Nested Rollout Policy for continuous variables' values choice
Here, the values sequence is represented as follows:
values_sequence = [departure_epoch, time_of_flight_0, ..., time_of_flight_n]

We would do the lambert leg values computation as well as the gravity assist impulses as a way to get heuristic values for GNRPA:
- At each step for action selection (except first step), generate candidates from the policy and store their weights
- Compute for each candidate TOF the Lambert leg and the flyby result
- Use these results to generate the biases values
- Add the biases to the candidates' weights and do the selection

Bias values are defined for each planet by the states which are represented as the normalized values of:
- Arrival relative velocity of the spacecraft
- The position of the planet
- The velocity of the planet
For each state, a division of the time of flight range, with a weight given to each value.
These weights change according to the ABGNRPA process of bias adaptation.
These biases are taken into consideration when the encountered state is close-by (speciefied radius).

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
    UNFEASIBILITY_VALUE,
    N_CANDIDATES,
    SAFE_RADIUS_FACTOR,
    DV_LAUNCHER,
    PLANETS_SEMI_MAJOR_AXIS_AND_EXCENTRICITY,
)
from utils.basic_functions import *
from utils.heuristic_functions import porkchop_scan
from utils.gaussian_kernel import (
    GaussianKernel,
)
from utils.orbital_elements import orbital_speeds
from solvers.continuous_variables_choice.cnrpa import adapt_policy


def adaptive_bias_policy_playout(
    policy: dict,
    biases_values: dict,
    bounds: list,
    planets_sequence: list = None,
    std_factor: float = 1,
    tau: float = 10,
    gamma: float = 0.1,
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

            epoch_list.append(epoch_list[-1] + chosen_value)
            planet_radius, planet_velocity = planet.eph(epoch_list[-1])
            planets_velocities_list.append(planet_velocity)
            current_state = states_sequence[:]
            planet_name = planets_sequence[advancement - 1].get_name()[:-8]

            if (advancement > 1) and (advancement < (len(planets_sequence) - 1)):
                # Bias values' based candidate
                ## Get candidates
                gaussian_kernel = GaussianKernel(current_state, 0.2)

                if planet_name in biases_values.keys():

                    bias_values, bias_candidates = list(), list()
                    for state, candidates in biases_values[planet_name].items():
                        weight = gaussian_kernel.pdf(state)
                        if weight > GAUSSIAN_KERNEL_THRESHOLD:
                            for candidate, candidate_weight in biases_values[
                                planet_name
                            ][state].items():
                                bias_candidates.append(candidate)
                                bias_values.append(weight * candidate_weight)
                    bias_values = np.array(bias_values)
                    relevant_nearby_states = (len(bias_values) > 0) and (
                        np.sum(bias_values) > 10 * GAUSSIAN_KERNEL_THRESHOLD
                    )

                else:
                    relevant_nearby_states = False
                if not relevant_nearby_states:
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
                    bias_values = 1000 / local_delta_v_list
                bias_values /= np.sum(bias_values)
                if list(bias_candidates):
                    bias_candidate, bias_sigma = fit_gaussian_from_density(
                        bias_candidates, bias_values
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
                    # Adding current_state to biases_values
                    if biases_values.get(planet_name, None) is None:
                        biases_values[planet_name] = dict()
                    adaptation_kernel = GaussianKernel(
                        chosen_value, (high_bound - low_bound) * std_factor
                    )
                    if code(current_state) not in biases_values[planet_name].keys():
                        biases_values[planet_name][code(current_state)] = {
                            code(tof): adaptation_kernel.pdf(tof)
                            for tof in range(
                                int(low_bound),
                                int(high_bound) + 1,
                                int((high_bound - low_bound) / N_CANDIDATES),
                            )
                        }
                    if relevant_nearby_states:
                        for state, candidates in biases_values[planet_name].items():
                            weight = gaussian_kernel.pdf(state)
                            for candidate, candidate_weight in candidates.items():
                                candidates[candidate] = (
                                    1
                                    + gamma * weight * adaptation_kernel.pdf(candidate)
                                ) * candidate_weight
                            # Normalization
                            noramlizing_factor = np.sum(list(candidates.values()))
                            for candidate in candidates.keys():
                                candidates[candidate] /= noramlizing_factor
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

            last_arrival_velocity = lambert_leg.v1[0]
            min_radius, max_radius, min_velocity, max_velocity = orbital_speeds(
                *PLANETS_SEMI_MAJOR_AXIS_AND_EXCENTRICITY[
                    planet_name
                ]  # 8 because get_name returns 'planet(jpl_elp)'
            )
            states_sequence.append(
                normalize(
                    chosen_value,
                    low_bound,
                    high_bound,
                )
            )

            planets_radii_list.append(planet_radius)
            planets_velocities_list.append(planet_velocity)

        values_sequence.append(chosen_value)
    return values_sequence, states_sequence


def run_cabgnrpa(
    evaluator: pk.trajopt.mga,
    policy: dict = dict(),
    biases_values: dict = dict(),
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    planets_sequence: list = None,
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
    gamma: float = 0.1,
    *args,
    **kwargs,
):
    assert not (planets_sequence is None), "planets_sequence is None"
    assert not (bounds is None), "bounds is None"
    current_time = time.time() - start_time
    if level == 0:

        values_sequence, states_sequence = adaptive_bias_policy_playout(
            policy=policy,
            biases_values=biases_values,
            bounds=bounds,
            planets_sequence=planets_sequence,
            std_factor=0.01 + 1 / np.sqrt(current_iteration + 1),
            tau=tau,
            gamma=gamma,
        )
        return (
            values_sequence,
            states_sequence,
            evaluator.fitness(values_sequence)[0],
        )
    else:
        current_policy = deepcopy(policy)
        current_biases_values = deepcopy(biases_values)
        for current_iteration in range(n_policies):
            values_sequence, states_sequence, total_delta_v = run_cabgnrpa(
                evaluator=evaluator,
                policy=current_policy,
                biases_values=current_biases_values,
                learning_rate=learning_rate,
                level=level - 1,
                n_policies=n_policies,
                planets_sequence=planets_sequence,
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
                gamma=gamma,
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


def cabgnrpa(
    evaluator: pk.trajopt.mga,
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    planets_sequence: list = None,
    learning_rate: float = 0.01,
    timeout: float = 10,
    tau: float = 10,
    gamma: float = 0.1,
    *args,
    **kwargs,
):
    start_time = time.time()
    best_values_list, time_list = list(), list()
    best_values_sequence, best_states_sequence, best_value = run_cabgnrpa(
        evaluator=evaluator,
        policy=dict(),
        biases_values=dict(),
        level=level,
        n_policies=n_policies,
        bounds=bounds,
        planets_sequence=planets_sequence,
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
        gamma=gamma,
    )
    return best_values_sequence, best_value, best_values_list, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = pk.trajopt.gym.cassini1
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
        "timeout": 300,
        "level": 2,
        "learning_rate": 0.1,
        "n_policies": 100,
        "tau": 5,
        "gamma": 0.5,
    }
    values__sequence, best_value, values_list, time_list = cabgnrpa(**inputs_values)
    print(f"Best Delta V: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
