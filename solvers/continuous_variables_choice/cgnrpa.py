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

from copy import deepcopy
import numpy as np
import pykep as pk
from utils.trajectory_evaluation import evaluate_mga_trajectory
from utils.constants import (
    GAUSSIAN_KERNEL_THRESHOLD,
    RANDOM_GENERATOR,
    VELOCITY_NORMALIZING_FACTOR,
    UNFEASIBILITY_VALUE,
    N_CANDIDATES,
    SAFE_RADIUS_FACTOR,
)
from utils.basic_functions import normalize, denormalize, code, truncate
from solvers.continuous_variables_choice.values_separators import (
    separate_values,
)
from solvers.continuous_variables_choice.cnrpa import adapt_policy
from utils.gaussian_kernel import (
    GaussianKernel,
)


def biased_policy_playout(
    policy: dict,
    bounds: list,
    multiple_values_policy: bool,
    planets_sequence: list = None,
    states_sequence: list = list(),
    std_factor: float = 1,
    tau: float = 10,
):
    planets_sequence = [
        pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence
    ]
    values_sequence = list()
    epoch_list, planets_radii_list, planets_velocities_list = list(), list(), list()
    last_arrival_velocity = None
    for advancement, planet in enumerate(planets_sequence):
        # Get bounds
        low_bound, high_bound = bounds[advancement]

        if multiple_values_policy:
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
                if policy:
                    current_policy = policy[advancement]
                    gaussian_kernel = GaussianKernel(
                        states_sequence[-1], sigma=std_factor
                    )
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
                        candidates_list = [
                            truncate(
                                denormalize(
                                    weights @ np.array(values).T, low_bound, high_bound
                                ),
                                low_bound,
                                high_bound,
                            )
                            for _ in range(N_CANDIDATES)
                        ]
                        weights_list = [
                            gaussian_kernel.pdf(candidate)
                            for candidate in candidates_list
                        ]
                    else:
                        candidates_list = RANDOM_GENERATOR.uniform(
                            low_bound, high_bound, size=N_CANDIDATES
                        )
                        weights_list = [1 / N_CANDIDATES] * N_CANDIDATES
                else:
                    candidates_list = RANDOM_GENERATOR.uniform(
                        low_bound, high_bound, size=N_CANDIDATES
                    )
                    weights_list = [1 / N_CANDIDATES] * N_CANDIDATES

                # Lambet legs
                epoch_list.append(epoch_list[-1] + chosen_value)
                planet_radius, planet_velocity = planet.eph(epoch_list[-1])
                planets_velocities_list.append(planet_velocity)
                local_delta_v_list = list()
                for candidate_number, candidate in enumerate(candidates_list):
                    lambert_leg = pk.lambert_problem(
                        tof=candidate * pk.DAY2SEC,
                        r0=planets_radii_list[-1],
                        r1=planet_radius,
                        mu=pk.MU_SUN,
                        cw=False,
                        multi_revs=0,
                    )
                    departure_velocity = np.array(lambert_leg.v0[0])

                    if advancement == 1:
                        # After departure, also valid if planet_sequence is of size 2
                        local_delta_v_list.append(
                            np.linalg.norm(
                                departure_velocity - planets_velocities_list[0]
                            )
                        )
                    else:
                        # Check feasibility and compute delta_v
                        try:
                            _, violation = pk.fb_con(
                                v_rel_in=(
                                    last_arrival_velocity
                                    - planets_velocities_list[advancement - 1]
                                ).tolist(),
                                v_rel_out=(
                                    departure_velocity
                                    - planets_velocities_list[advancement - 1]
                                ).tolist(),
                                mu=planets_sequence[advancement - 1].mu_self,
                                safe_radius=planets_sequence[advancement - 1].radius
                                * SAFE_RADIUS_FACTOR,
                                # pl        = planet  # pykep planet object knows its own mu, radius
                            )
                        except IndexError:
                            print("Advancement: ", advancement)
                            print("Planets: ", planets_sequence)
                            print("Velocities: ", planets_velocities_list)
                            raise IndexError
                        if (
                            violation > 0
                        ):  # value needs to be non-positive to be geometrically feasible
                            local_delta_v_list.append(UNFEASIBILITY_VALUE)
                        else:
                            delta_velocity_planet = pk.fb_dv(
                                v_rel_in=(
                                    last_arrival_velocity
                                    - planets_velocities_list[advancement - 1]
                                ).tolist(),
                                v_rel_out=(
                                    departure_velocity
                                    - planets_velocities_list[advancement - 1]
                                ).tolist(),
                                mu=planets_sequence[advancement - 1].mu_self,
                                safe_radius=planets_sequence[advancement - 1].radius
                                * SAFE_RADIUS_FACTOR,  # 5% margin above surface
                            )
                            local_delta_v_list.append(
                                np.linalg.norm(delta_velocity_planet)
                            )
                # Choosing the best candidate
                indices = np.argsort(local_delta_v_list)
                candidates_list = np.array(candidates_list)[indices]
                weights_list = np.array(weights_list)[indices]
                local_delta_v_list = np.array(local_delta_v_list)[indices]
                biases_list = [0] * N_CANDIDATES
                for index in range(N_CANDIDATES):
                    if local_delta_v_list[index] >= UNFEASIBILITY_VALUE:
                        weights_list[index] = 0
                    else:
                        biases_list[index] = 1 / (index + 1)
                probabilites = np.exp(
                    np.array(weights_list) / tau + np.array(biases_list)
                )
                probabilites /= np.sum(probabilites)
                chosen_value = RANDOM_GENERATOR.choice(candidates_list, p=probabilites)

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
                        arrival_velocity,
                        np.zeros(np.shape(arrival_velocity)),
                        VELOCITY_NORMALIZING_FACTOR
                        * np.ones(np.shape(arrival_velocity)),
                    )
                )
        else:
            if policy:
                chosen_value = RANDOM_GENERATOR.normal(
                    policy[advancement], std_factor * (high_bound - low_bound) / 10
                )
                chosen_value = truncate(chosen_value, low_bound, high_bound)
            else:
                chosen_value = RANDOM_GENERATOR.uniform(low_bound, high_bound)
        values_sequence.append(chosen_value)
    return values_sequence


def cgnrpa(
    values_sequence: list = list(),
    policy: dict = dict(),
    multiple_values_policy: bool = False,
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    planets_sequence: list = None,
    current_iteration: int = 0,
    states_sequence: list = list(),
    learning_rate: float = 0.01,
    tau: float = 10,
    *args,
    **kwargs,
):
    assert not (planets_sequence is None), "planets_sequence is None"
    assert not (bounds is None), "bounds is None"
    if level == 0:

        values_sequence = biased_policy_playout(
            policy=policy,
            bounds=bounds,
            multiple_values_policy=multiple_values_policy,
            planets_sequence=planets_sequence,
            states_sequence=states_sequence,
            std_factor=0.01 + 0.5 * np.exp(-10 * current_iteration / (n_policies)),
            tau=tau,
        )
        return (
            values_sequence,
            evaluate_mga_trajectory(
                planets_sequence,
                *separate_values(values_sequence),
            )[0],
        )
    else:
        best_delta_v = np.inf
        current_policy = deepcopy(policy)
        best_values_sequence = None
        best_states_sequence = None
        for current_iteration in range(n_policies):
            states_sequence = list()
            values_sequence, total_delta_v = cgnrpa(
                policy=current_policy,
                multiple_values_policy=multiple_values_policy,
                level=level - 1,
                n_policies=n_policies,
                planets_sequence=planets_sequence,
                bounds=bounds,
                current_iteration=current_iteration,
                states_sequence=states_sequence,
                tau=tau,
            )
            if total_delta_v < best_delta_v:
                best_delta_v = total_delta_v
                best_values_sequence = values_sequence
                best_states_sequence = states_sequence
            current_policy = adapt_policy(
                best_values_sequence=best_values_sequence,
                best_states_sequence=best_states_sequence,
                policy=current_policy,
                multiple_values_policy=multiple_values_policy,
                learning_rate=learning_rate,
                bounds=bounds,
            )
        return (
            best_values_sequence,
            evaluate_mga_trajectory(
                planets_sequence,
                *separate_values(best_values_sequence),
            )[0],
        )
