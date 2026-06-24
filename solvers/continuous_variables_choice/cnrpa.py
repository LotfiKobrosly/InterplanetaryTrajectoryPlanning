"""
Implements Continuous Nested Rollout Policy for continuous variables' values choice.
We will be implementing the Gaussian Kernel based variant
Here, the values sequence is represented as follows:
values_sequence = [departure_epoch, time_of_flight_0, ..., time_of_flight_n]

Values stored in the policy will be normalized between 0 and 1 to account for the huge differences in dimensions

Two variations are proposed:
- A state independent variation: we advance through the sequence, where for each step, we have a central value (of the corresponding
dimension) and we sample around it using a normal distribution (parameters are then the mean and the standard deviation)
- A state dependent variation:
    * In the first step we choose the departure epoch, which happend in a similar fashion to the other variation
    * In the second step, we take the departure epoch as our "state", and we compute the value of time_of_flight
    of the next encounter (second planet in the sequence) by using a gaussian kernel from the other states existing in this step, then we sample
    around it using a normal distribution
    * In the remaining steps, the state will become the departure velocity from the planet to the next one, computed from the lambert leg and
    subsequent gravity assist impulse. Simlarly to the previous step, we compute time_of_flight by sampling with a normal
    distribution around a mean computed with a gaussian kernel
"""

from copy import deepcopy
import numpy as np
import pykep as pk
from utils.trajectory_evaluation import evaluate_mga_trajectory
from utils.constants import (
    GAUSSIAN_KERNEL_THRESHOLD,
    RANDOM_GENERATOR,
    VELOCITY_NORMALIZING_FACTOR,
)
from utils.basic_functions import normalize, denormalize, code, truncate
from solvers.continuous_variables_choice.values_separators import (
    separate_values,
)
from utils.gaussian_kernel import (
    GaussianKernel,
)


def policy_playout(
    policy: dict,
    bounds: list,
    multiple_values_policy: bool,
    planets_sequence: list = None,
    states_sequence: list = list(),
    std_factor: float = 1,
):
    planets_sequence = [
        pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence
    ]
    values_sequence = list()
    epoch_list, planets_radii_list = list(), list()
    for advancement, planet in enumerate(planets_sequence):
        # Get bounds
        low_bound, high_bound = bounds[advancement]

        if multiple_values_policy:
            # Choose departure epoch
            if advancement == 0:
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
                radius, _ = planets_sequence[0].eph(chosen_value)
                planets_radii_list.append(radius)
                states_sequence.append(normalize(chosen_value, low_bound, high_bound))

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
                        chosen_value = truncate(
                            denormalize(
                                RANDOM_GENERATOR.normal(
                                    weights @ np.array(values).T, std_factor
                                ),
                                low_bound,
                                high_bound,
                            ),
                            low_bound,
                            high_bound,
                        )
                    else:
                        chosen_value = RANDOM_GENERATOR.uniform(low_bound, high_bound)
                else:
                    chosen_value = RANDOM_GENERATOR.uniform(low_bound, high_bound)

                # Computing Lambert leg
                epoch_list.append(epoch_list[-1] + chosen_value)
                planet_radius, planet_velocity = planet.eph(epoch_list[-1])
                lambert_leg = pk.lambert_problem(
                    tof=chosen_value * pk.DAY2SEC,
                    r0=planets_radii_list[-1],
                    r1=planet_radius,
                    mu=pk.MU_SUN,
                    cw=False,
                    multi_revs=0,
                )
                arrival_velocity = np.array(lambert_leg.v1[0])
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
                chosen_value = truncate(
                    denormalize(
                        RANDOM_GENERATOR.normal(policy[0], std_factor),
                        low_bound,
                        high_bound,
                    ),
                    low_bound,
                    high_bound,
                )
            else:
                chosen_value = RANDOM_GENERATOR.uniform(low_bound, high_bound)
        values_sequence.append(chosen_value)
    # assert (
    #     (len(values_sequence) == len(states_sequence) + 1) and multiple_values_policy
    # ) or not multiple_values_policy, (
    #     "Error in states_sequence and values_sequence generation: "
    #     + "values_sequence size: "
    #     + str(len(values_sequence))
    #     + ", VS states_sequence size: "
    #     + str(len(states_sequence))
    #     + ". Current lists: \n"
    #     + str(values_sequence)
    #     + "\n"
    #     + str(states_sequence)
    # )
    return values_sequence


def adapt_policy(
    best_values_sequence: list,
    best_states_sequence: list,
    policy: dict,
    multiple_values_policy: bool = False,
    learning_rate: float = 0.01,
    bounds: list = None,
):
    if policy:
        for advancement, element in enumerate(best_values_sequence):
            low_bound, high_bound = bounds[advancement]
            if advancement == 0:
                policy[advancement] += learning_rate * (
                    normalize(element, low_bound, high_bound) - policy[advancement]
                )
            else:
                if multiple_values_policy:
                    current_key = code(best_states_sequence[advancement - 1])
                    if current_key in policy[advancement].keys():
                        previous_value = policy[advancement][current_key]
                        policy[advancement][current_key] = (
                            previous_value
                            + learning_rate
                            * (
                                normalize(element, low_bound, high_bound)
                                - previous_value
                            )
                        )
                    else:
                        policy[advancement][current_key] = normalize(
                            element, low_bound, high_bound
                        )
                    gaussian_kernel = GaussianKernel(current_key, 0.2)
                    for key in policy[advancement].keys():
                        weight = gaussian_kernel.pdf(key)
                        if weight >= GAUSSIAN_KERNEL_THRESHOLD:
                            previous_value = policy[advancement][key]
                            policy[advancement][key] = (
                                previous_value
                                + learning_rate
                                * weight
                                * (
                                    normalize(element, low_bound, high_bound)
                                    - previous_value
                                )
                            )
                else:
                    previous_value = np.array(policy[advancement])
                    policy[advancement] = previous_value + learning_rate * (
                        normalize(element, low_bound, high_bound) - previous_value
                    )
    else:
        for advancement, element in enumerate(best_values_sequence):
            low_bound, high_bound = bounds[advancement]
            if multiple_values_policy:
                if advancement == 0:
                    policy[advancement] = normalize(element, low_bound, high_bound)
                else:
                    try:
                        policy[advancement] = {
                            code(best_states_sequence[advancement - 1]): normalize(
                                element, low_bound, high_bound
                            )
                        }
                    except IndexError:
                        print(
                            "Size values sequence: "
                            + str(len(best_values_sequence))
                            + ", vs states sequence: "
                            + str(len(best_states_sequence))
                        )
                        raise IndexError
            else:
                policy[advancement] = normalize(element, low_bound, high_bound)
    return policy


def cnrpa(
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
    *args,
    **kwargs,
):
    assert not (planets_sequence is None), "planets_sequence is None"
    assert not (bounds is None), "bounds is None"
    if level == 0:

        values_sequence = policy_playout(
            policy=policy,
            bounds=bounds,
            multiple_values_policy=multiple_values_policy,
            planets_sequence=planets_sequence,
            states_sequence=states_sequence,
            std_factor=0.01 + 0.5 * np.exp(-10 * current_iteration / (n_policies)),
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
            values_sequence, total_delta_v = cnrpa(
                policy=current_policy,
                multiple_values_policy=multiple_values_policy,
                level=level - 1,
                n_policies=n_policies,
                planets_sequence=planets_sequence,
                bounds=bounds,
                current_iteration=current_iteration,
                states_sequence=states_sequence,
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
