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

import time
from copy import deepcopy
import numpy as np
import pykep as pk
import matplotlib.pyplot as plt
from utils.constants import (
    GAUSSIAN_KERNEL_THRESHOLD,
    RANDOM_GENERATOR,
    VELOCITY_NORMALIZING_FACTOR,
    DV_LAUNCHER,
    UNFEASIBILITY_VALUE,
)
from utils.basic_functions import normalize, denormalize, code, truncate
from utils.gaussian_kernel import (
    GaussianKernel,
)
from utils.udp_wrapper import CountingEvaluator


def policy_playout(
    policy: dict,
    bounds: list,
    std_factor: float = 1,
    *args,
    **kwargs,
):
    values_sequence, states_sequence = list(), list()
    for advancement, bound in enumerate(bounds):
        # Get bounds
        low_bound, high_bound = bound

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

            states_sequence.append(normalize(chosen_value, low_bound, high_bound))

        else:
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

            states_sequence.append(normalize(chosen_value, low_bound, high_bound))

        values_sequence.append(chosen_value)
    return values_sequence

def adapt_policy(
    best_values_sequence: list,
    policy: dict,
    learning_rate: float = 0.01,
    bounds: list = None,
):
    best_values_sequence = [
        normalize(value, *bound)
        for (value, bound) in zip(best_values_sequence, bounds)
    ]
    if policy:
        for advancement, element in enumerate(best_values_sequence):
            low_bound, high_bound = bounds[advancement]
            if advancement == 0:
                policy[advancement] += learning_rate * (
                    element - policy[advancement]
                )
            else:
                current_key = code(best_values_sequence[:advancement])
                if current_key in policy[advancement].keys():
                    previous_value = policy[advancement][current_key]
                    policy[advancement][current_key] = (
                        previous_value
                        + learning_rate
                        * (element - previous_value)
                    )
                else:
                    policy[advancement][current_key] = element
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
                                element - previous_value
                            )
                        )

    else:
        for advancement, element in enumerate(best_values_sequence):
            low_bound, high_bound = bounds[advancement]
            if advancement == 0:
                policy[advancement] = element
            else:
                policy[advancement] = {
                    code(best_values_sequence[:advancement]): element
                }
    return policy


def run_cnrpa(
    evaluator: pk.trajopt.mga,
    policy: dict = dict(),
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    current_iteration: int = 0,
    learning_rate: float = 0.01,
    timeout: float = 10,
    start_time: float = 0,
    best_values_sequence: list = None,
    best_value: float = UNFEASIBILITY_VALUE,
    best_values_list: list = None,
    time_list: list = None,
    n_processes: int = 1,
    *args,
    **kwargs,
):
    assert not (bounds is None), "bounds is None"
    current_time = time.time() - start_time
    if level == 0:
        values_sequence = policy_playout(
            policy=policy,
            bounds=bounds,
            std_factor=0.01 + 1 / np.sqrt(current_iteration + 2),
        )

        return (
            values_sequence,
            evaluator.fitness(values_sequence)[0],
        )

    else:
        current_policy = deepcopy(policy)
        current_best_value = UNFEASIBILITY_VALUE
        current_best_sequence = None
        for current_iteration in range(n_policies):
            values_sequence, total_delta_v = run_cnrpa(
                evaluator=evaluator,
                policy=current_policy,
                level=level - 1,
                n_policies=n_policies,
                learning_rate=learning_rate,
                bounds=bounds,
                current_iteration=current_iteration,
                timeout=timeout,
                start_time=start_time,
                best_values_sequence=best_values_sequence,
                best_value=best_value,
                best_values_list=best_values_list,
                time_list=time_list,
            )
            if total_delta_v < current_best_value:
                current_best_value = total_delta_v
                current_best_sequence = values_sequence[:]
            
            if current_best_value < UNFEASIBILITY_VALUE:
                current_policy = adapt_policy(
                    best_values_sequence=current_best_sequence,
                    policy=current_policy,
                    learning_rate=learning_rate,
                    bounds=bounds,
                )
            if current_best_value < best_value:
                best_value = current_best_value
                best_values_sequence = current_best_sequence[:]
            current_time = time.time() - start_time
            if best_value < UNFEASIBILITY_VALUE:
                best_values_list.append(best_value)
                time_list.append(current_time)
            if current_time > timeout:
                break
        return (
            best_values_sequence,
            best_value,
        )


def cnrpa(
    evaluator: pk.trajopt.mga,
    level: int = 0,
    n_policies: int = 10,
    bounds: list = None,
    learning_rate: float = 0.01,
    timeout: float = 10,
    *args,
    **kwargs,
):
    start_time = time.time()
    best_values_list, time_list = list(), list()
    best_values_sequence, best_value = run_cnrpa(
        evaluator=evaluator,
        policy=dict(),
        level=level,
        n_policies=n_policies,
        bounds=bounds,
        current_iteration=0,
        learning_rate=learning_rate,
        timeout=timeout,
        start_time=start_time,
        best_values_sequence=None,
        best_value=UNFEASIBILITY_VALUE,
        best_values_list=best_values_list,
        time_list=time_list,
    )
    return best_values_sequence, best_value, best_values_list, time_list


if __name__ == "__main__":
    # Problem
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
        "timeout": 120,
        "level": 2,
        "learning_rate": 0.0030163,
        "n_policies": 5268,
    }
    values_sequence, best_value, values_list, time_list = cnrpa(**inputs_values)
    print(f"Best Delta V: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
    print(f"Total number of evaluations: {udp.count}")
    print("Details:")
    print(udp.pretty(values_sequence))

    figure = plt.figure(figsize=(10, 10))
    plt.plot(time_list, values_list)
    plt.title(f"Best value {best_value / 1000:.3f} km/s found first after {time_list[values_list.index(best_value)]:.3f} s")
    plt.show()

    axe = udp.plot(values_sequence, figsize=(20, 20))
    # figure = axe.figure
    axe.view_init(90, 0)
    axe.axis("off")
    axe.set_title(
        "GcNRPA"
        + r": $\Delta$V = "
        + f"{best_value / 1000:.3f} km/s"
    )
    plt.show()
