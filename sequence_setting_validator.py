"""
Validation script for MCTS-based algorithms on cassini1 MGA model
"""
from copy import deepcopy
import random
import numpy as np
import pykep as pk
from utils.constants import (
    PLANETS,
    VARIABLES_BOUNDS,
    SAMPLING_FUNCTIONS,
    UNFEASIBILITY_VALUE,
    RANDOM_GENERATOR,
)

def rollout(policy, planets_pool, t0_bounds, tof_bounds, values_vector):
    planets_sequence = ["earth"]
    # Do rollout
    for planet_id in range(4):
        # Select next planet
        choices = list()
        probabilities = list()
        for planet, weight in policy[planet_id].items():
            choices.append(planet)
            probabilities.append(np.exp(weight))
        probabilities = np.array(probabilities)
        probabilities /= np.sum(probabilities)
        next_planet = RANDOM_GENERATOR.choice(choices, p=probabilities)
        planets_sequence.append(next_planet)

    planets_sequence.append("saturn")

    # Evaluation
    trajectory = pk.trajopt.mga(
        seq=[pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence],
        t0=t0_bounds,
        tof=tof_bounds,
        vinf=3,
    )

    return planets_sequence, trajectory.fitness(values_vector)[0]

def adapt_policy(policy, planets_sequence, learning_rate):
    for id in range(4):
        normalizer = 0
        policy[id][planets_sequence[id + 1]] += learning_rate
        for weight in policy[id].values():
            normalizer += np.exp(weight)
        for planet, weight in policy[id].items():
            policy[id][planet] -= learning_rate * np.exp(weight) / normalizer


def nrpa_test(
    policy,
    level,
    n_policies,
    planets_pool,
    t0_bounds,
    tof_bounds,
    values_vector,
    learning_rate,
    global_best_sequence,
    global_best_dv,
):
    if level == 0:
        return rollout(policy, planets_pool, t0_bounds, tof_bounds, values_vector)
    else:
        local_best_dv = np.inf
        local_best_sequence = None
        for _ in range(n_policies):
            new_policy = deepcopy(policy)
            planets_sequence, total_dv = nrpa_test(
                new_policy,
                level - 1,
                n_policies,
                planets_pool,
                t0_bounds,
                tof_bounds,
                values_vector,
                learning_rate,
                local_best_sequence,
                local_best_dv,
            )
            if total_dv < local_best_dv:
                local_best_dv = total_dv
                local_best_sequence = planets_sequence
            adapt_policy(new_policy, local_best_sequence, learning_rate)

        if local_best_dv < global_best_dv:
            global_best_dv = local_best_dv
            global_best_sequence = local_best_sequence
        adapt_policy(policy, global_best_sequence, learning_rate)
        return global_best_sequence, global_best_dv


if __name__ == "__main__":
    udp = pk.trajopt.gym.cassini1
    bounds = [
        [low_bound, high_bound]
        for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    ]
    t0_bounds = bounds[0]
    tof_bounds = bounds[1:]
    values_vector = [-773.6954071581413, 180.49035683964848, 413.5167628675508, 52.96396625990729, 1005.2849447090892, 4543.903838074617]
    planets_pool = ["earth", "venus", "mercury", "mars", "jupiter"]
    policy = [
        {planet: 0 for planet in planets_pool}
        for _ in range(5)
    ]
    sequence, dv = nrpa_test(
        policy,
        2,
        300,
        planets_pool,
        t0_bounds,
        tof_bounds,
        values_vector,
        0.1,
        None,
        np.inf,
    )
    print("Best sequence:", sequence)
    print(f"Total DV: {dv / 1000:.3f} km/s")

    