"""
Verifying solvers capabilities
"""

import time
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pykep as pk
import pygmo as pg
from solvers.continuous_variables_choice.cnmcts import cnmcts
from solvers.continuous_variables_choice.cnrpa import cnrpa
from solvers.continuous_variables_choice.cgnrpa import cgnrpa
from solvers.continuous_variables_choice.cabgnrpa import cabgnrpa
from solvers.continuous_variables_choice.baselines import *
from utils.constants import UNFEASIBILITY_VALUE, DV_LAUNCHER


if __name__ == "__main__":
    print("Testing Cassini MGA:\n")

    SAMPLING_FUNCTIONS = {
        "cmaes": pygmo_baseline,
        "sade": pygmo_baseline,
        "sga": pygmo_baseline,
        "simulated_annealing": pygmo_baseline,
        # "pso": pygmo_baseline,
        "gaco": pygmo_baseline,
        "bee_colony": pygmo_baseline,
        # "uniform": uniform_variables_values_vector,
        # "gaussian_cma_es": gaussian_variables_values_vector,
        "cnmcts": cnmcts,
        "cnrpa": cnrpa,
        "cgnrpa": cgnrpa,
        "cabgnrpa": cabgnrpa,
        # "genetic": genetic_algorithm,
    }

    algorithms_list = list(SAMPLING_FUNCTIONS.keys())

    SAMPLING_FUNCTIONS_NAMES = {
        "cmaes": "Covariance Matrix Adaptation Evolution Strategy (pygmo)",
        "sade": "Self Adaptive Differential Evolution (pygmo)",
        "sga": "Simple Genetic Algorithm (pygmo)",
        "simulated_annealing": "Simulated Annealing (pygmo)",
        "pso": "Particle Swarm Optimization (pygmo)",
        "gaco": "Generalized Ant Colony Optimization (pygmo)",
        "bee_colony": "Artificial Bee Colony (pygmo)",
        "uniform": "Uniform sampling",
        "gaussian_cma_es": "Custom Gaussian CMA-ES",
        "cnmcts": "cNMCS",
        "cnrpa": "cNRPA",
        "cgnrpa": "cGNRPA",
        "cabgnrpa": "cABGNRPA",
        # "genetic": "Custom Genetic Algorithm",
    }

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
        "n_iterations": 2500,  # for uniform and gaussian sampling
        "n_generations": 1000,  # for Genetic Algorithm
        "population_size": 500,  # for Genetic Algorithm
        "mutation_probability": 0.15,  # for Genetic Algorithm
        "mutation_effect": 0.25,  # for Genetic Algorithm
    }

    # Timeouts list
    timeouts_list = [1, 5, 10, 30, 60]  # , 300, 600]

    # Specific parameters
    algorithm_parameters = {
        "cmaes": {"solver": "cmaes"},
        "sade": {"solver": "sade"},
        "sga": {"solver": "sga"},
        "pso": {"solver": "pso"},
        "gaco": {"solver": "gaco"},
        "bee_colony": {"solver": "bee_colony"},
        "simulated_annealing": {"solver": "simulated_annealing"},
        "uniform": {},
        "gaussian_cma_es": {},
        "cnmcts": {"level": 1, "bandwidth": int(5e6)},
        "cnrpa": {"level": 1, "n_policies": int(1e12), "learning_rate": 0.5},
        "cgnrpa": {
            "level": 1,
            "n_policies": int(1e12),
            "learning_rate": 0.5,
            "tau": 5,
        },
        "cabgnrpa": {
            "level": 1,
            "n_policies": int(1e12),
            "learning_rate": 0.5,
            "tau": 5,
            "gamma": 0.25,
        },
        # "genetic": {"n_generations": 5000, "mutation_probability": 0.15},
    }

    # Storing results
    results = np.zeros((len(SAMPLING_FUNCTIONS), len(timeouts_list)))

    for algorithm_id, algorithm in enumerate(algorithms_list):
        for timeout_id, timeout in enumerate(timeouts_list):
            print("\nFor " + algorithm + ":")
            specific_input_values = deepcopy(inputs_values)
            for arg, val in algorithm_parameters[algorithm].items():
                specific_input_values[arg] = val
            specific_input_values["timeout"] = timeout
            values_sequence, delta_v, _, _ = SAMPLING_FUNCTIONS[algorithm](
                **specific_input_values
            )
            results[algorithm_id, timeout_id] = min(delta_v, UNFEASIBILITY_VALUE)
            if delta_v < UNFEASIBILITY_VALUE:
                print(f"Delta V: {delta_v / 1000:.3f} km/s")

                # plt.plot(time_list, best_values_list)
                # plt.show()
                # axe = udp.plot(values_sequence, figsize=(20, 20))
                # figure = axe.figure
                # axe.view_init(90, 0)
                # axe.axis("off")
                # axe.set_title(
                #     algorithm.upper()
                #     + r": $\Delta$V = "
                #     + f"{udp.fitness(values_sequence)[0] / 1000:.3f} km/s"
                # )
                # figure.savefig(
                #     "./Cassini trials/" + SAMPLING_FUNCTIONS_NAMES[algorithm] + ".png"
                # )
                # plt.close(figure)
                # plt.show()
            # else:
            #     print("No valid solution produced")

    dataframe = pd.DataFrame(results, index=algorithms_list, columns=timeouts_list)
    dataframe.to_excel("Preliminary results Cassini.xlsx")
