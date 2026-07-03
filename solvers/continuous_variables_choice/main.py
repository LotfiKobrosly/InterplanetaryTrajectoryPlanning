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
from solvers.continuous_variables_choice.crbnmcts import crbnmcts
from solvers.continuous_variables_choice.cnrpa import cnrpa
from solvers.continuous_variables_choice.cgnrpa import cgnrpa
from solvers.continuous_variables_choice.cabgnrpa import cabgnrpa
from solvers.continuous_variables_choice.baselines import *
from classes.trajectory import Trajectory
from utils.constants import UNFEASIBILITY_VALUE, DV_LAUNCHER


if __name__ == "__main__":
    print("Testing Cassini MGA:\n")

    SAMPLING_FUNCTIONS = {
        # "cmaes": pygmo_baseline,
        # "sade": pygmo_baseline,
        # "sga": pygmo_baseline,
        # "simulated_annealing": pygmo_baseline,
        # "pso": pygmo_baseline,
        # "gaco": pygmo_baseline,
        # "bee_colony": pygmo_baseline,
        # "uniform": uniform_variables_values_vector,
        # "gaussian_cma_es": gaussian_variables_values_vector,
        "cnmcts": cnmcts,
        # "cnrpa": cnrpa,
        # "cgnrpa": cgnrpa,
        # "crbnmcts": crbnmcts,
        # "cabgnrpa": cabgnrpa,
        # "genetic": genetic_algorithm,
    }

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
        "crbnmcts": "cRbNMCS",
        "cabgnrpa": "cABGNRPA",
        "genetic": "Custom Genetic Algorithm",
    }

    trajectory = Trajectory()
    udp = pk.trajopt.gym.cassini1
    trajectory.instantiate("Earth", "Saturn")
    trajectory.variables["planets_sequence"] = [
        "Earth",
        "Venus",
        "Venus",
        "Earth",
        "Jupiter",
        "Saturn",
    ]
    trajectory.set_variables_bounds()

    # General input values
    inputs_values = {
        "evaluator": udp,
        "planets_sequence": trajectory.variables["planets_sequence"],
        "values_sequence": list(),
        "bounds": trajectory.bounds,
        "n_iterations": 2500,  # for uniform and gaussian sampling
        "level": 2,  # for cNMCTS, cNRPA and derivatives
        "bandwidth": 200,  # for cNMCTS
        "n_policies": 500,  # for cNRPA and derivatives
        "multiple_values_policy": True,  # for cNRPA and derivatives
        "learning_rate": 0.05,  # for cNRPA and derivatives
        "tau": 10,  # for cGNRPA and derivatives
        "gamma": 0.2,  # for cABGNRPA
        "n_generations": 1000,  # for Genetic Algorithm
        "population_size": 500,  # for Genetic Algorithm
        "mutation_probability": 0.15,  # for Genetic Algorithm
        "mutation_effect": 0.25,  # for Genetic Algorithm
    }

    # Timeouts list
    timeouts_list = [1, 5, 10, 30, 60, 300, 600]

    # Specific parameters
    algorithm_parameters = {
        "cmaes": {"solver": "cmaes"},
        "sade": {"solver": "sade"},
        "sga": {"solver": "sga"},
        "pso": {"solver": "pso"},
        "gaco": {"solver": "gaco"},
        "bee_colony": {"solver": "bee_colony"},
        "simulated_annealing": {"solver": "simulated_annealing"},
        "uniform": {"n_iterations": 10000},
        "gaussian_cma_es": {"n_iterations": 1500},
        "cnmcts": {"level": 2, "bandwidth": 200},
        "crbnmcts": {"level": 2, "bandwidth": 200},
        "cnrpa": {"level": 2, "n_policies": 1000, "multiple_values_policy": True},
        "cgnrpa": {
            "level": 2,
            "n_policies": 300,
            "multiple_values_policy": True,
            "tau": 10,
        },
        "cabgnrpa": {
            "level": 2,
            "n_policies": 100,
            "multiple_values_policy": True,
            "tau": 10,
            "gamma": 0.7,
        },
        "genetic": {"n_generations": 5000, "mutation_probability": 0.15},
    }

    for algorithm, function in SAMPLING_FUNCTIONS.items():
        print("\nFor " + algorithm + ":")
        specific_input_values = deepcopy(inputs_values)
        for arg, val in algorithm_parameters[algorithm].items():
            specific_input_values[arg] = val
        specific_input_values["timeout"] = timeouts_list[2]
        values_sequence, delta_v, best_values_list, time_list = function(**specific_input_values)
        if delta_v < UNFEASIBILITY_VALUE:
            print(f"Delta V: {np.linalg.norm(delta_v) / 1000:.3f} km/s")
            plt.plot(time_list, best_values_list)
            plt.show()
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
        else:
            print("No valid solution produced")
