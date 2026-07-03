"""
Verifying solvers capabilities
"""

import time
import matplotlib.pyplot as plt
from copy import deepcopy
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

SAMPLING_FUNCTIONS = {
    # "cmaes": cmaes_pygmo,
    # "sade": sade_pygmo,
    # "sga": genetic_pygmo,
    # "simulated_annealing": sa_pygmo,
    # "pso": pso_gen_pygmo,
    # "gaco": gaco_pygmo,
    # "bee_colony": bee_colony_pygmo,
    # "uniform": uniform_variables_values_vector,
    # "gaussian_cma_es": gaussian_variables_values_vector,
    # "cnmcts": cnmcts,
    # "cnrpa": cnrpa,
    # "cgnrpa": cgnrpa,
    # "crbnmcts": crbnmcts,
    "cabgnrpa": cabgnrpa,
    "genetic": genetic_algorithm,
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

if __name__ == "__main__":
    print("Testing Cassini MGA:\n")
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
    inputs_values = {
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
    algorithm_parameters = {
        "cmaes": {"n_iterations": 1500},
        "sade": {"n_iterations": 1500},
        "sga": {"n_generations": 5000, "mutation_probability": 0.15},
        "simulated_annealing": {},
        "pso": {"n_generations": 5000},
        "gaco": {"n_iterations": 1500, "kernel_size": 10, "learning_rate": 0.01},
        "bee_colony": {"n_generations": 5000},
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
        start_time = time.time()
        values_sequence, delta_v = function(**specific_input_values)
        end_time = time.time() - start_time
        if delta_v < UNFEASIBILITY_VALUE:
            # print(f"Delta V: {delta_v / 1000:.3f} km/s")
            print(f"Delta V: {udp.fitness(values_sequence)[0] / 1000:.3f} km/s")
            print(f"Total execution time: {end_time:.3f} s")
            axe = udp.plot(values_sequence, figsize=(20, 20))
            figure = axe.figure
            axe.view_init(90, 0)
            axe.axis("off")
            axe.set_title(
                algorithm.upper()
                + r": $\Delta$V = "
                + f"{udp.fitness(values_sequence)[0] / 1000:.3f} km/s"
            )
            figure.savefig(
                "./Cassini trials/" + SAMPLING_FUNCTIONS_NAMES[algorithm] + ".png"
            )
            plt.close(figure)
            # plt.show()
        else:
            print("No valid solution produced")
