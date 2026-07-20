import time
from copy import deepcopy
import click
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pykep as pk
import pygmo as pg
from solvers.continuous_variables_choice.cnmcts import cnmcts
from solvers.continuous_variables_choice.cnrpa import cnrpa
from solvers.continuous_variables_choice.cgnrpa import cgnrpa
from solvers.continuous_variables_choice.cabgnrpa import cabgnrpa
from solvers.continuous_variables_choice.gaco_cabgnrpa import gaco_cabgnrpa
from solvers.optimization_learning.optimizer_cnrpa import optimizer_cnrpa
from solvers.continuous_variables_choice.baselines import *
from utils.constants import UNFEASIBILITY_VALUE, DV_LAUNCHER


@click.command()
@click.option("--save_file", type=str)
@click.option("--best_sequence_file", type=str)
def main(save_file, best_sequence_file):

    print("Testing Cassini MGA:\n")

    SAMPLING_FUNCTIONS = {
        "cmaes": pygmo_baseline,
        "sade": pygmo_baseline,
        "sga": pygmo_baseline,
        "simulated_annealing": pygmo_baseline,
        "pso": pygmo_baseline,
        "gaco": pygmo_baseline,
        "bee_colony": pygmo_baseline,
        "uniform": uniform_variables_values_vector,
        "cnmcts": cnmcts,
        "cnrpa": cnrpa,
        "cgnrpa": cgnrpa,
        "cabgnrpa": cabgnrpa,
        "gaco_cabgnrpa": gaco_cabgnrpa,
        "optimizer_cnrpa": optimizer_cnrpa,
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
        "cnrpa": "GcNRPA",
        "cgnrpa": "GcGNRPA",
        "cabgnrpa": "GcABGNRPA",
        "gaco_cabgnrpa": "GACO GcABGNRA",
        "optimizer_cnrpa": "Optimizer GcNRPA",
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
        "planets_sequence": planets_sequence,
        "bounds": bounds,
        "n_iterations": 2500,  # for uniform and gaussian sampling
        "n_generations": 1000,  # for Genetic Algorithm
        "population_size": 500,  # for Genetic Algorithm
        "mutation_probability": 0.15,  # for Genetic Algorithm
        "mutation_effect": 0.25,  # for Genetic Algorithm
        "timeout": 300,  # 5min is maximum execution time
    }

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
        "cnmcts": {"level": 2, "bandwidth": 672},
        "cnrpa": {
            "level": 2,
            "n_policies": 33155,
            "learning_rate": 0.011047321079505495,
        },
        "cgnrpa": {
            "level": 3,
            "n_policies": 17908,
            "learning_rate": 0.154937318595788,
            "tau": 2.545153120408655,
        },
        "cabgnrpa": {
            "level": 1,
            "n_policies": 15930,
            "learning_rate": 0.034981412158836295,
            "tau": 1.3316892372714013,
            "gamma": 0.02938863839421489,
        },
        "gaco_cabgnrpa": {
            "level": 3,
            "n_policies": 29325,
            "learning_rate": 0.09447145922848464,
            "tau": 0.38210537834931463,
            "gamma": 0.13056072395524418,
            "zeta": 0.7086734675139809,
            "elitism_factor": 0.020372968801098613,
        },
        "optimizer_cnrpa": {
            "level": 3,
            "n_policies": 297,
            "learning_rate": 0.5744087681488131,
            "archive_size": 81,
            "max_steps": 50,
            "movement_range": 0.00064558763425047,
            "initial_state_strategy": "mixed",
            "score_type": "differences_sum",
        },
    }

    # Iterations
    n_iterations = 100

    # Storing results
    results = np.zeros((n_iterations, len(algorithms_list)))
    best_sequences = dict()

    for algorithm_id, algorithm in enumerate(algorithms_list):
        print("\nFor " + algorithm + ":")
        specific_input_values = deepcopy(inputs_values)
        best_delta_v = UNFEASIBILITY_VALUE
        best_sequence = None
        for arg, val in algorithm_parameters[algorithm].items():
            specific_input_values[arg] = val
        for iteration in range(n_iterations):
            specific_input_values["evaluator"] = deepcopy(udp)
            values_sequence, delta_v, _, _ = SAMPLING_FUNCTIONS[algorithm](
                **specific_input_values
            )
            results[iteration, algorithm_id] = min(delta_v, UNFEASIBILITY_VALUE)
            if delta_v < best_delta_v:
                best_delta_v = delta_v
                best_sequence = values_sequence
            if (iteration + 1) % 10 == 0:
                print(f"{iteration + 1} iterations done")
        best_sequences[algorithm] = {
            "best_sequence": list(best_sequence),
            "delta_v": best_delta_v / 1000,
        }

    # Storing all runs
    dataframe = pd.DataFrame(
        results, index=list(range(n_iterations)), columns=algorithms_list
    )
    dataframe.to_excel(save_file)

    # Storing best sequences
    with open(best_sequence_file, "w") as file:
        json.dump(best_sequences, file)


if __name__ == "__main__":
    main()
