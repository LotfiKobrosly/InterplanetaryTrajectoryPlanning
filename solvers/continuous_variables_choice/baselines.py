"""
Implements a one-shot choice of the control variables in a vector.

For the gaussian-based implementation, we update the means and covariance with the CMA-ES (Covariance Matrix Adaptation Evolution Strategy)

In this file, the values sequence is represented as follows:
[departure_epoch, time_of_flight_0,..., time_of_flight_n, (arrival_radius_0, arrival_angle_0),..., (arrival_radius_(n-1), arrival_angle(n-1))]
"""

import time
import warnings
from copy import deepcopy
import numpy as np
import pykep as pk
import pygmo as pg
import matplotlib.pyplot as plt
import cma
from utils.basic_functions import normalize, denormalize
from utils.constants import RANDOM_GENERATOR, DV_LAUNCHER, UNFEASIBILITY_VALUE

warnings.filterwarnings("ignore")

PYGMO_SOLVERS = {
    "bee_colony": {
        "function": pg.bee_colony,
        "solver_parameters": {"gen": 1500, "limit": 20},
    },
    "cmaes": {
        "function": pg.cmaes,
        "solver_parameters": {
            "gen": 1500,
            "force_bounds": True,
            "sigma0": 0.5,
            "ftol": 1e-4,
            "memory": True,
        },
    },
    "sade": {
        "function": pg.sade,
        "solver_parameters": {"gen": 1500, "ftol": 1e-4, "xtol": 1e-4, "memory": True},
    },
    "sga": {
        "function": pg.sga,
        "solver_parameters": {"gen": 200, "cr": 0.95, "m": 0.15},
    },
    "simulated_annealing": {
        "function": pg.simulated_annealing,
        "solver_parameters": {"Ts": 10.0, "Tf": 1e-5, "n_T_adj": 100},
    },
    "pso": {
        "function": pg.pso_gen,
        "solver_parameters": {"gen": 50, "memory": True},
    },
    "gaco": {
        "function": pg.gaco,
        "solver_parameters": {
            "gen": 100,
            "ker": 20,
            "q": 0.05,
            "oracle": 1e9,
            "memory": True,
        },
    },
}


def uniform_variables_values_vector(
    bounds: list, evaluator: pk.trajopt.mga, timeout: float = 10, *args, **kwargs
):
    """
    Returns a uniformly sampled vector from the given bounds
    """
    best_vector = None
    best_value = UNFEASIBILITY_VALUE
    best_values_list, time_list = list(), list()
    start_time = time.time()
    current_time = time.time() - start_time
    while current_time < timeout:
        # Imposing a strict while loop here without the max iterations condition might
        # result in an endless loop, given the pure random aspect of this algorithm
        vector = RANDOM_GENERATOR.uniform(
            np.array([bound[0] for bound in bounds]),
            np.array([bound[1] for bound in bounds]),
        )
        result = evaluator.fitness(vector)
        # print(result)

        if not (result is None) and (result[0] < best_value):
            best_value = result[0]
            best_vector = vector
        current_time = time.time() - start_time
        if best_value < UNFEASIBILITY_VALUE:
            best_values_list.append(best_value)
            time_list.append(current_time)
    if best_vector is None:
        best_vector = vector
    return (
        best_vector,
        best_value,
        best_values_list,
        time_list,
    )


def gaussian_variables_values_vector(
    bounds: list,
    evaluator: pk.trajopt.mga,
    timeout: float = 10,
    n_iterations: int = 5e10,
    *args,
    **kwargs,
):
    """
    Returns a sample of a fitted normal distribution through a CMA-ES.
    """
    lower_bounds = np.array([bound[0] for bound in bounds])
    upper_bounds = np.array([bound[1] for bound in bounds])

    def objective_function(normalized_vector: np.ndarray):
        result = evaluator.fitness(
            denormalize(normalized_vector, lower_bounds, upper_bounds)
        )
        if result is None:
            return 1e10  # penalty for invalid trajectory

        return result[0] / 1000

    # --- Initial mean: center of search space (normalized) ---
    initial_input = [0.95] * len(bounds)
    sigma0 = 10  # initial step size (in normalized space)

    # --- CMA-ES options ---
    options = cma.CMAOptions()
    options["bounds"] = [[0.0] * len(bounds), [1.0] * len(bounds)]  # normalized bounds
    options["maxiter"] = n_iterations
    options["popsize"] = 50  # λ : samples per iteration
    # options['CMA_diagonal'] = True
    options["tolx"] = 1e-6  # convergence on x
    options["tolfun"] = 1e-6  # convergence on f
    options["verbose"] = -9  # this value means no verboses

    # --- Run ---
    best_values_list, time_list = list(), list()
    estimator = cma.CMAEvolutionStrategy(initial_input, sigma0, options)
    start_time = time.time()
    current_time = time.time() - start_time
    while current_time < timeout:
        candidates = estimator.ask()  # sample λ candidates
        fitnesses = [objective_function(x) for x in candidates]
        estimator.tell(candidates, fitnesses)  # update means, covariance, sigma
        estimator.logger.add()
        estimator.disp()

        result_normalized = estimator.result.xbest
        best_vector = denormalize(result_normalized, lower_bounds, upper_bounds)
        best_value = evaluator.fitness(best_vector)[0]
        current_time = time.time() - start_time
        if best_value < UNFEASIBILITY_VALUE:
            best_values_list.append(best_value)
            time_list.append(current_time)

    return (
        best_vector,
        best_value,
        best_values_list,
        time_list,
    )


def pygmo_baseline(
    evaluator: pk.trajopt.mga,
    timeout: float = 10,
    solver: str = "sga",
    population_size: int = 100,
    *args,
    **kwargs,
):

    # Setting up
    problem = pg.problem(evaluator)
    solver = PYGMO_SOLVERS[solver]["function"](
        **PYGMO_SOLVERS[solver]["solver_parameters"],
    )
    algorithm = pg.algorithm(solver)
    time_list, best_values_list = list(), list()
    best_x = None
    best_value = UNFEASIBILITY_VALUE
    pop = pg.population(problem, population_size)
    start_time = time.time()
    current_time = time.time() - start_time

    # Running algorithm
    while current_time < timeout:
        pop = algorithm.evolve(pop)
        if np.linalg.norm(pop.champion_f) < best_value:
            best_x = deepcopy(pop.champion_x)
            best_value = np.linalg.norm(pop.champion_f)
        if best_value < UNFEASIBILITY_VALUE:
            best_values_list.append(best_value)
        current_time = time.time() - start_time
        time_list.append(current_time)

    return best_x, best_value, best_values_list, time_list


def genetic_algorithm(
    evaluator: pk.trajopt.mga,
    population: list = None,
    bounds: list = None,
    timeout: float = 10,
    population_size: int = 1000,
    mutation_probability: float = 0.1,
    mutation_effect: float = 0.2,
    *args,
    **kwargs,
):

    def crossover(parent_1: list, parent_2: list):
        child = list()
        for counter in range(len(parent_1)):
            if RANDOM_GENERATOR.uniform(0, 1) > 0.5:
                child.append(parent_1[counter])
            else:
                child.append(parent_2[counter])
        return child

    def mutation(chromosome: list, probability: float, perturbation: float):
        for element_index, element in enumerate(chromosome):
            if RANDOM_GENERATOR.uniform(0, 1) < probability:
                chromosome[element_index] = (
                    RANDOM_GENERATOR.uniform(1 - perturbation, 1 + perturbation)
                    * element
                )
        return chromosome

    chromosome_length = len(bounds)
    low_bound = [bound[0] for bound in bounds]
    high_bound = [bound[1] for bound in bounds]

    time_list, best_values_list = list(), list()
    best_x = None
    best_value = UNFEASIBILITY_VALUE
    start_time = time.time()
    current_time = time.time() - start_time

    # Initialization
    if population is None:
        population = RANDOM_GENERATOR.uniform(
            low_bound,
            high_bound,
            (population_size, chromosome_length),
        )

    while current_time < timeout:
        # Build new generation
        RANDOM_GENERATOR.shuffle(population)

        ## Crossovers
        descendants = [
            crossover(population[counter], population[counter + 1])
            for counter in range(0, population_size, 2)
        ]

        gene_pool = np.concatenate([population, descendants])

        ## Mutations
        for sequence_index, sequence in enumerate(gene_pool):
            if RANDOM_GENERATOR.uniform(0, 1) < mutation_probability:
                gene_pool[sequence_index] = mutation(
                    sequence, mutation_probability, mutation_effect
                )

        # Fitness evaluation
        fitness = [evaluator.fitness(sequence)[0] for sequence in gene_pool]

        # Selection
        indices = np.argsort(fitness)
        gene_pool = np.array(gene_pool)[indices]
        population = gene_pool[:population_size]

        # Store instantaneous results
        current_time = time.time() - start_time
        if evaluator.fitness(population[0])[0] < best_value:
            best_x = population[0]
            best_value = evaluator.fitness(best_x)[0]
        best_values_list.append(best_value)
        time_list.append(current_time)

    return best_x, best_value, best_values_list, time_list


if __name__ == "__main__":
    # Problem
    udp = pk.trajopt.gym.cassini2

    # Variables bounds
    bounds = [
        (low_bound, high_bound)
        for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    ]

    # General input values
    inputs_values = {
        "evaluator": deepcopy(udp),
        "bounds": bounds,
        "solver": "bee_colony",
        "timeout": 10,
        "population_size": 50,
    }

    # Baseline
    values_sequence, best_value, values_list, time_list = pygmo_baseline(
        **inputs_values
    )

    print(f"Best Delta V for {inputs_values["solver"]}: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
    # print(udp.pretty(values_sequence))

    axe = udp.plot(values_sequence, figsize=(20, 20))
    # figure = axe.figure
    axe.view_init(90, 0)
    axe.axis("off")
    axe.set_title(
        inputs_values["solver"].upper()
        + r": $\Delta$V = "
        + f"{best_value / 1000:.3f} km/s"
    )
    plt.show()
