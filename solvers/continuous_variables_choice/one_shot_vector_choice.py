"""
Implements a one-shot choice of the control variables in a vector.

For the gaussian-based implementation, we update the means and covariance with the CMA-ES (Covariance Matrix Adaptation Evolution Strategy)

In this file, the values sequence is represented as follows:
[departure_epoch, time_of_flight_0,..., time_of_flight_n, (arrival_radius_0, arrival_angle_0),..., (arrival_radius_(n-1), arrival_angle(n-1))]
"""

import warnings
import numpy as np
import pykep as pk
import pygmo as pg
import cma
from utils.basic_functions import normalize, denormalize
from utils.constants import RANDOM_GENERATOR, DV_LAUNCHER

warnings.filterwarnings("ignore")


def uniform_variables_values_vector(
    bounds: list, planets_sequence: list, n_iterations: int = 500, *args, **kwargs
):
    """
    Returns a uniformly sampled vector from the given bounds
    """
    best_vector = None
    best_value = 1e30
    planets_sequence = [pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence]
    evaluator = pk.trajopt.mga(
        planets_sequence,
        list(bounds[0]),
        [list(element) for element in bounds[1:]],
        vinf=DV_LAUNCHER
    )
    for iteration in range(n_iterations):
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
    if best_vector is None:
        best_vector = vector
    return (
        best_vector,
        evaluator.fitness(best_vector)[0],
    )


def gaussian_variables_values_vector(
    bounds: list, planets_sequence: list, n_iterations: int = 500, *args, **kwargs
):
    """
    Returns a sample of a fitted normal distribution through a CMA-ES.
    """
    lower_bounds = np.array([bound[0] for bound in bounds])
    upper_bounds = np.array([bound[1] for bound in bounds])
    planets_sequence = [pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence]
    evaluator = pk.trajopt.mga(
        planets_sequence,
        list(bounds[0]),
        [list(element) for element in bounds[1:]],
        vinf=DV_LAUNCHER
    )

    def objective_function(normalized_vector: np.ndarray):
        result = evaluator.fitness(denormalize(normalized_vector, lower_bounds, upper_bounds))
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
    estimator = cma.CMAEvolutionStrategy(initial_input, sigma0, options)

    while not estimator.stop():
        candidates = estimator.ask()  # sample λ candidates
        fitnesses = [objective_function(x) for x in candidates]
        estimator.tell(candidates, fitnesses)  # update means, covariance, sigma
        estimator.logger.add()
        estimator.disp()

    result_normalized = estimator.result.xbest
    best_vector = denormalize(result_normalized, lower_bounds, upper_bounds)
    return (
        best_vector,
        evaluator.fitness(best_vector)[0],
    )

def cmaes_pygmo(
    bounds: list, planets_sequence: list, n_iterations: int = 1500, *args, **kwargs
):
    """
    Use the PyGMO built-in Covariance Matrix Adaptation Evolution Strategy
    """
    planets_sequence = [pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence]
    evaluator = pk.trajopt.mga(
        planets_sequence,
        list(bounds[0]),
        [list(element) for element in bounds[1:]],
        vinf=DV_LAUNCHER
    )

    # Setting up
    problem = pg.problem(evaluator)
    solver = pg.cmaes(n_iterations, force_bounds=True, sigma0=0.5, ftol=1e-4)
    algorithm = pg.algorithm(solver)
    result = list()

    # Running algorithm
    for i in range(10):
        pop = pg.population(problem, 20)
        pop = algorithm.evolve(pop)
        result.append([pop.champion_f, pop.champion_x])
        # print(i, pop.champion_f[0], flush=True)
        
    best_value = sorted(result, key =  lambda x: x[0][0])[0][0][0]
    best_x = sorted(result, key =  lambda x: x[0][0])[0][1]

    return best_x, best_value

def sade_pygmo(
    bounds: list, planets_sequence: list, n_iterations: int = 2500, *args, **kwargs
):
    """
    Use the PyGMO built-in Self Adaptive Differential Evolution
    """
    planets_sequence = [pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence]
    evaluator = pk.trajopt.mga(
        planets_sequence,
        list(bounds[0]),
        [list(element) for element in bounds[1:]],
        vinf=DV_LAUNCHER
    )
    # Setting up
    problem = pg.problem(evaluator)
    solver = pg.sade(n_iterations, ftol=1e-4, xtol=1e-4)
    algorithm = pg.algorithm(solver)
    result = list()

    # Running algorithm
    for i in range(10):
        pop = pg.population(problem, 20)
        pop = algorithm.evolve(pop)
        result.append([pop.champion_f, pop.champion_x])
        # print(i, pop.champion_f[0], flush=True)
        
    best_value = sorted(result, key =  lambda x: x[0][0])[0][0][0]
    best_x = sorted(result, key =  lambda x: x[0][0])[0][1]

    return best_x, best_value


def genetic_algorithm(
    planets_sequence: list = None,
    population: list = None,
    bounds: list = None,
    n_generations: int = 1000,
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
                chromosome[element_index] = RANDOM_GENERATOR.uniform(1 - perturbation, 1 + perturbation) * element
        return chromosome

    planets_sequence = [pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence]
    evaluator = pk.trajopt.mga(
        planets_sequence,
        list(bounds[0]),
        [list(element) for element in bounds[1:]],
        vinf=DV_LAUNCHER
    )
    chromosome_length = len(planets_sequence)
    low_bound = [bound[0] for bound in bounds]
    high_bound = [bound[1] for bound in bounds]
    
    # Initialization
    if population is None:
        population = RANDOM_GENERATOR.uniform(
            low_bound,
            high_bound,
            (population_size, chromosome_length),
        )

    for _ in range(n_generations):
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
                gene_pool[sequence_index] = mutation(sequence, mutation_probability, mutation_effect)

        # Fitness evaluation
        fitness = [
            evaluator.fitness(sequence)[0]
            for sequence in gene_pool
        ]

        # Selection
        indices = np.argsort(fitness)
        gene_pool = np.array(gene_pool)[indices]
        population = gene_pool[:population_size]

    return (
        population[0],
        evaluator.fitness(population[0])[0]
    )
