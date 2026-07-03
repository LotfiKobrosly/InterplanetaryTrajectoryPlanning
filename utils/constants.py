# Author: Lotfi Kobrosly (2025)

import numpy as np

# Constants used accross the repository
RANDOM_SEED = 1
RANDOM_GENERATOR = np.random.default_rng(seed=RANDOM_SEED)
GAUSSIAN_KERNEL_THRESHOLD = 0.05

# Solar System Planets compatible with JPL data in pykep
PLANETS = [
    "Mercury",
    "Venus",
    "Earth",
    "Mars",
    "Jupiter",
    "Saturn",
    "Uranus",
    "Neptune",
]

ORBITAL_PERIODS = {
    "Mercury": 87.97,
    "Venus": 224.7,
    "Earth": 365.26,
    "Mars": 686.98,
    "Jupiter": 4332.82,
    "Saturn": 10755.7,
    "Uranus": 30687.15,
    "Neptune": 60190.03,
}

# Times of flight bounds in days from Claude estimations
TOF_BOUNDS = {
    ("Mercury", "Mercury"): (
        ORBITAL_PERIODS["Mercury"] // 4,
        ORBITAL_PERIODS["Mercury"] * 4,
    ),
    ("Venus", "Venus"): (ORBITAL_PERIODS["Venus"] // 4, ORBITAL_PERIODS["Venus"] * 4),
    ("Earth", "Earth"): (ORBITAL_PERIODS["Earth"] // 4, ORBITAL_PERIODS["Earth"] * 4),
    ("Mercury", "Venus"): (60, 250),
    ("Mercury", "Earth"): (90, 350),
    ("Mercury", "Mars"): (150, 500),
    ("Mercury", "Jupiter"): (350, 2000),
    ("Mercury", "Saturn"): (900, 4500),
    ("Venus", "Earth"): (60, 250),
    ("Venus", "Mars"): (100, 400),
    ("Venus", "Jupiter"): (400, 2500),
    ("Venus", "Saturn"): (900, 4500),
    ("Earth", "Mars"): (130, 500),
    ("Earth", "Jupiter"): (400, 2500),
    ("Earth", "Saturn"): (1000, 5000),
    ("Earth", "Uranus"): (2000, 12000),
    ("Earth", "Neptune"): (4000, 16000),
    ("Mars", "Jupiter"): (300, 2000),
    ("Mars", "Saturn"): (800, 4000),
    ("Mars", "Uranus"): (2000, 10000),
    ("Mars", "Neptune"): (3500, 14000),
    ("Jupiter", "Saturn"): (500, 2500),
    ("Jupiter", "Uranus"): (1000, 5000),
    ("Jupiter", "Neptune"): (2000, 8000),
    ("Saturn", "Uranus"): (1000, 6000),
    ("Saturn", "Neptune"): (1500, 6000),
    ("Uranus", "Neptune"): (2000, 8000),
}

# Heliocentric indices
INDICES = {
    "Sun": (0, 0),
    "Mercury": (0, 1),
    "Venus": (0, 2),
    "Earth": (0, 3),
    "Mars": (0, 4),
    "Jupiter": (0, 5),
    "Saturn": (0, 6),
    "Uranus": (0, 7),
    "Neptune": (0, 8),
    "Pluto": (0, 9),
}
PLANETS_SEMI_MAJOR_AXIS_AND_EXCENTRICITY = {
    "mercury": (0.387, 0.2056),
    "venus": (0.723, 0.0067),
    "earth": (1.000, 0.0167),
    "mars": (1.524, 0.0934),
    "jupiter": (5.203, 0.0489),
    "saturn": (9.537, 0.0565),
    "uranus": (19.19, 0.0463),
    "neptune": (30.07, 0.0097),
}


# Spacecraft constants
MAXIMUM_THRUST = 0.1
ISP = 3000

# Variables bounds
VARIABLES_BOUNDS = {
    "departure_epoch": (100, 2000),  # departure window (mjd2000)
}
DV_LAUNCHER = 10000  # m/s, generous launcher capability
SAFE_RADIUS_FACTOR = 1.05
UNFEASIBILITY_VALUE = 1e20
VELOCITY_NORMALIZING_FACTOR = 40000

# Algorithms
POLICY_ALGORITHMS = ["cnrpa", "cgnrpa", "nrpa", "gnrpa"]
SAMPLING_FUNCTIONS = [
    "sade",
    "cmaes",
    "sga",
    "simulatedd_annealing",
    "pso",
    "gaco",
    "bee_colony",
    "uniform",
    "gaussian_cma_es",
    "genetic",
    "cnmcts",
    "crbnmcts",
    "cnrpa",
    "cgnrpa",
    "cabgnrpa",
]

# cGNRPA
N_CANDIDATES = 25
