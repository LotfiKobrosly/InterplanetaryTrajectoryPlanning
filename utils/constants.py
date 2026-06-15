# Author: Lotfi Kobrosly (2025)

import numpy as np

# Constants used accross the repository

# Universal gravitational constant
GRAVITY_CONSTANT = 6.67430e-11
EARTH_GRAVITY_CONSTANT = 9.80665

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

# Times of flight bounds in days from Claude estimations
TOF_BOUNDS = {
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

# Bodies' masses in kilograms (Wikipedia)
MASSES = {
    "Sun": 1.9884e30,
    # Planets
    "Mercury": 3.301e23,
    "Venus": 4.867e24,
    "Earth": 5.972e24,
    "Mars": 6.417e23,
    "Jupiter": 1.899e27,
    "Saturn": 5.685e26,
    "Uranus": 8.682e25,
    "Neptune": 1.024e26,
    # Dwarf planets and asteroids
    "Pluto": 1.471e22,
    "Ceres": 9.3e20,
    "Vesta": 2.6e20,
    "Pallas": 2.0e20,
    # Planets' satellites
    "Moon": 7.348e22,
    "Io": 8.93e22,
    "Europa": 4.80e22,
    "Ganymede": 1.48e23,
    "Callisto": 1.08e23,
    "Titan": 1.35e23,
    "Titania": 3.52e21,
    "Oberon": 3.01e21,
    "Triton": 2.14e22,
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

# Spacecraft constants
MAXIMUM_THRUST = 0.1
ISP = 3000

# Variables bounds
VARIABLES_BOUNDS = {
    "departure_epoch": (1000, 2000),  # departure window (mjd2000)
    "departure_velocity": (0, 10000),  # m/s, generous launcher capability
    "planet_arrival_radius": (1.05, 10),  # must be > 1 (above planet surface)
    "planet_arrival_angle": (0, 2 * np.pi),
}

# Algorithms
POLICY_ALGORITHMS = ["cnrpa", "cgnrpa", "nrpa", "gnrpa"]
SAMPLING_FUNCTIONS = ["uniform", "gaussian_cma_es", "cnrpa", "cnmcts"]
SEQUENCE_FUNCTIONS = ["cnrpa", "cnmcts"]
VECTOR_FUNCTIONS = ["gaussian_cma_es", "uniform"]
