# Author: Lotfi Kobrosly (2025)

# Constants used accross the repository

# Universal gravitational constant
GRAVITY_CONSTANT = 6.67430e-11
EARTH_GRAVITY_CONSTANT = 9.80665;

# Solar System Planets compatible with JPL data in pykep
PLANETS = [
    "MERCURY"
    "VENUS",
    "EARTH",
    "MARS",
    "JUPITER",
    "SATURN",
    "URANUS",
    "NEPTUNE",
]

# Times of flight bounds in days from Claude estimations
TOF_BOUNDS = {
    ("mercury", "venus")   : (60,    250),
    ("mercury", "earth")   : (90,    350),
    ("mercury", "mars")    : (150,   500),
    ("mercury", "jupiter") : (350,   2000),
    ("mercury", "saturn")  : (900,   4500),
    ("venus",   "earth")   : (60,    250),
    ("venus",   "mars")    : (100,   400),
    ("venus",   "jupiter") : (400,   2500),
    ("venus",   "saturn")  : (900,   4500),
    ("earth",   "mars")    : (130,   500),
    ("earth",   "jupiter") : (400,   2500),
    ("earth",   "saturn")  : (1000,  5000),
    ("earth",   "uranus")  : (2000,  12000),
    ("earth",   "neptune") : (4000,  16000),
    ("mars",    "jupiter") : (300,   2000),
    ("mars",    "saturn")  : (800,   4000),
    ("mars",    "uranus")  : (2000,  10000),
    ("mars",    "neptune") : (3500,  14000),
    ("jupiter", "saturn")  : (500,   2500),
    ("jupiter", "uranus")  : (1000,  5000),
    ("jupiter", "neptune") : (2000,  8000),
    ("saturn",  "uranus")  : (1000,  6000),
    ("saturn",  "neptune") : (1500,  6000),
    ("uranus",  "neptune") : (2000,  8000),
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
