# Author: Lotfi Kobrosly (2025)

# Constants used accross the repository

# Bodies' masses in kilograms (Wikipedia)
MASSES = {

    "Sun" : 1.9884e30,

    # Planets
    "Mercury" : 3.301e23,
    "Venus" : 4.867e24,
    "Earth" : 5.972e24,
    "Mars" : 6.417e23,
    "Jupiter" : 1.899e27,
    "Saturn" : 5.685e26,
    "Uranus" : 8.682e25,
    "Neptune" : 1.024e26,

    # Dwarf planets and asteroids
    "Pluto" : 1.471e22,
    "Ceres" : 9.3e20,
    "Vesta" : 2.6e20,
    "Pallas" : 2.0e20,

    # Planets' satellites
    "Moon" : 7.348e22,
    "Io" : 8.93e22,
    "Europa" : 4.80e22,
    "Ganymede" : 1.48e23,
    "Callisto" : 1.08e23,
    "Titan" : 1.35e23,
    "Titania" : 3.52e21,
    "Oberon" : 3.01e21,
    "Triton" : 2.14e22,
}

# Heliocentric indices
INDICES = {
    "Sun": (0,0),
    "Mercury" : (0,1),
    "Venus" : (0,2),
    "Earth" : (0,3),
    "Mars" : (0,4),
    "Jupiter" : (0,5),
    "Saturn" : (0,6),
    "Uranus" : (0,7),
    "Neptune" : (0,8),
    "Pluto": (0, 9)
}

