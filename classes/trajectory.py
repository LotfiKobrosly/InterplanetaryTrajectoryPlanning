"""
Implements a class for trajectory building
"""

from copy import deepcopy
import numpy as np
import pykep as pk
from utils.constants import PLANETS, POLICY_ALGORITHMS, VARIABLES_BOUNDS, TOF_BOUNDS
from utils.trajectory_evaluation import evaluate_mga_trajectory


class Trajectory:

    def __init__(self):
        self.variables = {
            "planets_sequence": list(),  # Length msut be >= 2. If 2, a direct flight to goal is considered
            "departure_epoch": None,
            "time_of_flights_list": list(),  # Length must be len(plent_sequence) - 1
            "planets_flyby_parameters": list(),  # Length must be len(plent_sequence) - 2
        }
        self.start = None
        self.goal = None


    def set_start_planet(self, planet="Earth"):
        assert planet in PLANETS, "Given planet not in the provided list: " + str(
            PLANETS
        )
        self.start = planet
        self.variables["planets_sequence"].append(planet)

    def set_goal_planet(self, planet="Mars"):
        assert planet in PLANETS, "Given planet not in the provided list: " + str(
            PLANETS
        )
        self.goal = planet

    def define_planets_pool(self):
        indices = (PLANETS.index(self.start), PLANETS.index(self.goal))
        self.planets_pool = PLANETS[: max(indices) + 1]
        self.planets_pool.remove(self.start)

    def instantiate(self, start_planet: str, goal_planet: str):
        self.set_start_planet(start_planet)
        self.set_goal_planet(goal_planet)
        self.define_planets_pool()

    def reinitialize(self):
        self.variables = {
            "planets_sequence": [self.start],  # Length msut be >= 2. If 2, a direct flight to goal is considered
            "departure_epoch": None,
            "time_of_flights_list": list(),  # Length must be len(plent_sequence) - 1
            "planets_flyby_parameters": list(),  # Length must be len(planet_sequence) - 2
        }
        self.define_planets_pool()

    def planets_sequence_is_set(self):
        return self.variables["planets_sequence"] and (
            self.variables["planets_sequence"][-1] == self.goal
        )

    def is_terminal(self):
        return (
            self.planets_sequence_is_set()
            and not (self.variables["departure_epoch"] is None)
            and len(self.variables["time_of_flights_list"])
            == len(self.variables["planets_sequence"] - 1)
            and len(self.variables["planets_flyby_parameters"])
            == len(self.variables["planets_sequence"] - 2)
        )

    def evaluate_mga(self):
        self.mga_results = evaluate_mga_trajectory(**self.variables)
        if self.mga_results is None:
            return 1e30  # Invalid sequence
        return self.mga_results[0] / 1000  # Convert from m/s to km/s

    def add_planet(self, planet: str):
        self.variables["planets_sequence"].append(planet)
        self.planets_pool.remove(planet)

    def set_variables_bounds(self):
        assert self.planets_sequence_is_set(), "Planets' sequence not set yet"
        vector_bounds = [VARIABLES_BOUNDS["departure_epoch"]]
        sequence_bounds = [VARIABLES_BOUNDS["departure_epoch"]]
        for planet_id, planet in enumerate(self.variables["planets_sequence"][:-1]):
            if (
                planet,
                self.variables["planets_sequence"][planet_id + 1],
            ) in TOF_BOUNDS.keys():
                key = (planet, self.variables["planets_sequence"][planet_id + 1])
            elif (
                self.variables["planets_sequence"][planet_id + 1],
                planet,
            ) in TOF_BOUNDS.keys():
                key = (self.variables["planets_sequence"][planet_id + 1], planet)
            else:
                raise KeyError(
                    "Planet couple "
                    + str((self.variables["planets_sequence"][planet_id + 1], planet))
                    + " not found in bounds list"
                )
            vector_bounds.append(TOF_BOUNDS[key])
            if planet_id < len(self.variables["planets_sequence"]) - 2:
                sequence_element = [TOF_BOUNDS[key], VARIABLES_BOUNDS["planet_arrival_radius"], VARIABLES_BOUNDS["planet_arrival_angle"]]
            else:
                sequence_element = TOF_BOUNDS[key]
            sequence_bounds.append(sequence_element)

        for _ in range(len(self.variables["planets_sequence"]) - 2):
            vector_bounds.append(VARIABLES_BOUNDS["planet_arrival_radius"])
            vector_bounds.append(VARIABLES_BOUNDS["planet_arrival_angle"])

        self.vector_bounds = vector_bounds
        self.sequence_bounds = sequence_bounds
