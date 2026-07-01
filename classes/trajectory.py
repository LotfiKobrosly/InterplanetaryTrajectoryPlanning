"""
Implements a class for trajectory building
"""

from copy import deepcopy
import numpy as np
import pykep as pk
import matplotlib.pyplot as plt
from utils.constants import (
    PLANETS,
    POLICY_ALGORITHMS,
    VARIABLES_BOUNDS,
    TOF_BOUNDS,
    UNFEASIBILITY_VALUE,
    DV_LAUNCHER,
)


class Trajectory:

    def __init__(self):
        self.variables = {
            "planets_sequence": list(),  # Length msut be >= 2. If 2, a direct flight to goal is considered
            "values_sequence": list(),  # Length must be len(plent_sequence)
        }
        self.start = None
        self.goal = None
        self.evaluator = None

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
            "planets_sequence": [self.start],
            "values_sequence": list(),
        }
        self.define_planets_pool()
        self.evaluator = None

    def planets_sequence_is_set(self):
        return self.variables["planets_sequence"] and (
            self.variables["planets_sequence"][-1] == self.goal
        )

    def is_terminal(self):
        return self.planets_sequence_is_set() and len(
            self.variables["values_sequence"]
        ) == len(self.variables["planets_sequence"] - 1)

    def evaluate_mga(self):
        if self.evaluator is None:
            planets_sequence = [
                pk.planet(pk.udpla.jpl_lp(planet))
                for planet in self.variables["planets_sequence"]
            ]
            self.evaluator = pk.trajopt.mga(
                planets_sequence,
                list(self.bounds[0]),
                [list(element) for element in self.bounds[1:]],
                vinf=DV_LAUNCHER,
            )
        self.mga_results = self.evaluator.fitness(self.variables["values_sequence"])
        if self.mga_results is None:
            return UNFEASIBILITY_VALUE  # Invalid sequence
        return self.mga_results[0]

    def add_planet(self, planet: str):
        self.variables["planets_sequence"].append(planet)
        self.planets_pool.remove(planet)

    def set_variables_bounds(self):
        assert self.planets_sequence_is_set(), "Planets' sequence not set yet"
        bounds = [VARIABLES_BOUNDS["departure_epoch"]]
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
            bounds.append(TOF_BOUNDS[key])

        self.bounds = bounds

    def plot_trajectory(self):
        if self.evaluator is None:
            planets_sequence = [
                pk.planet(pk.udpla.jpl_lp(planet))
                for planet in self.variables["planets_sequence"]
            ]
            self.evaluator = pk.trajopt.mga(
                planets_sequence,
                list(self.bounds[0]),
                [list(element) for element in self.bounds[1:]],
                vinf=DV_LAUNCHER,
            )
        axe = self.evaluator.plot(self.variables["values_sequence"])
        axe.view_init(90, 0)
        axe.axis("off")
        axe.set_title(
            r"$\Delta$V = "
            + f"{self.evaluator.fitness(self.variables["values_sequence"])[0] / 1000:.3f} km/s"
        )
        plt.show()
