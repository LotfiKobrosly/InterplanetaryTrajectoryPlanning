"""
Implements a class for trajectory building
"""

import numpy as np
import pykep as pk
from utils.constants import PLANETS
from utils.trajectory_evaluation import evaluate_mga_trajectory

class Trajectory():

    def __init__(self):
        self.variables = {
            "planets_sequence": list(), # Length msut be >= 2. If 2, a direct flight to goal is considered
            "departure_epoch": None,
            "time_of_flights_list": list(), # Length must be len(plent_sequence) - 1
            "planets_flyby_parameters": list(), # Length must be len(plent_sequence) - 2
        }
        self.start = None
        self.goal = None

    def reinitialize(self):
        self.variables = {
            "planets_sequence": list(), # Length msut be >= 2. If 2, a direct flight to goal is considered
            "departure_epoch": None,
            "time_of_flights_list": list(), # Length must be len(plent_sequence) - 1
            "planets_flyby_parameters": list(), # Length must be len(plent_sequence) - 2
        }

    def set_start_planet(self, planet="Earth"):
        assert planet.upper() in PLANETS, "Given planet not in the provided list: " + str(PLANETS)
        self.start = pk.planet(pk.udpla.jpl_lp(planet))
        self.variables["planets_sequence"] = [self.start]

    def set_goal_planet(self, planet="Mars"):
        assert planet.upper() in PLANETS, "Given planet not in the provided list: " + str(PLANETS)
        self.goal = pk.planet(pk.udpla.jpl_lp(planet))

    def planets_sequence_is_set(self):
        return self.variables["planets_sequence"] and (self.variables["planets_sequence"][-1] == self.goal)

    def is_terminal(self):
        return (
            self.planets_sequence_is_set()
            and not (self.variables["departure_epoch"] is None)
            and len(self.variables["time_of_flights_list"]) == len(self.variables["planets_sequence"] - 1)
            and len(self.variables["planets_flyby_parameters"]) == len(self.variables["planets_sequence"] - 2)
        )

    def choose_solver(self, solver:str):
        self.solver = solver

    def evaluate(self):
        results = evaluate_mga_trajectory(**self.variables)
        return results[0] / 1000 # Convert from m/s to km/s

    def step(self, policy: dict=None):
        if not self.planets_sequence_is_set():
            pass
        else:
            pass
