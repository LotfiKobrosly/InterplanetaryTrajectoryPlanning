"""
Implements a Nested Monte Carlo Search
"""

from copy import deepcopy
import random
import numpy as np
import pykep as pk
from utils.constants import PLANETS, VARIABLES_BOUNDS, SAMPLING_FUNCTIONS
from classes.trajectory import Trajectory
from solvers.continuous_variables_choice.one_shot_vector_choice import get_variables_values


def nmcts(trajectory: Trajectory, level: int = 0, sampling_function: str = "uniform"):
    if level == 0:
        # Sequence setting
        while (not trajectory.planets_sequence_is_set()) and len(
            trajectory.variables["planets_sequence"]
        ) < len(PLANETS):
            chosen_planet = random.choice(trajectory.planets_pool)
            trajectory.add_planet(chosen_planet)

        # Values for continuous variables
        if len(trajectory.variables["planets_sequence"]) == len(PLANETS) and (
            not trajectory.planets_sequence_is_set()
        ):
            return trajectory, 1e30  # Invalid sequence
        else:
            # print(trajectory.variables["planets_sequence"])
            trajectory.set_variables_bounds()
            variables = trajectory.variables
            (
                variables["departure_epoch"],
                variables["time_of_flights_list"],
                variables["planets_flyby_parameters"],
            ) = get_variables_values(trajectory.bounds, variables["planets_sequence"], sampling_function)
            value = trajectory.evaluate_mga()
            return trajectory, value

    else:
        while not trajectory.planets_sequence_is_set() and len(
            trajectory.variables["planets_sequence"]
        ) < len(PLANETS):
            best_value = 1e30
            best_next_planet = None
            best_trajectory = None
            for planet in trajectory.planets_pool:
                if planet not in trajectory.variables["planets_sequence"]:
                    temporary_trajectory = deepcopy(trajectory)
                    temporary_trajectory.add_planet(planet)
                    temporary_trajectory, value = nmcts(temporary_trajectory, level - 1, sampling_function)
                    if value < best_value:
                        best_value = value
                        best_next_planet = planet
                        best_trajectory = temporary_trajectory
            if best_next_planet is None:
                best_next_planet = PLANETS[PLANETS.index(trajectory.variables["planets_sequence"][-1]) + 1]
            trajectory.add_planet(best_next_planet)
            # print("Current length:", len(trajectory.variables["planets_sequence"]))

        return best_trajectory, best_value


if __name__ == "__main__":
    trajectory = Trajectory()
    trajectory.instantiate("Earth", "Saturn")
    result, value = nmcts(trajectory, 1, "gaussian_cma_es")
    print(f"Best delta V: {value:.3f} km/s")
    print("Variables values:")
    print(result.variables)
    print(
        f"Departure velocity: {np.linalg.norm(result.mga_results[1][0]) / 1000:.3f} km/s"
    )
