"""
Implements a Nested Monte Carlo Search
"""

from copy import deepcopy
import random
import numpy as np
import pykep as pk
from utils.constants import PLANETS, VARIABLES_BOUNDS
from utils.trajectory_evaluation import evaluate_mga_trajectory
from classes.trajectory import Trajectory
from solvers.continuous_variables_choice.one_shot_vector_choice import (
    uniform_variables_values_vector,
    separate_values,
)


def nmcts(trajectory: Trajectory, level: int = 0):
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
            departure_velocity = 1e12
            n_iterations = 0

            while (departure_velocity > VARIABLES_BOUNDS["departure_velocity"][1]) and (n_iterations < 1000):
                n_iterations += 1
                (
                    variables["departure_epoch"],
                    variables["time_of_flights_list"],
                    variables["planets_flyby_parameters"],
                ) = separate_values(
                    uniform_variables_values_vector(trajectory.bounds),
                    len(variables["planets_sequence"]),
                )
                value = trajectory.evaluate()
                departure_velocity = np.linalg.norm(trajectory.mga_results[1][0])
                # if departure_velocity < VARIABLES_BOUNDS["departure_velocity"][1]:
                #     print(f"Valid departure velocity found at: {departure_velocity / 1000:.3f} km/s")
                #     print("For given sequence:", trajectory.variables["planets_sequence"])
                #     print(f"Giving delta V of {value:.3f} km/s")
            if (n_iterations >= 1000) and (departure_velocity > VARIABLES_BOUNDS["departure_velocity"][1]):
                value = 1e12
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
                    temporary_trajectory, value = nmcts(temporary_trajectory, level - 1)
                    if value < best_value:
                        best_value = value
                        best_next_planet = planet
                        best_trajectory = temporary_trajectory
            if best_next_planet is None:
                raise ValueError("No good trajectory found")
            trajectory.add_planet(best_next_planet)
            # print("Current length:", len(trajectory.variables["planets_sequence"]))

        return best_trajectory, best_value


if __name__ == "__main__":
    trajectory = Trajectory()
    trajectory.instantiate("Earth", "Saturn")
    result, value = nmcts(trajectory, 1)
    print(f"Best delta V: {value:.3f} km/s")
    print("Variables values:")
    print(result.variables)
    print(
        f"Departure velocity: {np.linalg.norm(result.mga_results[1][0]) / 1000:.3f} km/s"
    )
