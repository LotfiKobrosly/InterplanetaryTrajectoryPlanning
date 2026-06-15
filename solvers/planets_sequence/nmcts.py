"""
Implements a Nested Monte Carlo Search
"""

from copy import deepcopy
import random
import numpy as np
import pykep as pk
from utils.constants import PLANETS, VARIABLES_BOUNDS, SAMPLING_FUNCTIONS, SEQUENCE_FUNCTIONS, VECTOR_FUNCTIONS
from classes.trajectory import Trajectory
from solvers.continuous_variables_choice import (
    get_variables_values,
)



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
            input_values = {
                "planets_sequence": variables["planets_sequence"],
                "sampling_function": sampling_function,
                "n_iterations": 500,                               # for uniform and gaussian sampling
                "level": 2,                                        # for cNMCTS, cNRPA and derivatives
                "bandwidth": 15,                                   # for cNMCTS
                "values_sequence": list(),                         # for cNMCTS
            }
            if sampling_function in SEQUENCE_FUNCTIONS:
                input_values["bounds"] = trajectory.sequence_bounds
            elif sampling_function in VECTOR_FUNCTIONS:
                input_values["bounds"] = trajectory.vector_bounds
            else:
                raise ValueError("Specified sampling_function is unrecognizable")
            (
                variables["departure_epoch"],
                variables["time_of_flights_list"],
                variables["planets_flyby_parameters"],
                value
            ) = get_variables_values(input_values)
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
                    temporary_trajectory, value = nmcts(
                        temporary_trajectory, level - 1, sampling_function
                    )
                    if value < best_value:
                        best_value = value
                        best_next_planet = planet
                        best_trajectory = temporary_trajectory
            if best_next_planet is None:
                best_next_planet = PLANETS[
                    PLANETS.index(trajectory.variables["planets_sequence"][-1]) + 1
                ]
            trajectory.add_planet(best_next_planet)
            # print("Current length:", len(trajectory.variables["planets_sequence"]))
        best_trajectory.evaluate_mga()
        return best_trajectory, best_value


if __name__ == "__main__":
    trajectory = Trajectory()
    trajectory.instantiate("Earth", "Saturn")
    for sampling_function in ["uniform", "gaussian_cma_es", "cnmcts"]:
        print("\nSampling:", sampling_function)
        trajectory.reinitialize()
        trajectory, value = nmcts(trajectory, 1, sampling_function)
        print(f"Best delta V: {value:.3f} km/s")
        print("Planet sequence:", trajectory.variables["planets_sequence"])
        print("Departures velocities:")
        for velocity in trajectory.mga_results[1]:
            print(f"   {np.linalg.norm(velocity) / 1000:.3f} km/s ")
        print("Arrivals velocities:")
        for velocity in trajectory.mga_results[2]:
            print(f"   {np.linalg.norm(velocity) / 1000:.3f} km/s ")
        print("Planets velocities:")
        for velocity in trajectory.mga_results[-1]:
            print(f"   {np.linalg.norm(velocity) / 1000:.3f} km/s ")
