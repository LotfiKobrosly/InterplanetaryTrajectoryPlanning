"""
Implements a Nested Monte Carlo Search
"""

from copy import deepcopy
import random
import numpy as np
import pykep as pk
from utils.constants import (
    PLANETS,
    VARIABLES_BOUNDS,
    SAMPLING_FUNCTIONS,
    UNFEASIBILITY_VALUE,
)
from classes.trajectory import Trajectory
from solvers.continuous_variables_choice import (
    get_variables_values,
)


def nmcts(
    trajectory: Trajectory, level: int = 0, continuous_variables_parameters: dict = None
):
    assert (
        continuous_variables_parameters["sampling_function"] in SAMPLING_FUNCTIONS
    ), "Specified sampling_function is unrecognizable"
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
            return trajectory, UNFEASIBILITY_VALUE  # Invalid sequence
        else:
            # print(trajectory.variables["planets_sequence"])
            trajectory.set_variables_bounds()
            variables = trajectory.variables
            # print("Planets: ", variables["planets_sequence"])
            # print("Bounds: ", trajectory.bounds)
            input_values = deepcopy(continuous_variables_parameters)
            input_values["planets_sequence"] = variables["planets_sequence"]
            input_values["bounds"] = trajectory.bounds
            input_values["values_sequence"] = list()
            (
                variables["departure_epoch"],
                variables["time_of_flights_list"],
                delta_v,
            ) = get_variables_values(input_values)
            return trajectory, delta_v

    else:
        best_planet_sequence = None
        best_value = np.inf
        best_trajectory = None
        while not trajectory.planets_sequence_is_set() and len(
            trajectory.variables["planets_sequence"]
        ) < len(PLANETS):
            for planet in trajectory.planets_pool:
                if planet not in trajectory.variables["planets_sequence"]:
                    temporary_trajectory = deepcopy(trajectory)
                    temporary_trajectory.variables = deepcopy(trajectory.variables)
                    temporary_trajectory.add_planet(planet)
                    temporary_trajectory, value = nmcts(
                        temporary_trajectory, level - 1, continuous_variables_parameters
                    )
                    if value < best_value:
                        best_value = value
                        best_trajectory = deepcopy(temporary_trajectory)
            best_next_planet = best_trajectory.variables["planets_sequence"][
                len(trajectory.variables["planets_sequence"])
            ]
            trajectory.add_planet(best_next_planet)
            # print("Current length:", len(trajectory.variables["planets_sequence"]))
        # best_trajectory.evaluate_mga()
        return best_trajectory, best_value


if __name__ == "__main__":
    trajectory = Trajectory()
    trajectory.instantiate("Earth", "Jupiter")
    # for sampling_function in ["cnrpa"]: #["uniform", "gaussian_cma_es", "cnmcts"]:
    continuous_variables_parameters = {
        "sampling_function": "crbnmcts",
        "n_iterations": 500,  # for uniform and gaussian sampling
        "level": 2,  # for cNMCTS, cNRPA and derivatives
        "bandwidth": 100,  # for cNMCTS
        "n_policies": 4000,  # for cNRPA and derivatives
        "multiple_values_policy": True,  # for cNRPA and derivatives
        "learning_rate": 0.01,  # for cNRPA and derivatives
        "tau": 20,
    }
    print("\nSampling:", continuous_variables_parameters["sampling_function"])
    trajectory.reinitialize()
    trajectory, value = nmcts(trajectory, 2, continuous_variables_parameters)
    print(f"Value inside class instance: {trajectory.evaluate_mga() / 1000:.3f} km/s")
    print(f"Best delta V: {value / 1000:.3f} km/s")
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
