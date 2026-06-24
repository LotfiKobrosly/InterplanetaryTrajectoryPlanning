"""
Trajectory evaluation
"""

import pykep as pk
import numpy as np
from utils.constants import VARIABLES_BOUNDS, SAFE_RADIUS_FACTOR, UNFEASIBILITY_VALUE


def evaluate_mga_trajectory(
    planets_sequence: list = None,
    departure_epoch: float = None,
    time_of_flights_list: list = None,
    *args,
    **kwargs,
):  # with the help of Claude LLM from Anthropic
    assert len(planets_sequence) == len(time_of_flights_list) + 1, (
        "Mismatch between planets_sequence and time_of_flights_list: "
        + str(len(planets_sequence))
        + " planets VS "
        + str(len(time_of_flights_list))
        + " times of flight"
    )

    # Initializing varibales to store results
    total_delta_V = 0
    departure_velocities_list = list()
    arrival_velocities_list = list()
    planets_sequence = [
        pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence
    ]

    # Epochs
    epochs_list = [departure_epoch]
    for flight_time in time_of_flights_list:
        epochs_list.append(epochs_list[-1] + flight_time)
    epochs_list = [pk.epoch(element) for element in epochs_list]

    # Ephemeris
    planet_radii_list, planet_velocities_list = list(), list()
    for index, epoch in enumerate(epochs_list):
        radius, velocity = planets_sequence[index].eph(epoch)
        planet_radii_list.append(np.array(radius))
        planet_velocities_list.append(np.array(velocity))

    # Initializing departure planet velocity
    arrival_velocities_list.append(planet_velocities_list[0])

    # Lambert legs
    for planet_index in range(len(planets_sequence[:-1])):
        # Lambert leg between planet and planet_sequence[planet_index + 1]
        lambert_leg = pk.lambert_problem(
            r0=planet_radii_list[planet_index],
            r1=planet_radii_list[planet_index + 1],
            tof=time_of_flights_list[planet_index] * pk.DAY2SEC,
            mu=pk.MU_SUN,
            cw=False,
            multi_revs=0,
        )

        ## Storing departure and arrival velocities

        departure_velocities_list.append(np.array(lambert_leg.v0[0]))
        arrival_velocities_list.append(np.array(lambert_leg.v1[0]))

    # Verifying validity of departure velocity from departure planet
    if (
        np.linalg.norm(departure_velocities_list[0] - arrival_velocities_list[0])
        > VARIABLES_BOUNDS["departure_velocity"][1]
    ):
        total_delta_V += (
            UNFEASIBILITY_VALUE  # unfeasible departure: launch velocity too high
        )

    # Adding target planet velocity
    departure_velocities_list.append(planet_velocities_list[-1])

    # Computing velocity and total_delta_v at each gravity assist
    for planet_index, planet in enumerate(planets_sequence):
        # Computing velocity difference at launch or at arrival at target planet
        if (planet_index == 0) or (planet_index == len(planets_sequence) - 1):
            total_delta_V += np.linalg.norm(
                departure_velocities_list[planet_index]
                - arrival_velocities_list[planet_index]
            )

        else:
            # Gravity assist verification
            _, violation = pk.fb_con(
                v_rel_in=(
                    arrival_velocities_list[planet_index]
                    - planet_velocities_list[planet_index]
                ).tolist(),
                v_rel_out=(
                    departure_velocities_list[planet_index]
                    - planet_velocities_list[planet_index]
                ).tolist(),
                mu=planet.mu_self,
                safe_radius=planet.radius * SAFE_RADIUS_FACTOR,
                # pl        = planet  # pykep planet object knows its own mu, radius
            )
            if (
                violation > 0
            ):  # value needs to be non-positive to be geometrically feasible
                total_delta_V += UNFEASIBILITY_VALUE

            else:

                # Velocity difference between arrival and departure, what we want ot minimize
                delta_velocity_planet = pk.fb_dv(
                    v_rel_in=(
                        arrival_velocities_list[planet_index]
                        - planet_velocities_list[planet_index]
                    ).tolist(),
                    v_rel_out=(
                        departure_velocities_list[planet_index]
                        - planet_velocities_list[planet_index]
                    ).tolist(),
                    mu=planet.mu_self,
                    safe_radius=planet.radius
                    * SAFE_RADIUS_FACTOR,  # 5% margin above surface
                )
                total_delta_V += np.linalg.norm(delta_velocity_planet)

    return (
        total_delta_V,
        departure_velocities_list,
        arrival_velocities_list,
        planet_radii_list,
        planet_velocities_list,
    )
