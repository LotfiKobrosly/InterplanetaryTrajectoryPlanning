"""
Trajectory evaluation
"""

import pykep as pk
import numpy as np
from utils.constants import VARIABLES_BOUNDS


def evaluate_mga_trajectory(
    planets_sequence: list = None,
    departure_epoch: float = None,
    time_of_flights_list: list = None,
    planets_flyby_parameters: list = None,
    *args,
    **kwargs,
):

    # Initializing varibales to store results
    total_delta_V = 0
    departure_velocities = list()
    arrival_velocities = list()
    last_departure_velocity = None
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

    # Lambert legs and Gravity assists
    for planet_index, planet in enumerate(planets_sequence[1:]):
        # Lambert leg between planet and planet_sequence[planet_index + 1]
        lambert_leg = pk.lambert_problem(
            r0=planet_radii_list[planet_index],
            r1=planet_radii_list[planet_index + 1],
            tof=time_of_flights_list[planet_index] * pk.DAY2SEC,
            mu=pk.MU_SUN,
            cw=False,
            multi_revs=0,
        )

        ## Departure and arrival velocities

        departure_velocity = np.array(lambert_leg.v0[0])
        arrival_velocity = np.array(lambert_leg.v1[0])
        arrival_velocities.append(arrival_velocity)

        if planet_index == 0:
            if (
                np.linalg.norm(departure_velocity)
                > VARIABLES_BOUNDS["departure_velocity"][1]
            ):
                total_delta_V += 1e10 * (
                    np.linalg.norm(departure_velocity)
                    - VARIABLES_BOUNDS["departure_velocity"][1]
                )
            total_delta_V += np.linalg.norm(
                departure_velocity - planet_velocities_list[0]
            )
            departure_velocities.append(departure_velocity)
        if not (last_departure_velocity is None):
            total_delta_V += np.linalg.norm(
                departure_velocity - last_departure_velocity
            )  # Delta Patch velocity

        # Gravity assist at planet_index + 1 if not final
        if planet_index + 1 < len(planets_sequence) - 1:
            fly_by_radius, fly_by_angle = planets_flyby_parameters[planet_index]

            if fly_by_radius < 1:
                return None  # Spacecraft would crash into the planet

            # fb_vout takes absolute incoming velocity (not relative), planet velocity,
            # periapsis radius, beta angle, and planet mu
            last_departure_velocity = np.array(
                pk.fb_vout(
                    v_in=arrival_velocity.tolist(),
                    v_pla=planet_velocities_list[planet_index + 1].tolist(),
                    rp=fly_by_radius * planet.radius,
                    beta=fly_by_angle,
                    mu=planet.mu_self,
                )
            )
            departure_velocities.append(last_departure_velocity)
            # Delta-v at Venus: difference between Lambert arrival and flyby departure
            # (ideally 0 for a pure gravity assist)
            delta_velocity_planet = pk.fb_dv(
                v_rel_in=(
                    arrival_velocity - planet_velocities_list[planet_index + 1]
                ).tolist(),
                v_rel_out=(
                    last_departure_velocity - planet_velocities_list[planet_index + 1]
                ).tolist(),
                mu=planet.mu_self,
                safe_radius=planet.radius * 1.05,  # 5% margin above surface
            )
            total_delta_V += delta_velocity_planet

        else:
            total_delta_V += np.linalg.norm(
                arrival_velocity - planet_velocities_list[-1]
            )

    return (
        total_delta_V,
        departure_velocities,
        arrival_velocities,
        planet_radii_list,
        planet_velocities_list,
    )
