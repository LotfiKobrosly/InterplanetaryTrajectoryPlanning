"""
Implements a trajectory validation function to verify constraints are not violated
"""

import numpy as np
import pykep as pk


def trajectory_is_valid(
    planets_sequence: list = None,
    values_vector: list = None,
) -> bool:
    # Verify that at least one of the inputs was provided correctly
    assert not (
        (planets_sequence is None) or (values_vector is None)
    ), "Invalid arguments: either planets_sequence is None or values_vector is None"

    # Planets
    planets_sequence = [
        pk.planet(pk.udpla.jpl_lp(planet)) for planet in planets_sequence
    ]

    # Transfer legs computation
    lambert_legs_list = list()
    epoch = values_vector[0]
    last_planet = planets_sequence[0]
    planet_radius, planet_velocity = planets_sequence[0].eph(epoch)
    radii_list = [planet_radius]
    velocities_list = [planet_velocity]
    for planet_index, planet in enumerate(planets_sequence[1:]):

        epoch += values_vector[planet_index + 1]
        planet_radius, planet_velocity = planet.eph(epoch)

        lambert_legs_list.append(
            pk.lambert_problem(
                tof=values_vector[planet_index + 1] * pk.DAY2SEC,
                r0=radii_list[-1],
                r1=planet_radius,
                mu=pk.MU_SUN,
                cw=False,
                multi_revs=0,
            )
        )
        radii_list.append(planet_radius)
        velocities_list.append(planet_velocity)

    # Flyby constraint verification
    for planet_index, planet in enumerate(planets_sequence[1:-1]):
        dv, violation = pk.fb_con(
            v_rel_in=(
                np.array(lambert_legs_list[planet_index].v1[0])
                - np.array(velocities_list[planet_index + 1])
            ).tolist(),
            v_rel_out=(
                np.array(lambert_legs_list[planet_index + 1].v0[0])
                - np.array(velocities_list[planet_index + 1])
            ).tolist(),
            mu=planet.get_mu_self(),
            safe_radius=planet.get_safe_radius(),
        )
        if violation >= 0:
            return False

    return True
