"""
Implements a function to get range of orbital speeds of planets
"""

import pykep as pk

ASTRONOMICAL_UNIT = 149597870700  # Wikipedia


def orbital_speeds(a_au, e):
    a = a_au * ASTRONOMICAL_UNIT
    q = a * (1 - e)  # perihelion (m)
    Q = a * (1 + e)  # aphelion (m)
    # vis-viva: v = sqrt(mu * (2/r - 1/a))
    v_peri = (pk.MU_SUN * (2 / q - 1 / a)) ** 0.5  # max speed (at perihelion)
    v_aph = (pk.MU_SUN * (2 / Q - 1 / a)) ** 0.5  # min speed (at aphelion)
    return q, Q, v_peri, v_aph
