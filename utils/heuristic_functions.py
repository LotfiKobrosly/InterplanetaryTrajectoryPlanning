"""
Provided by Claude LLM
"""

import numpy as np
import pykep as pk
from utils.constants import SAFE_RADIUS_FACTOR, UNFEASIBILITY_VALUE


def synodic_period(T1: float, T2: float):
    """
    T1, T2: orbital periods of the two planets (days)
    """
    return abs(1 / (1 / T1 - 1 / T2))


def porkchop_scan(
    planet1: pk.udpla,
    planet2: pk.udpla,
    t0_range: list,
    tof_range: list,
    last_arrival_velocity: np.ndarray,
    n_t0: int = 50,
    n_tof: int = 50,
):
    """
    Returns a 2D grid of delta_v over (t0, tof).
    Use this to find good starting regions, then refine locally.
    """
    t0_vals = np.linspace(*t0_range, n_t0)
    tof_vals = np.linspace(*tof_range, n_tof)

    dv_grid = np.full((n_t0, n_tof), UNFEASIBILITY_VALUE)
    last_arrival_velocity = np.array(last_arrival_velocity)

    for i, t0 in enumerate(t0_vals):
        r1, v1 = planet1.eph(pk.epoch(t0))
        v1 = np.array(v1)
        for j, tof in enumerate(tof_vals):
            t1 = t0 + tof
            r2, v2 = planet2.eph(pk.epoch(t1))
            try:
                lambert_leg = pk.lambert_problem(
                    tof=tof * pk.DAY2SEC, r0=r1, r1=r2, mu=pk.MU_SUN
                )
                departure_velocity = np.array(lambert_leg.v0[0])
                _, violation = pk.fb_con(
                    v_rel_in=(last_arrival_velocity - v1).tolist(),
                    v_rel_out=(departure_velocity - v1).tolist(),
                    mu=planet1.mu_self,
                    safe_radius=planet1.radius * SAFE_RADIUS_FACTOR,
                )
                if violation > 0:
                    continue
                else:
                    dv_grid[i, j] = pk.fb_dv(
                        v_rel_in=(last_arrival_velocity - v1).tolist(),
                        v_rel_out=(departure_velocity - v1).tolist(),
                        mu=planet1.mu_self,
                        safe_radius=planet1.radius * SAFE_RADIUS_FACTOR,
                    )
            except Exception:
                continue

    return t0_vals, tof_vals, dv_grid
