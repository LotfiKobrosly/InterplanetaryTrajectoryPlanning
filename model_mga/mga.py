# Author: Lotfi Kobrosly (2025)

# This file contains the equations of the multiple gravity assist model

from collections.abc import Callable
import numpy as np
from jplephem.spk import SPK
from utils.constants import *

# ================================================================
# Utility Functions
# ================================================================


def norm(v):
    return np.linalg.norm(v)


def unit(v):
    return v / norm(v)


# ================================================================
# Two-body Propagation (universal variable formulation)
# ================================================================


def kepler_propagate(r0: np.ndarray, v0: np.ndarray, mu: float, dt: float):
    """
    Propagate a position/velocity vector under 2-body motion using
    universal variable formulation (works for elliptic & hyperbolic).
    """
    r0n = norm(r0)
    v0n = norm(v0)
    alpha = 2 / r0n - v0n**2 / mu

    # Initial guess for universal variable
    if alpha > 0:
        chi = np.sqrt(mu) * alpha * dt
    else:
        chi = (
            np.sign(dt)
            * np.sqrt(-1 / alpha)
            * np.log(
                (-2 * mu * alpha * dt)
                / (dot(r0, v0) + np.sign(dt) * np.sqrt(-mu / alpha) * (1 - r0n * alpha))
            )
        )

    # Iterative solve
    for _ in range(50):
        z = alpha * chi**2
        C = stumpff_C(z)
        S = stumpff_S(z)

        r = (
            chi**2 * C
            + dot(r0, v0) / np.sqrt(mu) * chi * (1 - z * S)
            + r0n * (1 - z * C)
        )
        f = r - mu * dt
        if abs(f) < 1e-8:
            break
        df = chi * (1 - z * S)
        chi = chi - f / df

    # f and g functions
    z = alpha * chi**2
    C = stumpff_C(z)
    S = stumpff_S(z)

    f = 1 - chi**2 / r0n * C
    g = dt - chi**3 / np.sqrt(mu) * S
    r = f * r0 + g * v0

    rnorm = norm(r)
    fdot = np.sqrt(mu) / (r0n * rnorm) * (z * S - 1) * chi
    gdot = 1 - chi**2 / rnorm * C

    v = fdot * r0 + gdot * v0

    return r, v


def stumpff_C(z: float):
    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        return 1 / 2


def stumpff_S(z: float):
    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (z**1.5)
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / ((-z) ** 1.5)
    else:
        return 1 / 6


# ================================================================
# Lambert Solver (universal-variable form)
# ================================================================


def lambert_universal(r1: np.ndarray, r2: np.ndarray, dt: float, mu: float):
    """
    Robust Lambert solver using universal variables.
    Returns v1, v2.
    """

    R1 = norm(r1)
    R2 = norm(r2)

    cos_dtheta = np.dot(r1, r2) / (R1 * R2)
    cos_dtheta = np.clip(cos_dtheta, -1.0, 1.0)
    dtheta = np.arccos(cos_dtheta)

    # Short-way transfer
    A = np.sin(dtheta) * np.sqrt(R1 * R2 / (1 - np.cos(dtheta)))

    if A == 0:
        raise ValueError("Lambert geometry degenerates (A=0).")

    # Start with z = 0 — THIS IS VALID
    z = 0.0

    def y(z):
        C = stumpff_C(z)
        S = stumpff_S(z)
        return R1 + R2 + A * (z * S - 1) / np.sqrt(C)

    def F(z):
        C = stumpff_C(z)
        S = stumpff_S(z)
        Y = y(z)
        return (Y / C) ** 1.5 * S + A * np.sqrt(Y) - np.sqrt(mu) * dt

    # Newton iteration on z
    for _ in range(100):
        C = stumpff_C(z)
        S = stumpff_S(z)
        Y = y(z)

        # Derivative dF/dz — written in numerically safe form
        if abs(z) < 1e-8:
            # L'Hôpital limit for z → 0
            dF = np.sqrt(2) * A / 40 * (Y**1.5)  # safe approximation
        else:
            dY = (A / 2) * (S / z + (C - 3 * S) / (2 * C))
            dF = (1.5 * np.sqrt(Y) / C - Y**1.5 * S / (2 * C**2)) * dY + A * dY / (
                2 * np.sqrt(Y)
            )

        z_next = z - F(z) / dF
        if abs(z_next - z) < 1e-10:
            z = z_next
            break
        z = z_next

    # Final values
    C = stumpff_C(z)
    S = stumpff_S(z)
    Y = y(z)

    f = 1 - Y / R1
    g = A * np.sqrt(Y / mu)
    gdot = 1 - Y / R2

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g

    return v1, v2


# ================================================================
# Gravity-Assist (Patched Conics)
# ================================================================


def gravity_assist(
    v_inf_in: np.ndarray, planet_velocity: np.ndarray, mu_planet: float, rp: float
):
    """
    Rotate v_infinity vector around planet using hyperbolic turn angle.
    v_inf_in: incoming heliocentric velocity minus planet velocity
    rp: periapsis radius of flyby
    """

    v_inf = norm(v_inf_in)
    e = 1 + rp * v_inf**2 / mu_planet  # hyperbolic eccentricity
    delta = 2 * np.arcsin(1 / e)  # turn angle

    # Choose arbitrary rotation axis orthogonal to v_inf_in
    axis = unit(np.cross(v_inf_in, np.array([0, 0, 1])))
    if norm(axis) < 1e-6:
        axis = unit(np.cross(v_inf_in, np.array([0, 1, 0])))

    # Rodrigues rotation
    v_inf_out = rotate(v_inf_in, axis, delta)

    # Convert back to heliocentric frame
    return v_inf_out + planet_velocity


def rotate(v: np.ndarray, axis: np.ndarray, angle: float):
    axis = unit(axis)
    return (
        v * np.cos(angle)
        + np.cross(axis, v) * np.sin(angle)
        + axis * np.dot(axis, v) * (1 - np.cos(angle))
    )


# ================================================================
# Multi-Gravity Assist Mission
# ================================================================


def multi_gravity_assist(
    kernel: SPK,
    sequence: list,
    times: list,
    periapses: list,
    ephemeris: Callable,
    masses: dict,
):
    """
    kernel   : the SPK kernel reading the ephemerids database
    sequence : list of planet names, e.g. ["Earth","Venus","Earth","Jupiter"]
    times    : list of encounter epochs
    periapses: list of periapsis radii for each flyby
    ephemeris: function body, t → (r, v)
    masses   : dict of planetary masses
    """
    trajectory = []

    # Initial state: spacecraft starts at first planet
    p0 = sequence[0]
    r0, v0 = ephemeris(kernel, p0, times[0])
    v_sc = v0  # assume spacecraft initially co-moving

    for i in range(1, len(sequence)):
        p_in = sequence[i - 1]
        p_out = sequence[i]

        t1 = times[i - 1]
        t2 = times[i]

        r1, v1_planet = ephemeris(kernel, p_in, t1)
        r2, v2_planet = ephemeris(kernel, p_out, t2)

        mu_sun = GRAVITY_CONSTANT * masses["Sun"]

        # Transfer arc
        v_depart, v_arrive = lambert_universal(r1, r2, t2 - t1, mu_sun)
        v_sc = v_depart  # heliocentric velocity leaving flyby planet

        # At next planet, compute v_inf
        v_inf_in = v_arrive - v2_planet

        # Gravity assist
        mu_planet = GRAVITY_CONSTANT * masses[p_out]
        rp = periapses[i]

        v_sc = gravity_assist(v_inf_in, v2_planet, mu_planet, rp)

        trajectory.append({"planet": p_out, "t": t2, "r": r2, "v": v_sc})

    return trajectory
