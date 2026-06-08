import pykep as pk
import numpy as np
import scipy as sp


def build_sf_leg(
    departure_planet_parameters,
    arrival_planet_parameters,
    time_of_flight,
    optimization_method: str="COBYLA", # Possible values that use the constraints argument: SLSQP, COBYQA, COBYLA, trust-constr
    mass=1000.0,
    final_mass=900.0,
    max_thrust=0.1,
    isp=3000.0,
    n_seg=20,
):
    """
    Build a Sims-Flanagan low-thrust leg
    using states from evaluate_trajectory().
    """
    departure_radius, departure_velocity = departure_planet_parameters
    arrival_radius, arrival_velocity = arrival_planet_parameters    

    # pykep 3.x: pk.leg.sims_flanagan(rvs, ms, throttles, rvf, mf, tof,
    #                                  max_thrust, veff, mu, cut)
    throttles = np.zeros(n_seg * 3)
    sf = pk.leg.sims_flanagan(
        rvs=list(departure_planet_parameters),  # [position, velocity] at start
        ms=mass,
        throttles=throttles,
        rvf=list(arrival_planet_parameters),  # [position, velocity] at end
        mf=final_mass,
        tof=time_of_flight * pk.DAY2SEC,
        max_thrust=max_thrust,
        veff=isp * pk.G0,  # effective exhaust velocity (m/s)
        mu=pk.MU_SUN,
        cut=0.5,  # midpoint of the leg
    )
    def compute_mismatch(throttles):
        sf.throttles = throttles[:]
        return np.linalg.norm(sf.compute_mismatch_constraints())

    first_mismatch = sf.compute_mismatch_constraints()
    # print(f"First mismatch: {np.linalg.norm(first_mismatch):.0f}")
    # print(type(first_mismatch))
    constraints = [
        {
            "type": "ineq",
            "fun": lambda x: 1 - np.linalg.norm(x[i:i+3])
        }
        for i in range(0, n_seg * 3, 3)
    ]

    # Option 1
    minimum_value = sp.optimize.minimize(
        compute_mismatch,
        throttles,
        method=optimization_method,
        bounds=tuple([(-1.0, 1.0) for _ in range(n_seg * 3)]),
        constraints=tuple(constraints),
    )
    sf.throttles = minimum_value.x

    return sf, first_mismatch
