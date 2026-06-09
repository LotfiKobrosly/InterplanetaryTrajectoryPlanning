"""
Implements a one-shot choice of the control variables in a vector.

For the gaussian-based implementation, we update the means and covariance with the CMA-ES (Covariance Matrix Adaptation Evolution Strategy)
"""

import numpy as np
import cma


def separate_values(input_vector: np.ndarray, planets_sequence_length: int):
    departure_epoch = input_vector[0]
    times_of_flight = input_vector[1:planets_sequence_length]
    flyby_parameters = input_vector[planets_sequence_length:]
    flyby_parameters = [
        (flyby_parameters[i], flyby_parameters[i + 1])
        for i in range(0, len(flyby_parameters), 2)
    ]
    return departure_epoch, times_of_flight, flyby_parameters


def uniform_variables_values_vector(bounds: list):
    return np.random.uniform(
        np.array([bound[0] for bound in bounds]),
        np.array([bound[1] for bound in bounds]),
    )


def gaussian_variables_values_vector(means: list, variances: list, bounds: list):
    pass


def adapt_gaussian_estimator():
    pass


if __name__ == "__main__":
    # --- Define your search vector ---
    # [t0, tof1, tof2, rp1_vrad, beta1]  for pure MGA Earth->Venus->Mars

    # Bounds per variable
    bounds_lo = [1000, 80, 100, 1.05, 0.0]
    bounds_hi = [2000, 300, 400, 10.0, 2 * np.pi]

    def normalize(x, lo, hi):
        return [(xi - l) / (h - l) for xi, l, h in zip(x, lo, hi)]

    def denormalize(x, lo, hi):
        return [xi * (h - l) + l for xi, l, h in zip(x, lo, hi)]

    def objective(x_norm):
        """
        CMA-ES minimizes this. Work in normalized [0,1] space.
        Returns total Δv, or a large penalty if invalid.
        """
        x = denormalize(x_norm, bounds_lo, bounds_hi)
        t0_mjd, tof1, tof2, rp1, beta1 = x

        result = evaluate_trajectory(
            t0_mjd=t0_mjd, tof1_days=tof1, tof2_days=tof2, rp1_vrad=rp1, beta1=beta1
        )

        if result is None:
            return 1e10  # penalty for invalid trajectory

        return result["total_dv_ms"]

    # --- Initial mean: center of search space (normalized) ---
    x0 = [0.5] * 5  # start at center
    sigma0 = 0.3  # initial step size (in normalized space)

    # --- CMA-ES options ---
    opts = cma.CMAOptions()
    opts["bounds"] = [[0.0] * 5, [1.0] * 5]  # normalized bounds
    opts["maxiter"] = 500
    opts["popsize"] = 20  # λ: samples per iteration
    opts["tolx"] = 1e-6  # convergence on x
    opts["tolfun"] = 1e-6  # convergence on f
    opts["verbose"] = 1

    # --- Run ---
    es = cma.CMAEvolutionStrategy(x0, sigma0, opts)

    while not es.stop():
        candidates = es.ask()  # sample λ candidates
        fitnesses = [objective(x) for x in candidates]
        es.tell(candidates, fitnesses)  # update μ, C, σ
        es.logger.add()
        es.disp()

    result_norm = es.result.xbest
    result_x = denormalize(result_norm, bounds_lo, bounds_hi)

    print("\nBest solution:")
    print(f"  t0       : {result_x[0]:.1f} MJD2000")
    print(f"  tof1     : {result_x[1]:.1f} days")
    print(f"  tof2     : {result_x[2]:.1f} days")
    print(f"  rp1      : {result_x[3]:.2f} Venus radii")
    print(f"  beta1    : {result_x[4]:.3f} rad")
    print(f"  Best Δv  : {es.result.fbest/1000:.3f} km/s")
