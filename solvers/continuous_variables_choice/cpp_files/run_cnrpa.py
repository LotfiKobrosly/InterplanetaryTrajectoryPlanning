"""
Example: run the C++ CNRPA (Claude code)
"""

import time
import pykep as pk
import pygmo as pg
import numpy as np
import cnrpa_cpp  # the compiled C++ module
from utils.constants import GAUSSIAN_KERNEL_THRESHOLD, UNFEASIBILITY_VALUE, RANDOM_SEED
from utils.udp_wrapper import CountingEvaluator

# --- Problem -----------------------------------------------------------------
udp = CountingEvaluator(pk.trajopt.gym.cassini1)
# udp = pk.trajopt.gym.cassini1

# IMPORTANT: always get bounds via pg.problem() wrapper,
# which guarantees clean numpy arrays regardless of pykep version
prob = pg.problem(udp)
lb, ub = prob.get_bounds()

# Convert to plain Python floats — pybind11 is strict about types
# [[lo0, hi0], [lo1, hi1], ...]
bounds = [[float(lo), float(hi)] for lo, hi in zip(lb, ub)]

print(f"Problem: {prob.get_name()}")
print(f"Number of variables: {prob.get_nx()}")
# print(f"Bounds: {bounds}")

# --- Run C++ CNRPA -----------------------------------------------------------
t0 = time.time()
result = cnrpa_cpp.cnrpa(
    evaluator=udp,
    level=2,
    n_policies=400,
    bounds=bounds,
    learning_rate=0.1,
    timeout=60.0,
    unfeasibility_value=UNFEASIBILITY_VALUE,
    gaussian_kernel_threshold=GAUSSIAN_KERNEL_THRESHOLD,
    seed=RANDOM_SEED,
)
elapsed = time.time() - t0

# --- Results -----------------------------------------------------------------
print(f"\nWall time    : {elapsed:.2f}s")
print(f"Best Δv      : {result['best_value']/1000:.3f} km/s")
print(f"Best x       : {result['best_values_sequence']}")
print(f"Improvements : {len(result['best_values_list'])}")
print(f"Evaluations: {udp.count}")

# if result["best_values_sequence"]:
#     udp.pretty(result["best_values_sequence"])

# --- Convergence plot --------------------------------------------------------
import matplotlib.pyplot as plt

if result["time_list"]:
    plt.figure()
    plt.semilogy(
        result["time_list"],
        [v / 1000 for v in result["best_values_list"]],
        marker="o",
        markersize=3,
    )
    plt.xlabel("Time (s)")
    plt.ylabel("Best Δv (km/s)")
    plt.title("CNRPA (C++) convergence on cassini2")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("cnrpa_convergence.png", dpi=150)
    plt.show()
