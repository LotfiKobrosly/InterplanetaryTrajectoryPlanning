"""
Made by Claude
Run GACO-CABGNRPA (C++ with native pagmo) on pk.trajopt.gym.cassini1/cassini2.
GACO now runs as pure C++ pagmo — no Python roundtrip for evolve().
Only fitness() callbacks re-enter Python (unavoidable: evaluator is pykep).
"""

import time
import pykep as pk
import pygmo as pg  # only needed to get bounds cleanly
import cnrpa_cpp
from utils.constants import UNFEASIBILITY_VALUE, GAUSSIAN_KERNEL_THRESHOLD, RANDOM_SEED

# --- Problem -----------------------------------------------------------------
udp = pk.trajopt.gym.cassini1  # swap to cassini2 for harder problem
prob = pg.problem(udp)
lb, ub = prob.get_bounds()
bounds = [[float(lo), float(hi)] for lo, hi in zip(lb, ub)]

print(f"Problem   : {prob.get_name()}")
print(f"Variables : {prob.get_nx()}")


# --- Run C++ GACO-CABGNRPA (pagmo native) ------------------------------------
t0 = time.time()
result = cnrpa_cpp.gaco_cabgnrpa(
    udp=udp,
    level=2,
    n_policies=100,
    zeta=0.2,
    bounds=bounds,
    learning_rate=0.25,
    timeout=60.0,
    tau=1.3,
    unfeasibility_value=UNFEASIBILITY_VALUE,
    gaussian_kernel_threshold=GAUSSIAN_KERNEL_THRESHOLD,
    kernel_size=63,
    n_generations=10,
    elitism_factor=0.2,
    random_seed=RANDOM_SEED,
)
elapsed = time.time() - t0

# --- Results -----------------------------------------------------------------
print(f"\nWall time    : {elapsed:.2f}s")
print(f"Best Δv      : {result['best_value']/1000:.3f} km/s")
print(f"Improvements : {len(result['best_values_list'])}")

if result["best_values_sequence"]:
    udp.pretty(result["best_values_sequence"])

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
    plt.title("GACO-CABGNRPA (C++ / native pagmo)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("gaco_cabgnrpa_convergence.png", dpi=150)
    plt.show()
