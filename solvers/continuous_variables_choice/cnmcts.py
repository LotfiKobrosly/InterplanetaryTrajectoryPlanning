"""
Implements a Continous Nested Monte Carlo Search for the continuous decision variables
Here, the values sequence is represented as follows:
values_sequence = [departure_epoch, time_of_flight_0, ..., time_of_flight_n]
"""

import time
import numpy as np
import pykep as pk
from utils.constants import RANDOM_GENERATOR, DV_LAUNCHER, UNFEASIBILITY_VALUE


def run_cnmcts(
    evaluator: pk.trajopt.mga,
    values_sequence: list = list(),
    bounds: list = None,
    level: int = 0,
    bandwidth: int = 10,
    timeout: float = 10.0,
    start_time: float = 0,
    best_sequence: list = None,
    best_value: float = UNFEASIBILITY_VALUE,
    best_values_list: list = None,
    time_list: list = None,
    *args,
    **kwargs,
):
    current_time = time.time() - start_time
    if level == 0:
        while len(values_sequence) < len(bounds):
            values_sequence.append(
                RANDOM_GENERATOR.uniform(*bounds[len(values_sequence)])
            )
        return values_sequence, evaluator.fitness(values_sequence)[0]

    else:
        while len(values_sequence) < len(bounds) and (current_time < timeout):
            temporary_values_sequences = [values_sequence[:] for _ in range(bandwidth)]
            for new_values_sequence in temporary_values_sequences:

                new_values_sequence.append(
                    RANDOM_GENERATOR.uniform(*bounds[len(values_sequence)])
                )
                new_values_sequence, result = run_cnmcts(
                    evaluator=evaluator,
                    values_sequence=new_values_sequence,
                    bounds=bounds,
                    level=level - 1,
                    bandwidth=bandwidth,
                    timeout=timeout,
                    start_time=start_time,
                    best_value=best_value,
                    best_sequence=best_sequence,
                    best_values_list=best_values_list,
                    time_list=time_list,
                )
                current_time = time.time() - start_time
                if result < best_value:
                    best_sequence = new_values_sequence[:]
                    best_value = result
                if best_value < UNFEASIBILITY_VALUE:
                    best_values_list.append(best_value)
                    time_list.append(current_time)
                if current_time > timeout:
                    break
            values_sequence.append(best_sequence[len(values_sequence)])
        return (best_sequence, best_value)


def cnmcts(
    evaluator: pk.trajopt.mga,
    bounds: list = None,
    level: int = 0,
    bandwidth: int = 10,
    timeout: float = 10.0,
    *args,
    **kwargs,
):
    start_time = time.time()
    best_values_list, time_list = list(), list()
    best_sequence, best_value = run_cnmcts(
        evaluator=evaluator,
        bounds=bounds,
        values_sequence=list(),
        level=level,
        bandwidth=bandwidth,
        timeout=timeout,
        start_time=start_time,
        best_sequence=None,
        best_value=UNFEASIBILITY_VALUE,
        best_values_list=best_values_list,
        time_list=time_list,
    )

    return best_sequence, best_value, best_values_list, time_list


if __name__ == "__main__":
    # Cassini problem
    udp = pk.trajopt.gym.cassini2

    # Variables bounds
    bounds = [
        (low_bound, high_bound)
        for (low_bound, high_bound) in zip(udp.get_bounds()[0], udp.get_bounds()[1])
    ]

    # General input values
    inputs_values = {
        "evaluator": udp,
        "bounds": bounds,
        "timeout": 10,
        "level": 2,
        "bandwidth": 200,
    }
    values__sequence, best_value, values_list, time_list = cnmcts(**inputs_values)
    print(f"Best Delta V: {best_value / 1000:.3f} km/s")
    print(f"Total time: {time_list[-1]:.2f} s")
