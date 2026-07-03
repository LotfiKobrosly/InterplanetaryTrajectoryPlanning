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
    planets_sequence: list = None,
    level: int = 0,
    bandwidth: int = 10,
    timeout: float = 10.0,
    start_time: float = 0,
    best_sequence: list = None,
    best_value: float = UNFEASIBILITY_VALUE,
    best_values_list: list = None,
    time_list: list = None,
    *args,
    **kwargs
):
    current_time = time.time() - start_time
    if level == 0:
        while len(values_sequence) < len(planets_sequence):
            values_sequence.append(
                RANDOM_GENERATOR.uniform(*bounds[len(values_sequence)])
            )
        return values_sequence, evaluator.fitness(values_sequence)[0]
        
    else:
        while len(values_sequence) < len(planets_sequence) and (current_time < timeout):
            temporary_values_sequences = [values_sequence[:] for _ in range(bandwidth)]
            for new_values_sequence in temporary_values_sequences:

                new_values_sequence.append(
                    RANDOM_GENERATOR.uniform(*bounds[len(values_sequence)])
                )
                new_values_sequence, result = run_cnmcts(
                    evaluator=evaluator,
                    values_sequence=new_values_sequence,
                    bounds=bounds,
                    planets_sequence=planets_sequence,
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
                    best_sequence = new_values_sequence
                    best_value = result
                if best_value < UNFEASIBILITY_VALUE:
                    best_values_list.append(best_value)
                    time_list.append(current_time)
            values_sequence.append(best_sequence[len(values_sequence)])
        return (best_sequence, best_value)

def cnmcts(
    evaluator: pk.trajopt.mga,
    bounds: list = None,
    planets_sequence: list = None,
    level: int = 0,
    bandwidth: int = 10,
    timeout: float = 10.0,
    *args,
    **kwargs
):
    start_time = time.time()
    best_values_list, time_list = list(), list()
    best_sequence, best_value = run_cnmcts(
        evaluator=evaluator,
        bounds=bounds,
        values_sequence=list(),
        planets_sequence=planets_sequence,
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
