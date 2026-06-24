"""
Values vector separation file
"""


def separate_values(input_vector: list):
    """
    Simply seperates the values vector into the desired variables compatible with Trajectory class from classes.trajectory
    """
    departure_epoch = input_vector[0]
    times_of_flight = input_vector[1:]
    return departure_epoch, times_of_flight
