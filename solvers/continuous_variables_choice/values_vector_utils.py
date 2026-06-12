"""
Values vector separation file
"""
def separate_values(input_vector: list, planets_sequence_length: int):
    """
    Simply seperates the values vector into the desired variables compatible with Trajectory class from classes.trajectory
    """
    departure_epoch = input_vector[0]
    times_of_flight = input_vector[1:planets_sequence_length]
    flyby_parameters = input_vector[planets_sequence_length:]
    flyby_parameters = [
        (flyby_parameters[i], flyby_parameters[i + 1])
        for i in range(0, len(flyby_parameters), 2)
    ]
    return departure_epoch, times_of_flight, flyby_parameters