"""
Implementing basic functions
"""

import numpy as np


def normalize(x, low, high):
    if isinstance(x, (float, int)):
        return (x - low) / (high - low)
    return [(xi - l) / (h - l) for xi, l, h in zip(x, low, high)]


def denormalize(x, low, high):
    if isinstance(x, (float, int)):
        return x * (high - low) + low
    return [xi * (h - l) + l for xi, l, h in zip(x, low, high)]


def code(input_value):
    if input_value is None:
        return None
    if isinstance(input_value, (float, int)):
        return round(input_value, 1)
    return tuple([round(component, 1) for component in input_value])


def cosine_similarity(vector_1: np.ndarray, vector_2: np.ndarray) -> float:
    return (vector_1 @ vector_2.T) / (
        np.linalg.norm(vector_1) * np.linalg.norm(vector_2)
    )


def truncate(value, min_value, max_value):
    if isinstance(value, float):
        if value < min_value:
            return min_value
        if value > max_value:
            return max_value
        return value
    else:
        value = np.where(value > min_value, value, min_value)
        value = np.where(value < max_value, value, max_value)
        return value
