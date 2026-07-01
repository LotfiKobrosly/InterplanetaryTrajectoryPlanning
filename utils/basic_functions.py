"""
Implementing basic functions
"""

import numpy as np
from utils.constants import RANDOM_GENERATOR


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


def fit_gaussian_from_density(x, density_values):
    """
    From Claude LLM
    x               : equally spaced points
    density_values  : unnormalized density/probability values at each point
    """
    x = np.array(x, dtype=float)
    p = np.array(density_values, dtype=float)

    # Normalize to a proper probability mass function
    p_norm = p / p.sum()

    mu = np.sum(p_norm * x)
    var = np.sum(p_norm * (x - mu) ** 2)
    sigma = np.sqrt(var)

    return mu, sigma


def sample_mixture_1d(n_samples, mu1, sigma1, mu2, sigma2, weight1=0.5):
    # From Claude LLM
    component = RANDOM_GENERATOR.uniform(0, 1, n_samples) < weight1
    samples = np.where(
        component,
        RANDOM_GENERATOR.normal(mu1, sigma1, n_samples),
        RANDOM_GENERATOR.normal(mu2, sigma2, n_samples),
    )
    return samples
