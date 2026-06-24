"""
Implements a gaussian kernel
"""

import numpy as np


class GaussianKernel:  # ChatGPT
    def __init__(self, center, sigma):
        """
        center: array-like -> mean of the Gaussian
        sigma: float -> standard deviation (same for both dimensions)
        """
        if isinstance(center, (float, int)):
            center = [center]
        self.mu = np.array(center, dtype=float)
        self.sigma = float(sigma)
        self.dim = np.shape(self.mu)[0]

        # Precompute constants
        self.norm_const = 1 / ((2 * np.pi * self.sigma) ** (self.dim / 2))

    def pdf(self, x):
        """
        Compute probability density at point x.
        x: array-like of shape (2,)
        """
        if isinstance(x, (float, int)):
            x = [x]
        x = np.array(x, dtype=float)
        diff = x - self.mu
        exponent = -0.5 * np.dot(diff, diff) / (self.sigma**2)
        return self.norm_const * np.exp(exponent)
