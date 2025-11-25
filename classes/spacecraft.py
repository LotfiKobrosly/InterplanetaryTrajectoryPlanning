# Author: Lotfi Kobrosly (2025)

# This file defines the class that represents the spacecraft

import numpy as np
from utils.constants import *


class SpaceCraft(object):

    def __init__(
        self,
        dry_mass: float,
        fuel_mass: float,
    ):
        self.dry_mass = dry_mass
        self.fuel_mass = fuel_mass

    def get_total_mass(self):
        return self.dry_mass + self.fuel_mass

    def propel(self, delta_velocity: float, isp: float):
        """
        Compute final spacecraft mass after a delta-velocity maneuver.

        Parameters
        ----------
        delta_velocity : float
            Required velocity change (m/s)
        isp : float
            Specific impulse of the engine (seconds)

        Returns
        -------
        None
        """
        effective_exhaust_velocity = isp * EARTH_GRAVITY_CONSTANT
        final_mass = self.get_total_mass() * np.exp(
            -delta_velocity / effective_exhaust_velocity
        )
        self.fuel_mass = final_mass - self.dry_mass
