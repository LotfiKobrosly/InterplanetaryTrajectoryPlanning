# Author: Lotfi Kobrosly (2025)

from jplephem.spk import SPK
from utils.constants import INDICES


def get_ephemeris(kernel: SPK, body: str, t: float):
    # Return (r, v) defining position and velocity of body at time t from the ephemerides source defined by kernel

    return kernel[*INDICES[body]].compute_and_differentiate(t)
