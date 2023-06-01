"""Module for general gravity force"""

import numpy as np
from ..constant import GRAVITY_CONSTANT


def general_gravity(mass_a, mass_b, position_a, position_b):
    """calculate the general gravity force

    Parameter
    --------
    mass_a: float
        the mass of object A (in Kg)
    mass_b: float
        the mass of object B (in Kg)
    position_a: array_like
        The coordinate position of object A
    position_b: array_like
        The coordinate position of object B

    Return
    ------
    F = -GMm/r^2

    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    dist = np.linalg.norm((position_b - position_a))
    dist_vec = (position_b - position_a) / dist
    force = -1 * GRAVITY_CONSTANT * mass_a * mass_b / (dist**2) * dist_vec
    return force
