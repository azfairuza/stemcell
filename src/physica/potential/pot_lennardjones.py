"""module for Lennard-Jones potential"""

import numpy as np


def lennardjones_6_12(position_a, position_b, epsilon=1, sigma=1):
    """calculate the Lennard-Jones 6-12 potential

    Parameter
    --------
    position_A: array_like
        The coordinate position of object A
    position_B: array_like
        The coordinate position of object B
    epsilon: float, default=1
        The depth of the potential creek/well
    sigma: float, default=1
        The distance when the energy value equal to zero

    Return
    ------
    the potential energy of lennard-jones potential with exponent of
    12-6.

    Notes
    -----
    - This function calculate the potential energy of object B due to
    force field of object A.
    """
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    dist_vec = position_b - position_a
    dist = np.linalg.norm(dist_vec)
    alpha = (sigma / dist) ** 6
    energy = 4 * epsilon * alpha * (alpha - 1)
    return energy
