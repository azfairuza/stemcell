"""module for coulombic force"""
import numpy as np
from ..constant import COULOMB_CONSTANT


def coulomb(charge_a, charge_b, position_a, position_b):
    """calculate the coulombic force

    Parameter
    --------
    charge_A: float
        the electronic charge of object A (in C)
    charge_B: float
        the electronic charge of object B (in C)
    position_A: array_like
        The coordinate position of object A
    position_B: array_like
        The coordinate position of object B

    Return
    ------
    F = Kq1q2/r^2

    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    dist = np.linalg.norm((position_b - position_a))
    dist_vec = (position_b - position_a) / dist
    force = (COULOMB_CONSTANT * charge_a * charge_b / (dist**2)) * dist_vec
    return force
