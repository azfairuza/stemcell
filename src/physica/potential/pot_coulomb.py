"""module for coulomb potential"""

import numpy as np
from ..constant import COULOMB_CONSTANT


def coulomb(charge_a, charge_b, position_a, position_b):
    """calculate the coulombic potential energy

    Parameter
    --------
    charge_a: float
        the electronic charge of object A (in C)
    charge_b: float
        the electronic charge of object B (in C)
    position_a: array_like
        The coordinate position of object A
    position_b: array_like
        The coordinate position of object B

    Return
    ------
    E = Kq1q2/r

    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    dist = np.linalg.norm((position_b - position_a))
    energy = COULOMB_CONSTANT * charge_a * charge_b / dist
    return energy
