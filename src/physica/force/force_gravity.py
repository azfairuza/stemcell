"""module for ordinary gravity force"""

import numpy as np


def gravity(mass: float, gravity_acceleration: np.ndarray):
    """calculate the gravitational force

    Parameter
    --------
    mass: float
        mass of the object
    gravity_acceleration: array_like
        the gravitational acceleration vector in the shape of
        array_like.

    Return
    ------
    m x g
    """
    gravity_acceleration = np.array(gravity_acceleration)
    return mass * gravity_acceleration
