"""module for spring force"""

import numpy as np


def spring(
    position_a,
    position_b,
    velocity_a,
    velocity_b,
    spring_constant,
    normal_length,
    damping_coefficient=0,
):
    """calculate the spring force

    Parameter
    --------
    position_a: array_like
        The coordinate position of object A
    position_b: array_like
        The coordinate position of object B
    velocity_a: array_like
        The velocity of object A
    velocity_b: array_like
        The velocity of object B
    spring_constant: float
        The spring constant
    normal_length: float
        The normal length of the spring (length when force is 0)
    damping_coefficient: float, default=0
        The damping coefficient related to the velocity of the object

    Return
    ------
    F = -kx-bv

    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    velocity_a = np.array(velocity_a)
    velocity_b = np.array(velocity_b)
    dist = np.linalg.norm((position_b - position_a))
    dist_vec = (position_b - position_a) / dist
    rel_velocity = velocity_b - velocity_a
    delta_dist = dist - normal_length
    force = -1 * (
        (spring_constant * delta_dist * dist_vec) + damping_coefficient * rel_velocity
    )
    return force
