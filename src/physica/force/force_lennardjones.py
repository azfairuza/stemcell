"""module for lennard jones forces"""

import numpy as np


def lj_6_12(position_a, position_b, epsilon=1, sigma=1):
    """calculate the Lennard-Jones 6-12 force

    Parameter
    --------
    position_a: array_like
        The coordinate position of object A
    position_b: array_like
        The coordinate position of object B
    epsilon: float, default=1
        The depth of the potential creek/well
    sigma: float, default=1
        The distance when the energy value equal to zero

    Return
    ------
    lennard-jones force with exponent of 12-6

    Notes
    -----
    - This function calculate the force acting on object B from
    object A.
    """
    position_a = np.array(position_a)
    position_b = np.array(position_b)
    dist_vec = position_b - position_a
    dist = np.linalg.norm(dist_vec)
    direction_vec = dist_vec / dist
    alpha = (sigma / dist) ** 6
    force = (48 / dist) * epsilon * alpha * (alpha - 0.5)
    force_vec = force * direction_vec
    return force_vec


def nearest_dist_LJ(epsilon, sigma, limit=10e-6, start_point=10e-3):
    """function to get the perimeter where the force greater than a limit

    parameter
    --------
    epsilon: float
        The depth of the potential creek/well.
    sigma: float
        The distance when the energy value equal to zero.
    limit: float, default 10E-6
        The minimum force.
    start_point: float, default 10E-3
        The start scanning point, 0 is impossible due to the value
        is limit to infinite.

    return
    ------
    the nearest distance from origin until the force become smaller
    than the limit.
    """
    nearest_dist = [start_point]
    searching = True
    while searching is True:
        force = lj_6_12([0.0], nearest_dist, epsilon, sigma)
        if abs(force[0]) < limit:
            searching = False
        else:
            exponent = np.log10(
                abs(force[0])
            )  # the increment based on the force exponent
            increment = min(10 ** (exponent + 3), 1)  # maximum increment is 1
            nearest_dist[0] = nearest_dist[0] + increment
    return nearest_dist[0]
