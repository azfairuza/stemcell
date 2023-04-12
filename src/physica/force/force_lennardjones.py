import numpy as np
import src.physica as psc

def force_LJ_6_12(position_A, position_B, epsilon=1, sigma=1):
    """calculate the Lennard-Jones 6-12 force
    
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
    $F = (\frac{48\epsilon}{r})*((\frac{\sigma}{r})^{12}-
        (\frac{\sigma}{r})^{6})$
    
    Notes
    -----
    - This function calculate the force acting on object B from
    object A.
    """
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    dist_vec = position_B - position_A
    dist = np.linalg.norm(dist_vec)
    direction_vec = dist_vec/dist
    alpha = (sigma/dist)**6
    force = (48/dist)*epsilon*alpha*(alpha - 0.5)
    force_vec = force*direction_vec
    return force_vec

def nearest_dist_LJ(epsilon, sigma, limit=10E-6, start_point=10E-3):
    x = [start_point]
    searching = True
    while searching is True:
        force = psc.force_LJ_6_12([0], x, epsilon, sigma)
        if abs(force[0]) < limit:
            searching = False
        else:
            exponent = np.log10(abs(force[0]))
            increment = min(10**(exponent+3), 1)
            x[0] = x[0] + increment
    return x[0]


