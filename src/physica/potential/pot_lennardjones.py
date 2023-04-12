import numpy as np

def pot_LJ_6_12(position_A, position_B, epsilon=1, sigma=1):
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
    $V = 4\epsilon((\frac{\sigma}{r})^{12}-
        (\frac{\sigma}{r})^{6})$
    
    Notes
    -----
    - This function calculate the potential energy of object B due to 
    force field of object A.
    """
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    dist_vec = position_B - position_A
    dist = np.linalg.norm(dist_vec)
    alpha = (sigma/dist)**6
    energy = 4*epsilon*alpha*(alpha - 1)
    return energy