import numpy as np
from physica.constant import COULOMB_CONSTANT

def force_coulomb(charge_A, charge_B, position_A, position_B):
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
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    dist = np.linalg.norm((position_B - position_A))
    dist_vec = (position_B - position_A)/dist
    force = (COULOMB_CONSTANT*charge_A*charge_B/(dist**2))*dist_vec
    return force