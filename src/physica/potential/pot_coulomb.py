import numpy as np
from physica.constant import COULOMB_CONSTANT

def pot_coulomb(charge_A, charge_B, position_A, position_B):
    """calculate the coulombic potential energy
    
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
    E = Kq1q2/r
    
    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    dist = np.linalg.norm((position_B - position_A))
    energy = COULOMB_CONSTANT*charge_A*charge_B/dist
    return energy