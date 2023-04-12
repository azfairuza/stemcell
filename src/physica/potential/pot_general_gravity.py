import numpy as np
from physica.constant import GRAVITY_CONSTANT

def pot_general_gravity(mass_A, mass_B, position_A, position_B):
    """calculate the general gravity potential energy
    
    Parameter
    --------
    mass_A: float
        the mass of object A (in Kg)
    mass_B: float
        the mass of object B (in Kg)
    position_A: array_like
        The coordinate position of object A
    position_B: array_like
        The coordinate position of object B
    
    Return
    ------
    EP = -GMm/r
    
    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    dist = np.linalg.norm((position_B - position_A))
    energy = -1*GRAVITY_CONSTANT*mass_A*mass_B/dist
    return energy
    