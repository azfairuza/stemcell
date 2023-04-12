import numpy as np

def pot_spring(position_A, position_B, spring_constant, normal_length):
    """calculate the spring potential
    
    Parameter
    --------
    position_A: array_like
        The coordinate position of object A
    position_B: array_like
        The coordinate position of object B
    spring_constant: float
        The spring constant
    normal_length: float
        The normal length of the spring (length when force is 0)
    
    Return
    ------
    F = (1/2)*k*x^2
    
    Notes
    -----
    This function calculate the force acting on object B from
    object A. So the main reference is object A
    """
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    dist = np.linalg.norm(position_B - position_A)
    dist_diff = dist - normal_length
    energy = 0.5*spring_constant*(dist_diff**2)
    return energy