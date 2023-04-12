import numpy as np

def force_spring(position_A, 
                 position_B,
                 velocity_A,
                 velocity_B, 
                 spring_constant,
                 normal_length,
                 damping_coefficient=0):
    """calculate the spring force
    
    Parameter
    --------
    position_A: array_like
        The coordinate position of object A
    position_B: array_like
        The coordinate position of object B
    velocity_A: array_like
        The velocity of object A
    velocity_B: array_like
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
    position_A = np.array(position_A)
    position_B = np.array(position_B)
    velocity_A = np.array(velocity_A)
    velocity_B = np.array(velocity_B)
    dist = np.linalg.norm((position_B - position_A))
    dist_vec = (position_B - position_A)/dist
    rel_velocity = velocity_B - velocity_A
    delta_dist = dist - normal_length
    force = -1*((spring_constant*delta_dist*dist_vec) + damping_coefficient*rel_velocity)
    return force

