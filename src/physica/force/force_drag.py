import numpy as np


def drag(obj_velocity, obj_radius, viscocity):
    """calculate the drag force

    Parameter
    --------
    obj_velocity: array_like
        the velocity vector of the object
    obj_radius: float
        the radius of the object
    viscocity: array_like
        The coordinate position of object A

    Return
    ------
    F = -6*pi*viscocity*radius*velocity

    """
    
    obj_velocity = np.array(obj_velocity)
    force = -6*np.pi*viscocity*obj_radius*obj_velocity
    return force