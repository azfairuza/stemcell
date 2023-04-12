import numpy as np

def pot_gravity(mass, gravity, position):
    """calculate the gravitational potential near surface
    
    Parameter
    --------
    mass: float
        mass of the object.
    gravity: array_like
        the gravitational acceleration vector in the shape of
        array_like. The direction should be in reverse with 
        surfave normal vector.
    position: array_like
        the position of the object relative to the surface.
    
    Return
    ------
    EP = mgh

    Notes
    -----
    Use this function to calculate only near surface potential
    energy.
    """
    position = np.array(position)
    gravity = np.array(gravity)
    energy = mass*np.dot(gravity, position)
    return abs(energy)