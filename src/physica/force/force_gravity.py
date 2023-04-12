import numpy as np

def force_gravity(m, g):
    """calculate the gravitational force
    
    Parameter
    --------
    m: float
        mass of the object 
    g: array_like
        the gravitational acceleration vector in the shape of
        array_like.
    
    Return
    ------
    m x g
    """
    g = np.array(g)
    return m*g