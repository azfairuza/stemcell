import numpy as np

def eom_euler(position, velocity, force, mass, dt):
    """a function to calculate equation of motion using eulerian 
    methods
    
    Parameters
    ----------
    position: array_like 
        the current position of object. it must be written in 
        array_like shape.
    velocity: array_like
        the current velocity of object. it must be written in 
        array_like shape.
    force: function
        the current force acting on the object. it must be written
        in as the function of (x, v).
    mass: float
        the mass of the object
    dt: float
        time step of integration
    
    Returns
    -------
    final_velocity: np.ndarray
        the calculated final velocity in the span of `dt`
    final_position: np.ndarray
        the calculated final position in the span of `dt` 
    """
    position = np.array(position)
    velocity = np.array(velocity)
    force_val = force(position, velocity)
    acceleration = force_val/mass
    final_velocity = velocity + acceleration*dt
    final_position = position + velocity*dt
    return final_position, final_velocity

