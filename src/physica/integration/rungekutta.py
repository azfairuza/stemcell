import numpy as np

def eom_rungekutta(position, velocity, force, mass, dt):
    """a function to calculate equation of motion using 4th order 
    Runge-Kutta methods
    
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
    # 1st order runge-kutta
    force_val = force(position, velocity)
    k1v = (force_val/mass)*dt
    k1x = velocity*dt
    # 2nd order runge-kutta
    temp_position = position + (k1x/2)
    temp_velocity = velocity + (k1v/2)
    force_val = force(temp_position, temp_velocity)
    k2v = (force_val/mass)*dt
    k2x = temp_velocity*dt
    # 3rd order runge-kutta
    temp_position = position + (k2x/2)
    temp_velocity = velocity + (k2v/2)
    force_val = force(temp_position, temp_velocity)
    k3v = (force_val/mass)*dt
    k3x = temp_velocity*dt
    # 4rd order runge-kutta
    temp_position = position + k3x
    temp_velocity = velocity + k3v
    force_val = force(temp_position, temp_velocity)
    k4v = (force_val/mass)*dt
    k4x = temp_velocity*dt
    #result
    final_velocity = velocity + (1/6)*(k1v + (2*k2v) + (2*k3v) + k4v)
    final_position = position + (1/6)*(k1x + (2*k2x) + (2*k3x) + k4x)
    return final_position, final_velocity