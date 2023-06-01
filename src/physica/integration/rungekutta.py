"""module for Runge-Kutta integration of equation of motion"""

import numpy as np


def eom_rungekutta(position, velocity, force, mass, timestep):
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
    timestep: float
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
    k1_vel = (force_val / mass) * timestep
    k1_pos = velocity * timestep
    # 2nd order runge-kutta
    temp_position = position + (k1_pos / 2)
    temp_velocity = velocity + (k1_vel / 2)
    force_val = force(temp_position, temp_velocity)
    k2_vel = (force_val / mass) * timestep
    k2_pos = temp_velocity * timestep
    # 3rd order runge-kutta
    temp_position = position + (k2_pos / 2)
    temp_velocity = velocity + (k2_vel / 2)
    force_val = force(temp_position, temp_velocity)
    k3_vel = (force_val / mass) * timestep
    k3_pos = temp_velocity * timestep
    # 4rd order runge-kutta
    temp_position = position + k3_pos
    temp_velocity = velocity + k3_vel
    force_val = force(temp_position, temp_velocity)
    k4_vel = (force_val / mass) * timestep
    k4_pos = temp_velocity * timestep
    # result
    final_velocity = velocity + (1 / 6) * (
        k1_vel + (2 * k2_vel) + (2 * k3_vel) + k4_vel
    )
    final_position = position + (1 / 6) * (
        k1_pos + (2 * k2_pos) + (2 * k3_pos) + k4_pos
    )
    return final_position, final_velocity
