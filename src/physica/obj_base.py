""" objbase module

This module contain the ObjBase class which can be used as a base to
construct an object.
"""

import numpy as np


class ObjBase:
    """Purpose of this class is as follow:

    it can be used for any kind of object base. For any other projects.

    """

    def __init__(self, x_poisition: float, y_position: float, id_: int) -> None:
        """initial function for Base class

        Parameters
        ----------
        x_position: float
            x position of the object
        y_position: float
            y position of the object
        id_: int
            the identity number of the object
        """
        self._position = np.array((x_poisition, y_position))
        self._velocity = np.array((0, 0), dtype=float)
        self._acceleration = np.array((0, 0), dtype=float)
        self._force = np.array((0, 0), dtype=float)
        self._id = id_
        self._mass = 1
        self._size = 1
        self._shape = "circle"

        # temporary
        self._temp_position = self._position
        self._temp_velocity = self._velocity
        self._temp_acceleration = self._acceleration
        self._temp_force = self._force

    # region <Method>
    def update(self):
        """Update the object state.

        This function is invoked in the end calculation to update the
        state of the object. Manual update can be used as well instead
        of using this function.
        """
        self.acceleration = self.temp_acceleration
        self.velocity = self.temp_velocity
        self.position = self.temp_position

    def get_distance(self, obj, bias=None):
        """Get the distance from two objects.

        Parameters
        ----------
        obj : :obj:`ObjBase`
            the object target point.
        bias : :obj:`ndarray`, default : `None`
            It should contain 2 dimension ndarray which corresponds to
            (x,y) bias.
        """
        if isinstance(obj, ObjBase):
            if isinstance(bias, np.ndarray):
                return np.linalg.norm(obj.position - (self.position + bias))
            return np.linalg.norm(obj.position - self.position)
        return None

    def get_x_distance(self, obj, bias=0.0):
        """Get the X distance from two objects.

        Parameters
        ----------
        obj : :obj:`ObjBase`
            the object target point.
        bias : float, default: 0.0
            bias value to alter the x position of first object.
        """
        if isinstance(obj, ObjBase):
            return obj.x_position - (self.x_position + bias)
        return None

    def get_y_distance(self, obj, bias=0.0):
        """Get the Y distance from two objects

        Parameters
        ----------
        obj : :obj:`ObjBase`
            the object target point.
        bias : float, default: 0.0
            bias value to alter the y position of first object.
        """
        if isinstance(obj, ObjBase):
            return obj.y_position - (self.y_position + bias)
        return None

    def get_angle(self, obj, bias=None):
        """Get the angle of two object in respect with X-axis

        Parameters
        ----------
        obj : :obj:`ObjBase`
            the object target point.
        bias : :obj: `ndarray`, default: none
            bias value to alter theposition of first object.

        Return
        ------
        The output is in radians
        """
        if isinstance(obj, ObjBase):
            if isinstance(bias, np.ndarray):
                return np.arctan2(
                    self.get_y_distance(obj, bias[1]), self.get_x_distance(obj, bias[0])
                )
            return np.arctan2(self.get_y_distance(obj), self.get_x_distance(obj))
        return None

    # endregion

    # region <Property>
    @property
    def id_(self):
        """The identity number of the object."""
        return self._id

    @property
    def x_position(self):
        """return the x axis position of the object."""
        return self._position[0]

    @x_position.setter
    def x_position(self, value: float):
        self._position[0] = value

    @property
    def y_position(self):
        """return the y axis position of the object."""
        return self._position[1]

    @y_position.setter
    def y_position(self, value: float):
        self._position[1] = value

    @property
    def position(self):
        """return the position (x,y) of the object."""
        return self._position

    @position.setter
    def position(self, coord):
        if isinstance(coord, ObjBase):
            self._position = coord.position
        elif isinstance(coord, (tuple, list, np.ndarray)):
            self.x_position = coord[0]
            self.y_position = coord[1]

    @property
    def x_velocity(self):
        """return the x axis speed of the object."""
        return self._velocity[0]

    @x_velocity.setter
    def x_velocity(self, value):
        self._velocity[0] = value

    @property
    def y_velocity(self):
        """return the y axis speed of the object."""
        return self._velocity[1]

    @y_velocity.setter
    def y_velocity(self, value):
        self._velocity[1] = value

    @property
    def velocity(self):
        """return the velocity (vx, vy) of the object."""
        return self._velocity

    @velocity.setter
    def velocity(self, value):
        if isinstance(value, ObjBase):
            self._velocity = value.velocity
        elif isinstance(value, (tuple, list, np.ndarray)):
            self.x_velocity = value[0]
            self.y_velocity = value[1]

    @property
    def x_acceleration(self):
        """return the x axis acceleration of the object."""
        return self._acceleration[0]

    @x_acceleration.setter
    def x_acceleration(self, value):
        self._acceleration[0] = value

    @property
    def y_acceleration(self):
        """return the y axis acceleration of the object."""
        return self._acceleration[1]

    @y_acceleration.setter
    def y_acceleration(self, value):
        self._acceleration[1] = value

    @property
    def acceleration(self):
        """return the acceleration (ax, ay) of the object."""
        return self._acceleration

    @acceleration.setter
    def acceleration(self, value):
        if isinstance(value, ObjBase):
            self._acceleration = value.acceleration
        elif isinstance(value, (tuple, list, np.ndarray)):
            self.x_acceleration = value[0]
            self.y_acceleration = value[1]

    @property
    def x_force(self):
        """return the x axis force acting on the object."""
        return self._force[0]

    @x_force.setter
    def x_force(self, value):
        self._force[0] = value

    @property
    def y_force(self):
        """return the y axis force acting on the object."""
        return self._force[1]

    @y_force.setter
    def y_force(self, value):
        self._force[1] = value

    @property
    def force(self):
        """return the force (fx, fy) acting on the object."""
        return self._force

    @force.setter
    def force(self, value):
        if isinstance(value, ObjBase):
            self._force = value.force
        elif isinstance(value, (tuple, list, np.ndarray)):
            self.x_force = value[0]
            self.y_force = value[1]

    @property
    def temp_x_position(self):
        """return the temporary x position of the object."""
        return self._temp_position[0]

    @temp_x_position.setter
    def temp_x_position(self, value: float):
        self._temp_position[0] = value

    @property
    def temp_y_position(self):
        """return the temporary y position of the object."""
        return self._temp_position[1]

    @temp_y_position.setter
    def temp_y_position(self, value: float):
        self._temp_position[1] = value

    @property
    def temp_position(self):
        """return the temporary position (x,y) of the object."""
        return self._temp_position

    @temp_position.setter
    def temp_position(self, coord):
        if isinstance(coord, (tuple, list, np.ndarray)):
            self.temp_x_position = coord[0]
            self.temp_y_position = coord[1]

    @property
    def temp_x_velocity(self):
        """return the temporary x component velocity of the object."""
        return self._temp_velocity[0]

    @temp_x_velocity.setter
    def temp_x_velocity(self, value):
        self._temp_velocity[0] = value

    @property
    def temp_y_velocity(self):
        """return the temporary y component velocity of the object."""
        return self._temp_velocity[1]

    @temp_y_velocity.setter
    def temp_y_velocity(self, value):
        self._temp_velocity[1] = value

    @property
    def temp_velocity(self):
        """return the temporary velocity of the object."""
        return self._temp_velocity

    @temp_velocity.setter
    def temp_velocity(self, value):
        if isinstance(value, ObjBase):
            self._temp_velocity = value.temp_velocity
        elif isinstance(value, (tuple, list, np.ndarray)):
            self.temp_x_velocity = value[0]
            self.temp_y_velocity = value[1]

    @property
    def temp_x_acceleration(self):
        """return the temporary x component acceleration of the object."""
        return self._temp_acceleration[0]

    @temp_x_acceleration.setter
    def temp_x_acceleration(self, value):
        self._temp_acceleration[0] = value

    @property
    def temp_y_acceleration(self):
        """return the temporary y component acceleration of the object."""
        return self._temp_acceleration[1]

    @temp_y_acceleration.setter
    def temp_y_acceleration(self, value):
        self._temp_acceleration[1] = value

    @property
    def temp_acceleration(self):
        """return the temporary acceleration (x, y) of the object."""
        return self._temp_acceleration

    @temp_acceleration.setter
    def temp_acceleration(self, value):
        if isinstance(value, ObjBase):
            self._temp_acceleration = value.temp_acceleration
        elif isinstance(value, (tuple, list, np.ndarray)):
            self.temp_x_acceleration = value[0]
            self.temp_y_acceleration = value[1]

    @property
    def temp_x_force(self):
        """return the temporary x component force acting on the
        object.
        """
        return self._temp_force[0]

    @temp_x_force.setter
    def temp_x_force(self, value):
        self._temp_force[0] = value

    @property
    def temp_y_force(self):
        """return the temporary y component force acting on the
        object.
        """
        return self._temp_force[1]

    @temp_y_force.setter
    def temp_y_force(self, value):
        self._temp_force[1] = value

    @property
    def temp_force(self):
        """return the temporary force (fx, fy) acting on the
        object.
        """
        return self._temp_force

    @temp_force.setter
    def temp_force(self, value):
        if isinstance(value, ObjBase):
            self._temp_force = value.temp_force
        elif isinstance(value, (tuple, list, np.ndarray)):
            self.temp_x_force = value[0]
            self.temp_y_force = value[1]

    @property
    def size(self):
        """return the size of the object"""
        return self._size

    @property
    def mass(self):
        """return the mass of the object"""
        return self._mass

    # endregion
