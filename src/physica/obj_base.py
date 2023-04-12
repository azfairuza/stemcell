""" objbase module

This module contain the ObjBase class which can be used as a base to
construct an object.
"""

import numpy as np

class ObjBase:
    """Purpose of this class is as follow:

    it can be used for any kind of object base. For any other projects.

    """

    def __init__(self, x: float, y: float, id: int) -> None:
        """initial function for Base class

        Parameters
        ----------
        x: float
            x position of the object
        y: float
            y position of the object
        id: int
            the identity number of the object
        """
        self._position = np.array((x, y))
        self._speed = np.array((0, 0), dtype=float)
        self._acceleration = np.array((0, 0), dtype=float)
        self._force = np.array((0, 0), dtype=float)
        self._id = id
        self._mass = 1
        self._size = 1
        self._shape = "circle"

        # temporary
        self._temp_position = self._position
        self._temp_speed = self._speed
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
        self.speed = self.temp_speed
        self.position = self.temp_position

    def getDistance(self, obj, bias=None):
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
            else:
                return np.linalg.norm(obj.position - self.position)

    def getXDistance(self, obj, bias=0.0):
        """Get the X distance from two objects.

        Parameters
        ----------
        obj : :obj:`ObjBase`
            the object target point.
        bias : float, default: 0.0
            bias value to alter the x position of first object.
        """
        if isinstance(obj, ObjBase):
            return obj.x - (self.x + bias)

    def getYDistance(self, obj, bias=0.0):
        """Get the Y distance from two objects

        Parameters
        ----------
        obj : :obj:`ObjBase`
            the object target point.
        bias : float, default: 0.0
            bias value to alter the y position of first object.
        """
        if isinstance(obj, ObjBase):
            return obj.y - (self.y + bias)

    def getAngle(self, obj, bias=None):
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
                    self.getYDistance(obj, bias[1]), self.getXDistance(obj, bias[0])
                )
            else:
                return np.arctan2(self.getYDistance(obj), self.getXDistance(obj))

    # endregion

    # region <Property>
    @property
    def id(self):
        """The identity number of the object."""
        return self._id

    @property
    def x(self):
        """return the x axis position of the object."""
        return self._position[0]

    @x.setter
    def x(self, value: float):
        self._position[0] = value

    @property
    def y(self):
        """return the y axis position of the object."""
        return self._position[1]

    @y.setter
    def y(self, value: float):
        self._position[1] = value

    @property
    def position(self):
        """return the position (x,y) of the object."""
        return self._position

    @position.setter
    def position(self, coord):
        if isinstance(coord, ObjBase):
            self._position = coord.position
        elif (
            isinstance(coord, tuple)
            or isinstance(coord, list)
            or isinstance(coord, np.ndarray)
        ):
            self.x = coord[0]
            self.y = coord[1]
        else:
            return NotImplemented

    @property
    def vx(self):
        """return the x axis velocity of the object."""
        return self._speed[0]

    @vx.setter
    def vx(self, value):
        self._speed[0] = value

    @property
    def vy(self):
        """return the y axis velocity of the object."""
        return self._speed[1]

    @vy.setter
    def vy(self, value):
        self._speed[1] = value

    @property
    def speed(self):
        """return the velocity (vx, vy) of the object."""
        return self._speed

    @speed.setter
    def speed(self, value):
        if isinstance(value, ObjBase):
            self._speed = value.speed
        elif (
            isinstance(value, tuple)
            or isinstance(value, list)
            or isinstance(value, np.ndarray)
        ):
            self.vx = value[0]
            self.vy = value[1]
        else:
            return NotImplemented

    @property
    def ax(self):
        """return the x axis acceleration of the object."""
        return self._acceleration[0]

    @ax.setter
    def ax(self, value):
        self._acceleration[0] = value

    @property
    def ay(self):
        """return the y axis acceleration of the object."""
        return self._acceleration[1]

    @ay.setter
    def ay(self, value):
        self._acceleration[1] = value

    @property
    def acceleration(self):
        """return the acceleration (ax, ay) of the object."""
        return self._acceleration

    @acceleration.setter
    def acceleration(self, value):
        if isinstance(value, ObjBase):
            self._acceleration = value.acceleration
        elif (
            isinstance(value, tuple)
            or isinstance(value, list)
            or isinstance(value, np.ndarray)
        ):
            self.ax = value[0]
            self.ay = value[1]
        else:
            return NotImplemented

    @property
    def fx(self):
        """return the x axis force acting on the object."""
        return self._force[0]

    @fx.setter
    def fx(self, value):
        self._force[0] = value

    @property
    def fy(self):
        """return the y axis force acting on the object."""
        return self._force[1]

    @fy.setter
    def fy(self, value):
        self._force[1] = value

    @property
    def force(self):
        """return the force (fx, fy) acting on the object."""
        return self._force

    @force.setter
    def force(self, value):
        if isinstance(value, ObjBase):
            self._force = value.force
        elif (
            isinstance(value, tuple)
            or isinstance(value, list)
            or isinstance(value, np.ndarray)
        ):
            self.fx = value[0]
            self.fy = value[1]
        else:
            return NotImplemented

    @property
    def temp_x(self):
        """return the temporary x position of the object."""
        return self._temp_position[0]

    @temp_x.setter
    def temp_x(self, value: float):
        self._temp_position[0] = value

    @property
    def temp_y(self):
        """return the temporary y position of the object."""
        return self._temp_position[1]

    @temp_y.setter
    def temp_y(self, value: float):
        self._temp_position[1] = value

    @property
    def temp_position(self):
        """return the temporary position (x,y) of the object."""
        return self._temp_position

    @temp_position.setter
    def temp_position(self, coord):
        if (
            isinstance(coord, tuple)
            or isinstance(coord, list)
            or isinstance(coord, np.ndarray)
        ):
            self.temp_x = coord[0]
            self.temp_y = coord[1]
        else:
            return NotImplemented

    @property
    def temp_vx(self):
        """return the temporary x component velocity of the object."""
        return self._temp_speed[0]

    @temp_vx.setter
    def temp_vx(self, value):
        self._temp_speed[0] = value

    @property
    def temp_vy(self):
        """return the temporary y component velocity of the object."""
        return self._temp_speed[1]

    @temp_vy.setter
    def temp_vy(self, value):
        self._temp_speed[1] = value

    @property
    def temp_speed(self):
        """return the temporary velocity of the object."""
        return self._temp_speed

    @temp_speed.setter
    def temp_speed(self, value):
        if isinstance(value, ObjBase):
            self._temp_speed = value.temp_speed
        elif (
            isinstance(value, tuple)
            or isinstance(value, list)
            or isinstance(value, np.ndarray)
        ):
            self.temp_vx = value[0]
            self.emp_vy = value[1]
        else:
            return NotImplemented

    @property
    def temp_ax(self):
        """return the temporary x component acceleration of the object."""
        return self._temp_acceleration[0]

    @temp_ax.setter
    def temp_ax(self, value):
        self._temp_acceleration[0] = value

    @property
    def temp_ay(self):
        """return the temporary y component acceleration of the object."""
        return self._temp_acceleration[1]

    @temp_ay.setter
    def temp_ay(self, value):
        self._temp_acceleration[1] = value

    @property
    def temp_acceleration(self):
        """return the temporary acceleration (x, y) of the object."""
        return self._temp_acceleration

    @temp_acceleration.setter
    def temp_acceleration(self, value):
        if isinstance(value, ObjBase):
            self._temp_acceleration = value.temp_acceleration
        elif (
            isinstance(value, tuple)
            or isinstance(value, list)
            or isinstance(value, np.ndarray)
        ):
            self.temp_ax = value[0]
            self.temp_ay = value[1]
        else:
            return NotImplemented

    @property
    def temp_fx(self):
        """return the temporary x component force acting on the
        object.
        """
        return self._temp_force[0]

    @temp_fx.setter
    def temp_fx(self, value):
        self._temp_force[0] = value

    @property
    def temp_fy(self):
        """return the temporary y component force acting on the
        object.
        """
        return self._temp_force[1]

    @temp_fy.setter
    def temp_fy(self, value):
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
        elif (
            isinstance(value, tuple)
            or isinstance(value, list)
            or isinstance(value, np.ndarray)
        ):
            self.temp_fx = value[0]
            self.temp_fy = value[1]
        else:
            return NotImplemented

    @property
    def size(self):
        """return the size of the object"""
        return self._size
    
    @property
    def mass(self):
        """return the mass of the object"""
        return self._mass

    # endregion
