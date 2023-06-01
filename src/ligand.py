""" ligand module

This module contain the ligand class which is main component of 
ligand-integrin base stem cell simulation.
"""

import physica as psc
import integrin as ign


class Ligand(psc.ObjBase):
    """Ligand description:

    Component which bounds to integrin protein in cell membrane.
    """

    count = 0

    def __init__(
        self, x_position: float, y_position: float, size: float, mass: float = 1
    ):
        """initial function for Ligand class

        Parameters
        ----------
        x_position: float
            x position of the ligand
        y_position: float
            y position of the ligand
        size: float
            the size of the ligand
        mass: float
            the mass of the ligand
        """
        self.__class__.count += 1
        super().__init__(x_position, y_position, self.__class__.count)
        self._size = size
        self._mass = mass
        self._bound = False
        self._target: ign.Integrin = None
        self._target_integrin_id = None
        self._target_cell_id = None

    @classmethod
    def reset_count(cls):
        """Reset the number of ligand created into 0."""
        cls.count = 0

    @property
    def bound(self):
        """ligand condition related to bonding status with other
        object.
        """
        return self._bound

    @bound.setter
    def bound(self, value):
        if isinstance(value, bool):
            self._bound = value

    @property
    def target(self):
        """the target object from the ligand"""
        return self._target

    @target.setter
    def target(self, obj_target):
        if isinstance(obj_target, psc.ObjBase):
            self._target = obj_target

    @property
    def target_cell_id(self):
        """the target cell id"""
        return self._target._cell.id_
    
    @property
    def target_integrin_id(self):
        """the target integrin id"""
        return self._target.id_