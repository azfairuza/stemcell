""" ligand module

This module contain the ligand class which is main component of 
ligand-integrin base stem cell simulation.

#FINISHED CHECKED
"""

import src.physica as psc


class Ligand(psc.ObjBase):
    """Ligand description:
    
    Component which bounds to integrin protein in cell membrane.
    """

    count = 0

    def __init__(self, x: float, y: float, size: float, mass: float = 1):
        """initial function for Ligand class

        Parameters
        ----------
        x: float
            x position of the ligand
        y: float
            y position of the ligand
        size: float
            the size of the ligand
        mass: float
            the mass of the ligand
        """
        self.__class__.count += 1
        super().__init__(x, y, self.__class__.count)
        self._size = size
        self._mass = mass
        self._bound = False
        self._target = None

    @classmethod
    def resetCount(cls):
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
    def target(self, objTarget):
        if isinstance(objTarget, psc.ObjBase):
            self._target = objTarget
