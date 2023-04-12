"""integrin module

This module contains the integrin class which is main component of
ligand-integrin base stem cell simulation.
"""

import numpy as np
import src.physica as psc
from .nanopattern import Nanopattern
from .cell import Cell, Cells
from .ligand import Ligand
from integrin import Integrin
from typing import Union

class Integrin(psc.ObjBase):
    """Integrin description:

    Integrin is protein which attached to a cell, especially in cell 
    membrane.
    """
    
    # class property
    count = 0
    
    def __init__(self, cell: Cell, x: float, y: float):
        """Initial function for Integrin class

        Parameters
        ----------
        cell: :obj: Cell
            Cell object which the integrin is located
        x: float
            x position of the Integrin
        y: float
            y position of the Integrin
        """
        self.__class__.count += 1
        super().__init__(x, y, self.__class__.count)
        self._cell = cell
        self._size = cell.integrin_size
        self._mass = cell._integrin_mass
        self._neighbors: list[Integrin] = []
        self._target = None
        self._bound = False
    
    def updateTargetBound(self, cells: Cells, substrate: Nanopattern):
        """procedure to update the target bound
        
        The updating process is based on the integrin position and the
        position of nearest unoccupied ligand/integrin.

        Parameters
        ----------
        cells: :obj: Cells
            the compilation of cells
        substrate: :obj: Nanopattern
            the nanopatterned substrate, consists of ligands
        """
        target_dist = -1
        # check whether it has previous target or not check the target_dist
        if self.target is not None:
            target_dist = self.getDistance(self.target)
        
        # if the integrin is surface integrin
        if self.isSurface:
            surface_integrins: list[Integrin] = cells.surfaceIntegrinsTarget(self._cell, 2*self.size)
            for integrin in surface_integrins:
                dist = self.getDistance(integrin)
                if dist <= (self.size + integrin.size):
                    if (self.target is None) or (dist < target_dist):
                        self.target = integrin
                        target_dist = self.getDistance(self.target)

        # if the integrin is not surface integrin
        else:   
            nearest_ligand: list[Ligand] = substrate.nearest(self.x, self.y, 2*substrate.ligand_size)
            for ligand in nearest_ligand:
                if ligand.bound is False:
                    dist = self.getDistance(ligand)
                    if dist <= (self.size + ligand.size):
                        if (self.target is None) or (dist < target_dist):
                            self.target = ligand
                            target_dist = self.getDistance(self.target)

    def bonding(self):
        """Procedure of integrin bonding on ligand or other integrin"""
        if (self.bound is False) and (self.target is not None) :
            if self.target.bound is False:
                self.bound = True
                self.target.bound = True
                self.target.target = self
            else:
                self.target = None


    @property
    def neighbors(self):
        """return the integrin neighbors"""
        return self._neighbors
    
    @property
    def x_distance_center(self):
        """return the x distance from the center of the host cell"""
        return self._cell.x - self.x
    
    @property
    def y_distance_center(self):
        """return the y distance from the center of the host cell"""
        return self._cell.y - self.y   

    @property
    def isSurface(self):
        """determine if the integrin is surface integrin based on
        neighbors number.

        Return
        ------
        True or false
        """
        if len(self.neighbors) < 6:
            return True
        else:
            return False  
    
    @property
    def target(self) -> Union[None, Ligand, Integrin]:
        """return the target object of integrin"""
        return self._target

    @property
    def bound(self):
        """Integrin condition related to bonding status with other 
        object.
        """
        return self._bound

    @bound.setter
    def bound(self, value):
        if isinstance(value, bool):
            self._bound = value
    
    @classmethod
    def resetCount(cls):
        """reset the number of integrin created into 0"""
        cls.count = 0