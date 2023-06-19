"""integrin module

This module contains the integrin class which is main component of
ligand-integrin base stem cell simulation.

"""

# bulit-in import
from __future__ import annotations
from typing import Union

# third party import
import numpy as np

# local import
import physica as psc
import nanopattern as npt
import ligand as lig
import cell as cel


class Integrin(psc.ObjBase):
    """Integrin description:

    Integrin is protein which attached to a cell, especially in cell
    membrane.
    """

    # class property
    count = 0

    def __init__(self, cell:cel.Cell, x_position: float, y_position: float):
        """Initial function for Integrin class

        Parameters
        ----------
        cell: Cell
            Cell object which the integrin is located
        x_position: float
            x position of the Integrin
        y_position: float
            y position of the Integrin
        """
        self.__class__.count += 1
        super().__init__(x_position, y_position, self.__class__.count)
        self._cell = cell
        self._size = cell.integrin_size
        self._mass = cell._integrin_mass
        self._kinetic_energy = 0.0
        self._potential_energy = 0.0
        self._bonding_energy = 0.0
        self._neighbors: list[Integrin] = []
        self._target = None
        self._bound = False
        self._nearest = []
        self._radar_radius = 0.0

    def update_target_bound(self, cells:cel.Cells, substrate: npt.Nanopattern):
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
            target_dist = self.get_distance(self.target)

        # if the integrin is surface integrin
        if self.issurface and cells.many:
            surface_integrins: list[Integrin] = cells.surface_integrins_target(
                self._cell, 2 * self.size
            )
            for integrin_ in surface_integrins:
                dist = self.get_distance(integrin_)
                if dist <= 1.2*(self.size + (2**(1/6))*integrin_.size):
                    if (self.target is None) or (dist < target_dist):
                        self.target = integrin_
                        target_dist = self.get_distance(self.target)

        # if the integrin is not surface integrin
        else:
            nearest_ligand: list[lig.Ligand] = substrate.nearest(
                self.x_position, self.y_position, 2 * substrate.ligand_size
            )
            for ligand_ in nearest_ligand:
                if ligand_.bound is False:
                    dist = self.get_distance(ligand_)
                    if dist <= (self.size + (2**(1/6))*ligand_.size):
                        if (self.target is None) or (dist < target_dist):
                            self._target = ligand_
                            target_dist = self.get_distance(self.target)

    def bonding(self):
        """Procedure of integrin bonding on ligand or other integrin"""
        if (self.bound is False) and (self.target is not None):
            if self.target.bound is False:
                # for self
                dist = self.get_distance(self.target)
                # change the position into lowest potential well
                if dist == 0:
                    self.x_position = self.x_position + (self.size + (2**(1/6))*self.target.size)
                else:
                    vector_dir = (self.position - self.target.position)/dist
                    self.position = self.target.position + vector_dir*(self.size + (2**(1/6))*self.target.size)
                self.bound = True
                self._bonding_energy = self.kinetic_energy
                self._velocity = np.array([0.0, 0.0])
                self._temp_velocity = self._velocity
                self._acceleration = np.array([0.0, 0.0])
                self._temp_acceleration = self._acceleration
                self._force = np.array([0.0, 0.0])
                self._temp_force = self._force
                self._temp_position = self._position

                # for target
                self.target.bound = True
                if isinstance(self.target, Integrin):
                    self.target._bonding_energy = self.target.kinetic_energy
                self.target._velocity = np.array([0.0, 0.0])
                self.target._temp_velocity = self.target._velocity
                self.target._acceleration = np.array([0.0, 0.0])
                self.target._temp_acceleration = self.target._acceleration
                self.target._force = np.array([0.0, 0.0])
                self.target._temp_force = self.target._force
                self.target._temp_position = self.target._position
                self.target.target = self
                return True
            self.target = None
        return False
                
    
    def calc_potential(self, 
                       nearest_obj, 
                       normal_length=1, 
                       spring_constant=1, 
                       epsilon=1, 
                       sigma=1 ):
        """Procedure to calculate the potential 
        
        Parameter
        --------
        nearest_obj: objBase
            The nearest objects which contribute to the 
            lenneard-jones potential.
        epsilon: float
            The epsilon value for Lennard-jones.
        sigma: float
            The sigma value for Lennard-jones.
        """
        energy = 0.0
        # Calculate the Lennard-Jones Potential
        for obj in nearest_obj:
            if isinstance(obj, lig.Ligand):
                lj_pot_energy = psc.potential.lennardjones_6_12(obj.position, self.position, epsilon, sigma)
            elif isinstance(obj, Integrin):
                lj_pot_energy = 0.5*psc.potential.lennardjones_6_12(obj.position, self.position, epsilon, sigma)
            energy = energy + lj_pot_energy
        if self.bound:
            if isinstance(self.target, Integrin):
                lj_pot_energy = 0.5*psc.potential.lennardjones_6_12(self.target.position, self.position, epsilon, sigma)
            if isinstance(self.target, lig.Ligand):
                lj_pot_energy = psc.potential.lennardjones_6_12(self.target.position, self.position, epsilon, sigma)
            energy = energy + lj_pot_energy
        # Calculate the spring potential
        for neighbor in self.neighbors:
            spring_pot_energy = 0.5*psc.potential.spring(neighbor.position, self.position, spring_constant, normal_length)
            energy = energy + spring_pot_energy
        self._potential_energy = energy
        return energy


    @property
    def neighbors(self):
        """return the integrin neighbors"""
        return self._neighbors
    
    @property
    def neighbors_id(self):
        """return the ids of the neighbors"""
        list_id = []
        for neighbors in self.neighbors:
            list_id.append(neighbors.id_)
        return list_id

    @property
    def x_distance_center(self):
        """return the x distance from the center of the host cell"""
        return self._cell.x_position - self.x_position

    @property
    def y_distance_center(self):
        """return the y distance from the center of the host cell"""
        return self._cell.y_position - self.y_position

    @property
    def issurface(self):
        """determine if the integrin is surface integrin based on
        neighbors number.

        Return
        ------
        True or false
        """
        return bool(len(self.neighbors) < 6)

    @property
    def target(self) -> Union[None, lig.Ligand, Integrin]:
        """return the target object of integrin"""
        return self._target
    
    @target.setter
    def target(self, value):
        if isinstance(value, lig.Ligand) or isinstance(value, Integrin):
            self._target = value

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

    @property
    def kinetic_energy(self):
        """return integrin kinetic energy"""
        return (0.5)*self.mass*np.dot(self.velocity, self.velocity)

    @property
    def bonding_energy(self):
        """return the bonding energy"""
        return self._bonding_energy
    
    @classmethod
    def reset_count(cls):
        """reset the number of integrin created into 0"""
        cls.count = 0
