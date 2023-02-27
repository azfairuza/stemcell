from __future__ import annotations
import numpy as np
from descartes import PolygonPatch
import alphashape
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from datetime import datetime
from pathlib import Path
import shutil
import os
import imageio
from PIL import Image
import warnings
import sys
from random import choice
from cv2 import GaussianBlur

class Base:
    """ Object base class.
    
    it can be used for any kind of object base. For any other projects.
  
    """
    def __init__(self, x: float, y: float, id: int) -> None:
        """ initial function for Base class
        
        Parameters
        ----------
        x: float
            x position of the object 
        y: float
            y position of the object
        id: int
            the identity number of the object
        """
        
        self._position = np.array((x, y))                       # Position of the object (x, y)
        self._speed = np.array((0, 0), dtype=float)             # Speed of the object (vx, vy)
        self._acceleration = np.array((0, 0), dtype=float)      # Acceleration of the object (ax, ay)
        self._force = np.array((0, 0), dtype=float)             # Force of the object (fx, fy)
        self._id = id                                           # Object identity number
        self._mass = 1                                          # Object mass (default is 1)
        self._bound = False                                     # Object bonding status with other object
        self._size = 1                                          # Object size (default is 1)
        self._target = []                                       # Object's target (list of another object)
        
        # temporary 
        self._temp_position = np.zeros(2)                       # Temporary position (0,0)
        self._temp_speed = np.zeros(2)                          # Temporary speed (0,0)
        self._temp_acceleration = np.zeros(2)                   # Temporary acceleration (0,0)
        self._temp_force = np.zeros(2)                          # Temporary force (0,0)
 
    # region <Method>
    def update(self):
        """Update the object state.

        This function is invoked in the end calculation to update the 
        state of the object. Manual update can be used as well instead
        of using this function.
        """
        if self.bound == False:
            self.acceleration = self.temp_acceleration
            self.speed = self.temp_speed
            self.position = self.temp_position
    
    def getDistance(self, obj, bias = None):
        """Get the distance from two objects.

        Parameters
        ----------
        obj : :obj:`base`
            the object target point.
        bias : :obj:`ndarray', default : None
            It should contain 2 dimension ndarray which corresponds to 
            (x,y) bias.
        """
        if isinstance(obj, Base) or isinstance(obj, Ligand) or isinstance(obj, Integrin):
            if isinstance(bias, np.ndarray):
                return np.linalg.norm(obj.position - (self.position + bias))
            else:
                return np.linalg.norm(obj.position - self.position )
        else:
            NotImplemented
    
    def getXDistance(self, obj, bias = 0.0):
        """ Get the X distance from two objects.

        Parameters
        ----------
        obj : :obj:`base`
            the object target point.
        bias : float, default: 0.0
            bias value to alter the x position of first object.
        """
        if isinstance(obj, Base):
            return obj.x - (self.x + bias)
        else:
            NotImplemented
    
    def getYDistance(self, obj, bias = 0.0):
        """Get the Y distance from two objects
        
        Parameters
        ----------
        obj : :obj:`base`
            the object target point.
        bias : float, default: 0.0
            bias value to alter the y position of first object.
        """
        if isinstance(obj, Base):
            return obj.y - (self.y + bias)
        else:
            NotImplemented            

    def getAngle(self, obj, bias = None):
        """Get the angle of two object in respect with X-axis
        
        Parameters
        ----------
        obj : :obj:`base`
            the object target point.
        bias : :obj: `ndarray`, default: none
            bias value to alter theposition of first object.

        Return
        ------
        The output is in radians
        """
        if isinstance(obj, Base):
            if isinstance(bias, np.ndarray):
                return np.arctan2(self.getYDistance(obj, bias[1]), self.getXDistance(obj, bias[0]))
            else:
                return np.arctan2(self.getYDistance(obj), self.getXDistance(obj))
        else:
            NotImplemented
    # endregion
    
    # region <Property>    
    @property
    def id(self):
        """The identity number of the object."""
        return self._id
    
    @property
    def number_target(self):
        """int: number of target objects."""
        return len(self._target)
   
    @property
    def distance_target(self):
        """float: the distance of the object from the target.
        
        Note
        ----
        This function can only be used if the object has single target.
        """
        if self.number_target == 1:
            return np.linalg.norm(self._target[0].position - self._position)
        else:
            NotImplemented
    
    @property
    def x_distance_target(self):
        """float: the x distance of the object from the target.
        
        Note
        ----
        This function can only be used if the object has single target.
        """
        if self.number_target == 1:
            return self._target[0].x - self.x
        else:
            NotImplemented
    
    @property
    def y_distance_target(self):
        """float: the y distance of the object from the target.
        
        Note
        ----
        This function can only be used if the object has single target.
        """
        if self.number_target == 1:
            return self._target[0].y - self.y
        else:
            NotImplemented
    
    @property
    def angle_target(self):
         """float: the angle of the object and the target with respect to x axis.
        
        Note
        ----
        This function can only be used if the object has single target.
        
        Return
        ------
        Return in radians value
        """
        if self.number_target == 1:
            return np.arctan2(self.y_distance_target, self.x_distance_target)
        else:
            NotImplemented

    @property
    def target(self):
        """return the list of targets"""
        return self._target
    
    @target.setter
    def target(self, value):
        if isinstance(value, Ligand) or isinstance(value, Integrin):
            self._target.append(value)
        else:
            NotImplemented

    @property
    def bound(self):
        """return the bonding status of the object."""
        return self._bound
    
    @bound.setter
    def bound(self, value: bool):
        self._bound = value

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
        if isinstance(coord, Base):
            self._position = coord.position
        elif isinstance(coord, tuple) or isinstance(coord, list) or isinstance(coord, np.ndarray):
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
        if isinstance(value, Base):
            self._speed = value.speed
        elif isinstance(value, tuple) or isinstance(value, list) or isinstance(value, np.ndarray):
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
        if isinstance(value, Base):
            self._acceleration = value.acceleration
        elif isinstance(value, tuple) or isinstance(value, list) or isinstance(value, np.ndarray):
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
        if isinstance(value, Base):
            self._force = value.force
        elif isinstance(value, tuple) or isinstance(value, list) or isinstance(value, np.ndarray):
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
        if isinstance(coord, tuple) or isinstance(coord, list) or isinstance(coord, np.ndarray):
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
        if isinstance(value, Base):
            self._temp_speed = value.temp_speed
        elif isinstance(value, tuple) or isinstance(value, list) or isinstance(value, np.ndarray):
            self.temp_vx = value[0]
            self.emp_vy = value[1]
        else:
            return NotImplemented

    @property
    def temp_ax(self):
        """return the temporary x component acceleration of the object.
        """
        return self._temp_acceleration[0]
    
    @temp_ax.setter
    def temp_ax(self, value):
        self._temp_acceleration[0] = value
    
    @property
    def temp_ay(self):
        """return the temporary y component acceleration of the object.
        """
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
        if isinstance(value, Base):
            self._temp_acceleration = value.temp_acceleration
        elif isinstance(value, tuple) or isinstance(value, list) or isinstance(value, np.ndarray):
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
        if isinstance(value, Base):
            self._temp_force = value.temp_force
        elif isinstance(value, tuple) or isinstance(value, list) or isinstance(value, np.ndarray):
            self.temp_fx = value[0]
            self.temp_fy = value[1]
        else:
            return NotImplemented
    
    @property
    def size(self):
        """return the size of the object"""
        return self._size

    # endregion

class Ligand(Base):
    
    # region <Class Variable>
    count = 0
    # endregion

    def __init__(self, x: float, y: float, size: float, mass: float=1) -> None:
        self.__class__.count += 1
        super().__init__(x, y, self.__class__.count)
        self._size = size
        self._mass = mass
    
    def decideBound(self):
        if self.bound == False  and len(self.target) > 0:
            target = choice(self.target)
            self.bound = True
            target.bound = True
        else:
            pass

    @classmethod
    def resetCount(cls):
        cls.count = 0

class Nanopattern:

    def __init__(self):
        '''
        Inital procedure when creating a nanopattern.
        '''
        Ligand.resetCount()
        patcon = readFile('PATCON')

        # get value
        substrate_size = getValue(patcon, 'size')
        scale = getValue(patcon, 'scale')
        gridsize = getValue(patcon, 'gridsize')

        # nanopattern properties
        self._height = substrate_size[0]*scale              # nanopattern height
        self._width = substrate_size[1]*scale                    # nanopattern width
        self._x_dist = getValue(patcon, 'xdist')             
        self._y_dist = getValue(patcon, 'ydist')
        self._ligand_size = getValue(patcon, 'ligandsize')*scale
        self._scale = scale
        self._ligands = [] 
        self._gridssize = gridsize   
        self._grid = [[[] for j in range(int(gridsize[1]))] for i in range(int(gridsize[0]))]         
        
        # build nanopattern
        self.build()


    def build(self):
        '''
        building the nanopattern from empty list of ligands
        '''
        x = 0
        y = 0
        x_iter = 0
        y_iter = 0
        while(y <=  self.height):
            y_dist_index = y_iter % len(self.y_dist)
            while( x <= self.width):
                x_dist_index = x_iter % len(self.x_dist)
                x_index = int(np.floor(x/self.x_gap))
                y_index = int(np.floor(y/self.y_gap))
                obj = Ligand(x, y, self._ligand_size)
                self._ligands.append(obj)
                try:
                    self._grid[y_index][x_index].append(obj)
                except:
                    print(f'in build nanopattern, xindex: {x_index}, yindex: {y_index}')
                    raise IndexError
                x += self.x_dist[x_dist_index]
                x_iter += 1
            y += self.y_dist[y_dist_index]
            y_iter += 1
            x_iter = 0
            x = 0
        print('SYSTEM: nanopattern has been created')

        # for row in range(self.row_number):
        #     for col in range(self.col_number):
        #         for pos in self.seed:
        #             x = (pos[0] + col)*self._unit_width
        #             y = (pos[1] + row)*self._unit_height
        #             x_index = int(np.floor(x/self.x_gap))
        #             y_index = int(np.floor(y/self.y_gap))
        #             obj = Ligand(x, y, self._ligand_size)
        #             self._ligands.append(obj)
        #             try:
        #                 self._grid[y_index][x_index].append(obj)
        #             except:
        #                 print(f'in build nanopattern, xindex: {x_index}, yindex: {y_index}')
        #                 raise IndexError
        # print('SYSTEM: nanopattern has been created')

    def getLigandById(self, id: int) -> Ligand:
        '''
        Procedure to get a ligand by id.
        '''
        for ligand in self.ligands:
            if ligand.id == id:
                return ligand
    
    def savePATMAP(self, time: datetime):
        namefolder = f'./output/{getTime(time)}-output/file'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/PATMAP.txt'
        head_text = 'ligand_id\tx_pos\ty_pos\n'
        with open(namefile, 'w') as output:
            output.write(head_text)
        for ligand in self.ligands:
            content = f'{ligand.id}\t{ligand.x/self._scale}\t{ligand.y/self._scale}\n'
            with open(namefile, 'a') as output:
                output.write(content)

    def nearest(self, x_index, y_index) -> list[Ligand]:
        ligandList = []
        regions = [True for i in range(9)]

        if x_index-1 < 0 or x_index-1 > (self._gridssize[0]-1):
            regions[0] = False
            regions[3] = False
            regions[6] = False
        if y_index-1 < 0 or y_index-1 > (self._gridssize[1]-1):
            regions[0] = False
            regions[1] = False
            regions[2] = False
        if x_index+1 >= (self._gridssize[0]-1) or x_index+1 < 0:
            regions[2] = False
            regions[5] = False
            regions[8] = False
        if y_index+1 >= (self._gridssize[1]-1) or y_index+1 < 0:
            regions[6] = False
            regions[7] = False
            regions[8] = False
        if x_index < 0 or x_index > (self._gridssize[0]-1) or y_index < 0 or y_index > (self._gridssize[1]-1):
            regions[4] = False
        try:
            if regions[0] == True:
                for ligand in self._grid[y_index-1][x_index-1]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[1] == True:
                for ligand in self._grid[y_index-1][x_index]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[2] == True:
                for ligand in self._grid[y_index-1][x_index+1]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[3] == True:
                for ligand in self._grid[y_index][x_index-1]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[4] == True:
                for ligand in self._grid[y_index][x_index]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[5] == True:
                for ligand in self._grid[y_index][x_index+1]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[6] == True:
                for ligand in self._grid[y_index+1][x_index-1]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[7] == True:
                for ligand in self._grid[y_index+1][x_index]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            if regions[8] == True:
                for ligand in self._grid[y_index+1][x_index+1]:
                    if ligand.bound == False:
                        ligandList.append(ligand)
            # for ligand in ligandList:
            #     if ligand.bound == True:
            #         ligandList.remove(ligand)
            return ligandList
        except:
            print(f'ERROR: {x_index}, {y_index} is exceding')
            return ligandList
    # region <Nanopattern property>
    
    @property
    def x_gap(self):
        return self._width/self._gridssize[0]

    @property
    def y_gap(self):
        return self._height/self._gridssize[1]

    @property
    def height(self):
        return self._height
    
    @property
    def height_scaled(self):
        return self._height/self._scale
    
    @property
    def width(self):
        return self._width
    
    @property
    def width_scaled(self):
        return self._width/self._scale

    @property
    def ligand_size(self):
        return self._ligand_size
    
    @property
    def ligand_size_scaled(self):
        return self._ligand_size/self._scale
    
    @property
    def pos_list(self):
        position = []
        for ligand in self.ligands:
            position.append((ligand.x, ligand.y))
        return position
    
    @property
    def pos_list_scaled(self):
        position = []
        for ligand in self.ligands:
            position.append((ligand.x/self._scale, ligand.y/self._scale))
        return position
    
    @property
    def x_pos_list(self):
        position = [[],[]]
        for ligand in self.ligands:
            if ligand.bound == False:
                position[0].append(ligand.x)
            else:
                position[1].append(ligand.x)
        return position
    
    @property
    def x_pos_list_scaled(self):
        position = []
        for ligand in self.ligands:
            position.append(ligand.x/self._scale)
        return position
    
    @property
    def y_pos_list(self):
        position = [[],[]]
        for ligand in self.ligands:
            if ligand.bound == False:
                position[0].append(ligand.y)
            else:
                position[1].append(ligand.y)
        return position
    
    @property
    def y_pos_list_scaled(self):
        position = []
        for ligand in self.ligands:
            position.append(ligand.y/self._scale)
        return position

    @property
    def bound_number(self):
        count = 0
        for ligand in self.ligands:
            if ligand.bound == True: 
                count += 1
        return count
    
    @property
    def dot_number(self):
        return Ligand.count
    
    @property
    def ligands(self) -> list[Ligand]:
        return self._ligands
    
    @property
    def x_dist(self):
        return self._x_dist

    @property
    def y_dist(self):
        return self._y_dist
    # endregion

class Integrin(Base):
    
    # class property
    count = 0
    
    def __init__(self, cell: Cell, x: float, y: float) -> None:
        self.__class__.count += 1
        super().__init__(x, y, self.__class__.count)
        self._cell = cell
        self._size = cell.integrin_size
        self._mass = cell._integrin_mass
        self._neighbors = []

    def rungeKutta(self, cells: Cells, substrate: Nanopattern, 
        springconstant, damper, epsilon, viscosity, dt):
        k1v = self.resultantAccel(cells, substrate, 
            springconstant, damper, epsilon, viscosity)*dt
        k1x = self.speed*dt
        k2v = self.resultantAccel(cells, substrate, 
            springconstant, damper, epsilon, viscosity, 
            (k1x/2), (k1v/2))*dt
        k2x = (self.speed + (k1v/2))*dt
        k3v = self.resultantAccel(cells, substrate, 
            springconstant, damper, epsilon, viscosity, 
            (k2x/2), (k2v/2))*dt
        k3x = (self.speed + (k2v/2))*dt
        k4v = self.resultantAccel(cells, substrate, 
            springconstant, damper, epsilon, viscosity, 
            k3x, k3v)*dt
        k4x = (self.speed + k3v)*dt
        self.temp_speed = self.speed + 1/6*(k1v + (2*k2v) + (2*k3v) + k4v)
        self.temp_position = self.position + 1/6*(k1x + (2*k2x) + (2*k3x) + k4x)  

    def euler(self, cells: Cells, substrate: Nanopattern, 
        springconstant, damper, epsilon, viscosity, dt):
        a = self.resultantAccel(cells, substrate,
        springconstant, damper, epsilon, viscosity)
        self.temp_speed =  self.speed + a*dt
        self.temp_position = self.position + self.speed*dt
       
    def resultantAccel(self, cells: Cells, substrate: Nanopattern, 
        springconstant, damper, epsilon, viscosity,
        xbias = np.array((0,0)), vbias = np.array((0,0))):
        total_force_x = 0
        total_force_y = 0
        if cells.number_cell > 1:
            if self.isSurface:
                #surface integrin interaction
                surface_integrin = cells.surfaceIntegrinsTarget(self._cell)
                surface_force = self.totalForceIntegrin(surface_integrin, epsilon, xbias)
                total_force_x += surface_force[0]
                total_force_y += surface_force[1]
            else:
                x_index = int(np.floor(self.x/substrate.x_gap))
                y_index = int(np.floor(self.y/substrate.y_gap))
                nearest_ligand = substrate.nearest(x_index, y_index)
                nearest_force = self.totalForceLigand(nearest_ligand, epsilon, xbias)
                total_force_x += nearest_force[0]
                total_force_y += nearest_force[1]
        else:
            if self.isSurface:
                pass
            else:
                x_index = int(np.floor(self.x/substrate.x_gap))
                y_index = int(np.floor(self.y/substrate.y_gap))
                nearest_ligand = substrate.nearest(x_index, y_index)
                nearest_force = self.totalForceLigand(nearest_ligand, epsilon, xbias)
                total_force_x += nearest_force[0]
                total_force_y += nearest_force[1]
        neighbors_force = self.totalForceSpring(springconstant, damper, xbias, vbias)
        total_force_x += neighbors_force[0]
        total_force_y += neighbors_force[1]
        drag_force = self.drag_force(viscosity, vbias)
        total_force_x += drag_force[0]
        total_force_y += drag_force[1]
        return np.array((total_force_x/self._mass, total_force_y/self._mass))
    
    def forceLJ(self, obj, epsilon, xbias = np.array((0,0))):
        dist = self.getDistance(obj, xbias)
        angle = self.getAngle(obj, xbias)
        # print(f'self ({self.x}, {self.y}) | obj ({obj.x}, {obj.y}) | angle {angle}')
        alpha = (1.78*self.size/dist)**6
        forceValue = (48/dist)*epsilon*alpha*(alpha-0.5)
        xForce = forceValue*np.cos(angle)
        yForce = forceValue*np.sin(angle)
        return (xForce, yForce)
    
    def totalForceIntegrin(self, surface_integrin, epsilon, xbias = np.array((0,0))):
        totalXForce = 0
        totalYForce = 0
        for integrin in surface_integrin:
            force = self.forceLJ(integrin, epsilon, xbias)
            totalXForce += force[0]
            totalYForce += force[1]
        return(totalXForce, totalYForce)
    
    def totalForceLigand(self, nearest_ligand, epsilon, xbias):
        totalXForce = 0
        totalYForce = 0
        for ligand in nearest_ligand:
            force = self.forceLJ(ligand, epsilon, xbias)
            totalXForce += force[0]
            totalYForce += force[1]
        return(totalXForce, totalYForce)
    
    def drag_force(self, viscosity, vbias = np.array((0,0))):
        drag_force_x = -6*np.pi*viscosity*self.size*(self.vx + vbias[0])
        drag_force_y = -6*np.pi*viscosity*self.size*(self.vy + vbias[1])
        return (drag_force_x, drag_force_y)


    def totalForceSpring(self, springconstant, damper, xbias = np.array((0,0)), vbias = np.array((0,0))):
        totalXForce = 0
        totalYForce = 0
        for obj in self.neighbors:
            dist = self.getDistance(obj, xbias)
            deltadist = dist - self._cell._min_dst
            angle = self.getAngle(obj, xbias)
            alpha = (1.78*self.size/dist)**12
            # print(f'self ({self.x}, {self.y}) | obj ({obj.x}, {obj.y}) | angle {angle}')
            if deltadist > 0:
                forceValue = (springconstant*deltadist)
            else:
                forceValue = (5*springconstant*deltadist)
            xForce = forceValue*np.cos(angle) - (damper*(self.vx+vbias[0]))
            yForce = forceValue*np.sin(angle) - (damper*(self.vy+vbias[1]))
            totalXForce += xForce
            totalYForce += yForce
        return (totalXForce, totalYForce)
    
    def move(self, timestep):
        self.temp_ax = self.fx/self._mass
        self.temp_ay = self.fy/self._mass
        self.temp_vx = self.vx + self.ax*timestep
        self.temp_vy = self.vy + self.ay*timestep
        self.temp_x = self.x + self.vx*timestep
        self.temp_y = self.y + self.vy*timestep
    
    def updateTargetBound(self, cells: Cells, substrate: Nanopattern):
        if self.isSurface:
            surface_integrins = cells.surfaceIntegrinsTarget(self._cell)
            for integrin in surface_integrins:
                if integrin.bound == False:
                    dist = self.getDistance(integrin)
                    if dist <= (self.size + integrin.size):
                        integrin.target = self
        else:
            x_index = int(np.floor(self.x/substrate.x_gap))
            y_index = int(np.floor(self.y/substrate.y_gap))
            # print(f'integrin {self.id}: Force: {self.fx},{self.fy}')
            # print(f'integrin {self.id}: ({self.x},{self.y}) -> index ({x_index},{y_index}) with gap ({substrate.x_gap},{substrate.y_gap})')
            nearest_ligand = substrate.nearest(x_index, y_index)
            for ligand in nearest_ligand:
                if ligand.bound == False:
                    dist = self.getDistance(ligand)
                    if dist <= (self.size + ligand.size):
                        ligand.target = self

    def decideBound(self):
        if self.isSurface and self.bound == False and len(self.target) > 0:
            target = choice(self.target)
            self.bound = True
            target.bound = True
        else:
            pass


    @property
    def neighbors(self):
        return self._neighbors
    
    @property
    def x_distance_center(self):
        return self._cell.x - self.x
    
    @property
    def y_distance_center(self):
        return self._cell.y - self.y
    
    @property
    def angle_center(self):
        return np.arctan2(self.y_distance_target, self.x_distance_target)   

    @property
    def isSurface(self):
        if len(self.neighbors) < 6:
            return True
        else:
            return False  

    @classmethod
    def resetCount(cls):
        cls.count = 0

class Cell:
    count = 0
    def __init__(self, x: float, y: float, size: float, 
    mass: float, radius: float, min_dst: float, scale: float) -> None:
        self.__class__.count += 1
        self._id = self.__class__.count
        self._mass = 0
        self._position = np.array((x*scale, y*scale))
        self._integrin_size = size*scale
        self._integrin_mass = mass
        self._min_dst = min_dst*scale
        self._scale = scale
        self._integrins = self.build(radius*scale, self.x, self.y, self._min_dst)
        self._alpha_shape = None

    def build(self, radius, x, y, a):
        # find the length total of parallelogram
        theta = (1/6)*np.pi
        diagonal = 2*radius/np.cos(theta)
        horizontal = 2*radius*(1+np.tan(theta))
        diag_index = int(np.floor(diagonal/a)) + 1
        hori_index = int(np.floor(horizontal/a)) + 1
        grid = [[(j,i) for j in range(hori_index)] for i in range(diag_index)]
        objs = []
        
        for i in range(diag_index):
            for j in range(hori_index):
                x_dot = (a*j - a*np.sin(theta)*i) + x - radius + a/2
                y_dot = (a*np.cos(theta)*i) + y - radius + a/2
                dist = np.sqrt((x_dot - x)**2 + (y_dot - y)**2)
                if dist > radius:
                    grid[i][j] = None
                else:
                    grid[i][j] = Integrin(self, x_dot, y_dot)
                    objs.append(grid[i][j])
        
        for i in range(diag_index):
            for j in range(hori_index):
                if grid[i][j] != None:
                    obj = grid[i][j]
                    neighbors = [True for i in range(6)]
                    #check index
                    if i-1 < 0:
                        neighbors[0] = False
                        neighbors[5] = False
                    if i+1 >= diag_index:
                        neighbors[2] = False
                        neighbors[3] = False
                    if j-1 < 0:
                        neighbors[0] = False
                        neighbors[1] = False
                    if j+1 >= hori_index:
                        neighbors[3] = False
                        neighbors[4] = False
                    # append the neighbors
                    if neighbors[0] == True:
                        if grid[i-1][j-1] != None:
                            obj.neighbors.append(grid[i-1][j-1])
                    if neighbors[1] == True:
                        if grid[i][j-1] != None:
                            obj.neighbors.append(grid[i][j-1])
                    if neighbors[2] == True:
                        if grid[i+1][j] != None:
                            obj.neighbors.append(grid[i+1][j])
                    if neighbors[3] == True:
                        if grid[i+1][j+1] != None:
                            obj.neighbors.append(grid[i+1][j+1])
                    if neighbors[4] == True:
                        if grid[i][j+1] != None:
                            obj.neighbors.append(grid[i][j+1])
                    if neighbors[5] == True:
                        if grid[i-1][j] != None:
                            obj.neighbors.append(grid[i-1][j]) 
        return objs
    
    def getIntegrinById(self, id):
        for integrin in self.integrins:
            if integrin.id == id:
                return integrin
    
    # def getIntegrinByPosition(self, x, y):
    #     origin = np.array((x,y))
    #     for integrin in self.integrins:
    #         if np.linalg.norm(integrin.position - origin) :
    #             return integrin
    
    def getSurfaceIntegrin(self):
        integrinList = []
        for integrin in self.integrins:
            if integrin.isSurface:
                integrinList.append(integrin)
        return integrinList

    def updatePosition(self):
        self.x = self.alpha_shape.centroid.x*self._scale
        self.y = self.alpha_shape.centroid.y*self._scale
    
    def updateAlphaShape(self, alphaValue=0):
        self.alpha_shape = alphashape.alphashape(self.pos_list, alphaValue)

    @property
    def integrin_size(self):
        return self._integrin_size
    
    @property
    def integrin_size_scaled(self):
        return self._integrin_size/self._scale

    @property
    def position(self):
        return self._position

    @property
    def x(self):
        return self._position[0]
    
    @x.setter
    def x(self, value):
        self._position[0] = value
    
    @property
    def x_scaled(self):
        return self._position[0]/self._scale
    
    @property
    def y(self):
        return self._position[1]
    
    @y.setter
    def y(self, value):
        self._position[1] = value
    
    @property
    def y_scaled(self):
        return self._position[1]/self._scale
    
    @property
    def pos_list(self):
        position = []
        for integrin in self.integrins:
            position.append((integrin.x, integrin.y))
        return position
    
    @property
    def pos_list_scaled(self):
        position =[]
        for integrin in self.integrins:
            position.append((integrin.x/self._scale, integrin.y/self._scale))
        return position
    
    @property
    def x_pos_list(self):
        position = []
        for integrin in self.integrins:
            position.append(integrin.x)
        return position
    
    @property
    def x_pos_list_scaled(self):
        position = []
        for integrin in self.integrins:
            position.append(integrin.x/self._scale)
        return position
    
    @property
    def y_pos_list(self):
        position = []
        for integrin in self.integrins:
            position.append(integrin.y)
        return position
    
    @property
    def y_pos_list_scaled(self):
        position = []
        for integrin in self.integrins:
            position.append(integrin.y/self._scale)
        return position

    @property
    def alpha_shape(self):
        return self._alpha_shape
    
    @alpha_shape.setter
    def alpha_shape(self, value):
        self._alpha_shape = value

    @property
    def area(self):
        if self.alpha_shape != None:
            return self.alpha_shape.area*(self._scale**2)
    
    @property
    def total_bound(self):
        count = 0
        for integrin in self.integrins:
            if integrin.bound == True:
                count += 1
        return count
    
    @property
    def integrins(self) -> list[Integrin]:
        return self._integrins

    @property
    def number_integrin(self):
        return len(self._integrins)
    
    @property
    def id(self):
        return self._id
    
    @property
    def surfaceIntegrin(self) -> list[Integrin]:
        listIntegrin = []
        for integrin in self.integrins:
            if integrin.isSurface:
                listIntegrin.append(integrin)
        return listIntegrin
    
    @property
    def kineticEnergy(self):
        energy = 0
        for integrin in self.integrins:
            if integrin.bound == False:
                kineticEn = 0.5*integrin._mass*(integrin.vx**2 + integrin.vy**2) 
                energy += kineticEn
        return energy

    @classmethod
    def resetCount(cls):
        cls.count = 0

class Cells:
    def __init__(self, cellsList=None) -> None:
        self._members = []
        
        if cellsList == None:
            # open file
            celcon = readFile('CELCON')

            # get value
            integrin_size = getValue(celcon, 'ressize')
            integrin_mass = getValue(celcon, 'resmass')
            cell_properties = getValue(celcon, 'CELL')
            scale = getValue(celcon, 'scale')
            
            
            for obj in cell_properties:
                member = Cell(obj[0], obj[1], integrin_size, integrin_mass, obj[2], obj[3], scale)
                self._members.append(member)
        
        elif isinstance(cellsList, list):
            for cell in cellsList:
                if isinstance(cell, Cell):
                    self._members.append(cell)
        elif isinstance(cellsList, Cell):
            self._members.append(cellsList)
            
        else:
            pass

        print(f'SYSTEM: {len(self.members)} cell(s) is created')

    def getCellbyId(self, id):
        '''
        Procedure to get cell by id.
        '''
        for cell in self.members:
            if cell.id == id:
                return cell
    
    def excludeCellbyId(self, id):
        '''
        Procedure to exclude a cell from list of cells by Id.
        '''
        new_cells = []
        for cell in self.members:
            if cell.id != id:
                new_cells.append(cell)
        return Cells(new_cells)
    
    def surfaceIntegrinsTarget(self, cell) -> list[Integrin]:
        surface_integrin = []
        for obj in self.members:
            if obj != cell:
                surface_integrin += obj.surfaceIntegrin
        for integrin in surface_integrin:
            if integrin.bound == True:
                surface_integrin.remove(integrin)
        return surface_integrin

    def saveCELLEN(self, num_iteration: int, time: datetime):
        namefolder = f'./output/{getTime(time)}-output/file'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/CELLEN.txt'
        if num_iteration <= 0:
            head_text = 't\t'
            for cell in self.members:
                head_cell = f'E{cell.id}\t'
                head_text += head_cell
            head_text += '\n'
            # save the data
            with open(namefile, 'w') as output:
                output.write(head_text)
        output_text = f'{num_iteration}\t'
        for cell in self.members:
            cell_output = f'{cell.kineticEnergy}\t'
            output_text += cell_output
        output_text += '\n'
        # save the data
        with open(namefile, 'a') as output:
            output.write(output_text)

    def saveCELMAP(self, num_iteration: int, time: datetime):
        namefolder = f'./output/{getTime(time)}-output/file/CELMAP'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/CELMAP{num_iteration:06}.txt'
        head_text = 'cell_id\tintegrin_id\tx_pos\ty_pos\n'
        # save the data
        with open(namefile, 'w') as output:
            output.write(head_text)
        
        for cell in self.members:
            for integrin in cell.integrins:
                content = f'{cell.id}\t{integrin.id}\t{integrin.x}\t{integrin.y}\n'
                with open(namefile, 'a') as output:
                    output.write(content)

    def saveCenterOfMass(self, num_iteration: int, time: datetime):
    
        namefolder = f'./output/{getTime(time)}-output/file'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/CELLCM.txt'
        if num_iteration <= 0:
            head_text = 't\t'
            for cell in self.members:
                head_cell = f'x{cell.id}\ty{cell.id}\tn{cell.id}\t'
                head_text += head_cell
            head_text += '\n'
            # save the data
            with open(namefile, 'w') as output:
                output.write(head_text)
        output_text = f'{num_iteration}\t'
        for cell in self.members:
            cell_output = f'{cell.x}\t{cell.y}\t{cell.total_bound}\t'
            output_text += cell_output
        output_text += '\n'
        # save the data
        with open(namefile, 'a') as output:
            output.write(output_text)
        print(f'SYSTEM: CELLCM updated on {namefolder}')
    
    def saveArea(self, num_iteration: int, time: datetime):
        namefolder = f'./output/{getTime(time)}-output/file'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/CELLAR.txt'
        if num_iteration <= 0:
            head_text = 't\t'
            for cell in self.members:
                head_cell = f'A{cell.id}\t'
                head_text += head_cell
            head_text += '\n'
            # save the data
            with open(namefile, 'w') as output:
                output.write(head_text)
        output_text = f'{num_iteration}\t'
        for cell in self.members:
            area = round(cell.alpha_shape.area, 2)
            cell_output = f'{area}\t'
            output_text += cell_output
        output_text += '\n'
        # save the data
        with open(namefile, 'a') as output:
            output.write(output_text)
        print(f'SYSTEM: CELLAR updated on {namefolder}')


    @property
    def members(self) -> list[Cell]:
        return self._members
    
    @property
    def number_cell(self):
        return len(self._members)
    
    @property
    def x_list(self):
        pos_list = []
        for cell in self.members:
            pos_list += cell.x_pos_list
        return pos_list
    
    @property
    def y_list(self):
        pos_list = []
        for cell in self.members:
            pos_list += cell.y_pos_list
        return pos_list
        
# region <general method>
def simulate(cond=None):
    warnings.filterwarnings("ignore", category=FutureWarning)
    # time = datetime.now()
    start_time =datetime.now()
    if cond == 'debug' or cond == -1:
        time = 'debug'
    elif cond == None:
        time = start_time
    else:
        return print('wrong condition')
    namefolder = f'./output/{getTime(time)}-output/file'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    namefile = f'{namefolder}/SIMLOG.txt'
    log = open(namefile, 'w')
    sys.stdout = log

    print('==================================================================')
    print(f'Program made by\t: Achmad Zacky Fairuza')
    print(f'email\t\t\t: fairuza.zacky1@gmail.com')
    print(f'this program is still under development')
    print('==================================================================')

    #current running time
    print(f'SYSTEM: Simulation is start \t: {time}')

    # read SIMCON.txt

    # open file
    simcon = readFile('SIMCON')
    metadata = getValue(simcon, 'METADATA')
    print(f'SYSTEM: simulation run by \t\t: {metadata["username"]}')
    print(f'SYSTEM: title of simulation\t\t: {metadata["title"]}')

    # get value
    n_iteration = int(getValue(simcon, 'iteration'))
    savefig = getValue(simcon, 'savefig')
    centerofmass = getValue(simcon, 'centerofmass')
    cellarea = getValue(simcon, 'cellarea')
    alphaValue = getValue(simcon, 'alpha')
    cellmaping = getValue(simcon, 'cellmaping')
    patternmaping = getValue(simcon, 'patternmaping')
    showintegrin = getValue(simcon, 'showintegrin')
    savegap = getValue(simcon, 'savegap')
    gif = getValue(simcon, 'gif')
    forcearrow = getValue(simcon, 'forcearrow')
    getcontour = getValue(simcon, 'contour')

    #physics
    epsilon = getValue(simcon, 'epsilon')*1000
    springconstant = getValue(simcon, 'springconstant')
    scale = getValue(simcon, 'scale')
    damper = getValue(simcon, 'damper')
    viscosity = getValue(simcon, 'viscosity')
    dt = getValue(simcon, 'timestep')

    Ligand.resetCount()
    Cell.resetCount()
    Integrin.resetCount()
    substrate = Nanopattern()
    cells = Cells()

    if patternmaping == 1:
        substrate.savePATMAP(time)
    if n_iteration == None:
        n_iteration = 1
    if savefig == None or savefig == 0:
        savefig = False
    else:
        savefig = True

    for cell in cells.members:
        cell.updateAlphaShape(alphaValue=alphaValue)
        for integrin in cell.integrins:
            integrin.updateTargetBound(cells, substrate)
    for cell in cells.members:
        for integrin in cell.integrins:
            integrin.decideBound()
    for ligand in substrate.ligands:
        ligand.decideBound()
    
    showAll(cells, substrate, time, dt,
            show_substrate=True,
            save=savefig,
            folder='newsimulate',
            number=0,
            forcearrow=forcearrow,
            showintegrin=showintegrin)
    
    if getcontour == 1:
        contourPlot(cells, substrate, time, dt,
                    number=0, 
                    folder='newsimulate')
    
    cells.saveCELLEN(0, time)
    cells.saveCELMAP(0, time)
    if cellarea == 1:
        cells.saveArea(0, time)
    if centerofmass == 1:
        cells.saveCenterOfMass(0, time)
    iter_simulation = 0
    # region <simulation>
    while(iter_simulation <= n_iteration):
        iter_simulation += 1
        print(f'SYSTEM: iteration number {iter_simulation}')
        for cell in cells.members:
            for integrin in cell.integrins:
                if integrin.bound == False:
                    integrin.rungeKutta(cells, substrate, 
                        springconstant, damper, epsilon, viscosity, dt)
                    # total_force_x = 0
                    # total_force_y = 0
                    # if cells.number_cell > 1:
                    #     if integrin.isSurface:
                    #         #surface integrin interaction
                    #         surface_integrin = cells.surfaceIntegrinsTarget(cell)
                    #         surface_force = integrin.totalForceIntegrin(surface_integrin, epsilon)
                    #         total_force_x += surface_force[0]
                    #         total_force_y += surface_force[1]
                    #     else:
                    #         x_index = int(np.floor(integrin.x/substrate.x_gap))
                    #         y_index = int(np.floor(integrin.y/substrate.y_gap))
                    #         nearest_ligand = substrate.nearest(x_index, y_index)
                    #         nearest_force = integrin.totalForceLigand(nearest_ligand, epsilon)
                    #         total_force_x += nearest_force[0]
                    #         total_force_y += nearest_force[1]
                    # else:
                    #     if integrin.isSurface:
                    #         pass
                    #     else:
                    #         x_index = int(np.floor(integrin.x/substrate.x_gap))
                    #         y_index = int(np.floor(integrin.y/substrate.y_gap))
                    #         nearest_ligand = substrate.nearest(x_index, y_index)
                    #         nearest_force = integrin.totalForceLigand(nearest_ligand, epsilon)
                    #         total_force_x += nearest_force[0]
                    #         total_force_y += nearest_force[1]
                    # neighbors_force = integrin.totalForceSpring(springconstant, damper)
                    # total_force_x += neighbors_force[0]
                    # total_force_y += neighbors_force[1]
                    # drag_force = integrin.drag_force(viscosity)
                    # total_force_x += drag_force[0]
                    # total_force_y += drag_force[1]
                    # integrin.fx = total_force_x
                    # integrin.fy = total_force_y
                    # integrin.move(dt)        

        if iter_simulation%savegap == 0 or iter_simulation > n_iteration:
            for cell in cells.members:
                cell.updateAlphaShape(alphaValue=alphaValue)
            showAll(cells, substrate, time, dt,
                show_substrate=True,
                save=savefig,
                folder='newsimulate',
                number=iter_simulation,
                showintegrin=showintegrin)

            if getcontour == 1:
                contourPlot(cells, substrate, time, dt,
                    number=iter_simulation, 
                    folder='newsimulate')
            
            if cellarea == 1:
                cells.saveArea(iter_simulation, time)
            if centerofmass == 1:
                cells.saveCenterOfMass(iter_simulation, time)
            if cellmaping == 1:
                cells.saveCELMAP(iter_simulation, time)   
        
        for cell in cells.members:
            for integrin in cell.integrins:
                integrin.update()
            try: 
                cell.updatePosition()
            except:
                print('ERROR: cell update')
        for cell in cells.members:
            for integrin in cell.integrins:
                integrin.updateTargetBound(cells, substrate)
        for cell in cells.members:
            for integrin in cell.integrins:
                integrin.decideBound()
        for ligand in substrate.ligands:
            ligand.decideBound()
    # endregion
        
        cells.saveCELLEN(iter_simulation, time)     

    if gif == 1:
        buildGIF(time)
        print(f'SYSTEM: GIF created!')
    
    saveInput(time)
    print(f'SYSTEM: simulation done!')
    elapse_time = datetime.now() - start_time
    print(f'SYSTEM: execution time: {elapse_time}')
    log.close()



def buildLogFolder(time: datetime):
    namefolder = f'./output/{getTime(time)}-ouput/file'
    Path(namefolder).mkdir(parents=True, exist_ok=True)

def starPatternConfig(angle: float, gap: float, size: tuple(float,float), rewrite = False):
    
    if angle == 90 or  angle == -90:
        y_gap = gap*np.sin(np.deg2rad(angle))
        num = abs(int(size[1]/y_gap))
        y = np.linspace(-0.5,0.5,num)
        x = np.zeros(num)
    elif angle == 0:
        x_gap = gap*np.cos(np.deg2rad(angle))
        num = abs(int(size[0]/x_gap))
        x = np.linspace(-0.5,0.5,num)
        y = np.zeros(num)
    elif angle <= 45 or angle >= -45:
        x_gap = gap*np.cos(np.deg2rad(angle))
        num = abs(int(size[0]/x_gap))
        m = np.tan(np.deg2rad(angle))
        x = np.linspace(-0.5,0.5,num)
        y = m*x
    elif angle > 45 or angle < -45:
        y_gap = gap*np.sin(np.deg2rad(angle))
        num = abs(int(size[1]/y_gap))
        m = np.tan(np.deg2rad(angle))
        y = np.linspace(-0.5,0.5,num)
        x = y/m
    
    if x.size != 0:        
        #save to file
        folder_dir = f'./output/patternmaker'
        Path(folder_dir).mkdir(parents=True, exist_ok=True)
        file_dir = f'{folder_dir}/patternmaker.txt'
        for i in range(len(x)-1):
            x_pos = round((x[i] + 0.5),4)
            y_pos = round((y[i] + 0.5),4)
            content = f'{x_pos}\t{y_pos}\n'
            if (i == 0 and rewrite == True):
                with open(file_dir, 'w') as output:
                    output.write(content)      
            else:
                with open(file_dir, 'a') as output:
                    output.write(content)      
    else:
        print('size is smaller than the gap')

def contourPlot(cells: Cells, substrate: Nanopattern, time: datetime, dt, 
    number=0, folder=None):
    '''
    Procedure to show contour plot of integrin in cells
    '''
    # determine the folder's name
    if folder is None:
        namefolder = f'./output/{getTime(time)}-output/figure'
    else:
        namefolder = f'./output/{getTime(time)}-output/figure/{folder}'
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    # determine the file name        
    namefile = f'{namefolder}/C{number:06}.jpg'
    # draw the figure
    fig = plt.figure(figsize=(20,20))
    fig.dpi = 100
    ax = plt.subplot(aspect='equal')
    

    x = cells.x_list
    y = cells.y_list
    xmin = 0
    xmax = substrate.width
    ymin = 0
    ymax = substrate.height

    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    f = np.zeros((100,100))
    x_gap = substrate.width/100
    y_gap = substrate.height/100
    for i in range(len(x)):
        x_index = int(np.floor(x[i]/x_gap))
        y_index = int(np.floor(y[i]/y_gap))
        f[x_index][y_index] += 1
    
    fBlur = GaussianBlur(f, (9,9), 0)
    
    # position = np.vstack([xx.ravel(), yy.ravel()])
    # values = np.vstack([x, y])
    # kernel = st.gaussian_kde(values)
    # f = np.reshape(kernel(position).T, xx.shape)

    fig.gca().set_xlim(xmin, xmax)
    fig.gca().set_ylim(ymin, ymax)
    cfset = fig.gca().contourf(xx, yy, fBlur, cmap='inferno', vmin=0, vmax=1)
    # cset = fig.gca().contour(xx, yy, f, colors='k')
    # fig.gca().clabel(cset, inline=1, fontsize=10)
    cbar = fig.colorbar(cfset, aspect=5, shrink=0.5)
    cbar.set_label('density (integrin/${\mu}m^2$)', rotation=270, fontsize=22, labelpad=40)
    cbar.ax.set_label('density (integrin/${\mu}m^2$)')
    cbar.ax.tick_params(labelsize=22)
    plt.xlabel('X (${\mu}m$)', fontsize=22)
    plt.ylabel('Y (${\mu}m$)', fontsize=22)
    plt.title(f"time = {number*dt}", fontsize=22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    # save if necessary
    fig.savefig(namefile, bbox_inches='tight', dpi=100)
    plt.close()

        



def showAll(cells : Cells, substrate: Nanopattern, time: datetime, dt,
    show_substrate=False, save=False, number=0, folder=None, showintegrin=1, forcearrow=0):
    '''
    Procedure to show all element of simulation including cells and nanopattern.
    '''
    # determine the folder's name
    if folder is None:
        namefolder = f'./output/{getTime(time)}-output/figure'
    else:
        namefolder = f'./output/{getTime(time)}-output/figure/{folder}'
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    # determine the file name        
    namefile = f'{namefolder}/{number:06}.jpg'
    # draw the figure
    fig = plt.figure(figsize=(20,20))
    fig.dpi = 100
    ax = plt.subplot(aspect='equal')
    # print the cells
    for cell in cells.members:
        if show_substrate is True:
            circles(substrate.x_pos_list[0], substrate.y_pos_list[0], substrate.ligand_size, 'green', alpha=1, ec = 'none')
            circles(substrate.x_pos_list[1], substrate.y_pos_list[1], substrate.ligand_size, 'yellow', alpha=1, ec = 'none')
        cm_position = (cell.x, cell.y)
        cm = plt.Circle(cm_position, 2, color='yellow', alpha=0.2)
        plt.gca().add_patch(cm)
        try:
            plt.gca().add_patch(PolygonPatch(cell.alpha_shape, alpha=0.2))
        except:
            pass
        if showintegrin == 1:
            circles(cell.x_pos_list, cell.y_pos_list, cell.integrin_size, 'red', alpha=1, ec='none')
        # draw line between integrins
        for integrin in cell.integrins:
            line = False
            if line == True:
                for neighbor in integrin.neighbors:
                    x_line = [integrin.x, neighbor.x]
                    y_line = [integrin.y, neighbor.y]
                    plt.plot(x_line, y_line, color="black", linewidth="1")
            if forcearrow == 1 and integrin.bound == False:
                plt.arrow(integrin.x, integrin.y, integrin.fx, integrin.fy, color='blue', width=0.3)
                # if integrin.id == 41:
                #     print(f'integrin {integrin.id} not bound ')
                #     x_index = int(np.floor(integrin.x/substrate.x_gap))
                #     y_index = int(np.floor(integrin.y/substrate.y_gap))
                #     nearest_ligand = substrate.nearest(x_index, y_index)
                #     for nearest in nearest_ligand:
                #         x_line = [integrin.x, nearest.x]
                #         y_line = [integrin.y, nearest.y]
                #         plt.plot(x_line, y_line, color="blue", linewidth="1")

    plt.xlim(0, substrate.width)
    plt.ylim(0, substrate.height)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.title(f"time = {number*dt}", fontsize=22)
    # save if necessary
    if save is True:
        print(f'SYSTEM: figure {number:06}.jpg saved on {namefolder}')
        fig.savefig(namefile, bbox_inches='tight', dpi=100)
        plt.close()

def saveInput(time: datetime):
    namefolder = f'./output/{getTime(time)}-output/input'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    shutil.copy2('./PATCON.txt', namefolder)
    shutil.copy2('./CELCON.txt', namefolder)
    shutil.copy2('./SIMCON.txt', namefolder)
    print(f'SYSTEM: input file has been copied on {namefolder}')

def getIntegrinDensity(num: int, a=1, b=1, c=0, width: float=1, gap: float=1, folder_dir=None):
    '''
    Using line with equation: `ax + by + c = 0.`
    
    Which `a` and `b` cannot both zero.

    '''
    distance = lambda a, b, c, x, y: round(abs(a*x + b*y + c)/np.sqrt(a**2 + b**2), 4)
    x_line = lambda a, b, c, x, y: round((b*(b*x - a*y)-a*c)/(a**2 + b**2), 4)
    y_line = lambda a, b, c, x, y: round((a*(-b*x + a*y)-b*c)/(a**2 + b**2), 4)
    
    #read PATCON get size
    patcon = readFile('PATCON')
    size = getValue(patcon, 'subsize')
    
    # check parameter
    if a == 0 and b == 0:
        return print('wrong line parameter [a = 0 and b = 0]')
    # if np.sign(a)*np.sign(b) == 1 and c == 0 :
    #     return print( 'wrong line parameter [a has same sign with b and c = 0]')
    # if np.sign(a)*np.sign(b) == 1 and np.sign(b)*np.sign(c) == 1:
    #     return print('wrong line parameter [a, b, and c have same sign')
    if width <= 0:
        return print('wrong width value')
    
    # read celmap  
    if folder_dir != None:
        file_dir = f'{folder_dir}/CELMAP/CELMAP{num:06}.txt'
    else:
        file_dir = f'./CELMAP/CELMAP{num:06}.txt'
    
    try:
        data_file = open(file_dir, 'r', encoding='utf-8')    
        lst_strng = data_file.readlines()
        data_file.close()
        stripped_strng = filterElement(lst_strng)
    except:
        print('File/Directory not found [file CELMAP must be inside folder CELMAP]')
    
    celmap = []
    for i in range(1,len(stripped_strng)):
        tmp = stripped_strng[i].split()
        tmp[0] = int(tmp[0])
        tmp[1] = int(tmp[1])
        tmp[2] = float(tmp[2])
        tmp[3] = float(tmp[3])
        celmap.append(tmp)
    
    #find start and end point of the line
    if a == 0 and c == 0:
        boundary = [0.0, 0.0, size[0], 0.0]
    elif b == 0 and c == 0:
        boundary = [0.0 , 0.0 , 0.0, size[1]]
    elif a == 0 and c != 0:
        boundary = [0.0, round(-c/b, 4), size[0], round(-c/b, 4)]
    elif b == 0 and c != 0:
        boundary = [round(-c/a, 4), 0.0, round(-c/a, 4), size[1]]
    else:
        boundary = []
        # if x = 0
        y1 = round(-c/b, 4)
        if y1 >= 0 and y1 <= size[1]:
            boundary.append(0.0)
            boundary.append(y1)
        # if y = 0
        x1 = round(-c/a, 4)
        if x1 >= 0 and x1 <= size[0]:
            boundary.append(x1)
            boundary.append(0.0)
        # if x = max
        y1 = round((-a/b)*size[0] + (-c/b), 4)
        if y1 >= 0  and y1 <= size[1]:
            boundary.append(size[0])
            boundary.append(y1)
        # if y = max
        x1 = round((-b/a)*size[1] + (-c/a), 4)
        if x1 >= 0  and x1 <= size[0]:
            boundary.append(x1)
            boundary.append(size[1])
    
    new_boundary = []
    for i in range(0, len(boundary), 2):
        if boundary[i] < 0 or boundary[i+1] < 0:
            pass
        elif boundary[i] > size[0] or boundary[i+1] > size[1]:
            pass
        else:
            new_boundary.append(boundary[i])
            new_boundary.append(boundary[i+1])
    boundary = new_boundary

    new_boundary = []
    valid_point = False
    for i in range(0,len(boundary), 2):
        if i == 0:
            valid_point = True
        else:
            for j in range(0, len(new_boundary), 2):
                if boundary[i] != new_boundary[j] or boundary[i+1] != new_boundary[j+1]:
                    valid_point = True
                else:
                    valid_point = False
        if valid_point == True:
            new_boundary.append(boundary[i])
            new_boundary.append(boundary[i+1])
    boundary = new_boundary
    print(boundary)

    if len(boundary) < 4:
        return print('the line lies outside the area')
    
    startPt = []
    endPt = []
    if boundary[0] == boundary[2]:
        if boundary[1] < boundary[3]:
            startPt.append(boundary[0])
            startPt.append(boundary[1])
            endPt.append(boundary[2])
            endPt.append(boundary[3])
        else:
            endPt.append(boundary[0])
            endPt.append(boundary[1])
            startPt.append(boundary[2])
            startPt.append(boundary[3])
    elif boundary[1] == boundary[3]:
        if boundary[0] < boundary[2]:
            startPt.append(boundary[0])
            startPt.append(boundary[1])
            endPt.append(boundary[2])
            endPt.append(boundary[3])
        else:
            endPt.append(boundary[0])
            endPt.append(boundary[1])
            startPt.append(boundary[2])
            startPt.append(boundary[3])
    elif boundary[0] < boundary[2]:
        startPt.append(boundary[0])
        startPt.append(boundary[1])
        endPt.append(boundary[2])
        endPt.append(boundary[3])
    else:
        endPt.append(boundary[0])
        endPt.append(boundary[1])
        startPt.append(boundary[2])
        startPt.append(boundary[3])
    
    # find the distance between start and end point
    line_length = np.sqrt((endPt[0]-startPt[0])**2 + (endPt[1]-startPt[1])**2)
    total_segment = int(line_length/gap)

    # create container
    integrin_distribution = np.zeros(total_segment)

    # count integrin
    for integrin in celmap:
        tangent_dst = distance(a, b, c, integrin[2], integrin[3])
        if width >= tangent_dst:
            x = x_line(a, b, c, integrin[2], integrin[3])
            y = y_line(a, b, c, integrin[2], integrin[3])
            line_distance = np.sqrt((x-startPt[0])**2 + (y-startPt[1])**2)
            segment_index = int(np.floor(line_distance/gap))
            integrin_distribution[segment_index] += 1
    
    # save file
    if folder_dir != None:
        file_dir = f'{folder_dir}/DOCMAP/DOCMAP{num:06}-{a}-{b}-{c}.txt'
        namefolder = f'{folder_dir}/DOCMAP'
    else:
        file_dir = f'./DOCMAP/DOCMAP{num:06}-{a}-{b}-{c}.txt'
        namefolder = './DOCMAP'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    head_text = 'num\tvalue\n'
    with open(file_dir, 'w') as output:
            output.write(head_text)
    counter = 0
    for segment in integrin_distribution:
        counter += 1
        content = f'{counter}\t{segment}\n'
        with open(file_dir, 'a') as output:
            output.write(content)

def getTime(time: datetime):
    if isinstance(time, datetime):
        year = time.strftime('%Y')
        month = time.strftime('%m')
        day = time.strftime('%d')
        hour = time.strftime('%H')
        minute = time.strftime('%M')
        return f'{year}-{month}-{day}-{hour}{minute}'
    else:
        return 'debug'

def buildGIF(time: datetime):
    image_dir = f'./output/{getTime(time)}-output/figure/newsimulate'
    gif_dir = f'./output/{getTime(time)}-output/figure/gif'
    Path(gif_dir).mkdir(parents=True, exist_ok=True)
    images = []
    for namefile in sorted(os.listdir(image_dir)):
        if namefile.endswith('.jpg'):
            pathfile = os.path.join(image_dir, namefile)
            images.append(Image.fromarray(imageio.imread(pathfile)).resize((1000,1000)))
    pathgif = os.path.join(gif_dir, 'simulation.gif')
    imageio.mimsave(pathgif, images, fps=5)

def readFile(filename: str):
    '''
    Procedure to read input file. There are three input file that can be read. 
    PATCON: contains the configuration for building nanopattern.
    CELCON: contains the configuration for building cells.
    SIMCON: contains the configuration for running the simulation.
    '''
    if filename == 'PATCON':
    # open PATCON file
        try:
            data_file = open('PATCON.txt', 'r', encoding='utf-8')    
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)           # filter out the empty string and '\n'              
        except:
            stripped_strng = 'Error in opening PATCON file'     # error flag if PATCON can not be opened
        finally:
            return stripped_strng

    elif filename == 'CELCON':
    # open CELCON file
        try:
            data_file = open('CELCON.txt', 'r', encoding = 'UTF-8')
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)           # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening CELCON file'     # error flag if CELCON can not be opened
        finally:
            return stripped_strng

    elif filename == 'SIMCON':
    # open SIMCON file
        try:
            data_file = open('SIMCON.txt', 'r', encoding='utf-8')
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)           # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening SIMCON file'     # error flag if SIMCON can not be opened
        finally:
            return stripped_strng

    else:
        print('ERROR: File name is not correct')
        return 0

def filterElement(input_lst: list):
    '''
    procedure to remove empty element in a list and remove the '\n' character.
    '''
    new_lst = []
    for element in input_lst:
        new_element = element.rstrip()          # remove the '\n'    
        if new_element:                         # check if the element is not empty
            new_lst.append(new_element)
    return new_lst

def getValue(lst_strng: list, property_name: str):
    '''
    Procedure to get the value of a property in XXXCON file. There are special
    property name to be used in this function.
    PATCON: property name = 'LIGAND'.
    CELCON: property name = 'CELL'.
    SIMCON: property name = 'METADATA'.
    '''
    # if property_name == 'LIGAND':                                       # get ligand position
    #     start_index = (lst_strng.index('#LIGAND') + 2)                  # start index to search
    #     end_index = lst_strng.index('#END')                             # last index to search
    #     ligand_position = []                                            # container for the ligands
    #     for i in range(start_index, end_index):
    #         dot_position = [float(j) for j in lst_strng[i].split()]     # make the value in float type
    #         ligand_position.append(dot_position)                        # also remove the spaces and tabs
    #     return ligand_position
    if property_name == 'xdist' or property_name == 'ydist':
        start_index = (lst_strng.index('#CONFIG') + 1)
        end_index =  lst_strng.index('#END')                                
        value = None 
        for i in range(start_index, end_index):             # search from starting index to the end index
            data = lst_strng[i].split()                     # split the string by space or tab
            raw_value = data[1:]
            value = [float(i) for i in raw_value]
        return value
    elif property_name == 'CELL':                                       # get cell properties
        start_index = (lst_strng.index('#CELL')+2)                      # start index to search
        end_index = lst_strng.index('#END')                             # last index to search
        cell_properties = []                                            # cell properties container
        for i in range(start_index, end_index):
            cell_property = [float(j) for j in lst_strng[i].split()]    # make the value in float type
            cell_properties.append(cell_property)                       # also remove the spaces and tabs
        return cell_properties
    elif property_name == 'METADATA':
        start_index = (lst_strng.index('#METADATA')+1)                  # start index to search
        end_index = lst_strng.index('#CONFIG')                          # last index to search
        cell_properties = {'username': None, 'title': None }                                            # cell properties container
        for i in range(start_index, end_index):
            stringdata = lst_strng[i].split()                           
            if stringdata[0] == 'username':                             # user of the simulation
                cell_properties['username'] = ' '.join(stringdata[1:])
            elif stringdata[0] == 'title':                              # the title of the simulation
                cell_properties['title'] = ' '.join(stringdata[1:])
        return cell_properties
    else: 
        start_index = (lst_strng.index('#CONFIG') + 1)
        end_index =  lst_strng.index('#END')                                
        value = None 
        for i in range(start_index, end_index):             # search from starting index to the end index
            data = lst_strng[i].split()                     # split the string by space or tab
            if data[0] == property_name:                    # if the first string is equal with argument
                if len(data) > 2:                           # if the string list contain more than 2 members
                    raw_value = data[1:]
                    value = [float(i) for i in raw_value]   # the value is from index 1 to end in float
                elif len(data) == 2:                        # if the string list contain 2 members
                    value = float(data[1])                  # the value is in the index 1
                break
        return value                     

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    '''
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    '''
  
    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection
# endregion




