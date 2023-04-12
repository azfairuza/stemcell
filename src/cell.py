"""cell module

This module contains the cell and cells classes which is the collection of 
integrins and cell, respectively.

#CHECK FINISHED
"""

import alphashape
import numpy as np
import src.physica as psc
from .integrin import Integrin
from .input_reader import readFile, getValue
from shapely.geometry import Polygon, Point
from misc import filter_by_dist

class Cell:
    count = 0
    def __init__(self, x: float, y: float, size: float, 
    mass: float, radius: float, min_dst: float) -> None:
        """ init function for cell class.

        This init function also call Cell._build function to generate
        the integrins inside a cell.

        parameters
        ----------
        x: float
            x-axis position of the cell center of mass
        y: float
            y-axis position of the cell center of mass
        size: float
            the original size of the cell
        mass: float
            the mass of the cell
        radius: float
            radius of the cell
        min_dst: float
            the minimum distance between the integrins. It is used as
            the neutral length of the spring model.

        """
        self.__class__.count += 1
        self._id = self.__class__.count
        self._mass = 0
        self._position = np.array((x, y,))
        self._integrin_size = size
        self._integrin_mass = mass
        self._min_dst = min_dst
        self._integrins: list[Integrin] = self._build(radius, self.x, self.y, self._min_dst)
        self._alpha_shape: Polygon = None

    def _build(self, radius, x, y, a):
        """ a function to build the cell.

        it is called in the init function only. 

        First, the system need to find the smallest parallelogram that
        can surround a circle (as the cell model is using circle). 
        There are 4 point coordinat in this paralellogram. Diameter of
        the circle is equal to the height of the parallelogram. From
        this, the diagonal and the horizontal side length can be found.

        After that, the system can create a grid system by dividing 
        those 2 sides with grid length (a). Then for every grid, there
        is one integrin located in the middle of each grid. The integrin
        will be valid if located inside the circle. 

        Parameters
        ----------
        x: float
            the x position of cell center of mass.
        y: float
            the y position of cell center of mass.
        a: float
            the dist between integrins. The integrins is spread in
            the hexagonal pattern.
        
        """
        # finding the diagonal and horizontal length
        theta = (1/6)*np.pi #angle: 30 degree
        diagonal = 2*radius/np.cos(theta) 
        horizontal = 2*radius*(1+np.tan(theta))

        # create grid base as the minimum lattice
        diag_index = int(np.floor(diagonal/a)) + 1
        hori_index = int(np.floor(horizontal/a)) + 1
        grid: list[list[Integrin]] = [[None for j in range(hori_index)] for i in range(diag_index)]
        objs = []
        
        # check the grid is valid to contain integrin
        # one grid one integrin

        for i in range(diag_index):
            for j in range(hori_index):
                # determine the position of integrin in the grid
                x_dot = (a*j - a*np.sin(theta)*i) + x - radius + a/2
                y_dot = (a*np.cos(theta)*i) + y - radius + a/2
                # get the distance from center of the cell (x, y)
                dist = np.sqrt((x_dot - x)**2 + (y_dot - y)**2)
                if dist > radius:
                    grid[i][j] = None # outside the cell
                else:
                    grid[i][j] = Integrin(self, x_dot, y_dot)
                    objs.append(grid[i][j]) # save the integrin
        
        # find the neighboring 
        for i in range(diag_index):
            for j in range(hori_index):
                if grid[i][j] is not None:
                    obj = grid[i][j]
                    # function to check if the index is valid
                    def checkGrid():
                        if (i-1<0) or (j-1<0) or (i+1>=diag_index) or (j+1>=hori_index):
                            return False
                        else:
                            return True
                    # input the neighboring integrin
                    if checkGrid() and grid[i-1][j-1] is not None:
                        obj.neighbors.append(grid[i-1][j-1])
                    if checkGrid() and grid[i-1][j] is not None:
                        obj.neighbors.append(grid[i-1][j])
                    if checkGrid() and grid[i][j-1] is not None:
                        obj.neighbors.append(grid[i][j-1])
                    if checkGrid() and grid[i][j+1] is not None:
                        obj.neighbors.append(grid[i][j+1])
                    if checkGrid() and grid[i+1][j] is not None:
                        obj.neighbors.append(grid[i+1][j])
                    if checkGrid() and grid[i+1][j+1] is not None:
                        obj.neighbors.append(grid[i+1][j+1])
        return objs
    
    def getIntegrinById(self, id):
        """procedure to get a specific integrin by its ID.

        Parameter
        --------
        id: int
            the id of the integrin.
        
        Return
        ------
        integrin with the id mentioned.
        """
        for integrin in self.integrins:
            if integrin.id == id:
                return integrin
    
    def getSurfaceIntegrin(self):
        """Procedure to get a list of surface integrins in the cell."""
        integrinList = []
        for integrin in self.integrins:
            if integrin.isSurface:
                integrinList.append(integrin)
        return integrinList

    def updatePosition(self):
        """Procedure to update the position of center of mass.
        """
        center_of_mass = np.array([0,0])
        for integrin in self.integrins:
            center_of_mass += integrin.position
        return center_of_mass/self.number_integrin


    def updateAlphaShape(self, alphaValue=0):
        """Alphashape is a method to create the surrounding area of 
        the cell, by updating the alphashape, we get the latest shape
        based on integrin position. 

        Parameters
        ----------
        alphavalue, default 0
            this value determines how tight the surrounding area
        """
        self.alpha_shape: Polygon = alphashape.alphashape(self.pos_list, alphaValue)

    @property
    def integrin_size(self):
        """return the size of the integrins in the cell"""
        return self._integrin_size

    @property
    def position(self):
        """return the position of the cell."""
        return self._position

    @property
    def x(self):
        """return the x-axis position of the cell."""
        return self._position[0]
    
    @x.setter
    def x(self, value):
        self._position[0] = value
    
    @property
    def y(self):
        """return the y-axis position of the cell."""
        return self._position[1]
    
    @y.setter
    def y(self, value):
        self._position[1] = value
    
    
    @property
    def pos_list(self):
        """return the list position (x,y) of integrins in the cell.
        
        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        position = [[], []]
        for integrin in self.integrins:
            if integrin.bound is False:
                position[0].append((integrin.x, integrin.y))
            else:
                position[1].append((integrin.x, integrin.y))

        return position
    
    @property
    def x_pos_list(self):
        """return the list of x-axis position of integrins in the cell.
        
        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        position = [[], []]
        for integrin in self.integrins:
            if integrin.bound is False:
                position[0].append(integrin.x)
            else:
                position[1].append(integrin.x)
        return position
    
    @property
    def y_pos_list(self):
        """return the list of y-axis position of integrins in the cell.
        
        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        position = [[], []]
        for integrin in self.integrins:
            if integrin.bound is False:
                position[0].append(integrin.y)
            else:
                position[1].append(integrin.y)
        return position

    @property
    def alpha_shape(self):
        """return the alphashape perimeter object of the cell."""
        return self._alpha_shape
    
    @alpha_shape.setter
    def alpha_shape(self, value):
        self._alpha_shape = value

    @property
    def area(self):
        """return the area size of the cell based on the alphashape 
        perimeter.
        """
        if self.alpha_shape != None:
            return self.alpha_shape.area*(self._scale**2)
    
    @property
    def total_bound(self):
        """return the total number of bonded integrin."""
        count = 0
        for integrin in self.integrins:
            if integrin.bound is True:
                count += 1
        return count
    
    @property
    def integrins(self):
        """return the list of integrins in the cell."""
        return self._integrins

    @property
    def number_integrin(self):
        """return the total number of integrin in the cell."""
        return len(self._integrins)
    
    @property
    def id(self):
        """return the identity number of the cell"""
        return self._id
    
    @property
    def surfaceIntegrin(self):
        """return the list of surface integrins"""
        listIntegrin: list[Integrin] = []
        for integrin in self.integrins:
            if integrin.isSurface:
                listIntegrin.append(integrin)
        return listIntegrin
    
    @property
    def mass(self):
        """return the mass of the cell"""
        return self._mass
    
    @property
    def kineticEnergy(self):
        """return the value of Kinetic energy of the integrin
        
        This kinetic energy is based on the integrin velocity in x-axis
        and y-axis.
        """
        energy = 0
        for integrin in self.integrins:
            if integrin.bound is False:
                kineticEn = 0.5*integrin.mass*np.dot(integrin.speed, integrin.speed)
                energy += kineticEn
        return energy

    @property
    def normal_length(self):
        """return the normal length of the spring connections between
        integrins
        """
        return self._min_dst    
    @classmethod
    def resetCount(cls):
        """reset the number of cell created into 0."""
        cls.count = 0

class Cells:
    def __init__(self, cellsList=None) -> None:
        """The initial procedure of creating cells object is as 
        follows:

        1. If there are no cell(s) created, read the celcon file 
        and get the integrin and cell properties then create the cell(s).
        2. If there is/are cell(s), then append the cell into cells' 
        member attribute.  
        """
        self._members: list[Cell] = []
        
        if cellsList == None:
            # open file
            celcon = readFile('CELCON')

            # get value
            integrin_size = getValue(celcon, 'ressize')
            integrin_mass = getValue(celcon, 'resmass')
            cell_properties = getValue(celcon, 'CELL')
            
            
            for obj in cell_properties:
                member = Cell(obj[0], obj[1], integrin_size, integrin_mass, obj[2], obj[3])
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
        """ Procedure to get cell by id. 
        
        Parameter
        ---------
        id: int
            the id of the cell.

        Return
        ------
        Cell with the id.
        """
        for cell in self.members:
            if cell.id == id:
                return cell
    
    def excludeCellbyId(self, id):
        """Procedure to exclude a cell from list of cells by Id.
        
        Parameter
        ---------
        id: int
            the id of the cell.

        Return
        ------
        Cells without the cell with specific id.
        """
        new_cells = []
        for cell in self.members:
            if cell.id != id:
                new_cells.append(cell)
        return Cells(new_cells)
    
    def surfaceIntegrinsTarget(self, cell: Cell, max_dist):
        """procedure to get the list of surface integrin of a 
        'cell'.

        Parameter
        ---------
        cell: :obj: Cell
            the cell as the reference. Other surfce integrins not from
            the cell will be treated as Integrin target.
        
        Return
        ------
        surface_integrin
            list of free integrins    

        Notes
        -----
        - this procedure has bound checking, so the result will be 
        free integrin list.
        """
        surface_integrin: list[Integrin] = []
        # find the average radius of the cell
        cell_bounding_rect= Polygon(cell.alpha_shape.minimum_rotated_rectangle)
        cell_diag_length = Point(cell_bounding_rect.exterior.coords[0]).distance(cell_bounding_rect.exterior.coords[2])
        cell_radius = cell_diag_length/2
        for obj in self.members:
            if obj != cell:
                dist = np.linalg.norm(obj.position - cell.position)
                # find the radius of the object
                min_bounding_rect = Polygon(obj.alpha_shape.minimum_rotated_rectangle)
                diag_length = Point(min_bounding_rect.exterior.coords[0]).distance(min_bounding_rect.exterior.coords[2])
                radius = diag_length/2
                # decide wether the cell is in the max_dist radius
                if (dist - (cell_radius + radius)) < max_dist:
                    surface_integrin += obj.surfaceIntegrin
        # filter the integrin
        for integrin in surface_integrin:
            if integrin.bound is True:
                surface_integrin.remove(integrin)
        return surface_integrin


    @property
    def members(self):
        """return list of cells from Cells object"""
        return self._members
    
    @property
    def number_cell(self):
        """return the total number of cells in the Cells object"""
        return len(self._members)
    
    @property
    def x_list(self):
        """return the list of x-axis position of all integrins from
        the member cells

        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin.
        """
        pos_list = [[], []]
        for cell in self.members:
            pos_list[0] += cell.x_pos_list[0]
            pos_list[1] += cell.x_pos_list[1]
        return pos_list
    
    @property
    def y_list(self):
        """return the list of y-axis position of all integrins from
        the member cells

        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        pos_list = [[], []]
        for cell in self.members:
            pos_list[0] += cell.y_pos_list[0]
            pos_list[1] += cell.y_pos_list[1]
        return pos_list
    
    @property
    def many(self):
        """return the condition whether it is single cell in the system
        or multicells. If multicells, it return True
        """
        if self.number_cell > 1:
            return True
        else:
            return False
    
    @property
    def integrin_size(self):
        """return the integrin size"""
        return self.members[0].integrin_size