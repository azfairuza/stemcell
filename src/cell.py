"""cell module

This module contains the cell and cells classes which is the collection of 
integrins and cell, respectively.

"""

# third party import
import alphashape
import numpy as np
import shapely.geometry as geo

# local import
import integrin as ign
import inputfile as ifile
import physica as psc


class Cell:
    """Cell represent collection of integrin as a unity."""

    count = 0

    def __init__(
        self,
        x_position: float,
        y_position: float,
        size: float,
        mass: float,
        radius: float,
        min_dst: float,
    ) -> None:
        """init function for cell class.

        This init function also call Cell._build function to generate
        the integrins inside a cell.

        parameters
        ----------
        x_position: float
            x-axis position of the cell center of mass
        y_position: float
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
        self._position = np.array(
            (
                x_position,
                y_position,
            )
        )
        self._integrin_size = size
        self._integrin_mass = mass
        self._min_dst = min_dst
        self._integrins: list[ign.Integrin] = self._build(
            radius, self.x_position, self.y_position, self._min_dst
        )
        self._alpha_shape = None
        self._energy_loss = 0.0
        # self._min_dst = 1.5*min_dst

    def _build(self, radius, x_position, y_position, grid_length):
        """a function to build the cell.

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
        x_position: float
            the x position of cell center of mass.
        y_position: float
            the y position of cell center of mass.
        grid_length: float
            the dist between integrins. The integrins is spread in
            the hexagonal pattern.

        """
        # finding the diagonal and horizontal length
        theta = (1 / 6) * np.pi  # angle: 30 degree
        diagonal = 2 * radius / np.cos(theta)
        horizontal = 2 * radius * (1 + np.tan(theta))

        # create grid base as the minimum lattice
        diag_index = int(np.floor(diagonal / grid_length)) + 1
        hori_index = int(np.floor(horizontal / grid_length)) + 1
        grid: list[list[ign.Integrin]] = [
            [None for j in range(hori_index)] for i in range(diag_index)
        ]
        objs = []

        # check the grid is valid to contain integrin
        # one grid one integrin

        for i in range(diag_index):
            for j in range(hori_index):
                # determine the position of integrin in the grid
                x_dot = (
                    (grid_length * j - grid_length * np.sin(theta) * i)
                    + x_position
                    - radius
                    + grid_length / 2
                )
                y_dot = (
                    (grid_length * np.cos(theta) * i)
                    + y_position
                    - radius
                    + grid_length / 2
                )
                # get the distance from center of the cell (x, y)
                dist = np.sqrt((x_dot - x_position) ** 2 + (y_dot - y_position) ** 2)
                if dist > radius:
                    grid[i][j] = None  # outside the cell
                else:
                    grid[i][j] = ign.Integrin(self, x_dot, y_dot)
                    objs.append(grid[i][j])  # save the integrin

        # find the neighboring
        for i in range(diag_index):
            for j in range(hori_index):
                if grid[i][j] is not None:
                    obj = grid[i][j]
                    # function to check if the index is valid
                    def check_grid(codenum):
                        if codenum==0:
                            if (i - 1 < 0):
                                return False
                            return True
                        if codenum==1:
                            if (j - 1 < 0):
                                return False
                            return True
                        if codenum==2:
                            if (i + 1 >= diag_index):
                                return False
                            return True                        
                        if codenum==3:
                            if (j + 1 >= diag_index):
                                return False
                            return True
                        

                    # input the neighboring integrin
                    if check_grid(0) and check_grid(1) and grid[i - 1][j - 1] is not None:
                        obj.neighbors.append(grid[i - 1][j - 1])
                    if check_grid(0) and grid[i - 1][j] is not None:
                        obj.neighbors.append(grid[i - 1][j])
                    if check_grid(1) and grid[i][j - 1] is not None:
                        obj.neighbors.append(grid[i][j - 1])
                    if check_grid(3) and grid[i][j + 1] is not None:
                        obj.neighbors.append(grid[i][j + 1])
                    if check_grid(2) and grid[i + 1][j] is not None:
                        obj.neighbors.append(grid[i + 1][j])
                    if check_grid(2) and check_grid(3) and grid[i + 1][j + 1] is not None:
                        obj.neighbors.append(grid[i + 1][j + 1])
        return objs

    def get_integrin_by_id(self, id_):
        """procedure to get a specific integrin by its ID.

        Parameter
        --------
        id_: int
            the id of the integrin.

        Return
        ------
        integrin with the id mentioned.
        """
        for integrin_ in self.integrins:
            if integrin_.id_ == id_:
                return integrin_
            return None

    def update_position(self):
        """Procedure to update the position of center of mass."""
        center_of_mass = np.array([0, 0], dtype=float)
        for integrin_ in self.integrins:
            center_of_mass += integrin_.position
        return center_of_mass / self.number_integrin

    def update_alphashape(self, alpha_value=0):
        """Alphashape is a method to create the surrounding area of
        the cell, by updating the alphashape, we get the latest shape
        based on integrin position.

        Parameters
        ----------
        alphavalue, default 0
            this value determines how tight the surrounding area
        """
        dot_position = self.pos_list[0] + self.pos_list[1]
        self.alpha_shape= alphashape.alphashape(dot_position, alpha_value)

    @property
    def integrin_size(self):
        """return the size of the integrins in the cell"""
        return self._integrin_size

    @property
    def position(self):
        """return the position of the cell."""
        return self._position

    @property
    def x_position(self):
        """return the x-axis position of the cell."""
        return self._position[0]

    @x_position.setter
    def x_position(self, value):
        self._position[0] = value

    @property
    def y_position(self):
        """return the y-axis position of the cell."""
        return self._position[1]

    @y_position.setter
    def y_position(self, value):
        self._position[1] = value

    @property
    def pos_list(self):
        """return the list position (x,y) of integrins in the cell.

        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        position = [[], []]
        for integrin_ in self.integrins:
            if integrin_.bound is False:
                position[0].append((integrin_.x_position, integrin_.y_position))
            else:
                position[1].append((integrin_.x_position, integrin_.y_position))

        return position

    @property
    def x_pos_list(self):
        """return the list of x-axis position of integrins in the cell.

        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        position = [[], []]
        for integrin_ in self.integrins:
            if integrin_.bound is False:
                position[0].append(integrin_.x_position)
            else:
                position[1].append(integrin_.x_position)
        return position

    @property
    def y_pos_list(self):
        """return the list of y-axis position of integrins in the cell.

        It is divided into two list, first is the position of unbounded
        integrin and second is the position of bounded integrin
        """
        position = [[], []]
        for integrin_ in self.integrins:
            if integrin_.bound is False:
                position[0].append(integrin_.y_position)
            else:
                position[1].append(integrin_.y_position)
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
        if self.alpha_shape is not None:
            return self.alpha_shape.area
        return -1

    @property
    def total_bound(self):
        """return the total number of bonded integrin."""
        count = 0
        for integrin_ in self.integrins:
            if integrin_.bound is True:
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
    def id_(self):
        """return the identity number of the cell"""
        return self._id

    @property
    def surface_integrin(self):
        """return the list of surface integrins"""
        integrin_list: list[ign.Integrin]  = []
        for integrin_ in self.integrins:
            if integrin_.issurface:
                integrin_list.append(integrin_)
        return integrin_list

    @property
    def mass(self):
        """return the mass of the cell"""
        return self._mass

    @property
    def kinetic_energy(self):
        """return the value of Kinetic energy of the integrin

        This kinetic energy is based on the integrin velocity in x-axis
        and y-axis.
        """
        energy = 0.0
        for integrin_ in self.integrins:
            # if not integrin_.bound:
            energy = energy + integrin_.kinetic_energy
        return energy
    
    @property
    def potential_energy(self):
        """return the value of Kinetic energy of the integrin

        This potential energy is based on the integrin: 
        1. spring connection.
        2. Lennard-jones potential
        """
        energy = 0.0
        for integrin_ in self.integrins:
           energy = energy + integrin_._potential_energy
        return energy - self._energy_loss

    @property
    def bonding_energy(self):
        """return the value of total Bonding energy of the integrin"""
        energy = 0.0
        for integrin_ in self.integrins:
            if integrin_.bound:
                energy = energy + integrin_.bonding_energy
        return energy

    @property
    def normal_length(self):
        """return the normal length of the spring connections between
        integrins
        """
        return self._min_dst

    @property
    def radius(self):
        """return the radius of the smallest circle which surroudnd the
        the cell
        """
        dist = -1
        for integrin in self.surface_integrin:
            calc_dist = np.linalg.norm(integrin.position-self.position)
            if dist < 0:
                dist = calc_dist
            elif dist < calc_dist:
                dist = calc_dist
        return dist
    @classmethod
    def reset_count(cls):
        """reset the number of cell created into 0."""
        cls.count = 0


class Cells:
    """Collection of cell objects"""

    def __init__(self, cells_list=None) -> None:
        """The initial procedure of creating cells object is as
        follows:

        1. If there are no cell(s) created, read the celcon file
        and get the integrin and cell properties then create the cell(s).
        2. If there is/are cell(s), then append the cell into cells'
        member attribute.
        """
        self._members: list[Cell] = []

        if cells_list is None:
            # open file
            celcon = ifile.Read("CELCON")

            # get value
            integrin_size = celcon.get("ressize")
            integrin_mass = celcon.get("resmass")
            cell_properties = celcon.get("CELL")

            for obj in cell_properties:
                member = Cell(
                    obj[0], obj[1], integrin_size, integrin_mass, obj[2], obj[3]
                )
                self._members.append(member)

        elif isinstance(cells_list, list):
            for cell in cells_list:
                if isinstance(cell, Cell):
                    self._members.append(cell)
        elif isinstance(cells_list, Cell):
            self._members.append(cells_list)

        else:
            pass

        print(f"SYSTEM: {len(self.members)} cell(s) is created")

    def get_cell_by_id(self, id_):
        """Procedure to get cell by id.

        Parameter
        ---------
        id_: int
            the id of the cell.

        Return
        ------
        Cell with the id.
        """
        for cell in self.members:
            if cell.id_ == id_:
                return cell
        return None
    
    def get_integrin_by_id(self, cell_id, integrin_id):
        """Procedure to get cell by id.

        Parameter
        ---------
        cell_id: int
            the id of the cell.
        integrin_id: int
            the id of the integrin

        Return
        ------
        integrin with the cell_id and integrin_id.
        """
        for cell in self.members:
            if cell.id_ == cell_id:
                for integrin_ in cell.integrins:
                    if integrin_.id_ == integrin_id:
                        return integrin_
        return None


    def exclude_cell_by_id(self, id_):
        """Procedure to exclude a cell from list of cells by Id.

        Parameter
        ---------
        id_: int
            the id of the cell.

        Return
        ------
        Cells without the cell with specific id.
        """
        new_cells = []
        for cell in self.members:
            if cell.id_ != id_:
                new_cells.append(cell)
        return Cells(new_cells)

    def surface_integrins_target(self, cell: Cell, max_dist, filter_bound: bool = True):
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
        surface_integrin: list[ign.Integrin] = []
        for obj in self.members:
            if obj != cell:
                dist = np.linalg.norm(obj.position - cell.position)
                # decide wether the cell is in the max_dist radius
                if (dist - (cell.radius + obj.radius)) < max_dist:
                    surface_integrin += obj.surface_integrin
        # filter the integrin
        if filter_bound is True:
            for integrin_ in surface_integrin:
                if integrin_.bound is True:
                    surface_integrin.remove(integrin_)
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
        return bool(self.number_cell > 1)

    @property
    def integrin_size(self):
        """return the integrin size"""
        return self.members[0].integrin_size
