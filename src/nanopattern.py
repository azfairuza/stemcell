""" nanopattern module

This module contain the nanopattern class which is the collection of 
ligands which resembles substrate in cell simulation.
"""

# third-party import
import numpy as np

# local import
import physica as psc
import misc
import inputfile as ifile
import ligand as lig


class Nanopattern:
    """The class of collective ligands in a specific pattern

    This object interact with cell through its ligand members.
    """

    def __init__(self, folder_code = None, timestamp = None):
        """Inital procedure when creating a nanopattern.

        The default procedure are as follows
        1. read the PATCON file as the base of nanopattern
        configuration
        2. get the properties of the nanopattern from the `PATCON` file
        3. build the nanopattern using build function

        The bulid from patmap are as follow:
        
        """
        lig.Ligand.reset_count()
        # default building
        if folder_code is None or timestamp is None:
            patcon = ifile.Read("PATCON")

            # get value
            substrate_size = patcon.get("size")
            gridnum = [int(i) for i in patcon.get("gridnum")]

            # nanopattern properties
            self._height = substrate_size[0]
            self._width = substrate_size[1]
            self._x_dist = patcon.get("xdist")
            self._y_dist = patcon.get("ydist")
            self._ligand_size = patcon.get("ligandsize")
            self._ligand_mass = patcon.get("ligandmass")
            self._ligands: list[lig.Ligand] = []
            self._gridnum = gridnum

            # build nanopattern
            self.build()
        # build from patmap
        else:
            #open input file
            try:
                with open(f"./output/{folder_code}-output/input/PATCON.txt", "r", encoding="utf-8") as data_file:
                    lst_strng = data_file.readlines()
                stripped_strng = ifile.filter_item(lst_strng)
            except:
                print("Error in opening PATCON file")
                stripped_strng = "Error in opening PATCON file"
            patcon = ifile.Read("OUTPUT", stripped_strng)

            # get value
            substrate_size = patcon.get("size")
            gridnum = [int(i) for i in patcon.get("gridnum")]

            # nanopattern properties
            self._height = substrate_size[0]
            self._width = substrate_size[1]
            self._x_dist = patcon.get("xdist")
            self._y_dist = patcon.get("ydist")
            self._ligand_size = patcon.get("ligandsize")
            self._ligand_mass = patcon.get("ligandmass")
            self._ligands: list[lig.Ligand] = []
            self._gridnum = gridnum
            self._grid: list[list[list[lig.Ligand]]] = [
                [[] for j in range(self.x_gridnum)] for i in range(self.y_gridnum)
            ]
            
            # open patmap file
            file_dir = f"./output/{folder_code}-output/file/PATMAP/PATMAP{int(timestamp):06}.txt"
            print(f"open: {file_dir}")
            with open(file_dir, "r", encoding="utf-8") as data_file:
                lst_strng = data_file.readlines()
            stripped_strng = ifile.filter_item(lst_strng)
            for i in range(2, len(stripped_strng)):
                obj_data = stripped_strng[i].split()
                obj = lig.Ligand(float(obj_data[4]), float(obj_data[5]), self._ligand_size, self._ligand_mass)
                obj._id = int(obj_data[0])
                obj._bound = bool(int(obj_data[1]))
                self._grid[int(obj_data[3])][int(obj_data[2])].append(obj)
                obj._target_cell_id = int(obj_data[6])
                obj._target_integrin_id = int(obj_data[7])
                self._ligands.append(obj)
            print("SYSTEM: nanopattern has been created")
    
    def sync(self, cells):
        """procedure to sync the data of the target"""
        for ligand_ in self.ligands:
            if ligand_.bound is True:
                if ligand_.target is not None:
                    ligand_._target_cell_id = ligand_.target_cell_id
                    ligand_._target_integrin_id = ligand_.target_integrin_id
                else:
                    ligand_.target = cells.get_integrin_by_id(ligand_.target_cell_id, ligand_.target_integrin_id)


    def build(self):
        """building the nanopattern from empty list of ligands

        The procedure of building are as follows:
        1. Since we can provide variation on the x and y distances
        between ligands. We need to determine which distance to use.
        2. Create the ligand and save it into ligands (list of ligands)
        3. Update the distance and iteration
        4. Create grid. We want to organize the ligands into grid-base
        region. Hence, it will be easier to find nearest neighbors as
        the ligand doesn't move.
        5. Arrange the ligand into appropriate grid base on its
        position

        """
        # create nanopattern
        x_position = 0
        y_position = 0
        x_iter = 0
        y_iter = 0
        while y_position <= self.height:
            y_dist_index = y_iter % len(self.y_dist)
            while x_position <= self.width:
                x_dist_index = x_iter % len(self.x_dist)
                obj = lig.Ligand(x_position, y_position, self._ligand_size, self._ligand_mass)
                self._ligands.append(obj)
                x_position += self.x_dist[x_dist_index]
                x_iter += 1
            y_position += self.y_dist[y_dist_index]
            y_iter += 1
            x_iter = 0
            x_position = 0

        # distribute nanopattern into grids
        self._grid: list[list[list[lig.Ligand]]] = [
            [[] for j in range(self.x_gridnum)] for i in range(self.y_gridnum)
        ]
        unlocated_ligand = self._ligands
        print(f"SYSTEM: There are {len(unlocated_ligand)} ligand(s) ungrouped.")
        for i in range(self.y_gridnum):
            for j in range(self.x_gridnum):
                located_ligand: list[lig.Ligand] = []
                for ligand in unlocated_ligand:
                    if (
                        self.x_grid_points[j]
                        <= ligand.x_position
                        < self.x_grid_points[j + 1]
                    ) and (
                        self.y_grid_points[i]
                        <= ligand.y_position
                        < self.y_grid_points[i + 1]
                    ):
                        located_ligand.append(ligand)
                self._grid[i][j] = located_ligand
                unlocated_ligand = [
                    ligand
                    for ligand in unlocated_ligand
                    if ligand not in located_ligand
                ]
                print(f"SYSTEM: There are {len(unlocated_ligand)} ligand(s) ungrouped.")
        print("SYSTEM: nanopattern has been created")

    def get_ligand_by_id(self, id_: int) -> lig.Ligand:
        """Procedure to get a ligand from ligand members of
        nanopattern by id.
        """
        for ligand in self._ligands:
            if ligand.id_ == id_:
                return ligand
        return None

    def nearest(self, x_position: float, y_position: float, radius):
        """Function to return list of nearest ligand from a position

        Parameters
        ----------
        x_position: float
            the x position of the object reference
        y_position: float
            the y position of the object reference
        radius: float
            the maximum distance of 'nearest' target
        """

        # Find the smallest kernel size needed
        kernel_size_x = int(2 * radius / self.x_gridsize) + 2
        kernel_size_y = int(2 * radius / self.y_gridsize) + 2
        # make sure the kernel size is odd
        if kernel_size_x % 2 == 0:
            kernel_size_x += 1
        if kernel_size_y % 2 == 0:
            kernel_size_y += 1
        # create the moore kernel
        moore_kernel = [
            [True for j in range(kernel_size_x)] for i in range(kernel_size_y)
        ]
        # find the index
        x_index = psc.get_index(x_position, self.x_grid_points)
        y_index = psc.get_index(y_position, self.y_grid_points)
        # put the index at the middle of the moore kernel
        mid_index_x = (kernel_size_x - 1) // 2
        mid_index_y = (kernel_size_y - 1) // 2
        # update kernel base on the position
        # check for the columns
        for k in range(kernel_size_x):
            x_grid_index = x_index + (k - mid_index_x)
            if x_grid_index < 0 or x_grid_index >= (self.x_gridnum - 1):
                for l in range(kernel_size_y):
                    moore_kernel[l][k] = False
        # check for the rows
        for k in range(kernel_size_y):
            y_grid_index = y_index + (k - mid_index_y)
            if y_grid_index < 0 or y_grid_index >= (self.y_gridnum - 1):
                for l in range(kernel_size_x):
                    moore_kernel[k][l] = False
        # get the nearest unbound ligand
        container: list[lig.Ligand] = []
        for i in range(kernel_size_y):
            for j in range(kernel_size_x):
                if moore_kernel[i][j] is True:
                    x_grid_index = x_index + (j - mid_index_x)
                    y_grid_index = y_index + (i - mid_index_y)
                    for ligand_ in self._grid[y_grid_index][x_grid_index]:
                        if ligand_.bound is False:
                            container.append(ligand_)
        # filter again to make circle
        nearest_ligand = misc.filter_by_dist(container, radius, (x_position, y_position))
        return nearest_ligand

    # region <Nanopattern property>

    @property
    def x_grid_points(self):
        """the boundary x point of the nanopattern grids"""
        return np.linspace(0, self.width, (self.x_gridnum + 1))

    @property
    def y_grid_points(self):
        """the boundary y point of the nanopattern grids"""
        return np.linspace(0, self.height, (self.y_gridnum + 1))

    @property
    def x_gridnum(self):
        """the number of grid columns"""
        return self._gridnum[0]

    @property
    def y_gridnum(self):
        """the number of grid rows"""
        return self._gridnum[1]

    @property
    def x_gridsize(self):
        """x size of one grid"""
        return self.width / self.x_gridnum

    @property
    def y_gridsize(self):
        """y size of one grid"""
        return self.height / self.y_gridnum

    @property
    def height(self):
        """height of the substrate (not the single nanopattern)"""
        return self._height

    @property
    def width(self):
        """width of the substrate (not the single nanopattern)"""
        return self._width

    @property
    def ligand_size(self):
        """size of every ligand"""
        return self._ligand_size

    @property
    def pos_list(self):
        """return list of ligands' position (x,y)

        It is divided into two list, first is the position of unbounded
        ligand and second is the position of bounded ligand.
        """
        position = [[], []]
        for ligand in self._ligands:
            if ligand.bound is False:
                position[0].append((ligand.x_position, ligand.y_position))
            else:
                position[1].append((ligand.x_position, ligand.y_position))
        return position

    @property
    def x_pos_list(self):
        """return the x position of ligands

        It is divided into two list, first is the position of unbounded
        ligand and second is the position of bounded ligand.
        """
        position = [[], []]
        for ligand in self._ligands:
            if ligand.bound is False:
                position[0].append(ligand.x_position)
            else:
                position[1].append(ligand.x_position)
        return position

    @property
    def y_pos_list(self):
        """return the y position of ligands

        It is divided into two list, first is the position of unbounded
        ligand and second is the position of bounded ligand.
        """
        position = [[], []]
        for ligand in self._ligands:
            if ligand.bound is False:
                position[0].append(ligand.y_position)
            else:
                position[1].append(ligand.y_position)
        return position

    @property
    def bound_number(self):
        """return the number of ligand with True bound status"""
        count = 0
        for ligand in self._ligands:
            if ligand.bound is True:
                count += 1
        return count

    @property
    def dot_number(self):
        """return the number of ligands"""
        return lig.Ligand.count

    @property
    def ligands(self):
        """return the list of ligands in MultiObjBase"""
        return self._ligands

    @property
    def x_dist(self):
        """return the list of distance between ligands in x-axis"""
        return self._x_dist

    @property
    def y_dist(self):
        """return the list of distance between ligands in y-axis"""
        return self._y_dist

    # endregion