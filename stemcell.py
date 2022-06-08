#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from cmath import sqrt
from encodings import utf_8
from random import uniform
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

class Ligand():
# class for simulate ligand in nanopattern

    # Class properties
    ligand_number = 0     #total ligand created
    
    def __init__(self, x_position: float, y_position: float):
    # initial procedure when creating a ligand
        
        Ligand.ligand_number += 1                # Updating ligand number

        # ligand properties
        self.x_position = x_position             # ligand x-position
        self.y_position = y_position             # ligand y-position
        self.bound_status = False                # False: ligand is not bounded to any integrin
        self.targeted_status = False             # False: ligand is not targeted by any integrin
        self.ligand_id = Ligand.ligand_number    # id number for ligand
        self.integrin_id = 0                     # id number for integrin connected to the ligand
    
    def resetLigandNumber():
        Ligand.ligand_number = 0
        

class Nanopattern():
# class for simulate nanopattern

    def __init__(self, height: float, width: float, grid_height: float, grid_width: float, 
        ligand_position: list , dot_size: float):
    # inital procedure when creating a nanopattern
        
        # nanopattern properties
        self.height = height                            # nanopattern height
        self.width = width                              # nanopattern width
        self.grid_height = grid_height                  # nanopattern grid height
        self.grid_width = grid_width                    # nanopattern grid width
        self.position_seed = ligand_position            # ligand position's list
        self.dot_size = dot_size                        # ligand size (dot size in simulation)
        row_number = int(np.ceil(height/grid_height))   # number of grid row created
        col_number = int(np.ceil(width/grid_width))     # number of grid column created

        # nanopattern ligand members
        self.ligand = []
        for row in range(row_number):
            for col in range(col_number):
                for position in ligand_position:
                    self.ligand.append(Ligand((position[0] + col)*grid_width, (position[1] + row)*grid_height))

        print('nanopattern has been created')
    
    def getNotTargetedLigand(self) -> list:
    # Procedure to get a list of ligand that is not targeted by integrin
        members = []
        for ligand in self.ligand:
            if ligand.bound_status == False:       
                members.append(ligand)
        return members
    
    def getXPositionLigand(self):
    # Procedure to get x position of ligands in a list
        position = []
        for ligand in self.ligand:
            position.append(ligand.x_position)
        return position
    
    def getYPositionLigand(self):
    # Procedure to get y position of ligand in a list
        position = []
        for ligand in self.ligand:
            position.append(ligand.y_position)
        return position

    def show(self):
    # Procedure to draw only nanopattern
        plt.figure(figsize=(20, 20))
        ax = plt.subplot(aspect='equal')
        x_position = self.getXPositionLigand()
        y_position = self.getYPositionLigand()

        out = circles(x_position, y_position, self.dot_size, 'green', alpha=0.5, ec='none')
        plt.xlim(0, self.width)
        plt.ylim(0, self.height)



class Integrin():
# class to simulate integrin

    # class properties
    number_integrin = 0
    
    def __init__(self, cell_id, x_center, y_center, max_radius):
        
        Integrin.number_integrin += 1         # Update integrin number

        # identity of integrin
        self.cell_id = cell_id                       # cell id number
        self.integrin_id = Integrin.number_integrin  # integrin id number 
        
        # position property
        radius = round(uniform(0, max_radius), 2)                      # distance from the cell's center
        theta = round(uniform(0,360), 2)                               # the angle from x-axis
        self.x_position = radius*np.cos(np.deg2rad(theta)) + x_center  # convert polar to x-position
        self.y_position = radius*np.cos(np.deg2rad(theta)) + y_center  # convert polat to y-position

        # condition properties
        if radius >= 0.9*max_radius:
            self.surface = True        # definition: integrin is on the cell surface
        else:
            self.surface = False       # definition: integrin is in the inside of the cell

        # bound properties
        self.bound_status = False      # False: integrin is not bounded to anyone
        self.object_type = 0           # 0: Empty ; 1: Ligand ; 2: Integrin
        self.cell_target_id = 0        # cell id number for integrin-integrin complex
        self.ligand_target_id = 0      # ligand id number for integrin-ligand complex 
        self.integrin_target_id = 0    # integrin id number for integrin-integrin complex
        self.x_target = 0              # x-position of the targeted object
        self.y_target = 0              # y-position of the targeted object
    
    def getInformation(self):
    # procedure to get integrin information 
        if self.bound_status == True:
            if self.object_type == 1:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: bound to ligand [ligand_id] at ([x_ligand],[y_ligand])
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', ' status: bound to ligand ', str(self.ligand_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
            
            elif self.object_type == 2:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: bound to integrin [cell_id*].[integrin_id*] at ([x_integrin*],[y_integrin*])
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', ' status: bound to integrin ', str(self.cell_target_id), '.', str(self.integrin_target_id),  ' at (', str(self.x_target), ',', str(self.y_target), ')')))
        
        else:
            if self.object_type == 0:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: not bound
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: not bound')))
            
            elif self.object_type == 1:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: targeting ligand [ligand_id] at ([x_ligand],[y_ligand])
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: targeting ligand ', str(self.ligand_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
            
            elif self.object_type == 2:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: targeting integrin [cell_id*].[integrin_id*]
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: targeting integrin ', str(self.cell_target_id), '.', str(self.integrin_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
    
    def getLigandDistance(self, ):
    # procedure to get distance to a ligand
        pass
    
    def getIntegrinDistance(self, self_cell, target_integrin, target_cell):
    # procedure to get distance to an integrin
        pass

    def searchNearestLigand(self, lst_of_ligands, distance):
    # procedure to get Nearest Ligand, returned in a list
    
        for ligand in lst_of_ligands:      # search every ligand
            pass
        pass

class Cell():
# class to simulate cell

    cell_number = 0       # number of cell created
       
    def __init__(self, x_center, y_center, max_radius, total_integrin, cell_mass):
        Cell.cell_number += 1   
        
        #cell properties
        self.mass = cell_mass                        # mass of the cell (in center of mass)
        self.x_center_of_mass = x_center        # x position of the center of mass of the cell
        self.y_center_of_mass = y_center        # y position of the center of mass of the cell
        
        #cell integrin members
        self.integrin = []                                                                      # empty list for integrin                                                                # 
        for i in range(total_integrin):                                                         
            self.integrin.append(Integrin(Cell.cell_number, x_center, y_center, max_radius))    # build the integrin

    def getIntegrinList(self):
    # get information of all integrins in the cell
        for reseptor in self.integrin:
            reseptor.getInformation()


def readFile(filename):
# procedure to read input file

    if filename == 'PATCON':
    # open PATCON file
        try:
            data_file = open('PATCON.txt', 'r', encoding='utf-8')    
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filter_element(lst_strng)       # filter out the empty string and '\n'              
        except:
            stripped_strng = 'Error in opening PATCON file'    # error flag if PATCON can not be opened
        finally:
            return stripped_strng

    elif filename == 'CELCON':
    # open CELCON file
        try:
            data_file = open('CELCON.txt', 'r', encoding = 'UTF-8')
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filter_element(lst_strng)       # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening CELCON file'    # error flag if CELCON can not be opened
        finally:
            return stripped_strng

    elif filename == 'SIMCON':
    # open SIMCON file
        try:
            data_file = open('SIMCON.txt', 'r', encoding='utf-8')
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filter_element(lst_strng)       # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening SIMCON file'    # error flag if SIMCON can not be opened
        finally:
            return stripped_strng

    else:
        print('File name is not correct')
        return 0

def get_value(lst_strng: list, property_name: string):
# procedure to get the value of a property in CON file
    if property_name == 'LIGAND':                                       # get ligand position
        start_index = (lst_strng.index('#LIGAND') + 2)                  # start index to search
        end_index = lst_strng.index('#END')                             # last index to search
        ligand_position = []                                            # container for the ligands
        for i in range(start_index, end_index):
            dot_position = [float(j) for j in lst_strng[i].split()]     # make the value in float type
            ligand_position.append(dot_position)                        # also remove the spaces and tabs
        return ligand_position
    elif property_name == 'CELL':                                       # get cell properties
        start_index = (lst_strng.index('#CELL')+2)                      # start index to search
        end_index = lst_strng.index('#END')                             # last index to search
        cell_properties = []                                            # cell properties container
        for i in range(start_index, end_index):
            cell_property = [float(j) for j in lst_strng[i].split()]    # make the value in float type
            cell_properties.append(cell_property)                       # also remove the spaces and tabs
        return cell_properties
    else: 
        start_index = (lst_strng.index('#CONFIG') + 1)
        end_index =  lst_strng.index('#END')                                
        value = -1 
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

def filter_element(input_lst):
# procedure to remove empty element in a list and remove the '\n' i
    new_lst = []
    for element in input_lst:
        new_element = element.rstrip()          # remove the '\n'    
        if new_element:                         # check if the element is not empty
            new_lst.append(new_element)
    return new_lst

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
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
    """
  
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

    
    

    


