#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from encodings import utf_8
from random import uniform
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from typing import List
from pathlib import Path

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
        self.object_type = 'ligand'
    
    def isTargeted(self):
    # procedure to get targeted status
        return self.targeted_status
    
    def getXpos(self):
    # procedure to get X position
        return self.x_position
    
    def getYpos(self):
    # procedure to get Y position
        return self.y_position
    
    def getId(self):
    # procedure to get Id
        return self.ligand_id

    @classmethod
    def resetNumber(cls):
    # procedure to reset ligand number
        cls.ligand_number = 0
        
class Nanopattern():
# class for simulate nanopattern

    def __init__(self, height: float, width: float, grid_height: float, grid_width: float, 
        ligand_position: list , dot_size: float):
    # inital procedure when creating a nanopattern
        
        # nanopattern properties
        self.height = height                                 # nanopattern height
        self.width = width                                   # nanopattern width
        self.grid_height = grid_height                       # nanopattern grid height
        self.grid_width = grid_width                         # nanopattern grid width
        self.position_seed = ligand_position                 # ligand position's list
        self.dot_size = dot_size                             # ligand size (dot size in simulation)
        self.row_number = int(np.ceil(height/grid_height))   # number of grid row created
        self.col_number = int(np.ceil(width/grid_width))     # number of grid column created
        self.object_type = 'nanopattern'

        # nanopattern ligand members
        self.ligand: List[Ligand] = []
        for row in range(self.row_number):
            for col in range(self.col_number):
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
    
    def getLigandById(self, id):
    # Procedure to get a ligand by id
        for ligand in self.ligand:
            if ligand.getId() == id:
                return ligand

    def show(self, save=False, number=0, folder=None):
    # Procedure to draw only nanopattern

        if folder is None:
            namefolder = './figure'
        else:
            namefolder = f'./figure/{folder}'
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        
        namefile = f'{namefolder}/{number:06}.jpg'
        
        fig = plt.figure(figsize=(20, 20))
        fig.dpi = 200
        ax = plt.subplot(aspect='equal')
        x_position = self.getXPositionLigand()
        y_position = self.getYPositionLigand()

        out = circles(x_position, y_position, self.dot_size, 'green', alpha=0.5, ec='none')
        plt.xlim(0, self.width)
        plt.ylim(0, self.height)

        if save is True:
            print(f'figure saved on {namefolder}')
            fig.savefig(namefile, bbox_inches='tight', dpi=600)
            plt.close()

    
class Integrin():
# class to simulate integrin

    # class properties
    integrin_number = 0
    
    def __init__(self, cell_id, x_center, y_center, max_radius, integrin_size, lst_integrin):
        
        Integrin.integrin_number += 1                   # Update integrin number

        # identity of integrin
        self.cell_id = cell_id                          # cell id number
        self.integrin_id = Integrin.integrin_number     # integrin id number 
        
        # position property
        while True:
            radius = round(uniform(0, max_radius), 2)                      # distance from the cell's center
            theta = round(uniform(0,360), 2)                               # the angle from x-axis
            self.x_position = radius*np.cos(np.deg2rad(theta)) + x_center  # convert polar to x-position
            self.y_position = radius*np.sin(np.deg2rad(theta)) + y_center  # convert polat to y-position
            if lst_integrin:
                short_distance = False
                for obj in lst_integrin:
                    distnce = np.sqrt((self.x_position - obj.x_position)**2 + (self.y_position - obj.y_position)**2)
                    if distnce < (integrin_size*4):
                        short_distance = True
                        break
                if short_distance == False:
                    break
            else:
                break

        # condition properties
        if radius >= 0.9*max_radius:
            self.surface = True        # definition: integrin is on the cell surface
        else:
            self.surface = False       # definition: integrin is in the inside of the cell

        # bound properties
        self.bound_status = False      # False: integrin is not bounded to anyone
        self.targeting_status = False  # False: integrin is not targeting to anyone
        self.object_type = 0           # 0: Empty ; 1: Ligand ; 2: Integrin
        self.cell_target_id = 0        # cell id number for integrin-integrin complex
        self.ligand_target_id = 0      # ligand id number for integrin-ligand complex 
        self.integrin_target_id = 0    # integrin id number for integrin-integrin complex
        self.x_target = 0              # x-position of the targeted object
        self.y_target = 0              # y-position of the targeted object
        self.mass = 0                  # integrin mass
    
    def getXpos(self):
    # procedure to get X position
        return self.x_position
    
    def getYpos(self):
    # procedure to get Y position
        return self.y_position
    
    def getId(self):
    # procedure to get integrin id
        return self.integrin_id

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
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: not targeting')))
            
            elif self.object_type == 1:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: targeting ligand [ligand_id] at ([x_ligand],[y_ligand])
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: targeting ligand ', str(self.ligand_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
            
            elif self.object_type == 2:
            # return: integrin [cell_id].[integrin_id]: ([x_integrin],[y_integrin]) status: targeting integrin [cell_id*].[integrin_id*]
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: targeting integrin ', str(self.cell_target_id), '.', str(self.integrin_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
    
    def getTargetDistance(self):
        return np.sqrt((self.x_position-self.x_target)**2 + (self.y_position-self.y_target)**2)

    def getLigandDistance(self, ligand: Ligand):
    # procedure to get distance to a ligand
        return np.sqrt((self.getXpos() - ligand.getXpos())**2 + (self.getYpos() - ligand.getYpos())**2)
    
    def getIntegrinDistance(self, target_integrin):
    # procedure to get distance to an integrin
        return np.sqrt((self.getXpos() - target_integrin.getXpos())**2 + (self.getYpos() - target_integrin.getYpos())**2)
    
    def isSurface(self):
        return self.surface
    
    def isTargeting(self):
        return self.targeting_status
    
    def isBound(self):
        return self.bound_status
    
    def move(self, movespeed):
        # find the angle
        # do trigonometri on the movespeed
        # update the x pos and y pos
        angle = np.arctan2((self.y_target - self.y_position), (self.x_target - self.x_position))
        dx = movespeed*np.cos(angle)
        dy = movespeed*np.sin(angle)
        self.x_position += dx
        self.y_position += dy
    
    def updateTarget(self, cells):
        if self.object_type == 2:
            cell_target = getCellbyId(cells, self.cell_target_id)
            target = cell_target.getIntegrinById(self.integrin_target_id)
            self.x_target = target.x_position
            self.y_target = target.y_position
        else: 
            pass
    
    @classmethod
    def resetNumber(cls):
    # procedure to reset Integrin number
        cls.integrin_number = 0

class Cell():
# class to simulate cell

    cell_number = 0       # number of cell created
       
    def __init__(self, x_center, y_center, max_radius, total_integrin, cell_mass, integrin_size):
        Cell.cell_number += 1   
        
        #cell properties
        self.mass = cell_mass                   # mass of the cell (in center of mass)
        self.x_center_of_mass = x_center        # x position of the center of mass of the cell
        self.y_center_of_mass = y_center        # y position of the center of mass of the cell
        self.integrin_size = integrin_size      # the size of integrin
        self.radius = max_radius                # the outer radius of the cell
        self.cell_id = Cell.cell_number         # cell_id
        self.object_type = 'cell'
        
        #cell integrin members
        self.integrin: List[Integrin] = []        
        max_number = int(max_radius**2/(8*(integrin_size**2)))    
        if total_integrin <= max_number:                                           # empty list for integrin                                                                # 
            for i in range(int(total_integrin)):      
                build_integrin = Integrin(Cell.cell_number, x_center, y_center, 
                    max_radius, integrin_size, self.integrin)                                                   
                self.integrin.append(build_integrin)                                 # build the integrin
        else:
            print(''.join(('maximum number integrin allowed is ', str(max_number))))

    def getIntegrinList(self):
    # get information of all integrins in the cell
        for reseptor in self.integrin:
            reseptor.getInformation()

    def getXPositionIntegrin(self):
    # procedure to get x position of integrins, contained in a list
        position = []
        for integrin in self.integrin:
            position.append(integrin.x_position)
        return position
    
    def getYPositionIntegrin(self):
    # procedure to get y position of integrins, contained in a list
        position = []
        for integrin in self.integrin:
            position.append(integrin.y_position)
        return position

    def getId(self):
    # procedure to get cell's Id
        return self.cell_id
    
    def getIntegrinById(self, id):
    # procedure to get integrin by Id
        for integrin in self.integrin:
            if integrin.getId() == id:
                return integrin
    
    def getUntargetedIntegrin(self):
        untargeted_integrin: List[Integrin] = []
        for integrin in self.integrin:
            if integrin.targeting_status == False:
                untargeted_integrin.append(integrin)
        return untargeted_integrin

    def show(self, substrate: Nanopattern, save=False, number=0, folder=None):
    # procedure to draw only nanopattern
        
        if folder is None:
            namefolder = './figure'
        else:
            namefolder = f'./figure/{folder}'
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        
        namefile = f'{namefolder}/{number:06}.jpg'
        
        fig = plt.figure(figsize=(20, 20))
        fig.dpi = 200
        ax = plt.subplot(aspect='equal')
        x_position = self.getXPositionIntegrin()                                        
        y_position = self.getYPositionIntegrin()
        surface_circle = plt.Circle((self.x_center_of_mass, self.y_center_of_mass), 
            self.radius, color='black', fill=False, ls='--')
        plt.gca().add_patch(surface_circle)

        out = circles(x_position, y_position, self.integrin_size, 'red', alpha=0.5, ec='none')
        plt.xlim(0, substrate.width)
        plt.ylim(0, substrate.height)

        if save is True:
            print(f'figure saved on {namefolder}')
            fig.savefig(namefile, bbox_inches='tight', dpi=600)
            plt.close()

    @classmethod
    def resetNumber(cls):
    # procedure to reset cell number
        cls.cell_number = 0

def readFile(filename):
# procedure to read input file

    if filename == 'PATCON':
    # open PATCON file
        try:
            data_file = open('PATCON.txt', 'r', encoding='utf-8')    
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)       # filter out the empty string and '\n'              
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
            stripped_strng = filterElement(lst_strng)       # filter out the empty string and '\n' 
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
            stripped_strng = filterElement(lst_strng)       # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening SIMCON file'    # error flag if SIMCON can not be opened
        finally:
            return stripped_strng

    else:
        print('File name is not correct')
        return 0

def getValue(lst_strng: list, property_name: string):
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

def filterElement(input_lst: list):
# procedure to remove empty element in a list and remove the '\n' i
    new_lst = []
    for element in input_lst:
        new_element = element.rstrip()          # remove the '\n'    
        if new_element:                         # check if the element is not empty
            new_lst.append(new_element)
    return new_lst

def getCellbyId(cells: List[Cell], id):
# procedure to get cell by id
    for cell in cells:
        if cell.getId() == id:
            return cell

def showAll(cells : List[Cell], substrate: Nanopattern, 
    show_substrate=False, save=False, number=0, folder=None, line=False):
# procedure to show all element of simulation including cells and nanopattern
    
    if folder is None:
        namefolder = './figure'
    else:
        namefolder = f'./figure/{folder}'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
        
    namefile = f'{namefolder}/{number:06}.jpg'
    
    fig = plt.figure(figsize=(20,20))
    fig.dpi = 200
    ax = plt.subplot(aspect='equal')

    # print the cells
    for cell in cells:
        if line == True:
            for integrin in cell.integrin:
                if integrin.targeting_status:
                    x = [integrin.x_position, integrin.x_target]
                    y = [integrin.y_position, integrin.y_target]
                    plt.plot(x, y, color="black", linewidth="1", linestyle="dashed")
        x_position = cell.getXPositionIntegrin()
        y_position = cell.getYPositionIntegrin()
        surface_circle = plt.Circle((cell.x_center_of_mass, cell.y_center_of_mass), 
            cell.radius, color='black', fill=False, ls='--')
        plt.gca().add_patch(surface_circle)
        out = circles(x_position, y_position, cell.integrin_size, 'red', alpha=0.5, ec='none')
    
    # print the nanopattern if the option is true
    if show_substrate is True:
        x_position = substrate.getXPositionLigand()
        y_position = substrate.getYPositionLigand()
        out = circles(x_position, y_position, substrate.dot_size, 'green', alpha=0.5, ec = 'none')
    
    
    plt.xlim(0, substrate.width)
    plt.ylim(0, substrate.height)

    if save is True:
        print(f'figure {number:06}.jpg saved on {namefolder}')
        fig.savefig(namefile, bbox_inches='tight', dpi=600)
        plt.close()

def excludeCellById(cells: List[Cell], cell_id: int):
    new_cells: List[Cell] = []
    for cell in cells:
        if cell.getId() != cell_id:
            new_cells.append(cell)
    return new_cells

def simulate1(cells: List[Cell], substrate: Nanopattern, dstlimit: float, 
    n_iteration: int=1, movespeed: float=1, line: bool=False):
# procedure to do simulation with single targeting    
    iter_simulation = 0

    # start the simulation
    # Targeting
    # pick a cell
    for cell in cells:
        #pick an integrin from the cell
        for integrin in cell.integrin:
            min_distance = 0
            target = None
            if integrin.targeting_status == False:
                #check wether the integrin is in surface and also there are more than 1 cell
                #print('integrin is not targeting')
                if integrin.surface and len(cells) > 1:
                    # get untargeted integrin from other cells
                    print(f'integrin {integrin.integrin_id} is in surface')
                    cell_targets = excludeCellById(cells, cell.getId())
                    for cell_target in cell_targets:
                        print(f'cell id: {cell_target.cell_id}')
                        for integrin_target in cell_target.integrin:
                            # check whether the integrin has target
                            if integrin_target.targeting_status == False:
                                if target == None:
                                    target = integrin_target
                                    # update min_distance
                                    min_distance = integrin.getIntegrinDistance(integrin_target)
                                    print(f'target is integrin with id {target.integrin_id}')
                                else:
                                    dstnce = integrin.getIntegrinDistance(integrin_target)
                                    if dstnce < min_distance:
                                        target = integrin_target
                                        # update min_distance
                                        min_distance = dstnce
                                        print(f'target is integrin with id {target.integrin_id}')
                    # get untargeted ligand from nanopattern
                    for ligand in substrate.ligand:
                        # check whether the ligand has targeted
                        #print('checking ligands for surface integrin')
                        if ligand.targeted_status == False:
                            if target == None:
                                target = ligand
                                # update min_distance
                                min_distance = integrin.getLigandDistance(ligand)
                                print(f'target is ligand with id {target.ligand_id}')
                            else:
                                dstnce = integrin.getLigandDistance(ligand)
                                if dstnce < min_distance:
                                    target = ligand
                                    # update min_distance
                                    min_distance = dstnce
                                    print(f'target is ligand with id {target.ligand_id}')
                else:
                    # get untargeted ligand from nanopattern
                    print(f'finding ligand for non surface integrin {integrin.integrin_id}')
                    for ligand in substrate.ligand:
                        # check whether the ligand has targeted
                        if ligand.targeted_status == False:
                            if target == None:
                                target = ligand
                                # update min_distance
                                min_distance = integrin.getLigandDistance(ligand)
                                print(f'target is ligand with id {target.ligand_id}')
                            else:
                                dstnce = integrin.getLigandDistance(ligand)
                                if dstnce < min_distance:
                                    target = ligand
                                    # update min_distance
                                    min_distance = dstnce
                                    print(f'target is ligand with id {target.ligand_id}')
                # update the attribute
                #print(target.__dict__)
                if target is not None:
                    if target.object_type == 'ligand':
                        print('update attribute due to ligand')
                        # integrin
                        integrin.targeting_status = True
                        integrin.object_type = 1                        # 0: Empty ; 1: Ligand ; 2: Integrin
                        integrin.ligand_target_id = target.getId()
                        integrin.x_target = target.getXpos()
                        integrin.y_target = target.getYpos()

                        # ligand target 
                        target.targeted_status = True
                        target.integrin_id = integrin.getId()

                    else:
                        print('update attribute due to integrin')
                        # integrin
                        integrin.targeting_status = True
                        integrin.object_type = 2
                        integrin.cell_target_id = target.cell_id
                        integrin.integrin_target_id = target.getId()
                        integrin.x_target = target.getXpos()
                        integrin.y_target = target.getYpos()

                        # integrin target
                        target.targeting_status = True
                        target.object_type = 2
                        target.cell_target_id = integrin.cell_id
                        target.integrin_target_id = integrin.getId()
                        target.x_target = integrin.getXpos()
                        target.y_target = integrin.getYpos()
    #print('targeting is finished')

    # Moving
    while(iter_simulation <= n_iteration): 
        for cell in cells:
            for integrin in cell.integrin:
                if integrin.targeting_status:
                    # check if the integrin is not bounded and the distance is below the limit
                    target_dst = integrin.getTargetDistance()
                    #print(target_dst)
                    if (integrin.bound_status == False) and (dstlimit > target_dst):
                        # check if the integrin can bounded
                        #print('integrin is not bound and below the distance limit')
                        if target_dst < 2*cell.integrin_size:
                            # update the bound status and mass of the integrin
                            integrin.bound_status = True
                            integrin.mass = 1
                            #print('integrin is bound')
                            if integrin.object_type == 1:
                                # if the integrin bound to ligand
                                target = substrate.getLigandById(integrin.ligand_target_id)
                                target.bound_status = True
                            elif integrin.object_type == 2:
                                # if the integrin bound to another integrin
                                cell_target = getCellbyId(cells, integrin.cell_target_id)
                                integrin_target = cell_target.getIntegrinById(integrin.integrin_target_id)
                                integrin_target.bound_status = True
                                integrin_target.mass = 1
                        else:
                            #print('integrin moving')
                            integrin.updateTarget(cells)
                            integrin.move(movespeed)
        showAll(cells, substrate, show_substrate=True, save=True, 
            folder='simulation_output', number=iter_simulation, line=line)
        iter_simulation += 1


                
    

        
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

    
    

    


