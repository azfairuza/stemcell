#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from cmath import sqrt
from encodings import utf_8
from random import uniform
from re import I
import string
import numpy as np

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
        

class Nanopattern():
# class for simulate nanopattern

    def __init__(self, height: float, width: float, grid_height: float, grid_width: float, 
        ligand_position: list , dot_size: float):
    # inital procedure when creating a nanopattern
        
        # nanopattern properties
        self.height = height                       # nanopattern height
        self.width = width                         # nanopattern width
        self.grid_height = grid_height             # nanopattern grid height
        self.grid_width = grid_width               # nanopattern grid width
        self.position_seed = ligand_position       # ligand position's list
        self.dot_size = dot_size                   # ligand size (dot size in simulation)
        row_number = np.ceil(height/grid_height)   # number of grid row created
        col_number = np.ceil(width/grid_width)     # number of grid column created

        # nanopattern ligand members
        self.ligand = []
        for row in range(row_number):
            for col in range(col_number):
                for position in ligand_position:
                    self.ligand.append(Ligand(position[0] + col*grid_width, position[1] + row*grid_height))
        
    
    
    def getNotTargetedLigand(self) -> list:
    # Procedure to get a list of ligand that is not targeted by integrin
        members = []
        for ligand in self.ligand:
            if ligand.bound_status == False:       
                members.append(ligand)
        return members


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
       
    def __init__(self, x_center, y_center, max_radius, total_integrin):
        Cell.cell_number += 1   
        
        #cell properties
        self.mass = 100                         # mass of the cell (in center of mass)
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
    print(filename)
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
    print(property_name)
    if property_name == 'LIGAND':
        start_index = (lst_strng.index('#LIGAND') + 1)
        end_index = lst_strng.index('#END')
        ligand_position = []
        for i in range(start_index, end_index):
            print(lst_strng[i])
            #dot_position = [float(j) for j in lst_strng[i]]
            ligand_position.append(lst_strng[i])
        return ligand_position
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


    
    

    


