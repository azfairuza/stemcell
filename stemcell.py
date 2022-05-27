#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from cmath import sqrt
from encodings import utf_8
from random import uniform
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
        position_list: list , dot_size: float):
    # inital procedure when creating a nanopattern
        
        # nanopattern properties
        self.height = height                       # nanopattern height
        self.width = width                         # nanopattern width
        self.grid_height = grid_height             # nanopattern grid height
        self.grid_width = grid_width               # nanopattern grid width
        self.position_seed = position_list         # ligand position's list
        self.dot_size = dot_size                   # ligand size (dot size in simulation)
        row_number = np.ceil(height/grid_height)   # number of grid row created
        col_number = np.ceil(width/grid_width)     # number of grid column created

        # nanopattern ligand members
        self.ligand = []
        for row in range(row_number):
            for col in range(col_number):
                for position in position_list:
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
        self.bound_status = False      # False: intgrin is not bounded to anyone
        self.object_type = 0           # 0: Empty ; 1: Ligand ; 2: Integrin
        self.cell_target_id = 0        # cell id number for integrin-integrin complex
        self.ligand_target_id = 0      # ligand id number for integrin-ligand complex 
        self.integrin_target_id = 0    # integrin id number for integrin-integrin complex
        self.x_target = 0              # x-position of the targeted object
        self.y_target = 0              # y-position of the 
    
    def getInformation(self):
        if self.bound_status == True:
            if self.object_type == 1:
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', ' status: bound to ligand ', str(self.ligand_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
            elif self.object_type == 2:
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', ' status: bound to integrin ', str(self.cell_target_id), '.', str(self.integrin_target_id),  ' at (', str(self.x_target), ',', str(self.y_target), ')')))
        else:
            if self.object_type == 0:
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: not bound')))
            elif self.object_type == 1:
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: targeting ligand ', str(self.ligand_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
            elif self.object_type == 2:
                return(''.join(('integrin ', str(self.cell_id), '.', str(self.integrin_id), ': (', str(self.x_position), ',', str(self.y_position), ')', 'status: targeting integrin ', str(self.cell_target_id), '.', str(self.integrin_target_id), ' at (', str(self.x_target), ',', str(self.y_target), ')')))
    
    def getLigandDistance(self, ):
        pass
    
    def getIntegrinDistance(self, self_cell, target_integrin, target_cell):
        pass

    def searchNearestLigand(self, list_of_ligands, distance):
        
        for ligand in list_of_ligands:
            pass
        pass

class Cell():
    cell_number = 0
       
    def __init__(self, x_center, y_center, max_radius, total_integrin):
        #cell properties
        self.mass = 100
        self.x_center_of_mass = x_center
        self.y_center_of_mass = y_center
        
        #cell integrin members
        self.integrin = []
        Cell.cell_number += 1
        for i in range(total_integrin):
            self.integrin.append(Integrin(Cell.cell_number, x_center, y_center, max_radius))

    def getIntegrinList(self):
        for reseptor in self.integrin:
            reseptor.getInformation()


def readFile(filename):
    print(filename)
    if filename == 'PATCON':
        try:
            datafile = open('PATCON.txt', 'r', encoding='utf-8')
            listofstring = datafile.readlines()
            datafile.close()
            strippedstring = []
            for line in listofstring:
                newline = line.rstrip('\n')
                strippedstring.append(newline)
        except:
            strippedstring = 'Error in opening PATCON file'
        finally:
            return strippedstring
    elif filename == 'CELCON':
        try:
            datafile = open('CELCON.txt', 'r', encoding = 'UTF-8')
            listofstring = datafile.readlines()
            datafile.close()
            strippedstring = []
            for line in listofstring:
                newline = line.rstrip('\n')
                strippedstring.append(newline)
        except:
            strippedstring = 'Error in opening CELCON file'
        finally:
            return strippedstring
    elif filename == 'SIMCON':
        try:
            datafile = open('SIMCON.txt', 'r', encoding='utf-8')
            listofstring = datafile.readlines()
            datafile.close()
            strippedstring = []
            for line in listofstring:
                newline = line.rstrip('\n')
                strippedstring.append(newline)
        except:
            strippedstring = 'Error in opening SIMCON file'
        finally:
            return strippedstring

    else:
        print('File name is not correct')
        return 0

def get_value(listofstring, startIndex, endIndex, arg):
    index = 0
    for i in range(startIndex, endIndex):
        data = listofstring[i].split()
        if data[0] == arg:
            index = i
            break
        else:
            index += 1
    
    if len(data) > 2:
        value = data[1:]
    elif len(data) == 2:
        value = data[1]
    else:
        value = -1
    return value


    
    

    


