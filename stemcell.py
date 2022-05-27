#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from cmath import sqrt
from encodings import utf_8
from random import uniform
import numpy as np

class Ligand():

    ligand_number = 0
    
    def __init__(self, x_position, y_position):
        Ligand.ligand_number += 1
        self.x_position = x_position
        self.y_position = y_position
        self.bound_status = False     # False status mean that ligand is not bounded to any integrin
        self.targeted_status = False    # False status mean that ligand is not targeted by any integrin
        self.ligand_id = Ligand.ligand_number
        self.integrin_id = 0
        

class Nanopattern():

    def __init__(self, height, width, grid_height, grid_width, position_list, dot_size):
        # nanopattern properties
        self.height = height
        self.width = width
        self.grid_height = grid_height
        self.grid_width = grid_width
        self.position_seed = position_list
        self.dot_size = dot_size
        row_number = np.ceil(height/grid_height)
        col_number = np.ceil(width/grid_width)

        # nanopattern ligand members
        self.ligand = []
        for row in range(row_number):
            for col in range(col_number):
                for position in position_list:
                    self.ligand.append(Ligand(position[0] + col*grid_width, position[2] + row*grid_height))
        
    
    
    def getNotTargetedLigand(self):
        members = []
        for ligand in self.ligand:
            if ligand.bound_status == False:
                members.append(ligand)
        return members


class Integrin():

    number_integrin = 0
    
    def __init__(self, cell_id, x_center, y_center, max_radius):
        # Update integrin number
        Integrin.number_integrin += 1

        # identity
        self.cell_id = cell_id
        self.integrin_id = Integrin.number_integrin
        
        # position property
        radius = round(uniform(0, max_radius), 2)
        theta = round(uniform(0,360), 2)
        self.x_position = radius*np.cos(np.deg2rad(theta)) + x_center
        self.y_position = radius*np.cos(np.deg2rad(theta)) + y_center

        # condition properties
        if radius >= 0.9*max_radius:
            self.surface = True
        else:
            self.surface = False

        # bound properties
        self.bound_status = False
        self.object_type = 0           # 0: Empty ; 1: Ligand ; 2: Integrin
        self.cell_target_id = 0
        self.ligand_target_id = 0
        self.integrin_target_id = 0
        self.x_target = 0
        self.y_target = 0
    
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


    
    

    


