#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from cmath import sqrt
from random import uniform
import numpy as np

class Ligand():
    
    def __init__(self, x_position, y_position):
        self.x_position = x_position
        self.y_position = y_position
        self.status = False
        #False status mean that ligand is not connected with any integrin


class Integrin():
    
    def __init__(self, cell_id, radius, theta, integrin_id):
        self.cell_id = cell_id
        self.radius = radius
        self.theta = theta
        self.integrin_id = integrin_id
        self.status = False
    
    def getInformation(self):
        print(''.join(("integrin ", str(self.cell_id),".", 
            str(self.integrin_id), ": (", str(self.radius), ",", str(self.theta), ")")))
    
    def getXYPosition(self, cell):
        x_position = self.radius*np.cos(np.deg2rad(self.theta)) + cell.x_center
        y_position = self.radius*np.cos(np.deg2rad(self.theta)) + cell.y_center
        return(x_position, y_position)
    
    def getLigandDistance(self, cell, target_ligand):
        integrin_x_position, integrin_y_position = self.getXYPosition(cell)
        x_distance = integrin_x_position - target_ligand.x_position
        y_distance = integrin_y_position - target_ligand.y_position
        return (sqrt(x_distance**2 + y_distance**2))
    
    def getIntegrinDistance(self, self_cell, target_integrin, target_cell):
        self_x_position, self_y_position = self.getXYPosition(self_cell)
        target_x_position, target_y_position = target_integrin.getXYPosition()
        x_distance = self_x_position - target_x_position
        y_distance = self_y_position - target_y_position
        return(sqrt(x_distance**2 + y_distance**2))

class Cell():
    cell_number = 0
       
    def __init__(self, x_center, y_center, max_radius, total_integrin):
        self.x_center = x_center
        self.y_center = y_center
        self.radius = max_radius
        self.integrin = []
        Cell.cell_number += 1
        for i in range(total_integrin):
            radius = round(uniform(0, max_radius), 2)
            theta = round(uniform(0,360), 2)
            self.integrin.append(Integrin(Cell.cell_number, radius, theta, i+1))
    
    def getIntegrinList(self):
        for reseptor in self.integrin:
            reseptor.getInformation()
    
    

    


