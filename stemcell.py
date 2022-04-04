#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from random import uniform
import numpy as np

class Ligand():
    
    def __init__(self):
        pass

class Integrin():
    
    def __init__(self, cell_id, radius, theta, integrin_id):
        self.cell_id = cell_id
        self.radius = radius
        self.theta = theta
        self.integrin_id = integrin_id
    
    def getIntegrinInformation(self):
        print(''.join(("integrin ", str(self.cell_id),".", 
            str(self.integrin_id), ": (", str(self.radius), ",", str(self.theta), ")")))

class Cell():
       
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
            reseptor.getIntegrinInformation()
    
    cell_number = 0

    


