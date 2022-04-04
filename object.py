# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:13:52 2022

@author: Asus
"""
import numpy as np

class Ligand():
    
    def __init__(self):
        pass

class Integrin():
    
    def __init__(self, cell_id, radius, theta, integrin_id):
        self.cellId = cell_id
        self.radius = radius
        self.theta = theta
        self.integrinId = integrin_id

class Cell():
    cell_number = 0
    
    def ___init__(self, x_center, y_center, max_radius, total_integrin):
        self.integrin = np.zeros(total_integrin)
        for i in range(total_integrin):
            
            self.integrin[i] = Integrin(Cell.cell_number, radius, theta, i)

