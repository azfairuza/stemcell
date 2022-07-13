#Author: azfairuza
#email: fairuza.zacky1@gmail.com

from __future__ import annotations
from encodings import utf_8
from random import uniform
from typing import Dict, List
from pathlib import Path
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from datetime import datetime
from descartes import PolygonPatch
import alphashape
import shutil

class Ligand():
    '''
    Class for simulate ligand in nanopattern.
    '''

    # Class properties
    ligand_number = 0     #total ligand created
    
    def __init__(self, x_position: float, y_position: float, dot_size):
        '''
        Initial procedure when creating a ligand.
        '''
        
        Ligand.ligand_number += 1                # Updating ligand number

        # ligand properties
        self.x_position = x_position            # ligand x-position
        self.y_position = y_position            # ligand y-position
        self.bound = False                      # False: ligand is not bounded to any integrin
        self.targeted = False                   # False: ligand is not targeted by any integrin
        self.id = Ligand.ligand_number          # id number for ligand
        self.target_id = []                     # id number for integrin connected to the ligand
        self.type_class = 'ligand'              # class type
        self.ligand_size = dot_size

    @classmethod
    def resetNumber(cls):
        '''
        Procedure to reset ligand number.
        '''
        cls.ligand_number = 0
        
class Nanopattern():
    '''
    Class for simulate nanopattern.
    '''

    def __init__(self, height: float, width: float, grid_height: float, grid_width: float, 
        ligand_position: list , dot_size: float):
        '''
        Inital procedure when creating a nanopattern.
        '''
        
        # nanopattern properties
        self.height = height                                 # nanopattern height
        self.width = width                                   # nanopattern width
        self.grid_height = grid_height                       # nanopattern grid height
        self.grid_width = grid_width                         # nanopattern grid width
        self.position_seed = ligand_position                 # ligand position's list
        self.dot_size = dot_size                             # ligand size (dot size in simulation)
        self.row_number = int(np.ceil(height/grid_height))   # number of grid row created
        self.col_number = int(np.ceil(width/grid_width))     # number of grid column created
        self.type_class = 'nanopattern'
        self.ligands: List[Ligand] = []

        # build nanopattern
        for row in range(self.row_number):
            for col in range(self.col_number):
                for position in ligand_position:
                    self.ligands.append(Ligand((position[0] + col)*grid_width, 
                        (position[1] + row)*grid_height, dot_size))

        print('nanopattern has been created')
    
    def getNotTargetedLigand(self) -> List[Ligand]:
        '''
        Procedure to get a list of ligand that is not targeted by integrin.
        '''
        members = []
        for ligand in self.ligands:
            if ligand.bound == False:       
                members.append(ligand)
        return members
    
    def getXPositionLigand(self) -> List[float]:
        '''
        Procedure to get x position of ligands in a list.
        '''
        position = []
        for ligand in self.ligands:
            position.append(ligand.x_position)
        return position
    
    def getYPositionLigand(self) -> List[float]:
        '''
        Procedure to get y position of ligand in a list.
        '''
        position = []
        for ligand in self.ligands:
            position.append(ligand.y_position)
        return position
    
    def getLigandById(self, id: int) -> Ligand:
        '''
        Procedure to get a ligand by id.
        '''
        for ligand in self.ligands:
            if ligand.id == id:
                return ligand

    def show(self, time: datetime, save=False, number=0, folder=None):
        '''
        Procedure to draw only nanopattern.
        '''
        # determine the name folder
        if folder is None:
            namefolder = f'./output/{getTime(time)}-output/figure'
        else:
            namefolder = f'./output/{getTime(time)}-output/figure/{folder}'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the file name
        namefile = f'{namefolder}/{number:06}.jpg'
        # plot the nanopattern
        fig = plt.figure(figsize=(20, 20))
        fig.dpi = 200
        ax = plt.subplot(aspect='equal')
        x_position = self.getXPositionLigand()
        y_position = self.getYPositionLigand()
        out = circles(x_position, y_position, self.dot_size, 'green', alpha=0.5, ec='none')
        plt.xlim(0, self.width)
        plt.ylim(0, self.height)
        # save the figure if necessary
        if save is True:
            print(f'figure saved on {namefolder}')
            fig.savefig(namefile, bbox_inches='tight', dpi=200)
            plt.close()
    
class Integrin():
    '''
    Class to simulate integrin.
    '''

    # class properties
    integrin_number = 0
    
    def __init__(self, cell_id: int, x_center: float, y_center: float, max_radius: float, 
        integrin_size: int, lst_integrin: List[Integrin]):
        '''
        Initial procedure when creating an integrin.
        '''
        
        Integrin.integrin_number += 1               # Update integrin number

        # identity of integrin
        self.cell_id = cell_id                      # cell id number
        self.id = Integrin.integrin_number          # integrin id number 
        self.type_class = 'integrin'                # object type: integrin
        # position property
        self.x_position = 0                         # x-position of the integrin
        self.y_position = 0                         # y-position of the integrin
        self.surface = False                        # Surface status of the integrin
        # bound properties
        self.bound = False                          # False: integrin is not bounded to anyone
        self.targeting = False                      # False: integrin is not targeting to anyone
        self.target_type = None                     # None ; Ligand ; Integrin
        self.cell_target_id = None                  # cell id number for integrin-integrin complex
        self.ligand_target_id = None                # ligand id number for integrin-ligand complex 
        self.integrin_target_id = None              # integrin id number for integrin-integrin complex
        self.x_target = None                        # x-position of the targeted object
        self.y_target = None                        # y-position of the targeted object
        self.mass = 0                               # integrin mass
        self.integrin_size = integrin_size          # integrin_size
        
        # generate position
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

        # determine surface status
        if radius >= 0*max_radius:
            self.surface = True        # definition: integrin is on the cell surface
        else:
            self.surface = False       # definition: integrin is in the inside of the cell

    def getInformation(self) -> string:
        '''
        Procedure to get integrin information.
        ''' 
        
        if self.bound == True:
            if self.target_type == 'ligand':
                return(f'integrin {self.cell_id}.{self.id}: ({self.x_position},{self.y_position}) status bound to ligand {self.ligand_target_id} at ({self.x_target},{self.y_target})')
            
            elif self.target_type == 'integrin':
                return(f'integrin {self.cell_id}.{self.id}: ({self.x_position},{self.y_position}) status bound to integrin {self.cell_target_id}.{self.integrin_target_id} at ({self.x_target},{self.y_target})')
        else:
            if self.target_type == None:
                return(f'integrin {self.cell_id}.{self.id}: ({self.x_position},{self.y_position}) status not targeting)')
            
            elif self.target_type == 'ligand':
                return(f'integrin {self.cell_id}.{self.id}: ({self.x_position},{self.y_position}) status targeting ligand {self.ligand_target_id} at ({self.x_target},{self.y_target})')
            
            elif self.target_type == 'integrin':
                return(f'integrin {self.cell_id}.{self.id}: ({self.x_position},{self.y_position}) status targeting integrin {self.cell_target_id}.{self.integrin_target_id} at ({self.x_target},{self.y_target})')
    
    def getTargetDistance(self) -> float:
        '''
        Procedure to get the distance between the integrin and the target.
        '''
        return np.sqrt((self.x_position-self.x_target)**2 + (self.y_position-self.y_target)**2)

    def getLigandDistance(self, ligand: Ligand) -> float:
        '''
        Procedure to get distance between the integrin and a ligand.
        '''
        return np.sqrt((self.x_position - ligand.x_position)**2 + (self.y_position - ligand.y_position)**2)
    
    def getIntegrinDistance(self, integrin: Integrin) -> float:
        '''
        Procedure to get distance between the integrin and an integrin.
        '''
        return np.sqrt((self.x_position - integrin.x_position)**2 + (self.y_position - integrin.y_position)**2)
    
    def targetingProcedure1(self, cells: List[Cell], substrate: Nanopattern) -> Integrin | Ligand | None:
        '''
        Procedure to get target object for the ligand using first method. The first method is a method
        that does not allow an object to be multi-targeted, or targeted by many objects. This method 
        will limit the number of target available. 
        '''
        min_distance = 0
        target = None
        if self.targeting == False:
                #check wether the integrin is in surface and also there are more than 1 cell
                if self.surface and len(cells) > 1:
                    # get untargeted integrin from other cells
                    cell_targets = excludeCellById(cells, self.cell_id)
                    for cell_target in cell_targets:
                        for integrin_target in cell_target.integrins:
                            # check whether the integrin has target
                            if integrin_target.targeting == False:
                                if target == None:
                                    target = integrin_target
                                    # update min_distance
                                    min_distance = self.getIntegrinDistance(integrin_target)
                                else:
                                    dstnce = self.getIntegrinDistance(integrin_target)
                                    if dstnce < min_distance:
                                        target = integrin_target
                                        # update min_distance
                                        min_distance = dstnce
                    # get untargeted ligand from nanopattern
                    for ligand in substrate.ligands:
                        # check whether the ligand has targeted
                        #print('checking ligands for surface integrin')
                        if ligand.targeted == False:
                            if target == None:
                                target = ligand
                                # update min_distance
                                min_distance = self.getLigandDistance(ligand)
                            else:
                                dstnce = self.getLigandDistance(ligand)
                                if dstnce < min_distance:
                                    target = ligand
                                    # update min_distance
                                    min_distance = dstnce
                else:
                    # get untargeted ligand from nanopattern
                    for ligand in substrate.ligands:
                        # check whether the ligand has targeted
                        if ligand.targeted == False:
                            if target == None:
                                target = ligand
                                # update min_distance
                                min_distance = self.getLigandDistance(ligand)
                            else:
                                dstnce = self.getLigandDistance(ligand)
                                if dstnce < min_distance:
                                    target = ligand
                                    # update min_distance
                                    min_distance = dstnce
        return target
    
    def updatingProcedure1(self, target: Ligand | Integrin):
        '''
        Procedure to update the attribute of integrin after it get the target. The update is
        happening on the target as well.
        '''
        if target is not None:
            if target.type_class == 'ligand':
                # integrin
                self.targeting = True
                self.target_type = 'ligand'                       
                self.ligand_target_id = target.id
                self.x_target = target.x_position
                self.y_target = target.y_position

                # ligand target 
                target.targeted = True
                target.target_id.append((self.cell_id, self.id))

            elif target.type_class == 'integrin':
                # integrin
                self.targeting = True
                self.target_type = 'integrin'
                self.cell_target_id = target.cell_id
                self.integrin_target_id = target.id
                self.x_target = target.x_position
                self.y_target = target.y_position

                # integrin target
                target.targeting = True
                target.target_type = 'integrin'
                target.cell_target_id = self.cell_id
                target.integrin_target_id = self.id
                target.x_target = self.x_position
                target.y_target = self.y_position
            else:
                pass
        else:
            pass
    
    def movingProcedure1(self, cells: List[Cell], substrate: Nanopattern, dstlimit: float, movespeed=1):
        '''
        Procedure to move the integrin towards the target. If the distance is out of limit
        the integrin can not be moved. Also, if the distance is close enough, binding can be 
        formed between the integrin and the object.
        '''
        if self.targeting:
            # check if the integrin is not bounded and the distance is below the limit
            target_dst = self.getTargetDistance()
            cell = getCellbyId(cells, self.cell_id)
            if (self.bound == False) and (dstlimit > target_dst):
                # check if the integrin can bounded
                #print('integrin is not bound and below the distance limit')
                if target_dst < 2*cell.integrin_size:
                    # update the bound status and mass of the integrin
                    self.bound = True
                    self.mass = 1
                    #print('integrin is bound')
                    if self.target_type == 'ligand':
                        # if the integrin bound to ligand
                        target = substrate.getLigandById(self.ligand_target_id)
                        target.bound = True
                        target.target_id = []
                        target.target_id.append((self.cell_id, self.id))
                    elif self.target_type == 'integrin':
                        # if the integrin bound to another integrin
                        cell_target = getCellbyId(cells, self.cell_target_id)
                        integrin_target = cell_target.getIntegrinById(self.integrin_target_id)
                        integrin_target.bound = True
                        integrin_target.mass = 1
                else:
                    #print('integrin moving')
                    self.updateTarget(cells)
                    self.move(movespeed)
            else:
                pass
        else:
            pass
    
    def targetingProcedure2(self, cells: List[Cell], substrate: Nanopattern) -> Integrin | Ligand | None:
        '''
        Procedure to get target object for the ligand using second method. The second method is a method
        that does allow an object to be multi-targeted, or targeted by many objects. This method 
        will ensure that the target is the closest object possible. 
        '''
        min_distance = 0
        target = None
        if self.bound == False:
                #check wether the integrin is in surface and also there are more than 1 cell
                if self.surface and len(cells) > 1:
                    # get untargeted integrin from other cells
                    cell_targets = excludeCellById(cells, self.cell_id)
                    for cell_target in cell_targets:
                        for integrin_target in cell_target.integrins:
                            # check whether the integrin has target
                            if integrin_target.bound == False:
                                if target == None:
                                    target = integrin_target
                                    # update min_distance
                                    min_distance = self.getIntegrinDistance(integrin_target)
                                else:
                                    dstnce = self.getIntegrinDistance(integrin_target)
                                    if dstnce < min_distance:
                                        target = integrin_target
                                        # update min_distance
                                        min_distance = dstnce
                    # get untargeted ligand from nanopattern
                    for ligand in substrate.ligands:
                        # check whether the ligand has targeted
                        #print('checking ligands for surface integrin')
                        if ligand.bound == False:
                            if target == None:
                                target = ligand
                                # update min_distance
                                min_distance = self.getLigandDistance(ligand)
                            else:
                                dstnce = self.getLigandDistance(ligand)
                                if dstnce < min_distance:
                                    target = ligand
                                    # update min_distance
                                    min_distance = dstnce
                else:
                    # get untargeted ligand from nanopattern
                    for ligand in substrate.ligands:
                        # check whether the ligand has targeted
                        if ligand.bound == False:
                            if target == None:
                                target = ligand
                                # update min_distance
                                min_distance = self.getLigandDistance(ligand)
                            else:
                                dstnce = self.getLigandDistance(ligand)
                                if dstnce < min_distance:
                                    target = ligand
                                    # update min_distance
                                    min_distance = dstnce
        return target
    
    def updatingProcedure2(self, target: Ligand | Integrin):
        '''
        Procedure to update the attribute of integrin after it get the target. The update is
        happening on the target as well.
        '''
        if target is not None:
            if target.type_class == 'ligand':
                # integrin
                self.targeting = True
                self.target_type = 'ligand'                       
                self.ligand_target_id = target.id
                self.x_target = target.x_position
                self.y_target = target.y_position

                # ligand target 
                target.targeted = True
                target.target_id.append((self.cell_id, self.id))

                #check bound
                target_dst = self.getTargetDistance()
                if target_dst < (self.integrin_size + target.ligand_size):
                    self.bound = True
                    target.bound = True
                    target.target_id = []
                    target.target_id.append((self.cell_id, self.id))


            elif target.type_class == 'integrin':
                # integrin
                self.targeting = True
                self.target_type = 'integrin'
                self.cell_target_id = target.cell_id
                self.integrin_target_id = target.id
                self.x_target = target.x_position
                self.y_target = target.y_position

                #check bound
                target_dst = self.getTargetDistance()
                if target_dst < (self.integrin_size + target.integrin_size):
                    self.bound = True
                    target.bound = True
                    target.target_type = 'integrin'
                    target.cell_target_id = self.cell_id
                    target.integrin_target_id = self.id
                    target.x_target = self.x_position
                    target.y_target = self.y_position
                    self.mass = 1

            else:
                pass
        else:
            pass

    
    def movingProcedure2(self, cells: List[Cell], substrate: Nanopattern, dstlimit: float, movespeed=1):
        '''
        Procedure to move the integrin towards the target. If the distance is out of limit
        the integrin can not be moved. Also, if the distance is close enough, binding can be 
        formed between the integrin and the object. This moving object using mark II method
        '''
        
        # check if the integrin is not bounded and the distance is below the limit
        target_dst = self.getTargetDistance()
        cell = getCellbyId(cells, self.cell_id)
        if dstlimit > target_dst:
            # check if the integrin can bounded
            #print('integrin is not bound and below the distance limit')
            if target_dst < 2*cell.integrin_size:
                # update the bound status and mass of the integrin
                self.bound = True
                #print('integrin is bound')
                if self.target_type == 'ligand':
                    # if the integrin bound to ligand
                    target = substrate.getLigandById(self.ligand_target_id)
                    target.bound = True
                    target.target_id = []
                    target.target_id.append((self.cell_id, self.id))
                elif self.target_type == 'integrin':
                    # if the integrin bound to another integrin
                    cell_target = getCellbyId(cells, self.cell_target_id)
                    integrin_target = cell_target.getIntegrinById(self.integrin_target_id)
                    integrin_target.bound = True
                    integrin_target.mass = 1
                    if integrin_target.target_type == 'ligand':
                        integrin_target.ligand_target_id = None
                    integrin_target.cell_target_id = self.cell_id
                    integrin_target.integrin_target_id = self.id
                    self.mass = 1
            else:
                #print('integrin moving')
                self.updateTarget(cells)
                self.move(movespeed)
        else:
            pass
    def validateTarget2(self, cells: List[Cell], substrate: Nanopattern):
        if self.target_type == 'ligand':
            ligand_target = substrate.getLigandById(self.ligand_target_id)
            if ligand_target.bound == True:
                self.targeting = False
                self.target_type = None
                self.ligand_target_id = None
        elif self.target_type == 'integrin':
            cell_target = getCellbyId(cells, self.cell_target_id)
            integrin_target = cell_target.getIntegrinById(self.integrin_target_id)
            if integrin_target.bound == True:
                self.targeting = False 
                self.target_type = None
                self.cell_target_id = None
                self.integrin_target_id = None

    def move(self, movespeed=1):
        '''
        Procedure to update the position of integrin towards target.
        '''
        # find the angle
        # do trigonometri on the movespeed
        # update the x pos and y pos
        angle = np.arctan2((self.y_target - self.y_position), (self.x_target - self.x_position))
        dx = movespeed*np.cos(angle)
        dy = movespeed*np.sin(angle)
        self.x_position += dx
        self.y_position += dy
    
    def updateTarget(self, cells: List[Cell]):
        '''
        Procedure update the target position
        '''
        if self.target_type == 'integrin':
            cell_target = getCellbyId(cells, self.cell_target_id)
            target = cell_target.getIntegrinById(self.integrin_target_id)
            self.x_target = target.x_position
            self.y_target = target.y_position
        else: 
            pass
    
    @classmethod
    def resetNumber(cls):
        '''
        Procedure to reset Integrin number.
        '''
        cls.integrin_number = 0

class Cell():
    '''
    class to simulate cell.
    '''
    cell_number = 0       # number of cell created
       
    def __init__(self, x_center: float, y_center: float, max_radius: float, total_integrin: int, 
        cell_mass: float, integrin_size: int):
        Cell.cell_number += 1   
        
        #cell properties
        self.mass = cell_mass                   # mass of the cell (in center of mass)
        self.x_center = x_center                # x position of the center of mass of the cell
        self.y_center = y_center                # y position of the center of mass of the cell
        self.integrin_size = integrin_size      # the size of integrin
        self.radius = max_radius                # the outer radius of the cell
        self.id = Cell.cell_number              # cell_id
        self.type_class = 'cell'                # class type
        self.integrins: List[Integrin] = []      # list of integrin
        self.alpha_shape = None


        #cell integrin members       
        max_number = int(max_radius**2/(8*(integrin_size**2)))    
        if int(total_integrin) <= max_number:                                                                                                           # 
            for i in range(int(total_integrin)):      
                build_integrin = Integrin(Cell.cell_number, x_center, y_center, 
                    max_radius, integrin_size, self.integrins)                                                   
                self.integrins.append(build_integrin)                                 # build the integrin
        else:
            print(f'maximm number integrin allowed is {max_number}')

    def getIntegrinInfo(self):
        '''
        Get information of all integrins in the cell.
        '''
        for integrin in self.integrins:
            integrin.getInformation()

    def getXPositionIntegrin(self) -> List[float]:
        '''
        Procedure to get x position of integrins, contained in a list.
        '''
        position = []
        for integrin in self.integrins:
            position.append(integrin.x_position)
        return position
    
    def getYPositionIntegrin(self) -> List[float]:
        '''
        Procedure to get y position of integrins, contained in a list.
        '''
        position = []
        for integrin in self.integrins:
            position.append(integrin.y_position)
        return position
    
    def getIntegrinById(self, id: int) -> Integrin:
        '''
        Procedure to get integrin by Id.
        '''
        for integrin in self.integrins:
            if integrin.id == id:
                return integrin
    
    def getAvailableIntegrin(self) -> List[Integrin]:
        '''
        Procedure to get integrin that has not found a target
        '''
        available_integrin: List[Integrin] = []
        for integrin in self.integrins:
            if integrin.targeting == False:
                available_integrin.append(integrin)
        return available_integrin

    def show(self, time: datetime, substrate: Nanopattern, alphaValue=0, save=False, number=0, folder=None):
        '''
        Procedure to draw only nanopattern
        '''        
        # determine the folder name
        if folder is None:
            namefolder = f'./output/{getTime(time)}-output/figure'
        else:
            namefolder = f'./output/{getTime(time)}-output/figure/{folder}'
        # build the folder
        Path(namefolder).mkdir(parents=True, exist_ok=True)
        # determine the name of the file
        namefile = f'{namefolder}/{number:06}.jpg'
        # draw the nanopattern
        fig = plt.figure(figsize=(20, 20))
        fig.dpi = 200
        ax = plt.subplot(aspect='equal')
        x_position = self.getXPositionIntegrin()                                        
        y_position = self.getYPositionIntegrin()
        listpoint =[]
        for i in range(len(x_position)):
            point = [round(x_position[i], 2), round(y_position[i],2)]
            listpoint.append(point)
        self.alpha_shape = alphashape.alphashape(listpoint, alphaValue)
        plt.gca().add_patch(PolygonPatch(self.alpha_shape, alpha=0.2))
        # points =np.array(listpoint)
        # hull = ConvexHull(points)
        # for simplex in hull.simplices:
        #     plt.plot(points[simplex, 0], points[simplex, 1], linestyle='--', color='k', linewidth='2')
        # surface_circle = plt.Circle((self.x_center, self.y_center), 
        #     self.radius, color='black', fill=False, ls='--')
        # plt.gca().add_patch(surface_circle)
        out = circles(x_position, y_position, self.integrin_size, 'red', alpha=0.5, ec='none')
        plt.xlim(0, substrate.width)
        plt.ylim(0, substrate.height)
        # save figure if necessary
        if save is True:
            print(f'figure saved on {namefolder}')
            fig.savefig(namefile, bbox_inches='tight', dpi=200)
            plt.close()
    
    def getCenterofMass(self):
        '''
        procedure to get the center of mass of the cell
        '''
        # initiation
        total_mass = 0
        total_x_mass = 0
        total_y_mass = 0
        # calculate the mass of integrins
        for integrin in self.integrins:
            if integrin.bound == True:
                total_mass += integrin.mass
                total_x_mass += (integrin.mass*integrin.x_position)
                total_y_mass += (integrin.mass*integrin.y_position)
            else:
                pass
        # calculate the mass of cell
        total_mass += self.mass
        total_x_mass += (self.mass*self.x_center)
        total_y_mass += (self.mass*self.y_center)
        # get the center of mass
        x_cm = round(total_x_mass/total_mass, 2)
        y_cm = round(total_y_mass/total_mass, 2)
        return (x_cm, y_cm)
    
    def totalIntegrinBound(self):
        '''
        Procedure to get the total of integrin bound.
        '''
        total_integrin = 0
        for integrin in self.integrins:
            if integrin.bound == True:
                total_integrin += 1
        return total_integrin
    
    def totalMass(self):
        '''
        Procedure to get total mass of the cell.
        '''
        total_mass = self.mass
        for integrin in self.integrins:
            total_mass += integrin.mass
        return total_mass

    @classmethod
    def resetNumber(cls):
        '''
        Procedure to reset cell number.
        '''
        cls.cell_number = 0

def readFile(filename: string):
    '''
    Procedure to read input file. There are three input file that can be read. 
    PATCON: contains the configuration for building nanopattern.
    CELCON: contains the configuration for building cells.
    SIMCON: contains the configuration for running the simulation.
    '''
    if filename == 'PATCON':
    # open PATCON file
        try:
            data_file = open('PATCON.txt', 'r', encoding='utf-8')    
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)           # filter out the empty string and '\n'              
        except:
            stripped_strng = 'Error in opening PATCON file'     # error flag if PATCON can not be opened
        finally:
            return stripped_strng

    elif filename == 'CELCON':
    # open CELCON file
        try:
            data_file = open('CELCON.txt', 'r', encoding = 'UTF-8')
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)           # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening CELCON file'     # error flag if CELCON can not be opened
        finally:
            return stripped_strng

    elif filename == 'SIMCON':
    # open SIMCON file
        try:
            data_file = open('SIMCON.txt', 'r', encoding='utf-8')
            lst_strng = data_file.readlines()
            data_file.close()
            stripped_strng = filterElement(lst_strng)           # filter out the empty string and '\n' 
        except:
            stripped_strng = 'Error in opening SIMCON file'     # error flag if SIMCON can not be opened
        finally:
            return stripped_strng

    else:
        print('File name is not correct')
        return 0

def getValue(lst_strng: list, property_name: string):
    '''
    Procedure to get the value of a property in XXXCON file. There are special
    property name to be used in this function.
    PATCON: property name = 'LIGAND'.
    CELCON: property name = 'CELL'.
    SIMCON: property name = 'METADATA'.
    '''
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
    elif property_name == 'METADATA':
        start_index = (lst_strng.index('#METADATA')+1)                  # start index to search
        end_index = lst_strng.index('#CONFIG')                          # last index to search
        cell_properties = {'username': None, 'title': None, 'date': None }                                            # cell properties container
        for i in range(start_index, end_index):
            stringdata = lst_strng[i].split()                           
            if stringdata[0] == 'username':                             # user of the simulation
                cell_properties['username'] = ' '.join(stringdata[1:])
            elif stringdata[0] == 'title':                              # the title of the simulation
                cell_properties['title'] = ' '.join(stringdata[1:])
            elif stringdata[0] == 'date':                               # simulation date
                cell_properties['date'] = stringdata[1]
        return cell_properties
    else: 
        start_index = (lst_strng.index('#CONFIG') + 1)
        end_index =  lst_strng.index('#END')                                
        value = None 
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
    '''
    procedure to remove empty element in a list and remove the '\n' character.
    '''
    new_lst = []
    for element in input_lst:
        new_element = element.rstrip()          # remove the '\n'    
        if new_element:                         # check if the element is not empty
            new_lst.append(new_element)
    return new_lst

def saveInput(time: datetime):
    namefolder = f'./output/{getTime(time)}-output/input'
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    shutil.copy2('./PATCON.txt', namefolder)
    shutil.copy2('./CELCON.txt', namefolder)
    shutil.copy2('./SIMCON.txt', namefolder)
    print('input file has been copied!')

def getCellbyId(cells: List[Cell], id):
    '''
    Procedure to get cell by id.
    '''
    for cell in cells:
        if cell.id == id:
            return cell

def showAll(cells : List[Cell], substrate: Nanopattern, time: datetime,
    show_substrate=False, save=False, number=0, folder=None, line=False, alphaValue=0):
    '''
    Procedure to show all element of simulation including cells and nanopattern.
    '''
    # determine the folder's name
    if folder is None:
        namefolder = f'./output/{getTime(time)}-output/figure'
    else:
        namefolder = f'./output/{getTime(time)}-output/figure/{folder}'
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    # determine the file name        
    namefile = f'{namefolder}/{number:06}.jpg'
    # draw the figure
    fig = plt.figure(figsize=(20,20))
    fig.dpi = 200
    ax = plt.subplot(aspect='equal')
    # print the cells
    for cell in cells:
        # draw line between cells
        if line == True:
            for integrin in cell.integrins:
                if integrin.targeting:
                    x = [integrin.x_position, integrin.x_target]
                    y = [integrin.y_position, integrin.y_target]
                    plt.plot(x, y, color="black", linewidth="1", linestyle="dashed")
        x_position = cell.getXPositionIntegrin()
        y_position = cell.getYPositionIntegrin()
        listpoint =[]
        for i in range(len(x_position)):
            point = [round(x_position[i], 2), round(y_position[i],2)]
            listpoint.append(point)
        cell.alpha_shape = alphashape.alphashape(listpoint, alphaValue)
        # print(f"cell {cell.id} area: {cell.alpha_shape.area}")
        plt.gca().add_patch(PolygonPatch(cell.alpha_shape, alpha=0.2))
        out = circles(x_position, y_position, cell.integrin_size, 'red', alpha=0.5, ec='none')
        # if edge == 'ConvexHull':
        #     points =np.array(listpoint)
        #     hull = ConvexHull(points)
        #     for simplex in hull.simplices:
        #         plt.plot(points[simplex, 0], points[simplex, 1], linestyle='--', color='k', linewidth='2')
        # elif edge == 'Circle':
        #     surface_circle = plt.Circle((cell.x_center, cell.y_center), 
        #         cell.radius, color='black', fill=False, ls='--')
        #     plt.gca().add_patch(surface_circle)
        # elif edge == 'AlphaShape':
        #     alpha_shape = alphashape.alphashape(listpoint, 0.01)
        #     plt.gca().add_patch(PolygonPatch(alpha_shape, alpha=0.2))
        # out = circles(x_position, y_position, cell.integrin_size, 'red', alpha=0.5, ec='none')
    # print the nanopattern if the option is true
    if show_substrate is True:
        x_position = substrate.getXPositionLigand()
        y_position = substrate.getYPositionLigand()
        out = circles(x_position, y_position, substrate.dot_size, 'green', alpha=0.5, ec = 'none')
    plt.xlim(0, substrate.width)
    plt.ylim(0, substrate.height)
    plt.title(f"iteration = {number:06}")
    # save if necessary
    if save is True:
        print(f'figure {number:06}.jpg saved on {namefolder}')
        fig.savefig(namefile, bbox_inches='tight', dpi=200)
        plt.close()

def excludeCellById(cells: List[Cell], cell_id: int):
    '''
    Procedure to exclude a cell from list of cells by Id.
    '''
    new_cells: List[Cell] = []
    for cell in cells:
        if cell.id != cell_id:
            new_cells.append(cell)
    return new_cells

def saveCenterOfMass(cells: List[Cell], num_iteration, time: datetime):
    
    namefolder = f'./output/{getTime(time)}-output/file'
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    # determine the name of the file
    namefile = f'{namefolder}/CELLCM.txt'
    if num_iteration <= 0:
        head_text = 't\t'
        for cell in cells:
            head_cell = f'x{cell.id}\ty{cell.id}\tm{cell.id}\tn{cell.id}\t'
            head_text += head_cell
        head_text += '\n'
        # save the data
        with open(namefile, 'w') as output:
            output.write(head_text)
    output_text = f'{num_iteration}\t'
    for cell in cells:
        cm_position = cell.getCenterofMass()
        mass = cell.totalMass()
        integrin_bound = cell.totalIntegrinBound()
        cell_output = f'{cm_position[0]}\t{cm_position[1]}\t{mass}\t{integrin_bound}\t'
        output_text += cell_output
    output_text += '\n'
    # save the data
    with open(namefile, 'a') as output:
        output.write(output_text)

def simulate1(cells: List[Cell], substrate: Nanopattern, time: datetime):
    '''
    Procedure to do simulation with single targeting or first method.
    '''   

    simcon = readFile('SIMCON')

    # get value
    n_iteration = int(getValue(simcon, 'iteration'))
    print(f'iter {n_iteration}')
    dstlimit = getValue(simcon, 'dstlimit')
    print(f'dstlimit {dstlimit}')
    movespeed = getValue(simcon, 'movespeed')
    line = getValue(simcon, 'line')
    savefig = getValue(simcon, 'savefig')
    centerofmass = getValue(simcon, 'centerofmass')
    cellarea = getValue(simcon, 'cellarea')
    alphaValue = getValue(simcon, 'alpha')

    if n_iteration == None:
        n_iteration = 1
    if dstlimit == None:
        dstlimit = 1
    if movespeed == None:
        movespeed = 1
    if line == None:
        line = False
    if savefig == None or savefig == 0:
        savefig = False
    else:
        savefig = True

    iter_simulation = 0
    # start the simulation
    # Targeting
    # pick a cell
    for cell in cells:
        #pick an integrin from the cell
        for integrin in cell.integrins:
            # get the target
            target = integrin.targetingProcedure1(cells, substrate)
            # update the data
            integrin.updatingProcedure1(target)
    # Moving
    while(iter_simulation <= n_iteration): 
        # pick a cell
        for cell in cells:
            # pick an integrin from the cell
            for integrin in cell.integrins:
                # move based on dstlimit, bound, and target
                integrin.movingProcedure1(cells, substrate, dstlimit, movespeed)
        showAll(cells, substrate, time,
            show_substrate=True, 
            save=savefig, 
            folder='simulation_output', 
            number=iter_simulation, 
            line=line, 
            alphaValue=alphaValue)
        
        if centerofmass == 1:
            saveCenterOfMass(cells, iter_simulation, time)
        if cellarea == 1:
            saveArea(cells, iter_simulation, time)
    
        iter_simulation += 1

def simulate2(cells: List[Cell], substrate: Nanopattern, time: datetime):
    '''
    Procedure to do simulation with multi targeting method
    '''   

    simcon = readFile('SIMCON')

    # get value
    n_iteration = int(getValue(simcon, 'iteration'))
    print(f'iter {n_iteration}')
    dstlimit = getValue(simcon, 'dstlimit')
    print(f'dstlimit {dstlimit}')
    movespeed = getValue(simcon, 'movespeed')
    line = getValue(simcon, 'line')
    savefig = getValue(simcon, 'savefig')
    centerofmass = getValue(simcon, 'centerofmass')
    cellarea = getValue(simcon, 'cellarea')
    alphaValue = getValue(simcon, 'alpha')

    if n_iteration == None:
        n_iteration = 1
    if dstlimit == None:
        dstlimit = 1
    if movespeed == None:
        movespeed = 1
    if line == None:
        line = False
    if savefig == None or savefig == 0:
        savefig = False
    else:
        savefig = True

    iter_simulation = 0
    # start the simulation
    # Targeting
    # pick a cell
    for cell in cells:
        #pick an integrin from the cell
        for integrin in cell.integrins:
            # get the target
            target = integrin.targetingProcedure2(cells, substrate)
            # update the data
            integrin.updatingProcedure2(target)
    # Moving
    while(iter_simulation <= n_iteration): 
        # pick a cell
        for cell in cells:
            # pick an integrin from the cell
            for integrin in cell.integrins:
                # move based on dstlimit, bound, and target
                if integrin.bound == False:
                    if integrin.targeting:
                        integrin.validateTarget2(cells, substrate)
                        if integrin.targeting == False:
                            target = integrin.targetingProcedure2(cells, substrate)
                            integrin.updatingProcedure2(target)
                        else:
                            integrin.movingProcedure2(cells, substrate, dstlimit, movespeed)
                    else:
                        target = integrin.targetingProcedure2(cells, substrate)
                        integrin.updatingProcedure2(target)
                else:
                    pass
        showAll(cells, substrate, time,
            show_substrate=True, 
            save=savefig, 
            folder='simulation_output', 
            number=iter_simulation, 
            line=line, 
            alphaValue=alphaValue)
        
        if centerofmass == 1:
            saveCenterOfMass(cells, iter_simulation, time)
        if cellarea == 1:
            saveArea(cells, iter_simulation, time)
    
        iter_simulation += 1

def getTime(time: datetime):
    year = time.strftime('%Y')
    month = time.strftime('%m')
    day = time.strftime('%d')
    hour = time.strftime('%H')
    minute = time.strftime('%M')
    return f'{year}-{month}-{day}-{hour}{minute}'

def saveArea(cells: List[Cell], num_iteration, time):
    namefolder = f'./output/{getTime(time)}-output/file'
    # build the folder
    Path(namefolder).mkdir(parents=True, exist_ok=True)
    # determine the name of the file
    namefile = f'{namefolder}/CELLAR.txt'
    if num_iteration <= 0:
        head_text = 't\t'
        for cell in cells:
            head_cell = f'A{cell.id}\t'
            head_text += head_cell
        head_text += '\n'
        # save the data
        with open(namefile, 'w') as output:
            output.write(head_text)
    output_text = f'{num_iteration}\t'
    for cell in cells:
        area = round(cell.alpha_shape.area, 2)
        cell_output = f'{area}\t'
        output_text += cell_output
    output_text += '\n'
    # save the data
    with open(namefile, 'a') as output:
        output.write(output_text)

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    '''
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
    '''
  
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

    
    

    


