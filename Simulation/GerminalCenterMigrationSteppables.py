from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np
import random as rand
import math

#VOCABULARY:
#Ag = antigen

NUM_AG_KEY = "ag" #a dictionary key for the number of antigen a cell holds
LAST_DIV_TIME_KEY = "div_time"

ANTIGEN_VALUE = "ACTG"

DIVIDE_AG_ASYMMETRIC = 0.72
POLARITY_INDEX = 0.88
FDC_ANTIGEN = 3000 #from MerinoTejero2021

#also from MerinoTejero2021
# B-Cell Speed (um / hr)
# 0.18
# B-Cell Persistent Time average (hr.)
# 0.08

#Duration of CC collection of Antigen by serial encounters with FDC (hr.)
#i.e. how long a b cell touches an FDC to receive ag
# 0.7

#one idea: let FDCs (receptor length 8) give antigen from far away with cell links



class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        # self.create_scalar_field_cell_level_py("AFFINITY_FIELD")
        self.num_ag_field = self.create_scalar_field_cell_level_py("NUM_AG_FIELD")
        
        

    def setup_cells(self):
        TFH_CELL_SIZE = 5
        DENDRITIC_CELL_SIZE = 10
        B_CELL_SIZE = 3
        STROMAL_CELL_SIZE = B_CELL_SIZE
        
        for x in range(136, 256-TFH_CELL_SIZE, TFH_CELL_SIZE*3):
            for y in range(TFH_CELL_SIZE, 256-TFH_CELL_SIZE, TFH_CELL_SIZE*3):
                if rand.random() < 0.33:
                    self.cell_field[x:x+TFH_CELL_SIZE, y:y+TFH_CELL_SIZE, 0] = self.new_cell(self.TFH)
        
        for x in range(0, 120-STROMAL_CELL_SIZE, STROMAL_CELL_SIZE*3):
            for y in range(STROMAL_CELL_SIZE, 256-STROMAL_CELL_SIZE, STROMAL_CELL_SIZE*3):
                if rand.random() < 0.4:
                    self.cell_field[x:x+STROMAL_CELL_SIZE, y:y+STROMAL_CELL_SIZE, 0] = self.new_cell(self.STROMAL)
                    
        for x in range(0, 120-B_CELL_SIZE, B_CELL_SIZE*5):
            for y in range(B_CELL_SIZE, 256-B_CELL_SIZE, B_CELL_SIZE*5):
                if rand.random() < 0.5:
                    b_cell = self.new_cell(self.CENTROBLAST)
                    b_cell.dict[NUM_AG_KEY] = FDC_ANTIGEN
                    b_cell.dict[LAST_DIV_TIME_KEY] = -10000
                    self.cell_field[x:x+B_CELL_SIZE, y:y+B_CELL_SIZE, 0] = b_cell

    def start(self):
        rand.seed(667)
        self.setup_cells()

        # for cell in self.cell_list:

            # cell.targetVolume = 25
            # cell.lambdaVolume = 2.0
        
        
            
    def step(self, mcs):
        for cell in self.cell_list_by_type(self.CENTROCYTE):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor and neighbor.type == self.DENDRITIC:
                    #Robert2017: bindprobability ← affinity(BCR, antigen) . ( F.antigenAmount[f] / antigenSaturation)
                    #TODO
                    cell.dict[NUM_AG_KEY] = FDC_ANTIGEN
            
            self.num_ag_field[cell] = rand.random()
        
        
            
        
        
# class GrowthSteppable(SteppableBasePy):
    # def __init__(self,frequency=1):
        # SteppableBasePy.__init__(self, frequency)

    # def step(self, mcs):
    
        # for cell in self.cell_list_by_type(self.CENTROBLAST):
            # cell.targetVolume += 1        

        # # # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        

        # # field = self.field.CHEMICAL_FIELD_NAME
        
        # # for cell in self.cell_list:
            # # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            # # # you can use here any fcn of concentrationAtCOM
            # # cell.targetVolume += 0.01 * concentrationAtCOM       

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        
        for cell in self.cell_list_by_type(self.CENTROBLAST):
            last_div_time = cell.dict[LAST_DIV_TIME_KEY]
            if mcs - last_div_time >= 300:
                cell.dict[LAST_DIV_TIME_KEY] = mcs
                cells_to_divide.append(cell)
        
        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)
    

    def update_attributes(self):
        # reducing parent target volume
        # self.parent_cell.targetVolume /= 2.0                  

        self.clone_parent_2_child()

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        
        total_ag = self.parent_cell.dict[NUM_AG_KEY]
        if rand.random() < DIVIDE_AG_ASYMMETRIC:
            #Asymmetric division: one centroblast loses Ag
            sep = rand.uniform(0, POLARITY_INDEX)
            self.parent_cell.dict[NUM_AG_KEY] = math.floor(total_ag * sep)
            self.child_cell.dict[NUM_AG_KEY] = math.floor(total_ag * (1-sep))
        else:
            #Even split of Ag
            self.parent_cell.dict[NUM_AG_KEY] = total_ag // 2
            self.child_cell.dict[NUM_AG_KEY] = total_ag // 2
            
        
        if self.parent_cell.type==1:
            self.child_cell.type=2
        else:
            self.child_cell.type=1
            
        for cell in [self.parent_cell, self.child_cell]:
            if cell.dict[NUM_AG_KEY] < 100: #about 5 splits from 3000 antigen
                cell.type = self.CENTROCYTE

        
# class DeathSteppable(SteppableBasePy):
    # def __init__(self, frequency=1):
        # SteppableBasePy.__init__(self, frequency)

    # def step(self, mcs):
        # pass
        # # if mcs == 1000:
            # # for cell in self.cell_list:
                # # if cell.type==1:
                    # # cell.targetVolume=0
                    # # cell.lambdaVolume=100

        