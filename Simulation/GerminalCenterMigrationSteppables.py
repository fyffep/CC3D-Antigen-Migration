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
DNA_KEY = "dna"
AFFINITY_KEY = "aff"

ANTIGEN_VALUE = ["A","C","T"]
INITIAL_DNA = ["T","A","T"]
ANTIGEN_SYMBOLS = ["A","C","T","G"]
ANTIGEN_LEN = len(ANTIGEN_VALUE)

DIVIDE_AG_ASYMMETRIC = 0.72
POLARITY_INDEX = 0.88
FDC_ANTIGEN = 3000 #from MerinoTejero2021

affinity_field = None
num_ag_field = None

num_plasma = 0

#also from MerinoTejero2021
# B-Cell Speed (um / hr)
# 0.18
# B-Cell Persistent Time average (hr.)
# 0.08

#Duration of CC collection of Antigen by serial encounters with FDC (hr.)
#i.e. how long a b cell touches an FDC to receive ag
# 0.7

#one idea: let FDCs (receptor length 8) give antigen from far away with cell links

def mutate(dna):
    # print("mutated ",dna,end="  ")
    i = rand.randint(0, ANTIGEN_LEN-1)
    symbol = ANTIGEN_SYMBOLS[rand.randint(0, len(ANTIGEN_SYMBOLS)-1)]
    dna[i] = symbol
    # print("to",dna)
    return dna
    
def judge_affinity(dna):
    ham = 0
    for i in range(ANTIGEN_LEN):
        if dna[i] == ANTIGEN_VALUE[i]:
            ham += 1
    return ham / ANTIGEN_LEN


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        rand.seed(667)
        
        global affinity_field
        global num_ag_field
        affinity_field = self.create_scalar_field_cell_level_py("AFFINITY_FIELD")
        num_ag_field = self.create_scalar_field_cell_level_py("NUM_AG_FIELD")
        
        
        

    def setup_cells(self):
        TFH_CELL_SIZE = 5
        DENDRITIC_CELL_SIZE = 10
        B_CELL_SIZE = 3
        STROMAL_CELL_SIZE = B_CELL_SIZE
        
        for x in range(136, 256-TFH_CELL_SIZE, TFH_CELL_SIZE*3):
            for y in range(TFH_CELL_SIZE, 256-TFH_CELL_SIZE, TFH_CELL_SIZE*3):
                if rand.random() < 0.33:
                    self.cell_field[x:x+TFH_CELL_SIZE, y:y+TFH_CELL_SIZE, 0] = self.new_cell(self.TFH)
                    
        for x in range(136, 256-DENDRITIC_CELL_SIZE, DENDRITIC_CELL_SIZE*2):
            for y in range(DENDRITIC_CELL_SIZE, 256-DENDRITIC_CELL_SIZE, DENDRITIC_CELL_SIZE*2):
                if rand.random() < 0.5:
                    self.cell_field[x:x+DENDRITIC_CELL_SIZE, y:y+DENDRITIC_CELL_SIZE, 0] = self.new_cell(self.DENDRITIC)
        
        # for x in range(0, 120-STROMAL_CELL_SIZE, STROMAL_CELL_SIZE*3):
            # for y in range(STROMAL_CELL_SIZE, 256-STROMAL_CELL_SIZE, STROMAL_CELL_SIZE*3):
                # if rand.random() < 0.4:
                    # self.cell_field[x:x+STROMAL_CELL_SIZE, y:y+STROMAL_CELL_SIZE, 0] = self.new_cell(self.STROMAL)
        
        original_affinity = judge_affinity(INITIAL_DNA)
        for x in range(0, 120-B_CELL_SIZE, B_CELL_SIZE*5):
            for y in range(B_CELL_SIZE, 256-B_CELL_SIZE, B_CELL_SIZE*5):
                if rand.random() < 0.5:
                    b_cell = self.new_cell(self.CENTROBLAST)
                    b_cell.dict[NUM_AG_KEY] = FDC_ANTIGEN
                    b_cell.dict[LAST_DIV_TIME_KEY] = -10000
                    b_cell.dict[DNA_KEY] = INITIAL_DNA
                    b_cell.dict[AFFINITY_KEY] = original_affinity
                    self.cell_field[x:x+B_CELL_SIZE, y:y+B_CELL_SIZE, 0] = b_cell

    def start(self):
        self.setup_cells()

        # # for cell in self.cell_list_by_type(self.CENTROBLAST):
            # # cell.targetVolume = 25
            # # cell.lambdaVolume = 2.0
            
        # # iterating over all cells in simulation        
        # for cell in self.cell_list:
            # # you can access/manipulate cell properties here
            # print("id=", cell.id, " type=", cell.type, " volume=", cell.volume)
        
        
        
        
        
        
        
            
    def step(self, mcs):
        global num_ag_field
        global num_plasma
        # for cell in self.cell_list_by_type(self.CENTROBLAST):
            # num_ag_field[cell] = cell.dict[NUM_AG_KEY]
            
        for cell in self.cell_list_by_type(self.CENTROCYTE):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.DENDRITIC:
                        #Robert2017: bindprobability ← affinity(BCR, antigen) . ( F.antigenAmount[f] / antigenSaturation)
                        #TODO
                        cell.dict[NUM_AG_KEY] = FDC_ANTIGEN
                        # cell.type = 
            
                    if neighbor.type == self.TFH:
                        affinity = cell.dict[AFFINITY_KEY]
                        if affinity >= 1.0:
                            cell.type = self.PLASMA
                            global num_plasma
                            num_plasma += 1
                            print("num_plasma",num_plasma)
                        else:
                            cell.type = self.CENTROBLAST
                            
            
            
        
        
        
        
        
        
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
            if mcs - last_div_time >= 200:
                cell.dict[LAST_DIV_TIME_KEY] = mcs
                if len(self.get_cell_neighbor_data_list(cell)) < 4:
                    cells_to_divide.append(cell)
        
        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)
    

    def update_attributes(self):
        global affinity_field
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
                    
        for cell in [self.parent_cell, self.child_cell]:
            
            #Mutate and update affinity
            new_dna = mutate(cell.dict[DNA_KEY])
            cell.dict[DNA_KEY] = new_dna
            
            affinity = judge_affinity(new_dna)
            cell.dict[AFFINITY_KEY] = affinity
            affinity_field[cell] = affinity
            
            if cell.dict[NUM_AG_KEY] < 100: #about 5 splits from 3000 antigen
                cell.type = self.CENTROCYTE
                cell.dict[NUM_AG_KEY] = 0

        
class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        for cell in self.cell_list:
            if cell.volume < 2:
                # if cell.type == self.CENTROBLAST or cell.type == self.CENTROCYTE:
                    # global num_plasma
                    # num_plasma += 1
                self.delete_cell(cell)
                

                
class UpdatePlotsSteppable(SteppableBasePy):

    def start(self):
        # initialize setting for Histogram
        self.plot_win = self.add_new_plot_window(title='Histogram of Cell Volumes', x_axis_title='Number of Cells',
                                                 y_axis_title='Volume Size in Pixels')
        # _alpha is transparency 0 is transparent, 255 is opaque
        self.plot_win.add_histogram_plot(plot_name='Current Affinity Levels', color='green', alpha=200)

    def step(self, mcs):
        global num_plasma
        hist_list = [1] * (len(self.cell_list_by_type(self.CENTROBLAST, self.CENTROCYTE)) + num_plasma)
        i = 0
        for cell in self.cell_list_by_type(self.CENTROBLAST, self.CENTROCYTE):
            num_ag_field[cell] = cell.dict[NUM_AG_KEY]
            hist_list[i] = cell.dict[AFFINITY_KEY]
            i += 1
        self.plot_win.add_histogram(plot_name='Current Affinity Levels', value_array=hist_list, number_of_bins=10)

    def on_stop(self):
        '''
        this gets called each time user stops simulation
        '''        
        # PLACE YOUR CODE BELOW THIS LINE
        
        return
