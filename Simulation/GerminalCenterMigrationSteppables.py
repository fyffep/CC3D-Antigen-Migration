from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np
import random as rand
import math
from copy import deepcopy
import os

#VOCABULARY:
#Ag = antigen
#100 mcs = one hour

#division time of centroblasts = ~7 hours = 700 mcs
DIVISION_TIME = 140 #FIXME...700 or 140?

NUM_AG_KEY = "ag" #a dictionary key for the number of antigen a cell holds
TOTAL_AG_COLLECTED_KEY = "tag"
LAST_DIV_TIME_KEY = "div_time"
DNA_KEY = "dna"
AFFINITY_KEY = "aff"

ANTIGEN_VALUE = ["A","C","T"]
INITIAL_DNA = ["T","A","T"]
ANTIGEN_SYMBOLS = ["A","C","T","G"]
ANTIGEN_LEN = len(ANTIGEN_VALUE)

DIVIDE_AG_ASYMMETRIC = 0.72 #from Robert2017
POLARITY_INDEX = 0.88 #from Robert2017
FDC_ANTIGEN = 3000 #from MerinoTejero2021

COLLECT_FDC_TIME_KEY = "collect_period"
COLLECT_FDC_TIME = 70 #from Robert2017

TFH_TOUCH_TIME_KEY = "tfh_time"
#below times are measured in mcs
TFH_RESCUE_TIME = 50 #from Robert2017
TFH_TOUCH_TIME = 60 #from Robert2017

affinity_field = None
num_ag_field = None

num_plasma = 0



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
        
        for x in range(156, 256-TFH_CELL_SIZE, TFH_CELL_SIZE*3):
            for y in range(TFH_CELL_SIZE, 256-TFH_CELL_SIZE, TFH_CELL_SIZE*2):
                if rand.random() < 0.4:
                    self.cell_field[x:x+TFH_CELL_SIZE, y:y+TFH_CELL_SIZE, 0] = self.new_cell(self.TFH)
                    
        for x in range(156, 256-DENDRITIC_CELL_SIZE, DENDRITIC_CELL_SIZE*3):
            for y in range(DENDRITIC_CELL_SIZE, 256-DENDRITIC_CELL_SIZE, DENDRITIC_CELL_SIZE*2):
                if rand.random() < 0.6:
                    self.cell_field[x:x+DENDRITIC_CELL_SIZE, y:y+DENDRITIC_CELL_SIZE, 0] = self.new_cell(self.DENDRITIC)
        
        original_affinity = judge_affinity(INITIAL_DNA)
        for x in range(0, 120-B_CELL_SIZE, B_CELL_SIZE*5):
            for y in range(B_CELL_SIZE, 256-B_CELL_SIZE, B_CELL_SIZE*5):
                if rand.random() < 0.5:
                    b_cell = self.new_cell(self.CENTROBLAST)
                    b_cell.dict[NUM_AG_KEY] = FDC_ANTIGEN
                    b_cell.dict[TOTAL_AG_COLLECTED_KEY] = 0
                    b_cell.dict[LAST_DIV_TIME_KEY] = -10000
                    b_cell.dict[DNA_KEY] = deepcopy(INITIAL_DNA) #shallow copy of DNA list
                    b_cell.dict[AFFINITY_KEY] = original_affinity
                    b_cell.dict[COLLECT_FDC_TIME_KEY] = 0
                    self.cell_field[x:x+B_CELL_SIZE, y:y+B_CELL_SIZE, 0] = b_cell
        
        print("Initial no. TFH cells:", len(self.cell_list_by_type(self.TFH)))
        print("Initial no. DENDRITIC cells:", len(self.cell_list_by_type(self.DENDRITIC)))
        print("Initial no. CENTROBLAST cells:", len(self.cell_list_by_type(self.CENTROBLAST)))


    def start(self):
        self.setup_cells()
        
        self.ag_plot_win = self.add_new_plot_window(title='Plasmablast Exit Time v. Antigen Collected',
                                                 x_axis_title='Time Plasmablast Formed (log2 of MCS)',
                                                 y_axis_title='Total Antigen Received from FDCs', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        self.ag_plot_win.add_plot("plot_ag_collected", style='Lines', color='red', size=5)
        
        self.death_plot_win = self.add_new_plot_window(title='Affinity Level Upon Apoptosis',
                                                 x_axis_title='Time of B Cell Death (log2 of MCS)',
                                                 y_axis_title='Affinity Level of Cell', x_scale_type='linear', y_scale_type='linear',
                                                 grid=True)
        self.death_plot_win.add_plot("plot_death", style='Lines', color='green', size=5)
        
        
        
        
        
        
            
    def step(self, mcs):
        global num_ag_field
        global num_plasma
        
        for cell in self.cell_list_by_type(self.CENTROCYTE):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    interacting = False
                    
                    if neighbor.type == self.DENDRITIC:
                        interacting = True
                        #Robert2017: bindprobability ← affinity(BCR, antigen) . ( F.antigenAmount[f] / antigenSaturation)
                        #Cell always picks up 400 antigen TODO improve
                        # ag_count = 400
                        # cell.dict[TOTAL_AG_COLLECTED_KEY] += ag_count
                        # cell.dict[NUM_AG_KEY] += ag_count
                        cell.dict[TOTAL_AG_COLLECTED_KEY] = 3000
                        cell.dict[NUM_AG_KEY] = 3000
                        
                        #add 1 for the elapsed mcs
                        # collect_period = 1 + cell.dict[COLLECT_FDC_TIME_KEY]
                        # if collect_period >= COLLECT_FDC_TIME:
                            # collect_period = 0
                            # if cell.dict[NUM_AG_KEY] < 2000:
                                # #kill cell that didn't pick up enough antigen
                                # cell.targetVolume = 0
                        # cell.dict[COLLECT_FDC_TIME_KEY] = collect_period
            
        for tfh_cell in self.cell_list_by_type(self.TFH):
            for b_cell, common_surface_area in self.get_cell_neighbor_data_list(tfh_cell):
                    if b_cell and b_cell.type == self.CENTROCYTE:
                        if b_cell.dict[NUM_AG_KEY] >= 2000: #If FDC selected
                            # interacting = True
                            affinity = b_cell.dict[AFFINITY_KEY]
                            if affinity >= 1.0:
                                b_cell.type = self.PLASMA
                                global num_plasma
                                num_plasma += 1
                                print("num_plasma",num_plasma)
                                
                                self.ag_plot_win.add_data_point("plot_ag_collected", round(math.log2(mcs), 2), cell.dict[TOTAL_AG_COLLECTED_KEY])
                            else:
                                b_cell.type = self.CENTROBLAST
                            
                    # if interacting:
                        # cell.
                        # # Make sure Chemotaxis Plugin is loaded
                        # # modifying chemotaxis properties of individual cell 'cell'
                        # cd = self.chemotaxisPlugin.getChemotaxisData(cell, "CXCL13")
                        # if cd:
                            # l = cd.getLambda() - 3
                            # cd.setLambda(l)
            
            
    # def on_stop(self):
        # print('stop')
        
        
        
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        
        for cell in self.cell_list_by_type(self.CENTROBLAST):
            last_div_time = cell.dict[LAST_DIV_TIME_KEY]
            #Here, there are some arbitrary rules to prevent crowding
            if mcs - last_div_time >= DIVISION_TIME and cell.xCOM < 140:
                cell.dict[LAST_DIV_TIME_KEY] = mcs
                if len(self.get_cell_neighbor_data_list(cell)) < 5:
                    cells_to_divide.append(cell)
        
        for cell in cells_to_divide:
            self.divide_cell_random_orientation(cell)
    

    def update_attributes(self):
        global affinity_field             

        self.clone_parent_2_child()

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        
        total_ag = self.parent_cell.dict[NUM_AG_KEY]
        #random chance to perform asymmetric division 
        #and a centroblast must divide at least 2-6 times before becoming a centrobyte
        if rand.random() < DIVIDE_AG_ASYMMETRIC and not total_ag > FDC_ANTIGEN//4:
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
        self.ag_plot_win = self.add_new_plot_window(title='Number of Cells with Affinity', x_axis_title='Affinity Level',
                                                 y_axis_title='Number of Cells')
        # _alpha is transparency 0 is transparent, 255 is opaque
        self.ag_plot_win.add_histogram_plot(plot_name='Current Affinity Levels', color='green', alpha=200)

    def step(self, mcs):
        global num_plasma
        all_b_cells = self.cell_list_by_type(self.CENTROBLAST, self.CENTROCYTE)
        hist_list = [1] * (len(all_b_cells) + num_plasma)
        i = 0
        for cell in all_b_cells:
            num_ag_field[cell] = cell.dict[NUM_AG_KEY]
            hist_list[i] = cell.dict[AFFINITY_KEY]
            i += 1
        self.ag_plot_win.add_histogram(plot_name='Current Affinity Levels', value_array=hist_list, number_of_bins=10)
        
        
        

    def on_stop(self):
        '''
        this gets called each time user stops simulation
        '''        
        
        # file_name = "C:\\Users\\Pete\\Desktop\\b-cell-project\\output\\hist"
        # i = 0
        # while os.path.exists(file_name + str(i) + ".png"):
            # print(file_name + str(i) + ".png exists")
            # i += 1
        # file_name = file_name + str(i) + ".png"
        
        # self.ag_plot_win.save_plot_as_png(file_name)
        
        
        
        