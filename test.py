from pcraster import *
from pcraster.framework import *
import numpy as np
import random 

shrubdevelopment = 0

class shrubmanagement(DynamicModel):

  def __init__(self):
    DynamicModel.__init__(self)
    setclone('mask.map')
    #self.burn_rate = 'brnR'
      
  def initial(self):
    # initial distribution of shrubs 10%, grass 80%, empty 10%
    map = uniform(1)
    initialdistr = ifthenelse(map < 0.1, nominal(3), ifthenelse(map>0.9, nominal(0), nominal(2))) 
    self.report(initialdistr, "initialdistr")
    self.resultMap = self.readmap("initialdistr")
    
    #self.report(self.resultMap, "result_map")
    self.adultplants = 40
    self.totalCells = 200 * 200
    self.currentTimestep = 0

  def dynamic(self): # 0 empty, 1 fire, 2 grass, 3 shrub
    self.currentTimestep = self.currentTimestep +1 
    global shrubdevelopment
      
    #establish the parameters, detailed in the paper Table 1
    shrub_empty_rate = 6.8+mapnormal()*2.323
    shrub_grass_rate = 0.387+mapnormal()*0.082
    shrub_burned_rate = 38.8+mapnormal()*9.768
    shrub_empty_reprd = 0.109+mapnormal()*0.01
    shrub_grass_reprd = 0.061+mapnormal()*0.016
    grass_empty_clonal_rate = 0.5
    shrub_compete_rate = 0.00015+mapnormal()*0.001
    shrub_mortality_rate = 0.028+mapnormal()*0.0005
    grass_mortality_rate = 0.125
    grass_seeds_empty_rate = 0.8


    #defining local probabilities for the transition functions (q = density), ie how the von neumann neighborhood effects the state of the center cell (probability that if I am a shrub, for example, that I become a grass cell)
    prob_qSE=scalar(ifthen(self.resultMap == 3, window4total(scalar(self.resultMap == 0))/4)) #probability that a shrub cell becomes empty 
    prob_qSG=scalar(ifthen(self.resultMap == 3, window4total(scalar(self.resultMap == 2))/4))
    prob_qSB=scalar(ifthen(self.resultMap == 3, window4total(scalar(self.resultMap == 1))/4))
    prob_qSS=scalar(ifthen(self.resultMap == 3, window4total(scalar(self.resultMap == 3))/4))
    
    prob_qGE=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 0))/4))
    prob_qGS=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 3))/4))
    prob_qGB=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 1))/4))

    prob_qES=scalar(ifthen(self.resultMap == 0, window4total(scalar(self.resultMap == 3))/4))
    prob_qEG=scalar(ifthen(self.resultMap == 0, window4total(scalar(self.resultMap == 2))/4))

    prob_qBS=scalar(ifthen(self.resultMap == 1, window4total(scalar(self.resultMap == 3))/4))
    prob_qBG=scalar(ifthen(self.resultMap == 1, window4total(scalar(self.resultMap == 2))/4))
    prob_qBE=scalar(ifthen(self.resultMap == 1, window4total(scalar(self.resultMap == 0))/4))

    #write to disk
    self.report(prob_qSE, "Results/SE")
    self.report(prob_qSG, "Results/SG")
    self.report(prob_qSB, "Results/SB")
    self.report(prob_qSS, "Results/SS")
    self.report(prob_qGE, "Results/GE")
    self.report(prob_qGS, "Results/GS")
    self.report(prob_qGB, "Results/GB")
    self.report(prob_qES, "Results/ES")
    self.report(prob_qEG, "Results/EG")
    self.report(prob_qBS, "Results/BS")
    self.report(prob_qBG, "Results/BG")
    self.report(prob_qBE, "Results/BE")

    #t is the total number of cells for empty, burned, grass, and shrub
    self.resultMap = ifthenelse(self.resultMap == 1, 0, self.resultMap)
    t_ec = maptotal(scalar(self.resultMap == 0))
    t_bc = maptotal(scalar(self.resultMap == 1))
    t_gc= maptotal(scalar(self.resultMap == 2))
    t_sc = maptotal(scalar(self.resultMap == 3))

    #gp is global probability
    gp_e = t_ec / self.totalCells #10 percent empty in intial time step
    gp_g = t_gc / self.totalCells #80 percent grass in initial time step
    gp_s = t_sc / self.totalCells #10 percent shrub in initial time step

    self.report(t_ec, "Results/t_ec")
    self.report(t_bc, "Results/t_bc")
    self.report(gp_e, "Results/gp_e")
    self.report(gp_g, "Results/gp_g")
    self.report(gp_s, "Results/gp_s")


    #writing the transition rules
    #Colonization of an empty cell by a shrub
    transition_ES = ((1-prob_qSE) * shrub_empty_rate + shrub_empty_reprd) * prob_qSE
    
    #Colonization of an empty cell by grass
    transition_EG = grass_empty_clonal_rate * prob_qGE + grass_seeds_empty_rate * gp_g
    
    #Colonization of grass cells by shrub
    transition_GS = ((1-prob_qSG) *shrub_grass_rate + shrub_grass_reprd) * prob_qSG
    
    #Mortality rates
    shrub_death = shrub_mortality_rate + shrub_compete_rate * prob_qSS
    grass_death = grass_mortality_rate
    
    #Colonization of burned cells by shrub
    transition_BS = ((1-prob_qSB) * shrub_burned_rate + shrub_empty_reprd) * prob_qSB
    
    #Colonization of burned cells by grass
    transition_BG = grass_empty_clonal_rate * prob_qGB + grass_seeds_empty_rate * gp_g


    self.report(transition_ES, "Results/transES")
    self.report(transition_EG, "Results/transEG")
    self.report(transition_GS, "Results/transGS")
    self.report(shrub_death, "Results/s_death")
    self.report(grass_death, "Results/g_death")
    self.report(transition_BS, "Results/transBS")
    self.report(transition_BG, "Results/transBG")

    #fire management strategy

    N = 100 
    ON = 1
    OFF = 0
    random.seed = (1379)
    global cells
    cells = np.zeros((N, N)).reshape(N, N)
    
    for x in range (1, 200):
      cells.itemset((x, x), ON)
      cells.itemset((2 * x, int(x / 2)), ON)
      cells.itemset((3 * x, int(x / 2)), ON)
    
    newGrid = cells.copy()
    for x in range(1, N - 1):
        for y in range(1, N - 1):
            cell_state = window4total(scalar(self.resultMap == 3))
            if cell_state == 1: 
                for neighbor in window4total:
                    if newGrid[neighbor] == 3 and random.random() < prob_qSB:
                        newGrid[neighbor] = ON
                if random.random() < prob_qBE:
                    newGrid[x][y] = 1
            elif cell_state == 1 and random.random() < prob_qSE:
                newGrid[x][y] = 3
            elif curr_state == 2:
                for neighbor in window4total[(x, y)]:
                    if newGrid[neighbor] == 1 and random.random() < prob_qGE:
                        newGrid[neighbor] = 2


nrOfTimeSteps=20
myModel = shrubmanagement() 
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run() 
