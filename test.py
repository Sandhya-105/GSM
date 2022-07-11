from pcraster import *
from pcraster.framework import *
import random 


class shrubmanagement(DynamicModel):

  def __init__(self):
    DynamicModel.__init__(self)
    setclone('mask.map')
      
  def initial(self):
    # initial distribution of shrubs 10%, grass 80%, empty 10%
    #map = uniform(1)
    #initialdistr = ifthenelse(map < 0.1, nominal(3), ifthenelse(map>0.9, nominal(0), nominal(2))) 
    #self.report(initialdistr, "initialdistr")
    #self.resultMap = self.readmap("initialdistr")
    self.resultMap = self.readmap("initial3")
    self.totalCells = 200 * 200

  def dynamic(self): # 0 empty, 2 grass, 3 shrub
      
    #establish the parameters, detailed in the paper Table 1
    shrub_empty_rate = 6.8+mapnormal()*2.323
    shrub_grass_rate = 0.387+mapnormal()*0.082
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
    prob_qSS=scalar(ifthen(self.resultMap == 3, window4total(scalar(self.resultMap == 3))/4))
    
    prob_qGE=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 0))/4))
    prob_qGS=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 3))/4))

    prob_qES=scalar(ifthen(self.resultMap == 0, window4total(scalar(self.resultMap == 3))/4))
    prob_qEG=scalar(ifthen(self.resultMap == 0, window4total(scalar(self.resultMap == 2))/4))


    #write to disk
    self.report(prob_qSE, "Results/SE")
    self.report(prob_qSG, "Results/SG")
    self.report(prob_qSS, "Results/SS")
    self.report(prob_qGE, "Results/GE")
    self.report(prob_qGS, "Results/GS")
    self.report(prob_qES, "Results/ES")
    self.report(prob_qEG, "Results/EG")

    #t is the total number of cells for empty, grass, and shrub
    t_ec = maptotal(scalar(self.resultMap == 0))
    t_gc= maptotal(scalar(self.resultMap == 2))
    t_sc = maptotal(scalar(self.resultMap == 3))

    #gp is global probability
    gp_e = t_ec / self.totalCells #10 percent empty in intial time step
    gp_g = t_gc / self.totalCells #80 percent grass in initial time step
    gp_s = t_sc / self.totalCells #10 percent shrub in initial time step

    self.report(t_ec, "Results/t_ec")
    self.report(gp_e, "Results/gp_e")
    self.report(gp_g, "Results/gp_g")
    self.report(gp_s, "Results/gp_s")


    #writing the transition rules
    #Colonization of an empty cell by a shrub
    transition_ES = ((1-prob_qES) * shrub_empty_rate + shrub_empty_reprd) * prob_qES
    
    #Colonization of an empty cell by grass
    transition_EG = grass_empty_clonal_rate * prob_qEG + grass_seeds_empty_rate * gp_g
    
    #Colonization of grass cells by shrub
    transition_GS = ((1-prob_qGS) *shrub_grass_rate + shrub_grass_reprd) * prob_qGS
    
    #Mortality rates
    shrub_death = shrub_mortality_rate + shrub_compete_rate * prob_qSS
    grass_death = grass_mortality_rate

    self.report(transition_ES, "Results/transES")
    self.report(transition_EG, "Results/transEG")
    self.report(transition_GS, "Results/transGS")
    self.report(shrub_death, "Results/s_death")

    #final distribution map 

    randomNumber = uniform(1)
    self.report(randomNumber, "Results/f_distr")

    fdistr_ES = ifthenelse((self.resultMap == 0) & (randomNumber<transition_ES), boolean(1), boolean(0))
    fdistr_EG = ifthenelse((self.resultMap == 0) & (randomNumber<transition_EG), boolean(1), boolean(0))
    fdistr_GS = ifthenelse((self.resultMap == 2) & (randomNumber<transition_GS), boolean(1), boolean(0))
    fdistr_shrub_death = ifthenelse((self.resultMap == 3) & (randomNumber<shrub_death), boolean(1), boolean(0))
    
    grass = cover(fdistr_ES, boolean(0))
    grassclass = ifthen(grass, nominal(2))
    self.report(fdistr_ES, "Results/f_ES")

    shrub = cover(fdistr_EG, fdistr_GS, boolean(0))
    shrubclass = ifthen(shrub, nominal(3))
    self.report(fdistr_EG, "Results/f_EG")
    self.report(fdistr_GS, "Results/f_GS")

    empty = cover(fdistr_shrub_death, boolean(0))
    emptyclass = ifthen(empty, nominal(0))
    self.report(fdistr_shrub_death, "Results/f_sdeath")
    self.report(emptyclass, "Results/empty")

    self.resultMap = cover(grassclass, shrubclass, emptyclass, self.resultMap)
    self.report(self.resultMap, "Results/f_result")
   

nrOfTimeSteps=20
myModel = shrubmanagement() 
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run() 
