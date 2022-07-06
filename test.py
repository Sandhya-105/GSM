from pcraster import *
from pcraster.framework import *

shrubdevelopment = 0

class ShrubEncroachment(DynamicModel):

  def __init__(self, brnR, grzR):
      DynamicModel.__init__(self)
      setclone('clone.map')
      self.burn_rate = brnR
      self.graze_rate = grzR
      
  def initial(self):
      # initial distribution of shrubs 10%, grass 80%, empty 10%
      self.resultMap = self.readmap("initial_distribution")
      self.report(self.resultMap, "Result_map")
      self.adultplants = 40
      self.totalCells = 200 * 200
      self.currentTimestep = 0

  def dynamic(self): # 0 empty, 1 fire, 2 grass, 3 shrub
    self.currentTimestep = self.currentTimestep +1 
      global shrubDevelopment
      
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
      
      prob_qGE=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 0))/4))
      prob_qGS=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 3))/4))
      prob_qGB=scalar(ifthen(self.resultMap == 2, window4total(scalar(self.resultMap == 1))/4))

      prob_qES=scalar(ifthen(self.resultMap == 0, window4total(scalar(self.resultMap == 3))/4))
      prob_qEG=scalar(ifthen(self.resultMap == 0, window4total(scalar(self.resultMap == 2))/4))

      prob_qBS=scalar(ifthen(self.resultMap == 1, window4total(scalar(self.resultMap == 3))/4))
      prob_qBG=scalar(ifthen(self.resultMap == 1, window4total(scalar(self.resultMap == 2))/4))

      #write to disk
      self.report(prob_qSE, "Results/SE")
      self.report(prob_qSG, "Results/SG")
      self.report(prob_qSB, "Results/SB")
      self.report(prob_qGE, "Results/GE")
      self.report(prob_qGS, "Results/GS")
      self.report(prob_qGB, "Results/GB")
      self.report(prob_qBS, "Results/BS")
      self.report(prob_qBG, "Results/BG")


      self.resultMap = ifthenelse(self.resultMap == 1, 0, self.resultMap)
      total_empty_cells = maptotal(scalar(self.resultMap == 0))
      total_burn_cells = maptotal(scalar(self.resultMap == 1))
      total_grass_cells= maptotal(scalar(self.resultMap == 2))
      total_shrub_cells = maptotal(scalar(self.resultMap == 3))

      global_prob_empty = total_empty_cells / self.totalCells #10 percent empty in intial time step
      global_prob_grass = total_grass_cells / self.totalCells #80 percent grass in initial time step
      global_prob_shrub = total_shrub_cells / self.totalCells #10 percent shrub in initial time step

      self.report(total_empty_cells, "Results/total_empty_cells")
      self.report(global_prob_empty, "Results/prob_empty_cells")
      self.report(global_prob_grass, "Results/prob_grass_cells")
      self.report(global_prob_shrub, "Results/prob_shrub_cells")


      #writing the transition rules
      # 1 Colonization of an empty cell by a shrub
      transition_ES = ((1-prob_qES) * shrub_empty_rate + shrub_empty_reprd) * prob_qES*(1-self.graze_rate)
      
      # 2 Colonization of an empty cell by grass
      transition_EG = grass_empty_clonal_rate * prob_qEG + grass_seeds_empty_rate * global_prob_grass
      
      # 3 Colonization of grass cells by shrub
      transition_GS = ((1-local_prob_qGS) *estab_rate_shrub_on_grass + repro_rate_shrub_on_grass) * local_prob_qGS*(1-self.grazingRate)
      
      # 4 Mortality
      shrub_death = mortality_rate_shrub + competition_shrub_center * local_prob_qSS
      grass_death = mortality_rate_grass
      
      # 5 Colonization of burned cells by shrub
      transition_BS = ((1-local_prob_qBS) * estab_rate_shrub_on_burned + repro_rate_shrub_on_empty) * local_prob_qBS*(1-self.grazingRate)
      
      # 6 Colonization of burned cells by grass
      transition_BG = clonal_rate_gras_on_empty * local_prob_qBG + estab_rate_grass_on_empty * global_prob_grass
      #transition_burned_grass = 0.5


      # 9 shrub removed -> empty
      # 10 grass removed -> empty




      self.report(transition_empty_shrub, "Report/tes")
      self.report(transition_empty_grass, "Report/teg")
      self.report(transition_grass_shrub, "Report/tgs")
      self.report(mortality_shrub, "Report/ms")

nrOfTimeSteps=1
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()
