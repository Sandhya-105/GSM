from pcraster import *
from pcraster.framework import *

class MyFirstModel(DynamicModel):
  def __init__(self):
    DynamicModel.__init__(self)
    setclone('clone.map')

  def initial(self):
    aUniformMap = uniform(1)
    self.alive < .10
    alive = uniform(1)<.10

  def dynamic(self):
    pass

nrOfTimeSteps=1
myModel = MyFirstModel()
dynamicModel = DynamicFramework(myModel,nrOfTimeSteps)
dynamicModel.run()
