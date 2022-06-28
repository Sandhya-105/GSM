from pcraster import *
from pcraster.framework import *


class MyFirstModel(DynamicModel, MonteCarloModel):
  def __init__(self):
    DynamicModel.__init__(self)
    MonteCarloModel.__init__(self)
    setclone('clone.map')
