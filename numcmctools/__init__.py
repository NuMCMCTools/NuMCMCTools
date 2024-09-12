import logging
from numcmctools.mcmcsamples import MCMCSamples
from numcmctools.variable import Variable
from numcmctools.plotstack import PlotStack
from numcmctools.plot import Plot 
from numcmctools.jacobiangraph import JacobianGraph

logging.basicConfig(level=logging.INFO, format='%(asctime)s::%(levelname)s::%(name)s: - %(message)s')

def set_logging_level(level):
    logging.getLogger().setLevel(level)
