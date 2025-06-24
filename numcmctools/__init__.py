import logging
from numcmctools.mcmcsamples import MCMCSamples
from numcmctools.variable import Variable
from numcmctools.plotstack import PlotStack
from numcmctools.plot import Plot 
from numcmctools.jacobiangraph import JacobianGraph

logging.basicConfig(level=logging.INFO, format='%(asctime)s::%(levelname)s::%(name)s::%(lineno)s: - %(message)s')

def set_external_logging_level(level):
    """
    Sets the log level for all the external packages

    :param level: logging level
    """
    logging.getLogger("asyncio").setLevel(level)
    logging.getLogger("fsspec").setLevel(level)
    logging.getLogger("matplotlib").setLevel(level)

def set_logging_level(level):
    logging.getLogger().setLevel(level)

# Default log level for external packages is "WARNING"
set_external_logging_level(logging.WARNING)
