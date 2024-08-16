from .plot import Plot
from .mcmcsamples import MCMCSamples
import numpy as np

class PlotStack:
    def __init__(self, chain: MCMCSamples):

        # Crash if passed object is not an instance of MCMCSamples
        if not isinstance(chain, MCMCSamples):
            raise TypeError(f"Expected MCMCSamples instance, got {type(chain).__name__}")

        self.chain = chain
        self.plots = []

    def add_plot(self,variables, priors, bins, axrange=None, mo_option=False):
        self.plots.append(Plot(variables, priors, bins, axrange, mo_option))

    def fill_plots(self,n_steps=None, batchsize=100000):
        for batch in self.chain.tree.iterate(step_size=batchsize, library="np", entry_stop=n_steps):
            for plot in self.plots:
                if (plot.nvar==1):
                    plot.fill_plot(batch[plot.variables[0]])
                if(plot.nvar==2):
                    data = np.stack((batch[plot.variables[0]],batch[plot.variables[1]]))
                    plot.fill_plot(data)

        for plot in self.plots:
            plot.finalize_histogram()

        
 
