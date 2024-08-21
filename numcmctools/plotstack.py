from .plot import Plot
from .mcmcsamples import MCMCSamples
import numpy as np

class PlotStack:
    def __init__(self, chain: MCMCSamples):

        # Crash if passed object is not an instance of MCMCSamples
        if not isinstance(chain, MCMCSamples):
            raise TypeError(f"Expected MCMCSamples instance, got {type(chain).__name__}")

        self.chain = chain
        self.plotted_variables = []
        self.plots = []

    def add_plot(self,variables, priors, bins, axrange=None, mo_option=False):

        # Crash if user supplied a non-existant variable
        for var in variables:
            if var not in self.chain.variables:
                raise TypeError(f"The variable you supplied in plot, {var}, is not defined in the MCMCSamples chain")

        # Add variable to the internal list
        for var in variables:
            if var in self.plotted_variables:
                continue
            self.plotted_variables.append(var)

        self.plots.append(Plot(variables, priors, bins, axrange, mo_option))

    def fill_plots(self,n_steps=None, batchsize=100000):
        for batch in self.chain.tree.iterate(step_size=batchsize, library="np", entry_stop=n_steps):

            for var in self.plotted_variables:
                if var in batch:
                    continue
                batch[var] = self.chain.variables[var].evaluate(batch)

            for plot in self.plots:
                plot.fill_plot(batch)
                
        for plot in self.plots:
            plot.finalize_histogram()

        
 
