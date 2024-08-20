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

            evaluated_data = {}
            for variable in self.plotted_variables:
                evaluated_data[variable] = self.chain.variables[variable].evaluate(batch)

            for plot in self.plots:
                if (plot.nvar==1):
                    plot.fill_plot(evaluated_data[plot.variables[0]])
                if(plot.nvar==2):
                    data = np.stack((evaluated_data[plot.variables[0]],evaluated_data[plot.variables[1]]))
                    plot.fill_plot(data)

        for plot in self.plots:
            plot.finalize_histogram()

        
 
