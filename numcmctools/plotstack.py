from .plot import Plot
from .mcmcsamples import MCMCSamples
import numpy as np

class PlotStack:
    def __init__(self, chain: MCMCSamples):
       """
        Initialise a PlotStack instance.

        :chain: An instance of MCMCSamples that plots will be made from
        """
        # Crash if passed object is not an instance of MCMCSamples
        if not isinstance(chain, MCMCSamples):
            raise TypeError(f"Expected MCMCSamples instance, got {type(chain).__name__}")

        self.chain = chain
        self.plotted_variables = []
        self.plots = []

    def add_plot(self,variables, priors, bins, axrange=None, mo_option=False):
        """
        Add a plot to the stack
  
        :variables: Array of strings indicating the variables to be plotted.
                    Only 1D and 2D are currently supported. Custom variables can be declared
                    when through the mcmcsamples class
        :priors: TBD
        :bins: number of bins or bin edges, formatted for either 1 or 2D as in the numpy
               documentation for histogram (https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)
               or histogram2D (https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html)
        :axrange: See as for bins
        :mo_option: When creating intervals, calculate intervals jointly over the mass hierarchies (True) or 
                    marginalized over hierarchies (False). Default is False; to be implemented.
        """
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
        """
        Fill the plots from the chain

        :n_steps: Fill a maximum number of steps from the chain. None means to fill all steps
        :batchsize: number of steps to draw simultaneously from the chain, to manage memory 
                    requirements.
        """
        for batch in self.chain.tree.iterate(step_size=batchsize, library="np", entry_stop=n_steps):

            for var in self.plotted_variables:
                if var in batch:
                    continue
                batch[var] = self.chain.variables[var].evaluate(batch)

            for plot in self.plots:
                plot.fill_plot(batch)
                
        for plot in self.plots:
            plot.finalize_histogram()

     def make_intervals(self,levels):
        """
        Make intervals for all the plots in the stack at the same credible interval levels.

        :levels: an array of numbers between 0 and 1 representing the credible interval
                 levels required. 
        """    
        for plot in self.plots:
            plot.make_intervals(levels)       
 
