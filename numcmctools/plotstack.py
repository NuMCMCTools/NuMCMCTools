from .plot import Plot
from .mcmcsamples import MCMCSamples
from .jacobiangraph import JacobianGraph
import numpy as np
from tqdm import tqdm
import logging
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

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
        self.jacobian_graph = JacobianGraph()
        self.plotted_variables = []
        self.plots = []

    def add_plot(self,variables, priors, bins, axrange=None, mo_option=False):
        """
        Add a plot to the stack
  
        :variables: Array of strings indicating the variables to be plotted.
                    Only 1D and 2D are currently supported. Custom variables can be declared
                    when through the mcmcsamples class
        :priors: a list of priors in the Funtion:Variable format
        :bins: number of bins or bin edges, formatted for either 1 or 2D as in the numpy
               documentation for histogram (https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)
               or histogram2D (https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html)
        :axrange: See as for bins
        :mo_option: When creating intervals, calculate intervals jointly over the mass hierarchies (True) or 
                    marginalized over hierarchies (False). Default is False; to be implemented.
        """

        # Add the priors to the plotting function
        # TODO: Would be nice to have some "verbose" option.
        plot_jacobians = {}
        if not priors:
            logger.debug(f"No priors supplied for plot with variables: {variables}, will be uniform in whatever the supplied chain is in.")
            for var in self.chain.compulsory_variables:
                plot_jacobians[var] = None
        else:
            parsed_priors = self.jacobian_graph.parse_priors(priors, self.chain.compulsory_variables)
            for var in self.chain.compulsory_variables:
                if var not in parsed_priors:
                    logger.debug(f"No prior for variable {var} supplied in plot: {variables}, will be uniform in whatever the supplied chain is in.")
                    plot_jacobians[var] = None
                else:
                    logger.debug(f"Prior for variable {var} supplied in plot: {variables}: {parsed_priors[var]}")
                    plot_jacobians[var] = self.jacobian_graph.get_jacobian_func(self.chain.variable_priors[var], parsed_priors[var])

        # Crash if user supplied a non-existant variable
        for var in variables:
            if var not in self.chain.variables:
                raise TypeError(f"The variable you supplied in plot, {var}, is not defined in the MCMCSamples chain")

        # Add variable to the internal list
        for var in variables:
            if var in self.plotted_variables:
                continue
            self.plotted_variables.append(var)

        self.plots.append(Plot(variables, plot_jacobians, bins, axrange, mo_option))

    def fill_plots(self,n_steps=None, batchsize=100000):
        """
        Fill the plots from the chain

        :n_steps: Fill a maximum number of steps from the chain. None means to fill all steps
        :batchsize: number of steps to draw simultaneously from the chain, to manage memory 
                    requirements.
        """
        n_batches = 0
        if n_steps == None:
            n_batches = self.chain.tree.num_entries / batchsize
        else:
            n_batches = n_steps / batchsize

        for batch in tqdm(self.chain.tree.iterate(step_size=batchsize, library="np", entry_stop=n_steps), total=n_batches):

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
                 levels required. """    

        for plot in self.plots:
            plot.make_intervals(levels)       
 
    def draw_plots(self, plot_array_dim = [],mo_separate=True):
        """
        Draw all the plots on a freshly created figure with
        automatic alloction of the subplot array dimensions. Returns
        the figure and the axes array.

        :plot_array_dim: an array of two numbers to override the 
        automatic dimensioning

        returns the figure, subfigures, and subplots

        """
        xplt=1
        yplt=1
        if len(plot_array_dim)==2:
            xplt = plot_array_dim[0]
            yplt = plot_array_dim[1]
        else:
            xplt, yplt = self.__determine_plot_array()

        self.figplt = plt.figure()
        self.sfigplt = self.figplt.subfigures(xplt, yplt)
        self.saxesplt = []

        for index, (plot, subfig) in enumerate(zip(self.plots, self.sfigplt.flat)):
            if plot.mo_option:
                ax = subfig.subplots(1,2, sharey='row', squeeze=True)
                subfig.subplots_adjust(wspace=0)
                self.saxesplt.append(ax)
            else:
                ax = subfig.add_subplot()
                self.saxesplt.append([ax])

            plot.draw_plot(ax)

            for ax in subfig.get_axes():
                ax.label_outer()
        return self.figplt, self.sfigplt, self.saxesplt

    def draw_intervals(self, plot_array_dim = [], mo_separate=True):
        """
        Draw all the intervals plots on a freshly created figure with
        automatic alloction of the subplot array dimensions. Returns
        the figure and the axes array.

        :plot_array_dim: an array of two numbers to override the 
        automatic dimensioning

        """
        xplt=1
        yplt=1
        if len(plot_array_dim)==2:
            xplt = plot_array_dim[0]
            yplt = plot_array_dim[1]
        else:
            xplt, yplt = self.__determine_plot_array()

        self.figint = plt.figure()
        self.sfigint = self.figint.subfigures(xplt, yplt)
        self.axesint = []
        
        for index, (plot, subfig) in enumerate(zip(self.plots, self.sfigint.flat)):
            ax = None
            if plot.mo_option:
                ax = subfig.subplots(1,2, sharey='row', squeeze=True)
                self.axesint.append(ax)
                subfig.subplots_adjust(wspace=0)
            else:
                ax = subfig.add_subplot()
                self.axesint.append([ax])
            plot.draw_interval(ax)

            for ax in subfig.get_axes():
                ax.label_outer()

        return self.figint, self.sfigint, self.axesint

    def __determine_plot_array(self):
        """
        Determine the optimal array for plotting all the plots
        Original code, lightly modified from https://github.com/matplotlib/grid-strategy
        
        """
        # this is used from matplotlib/gridstrategy
        n = len(self.plots)
        special_cases = {3: (2, 2), 5: (2, 3)}

        if n in special_cases:
            return special_cases[n]

        # May not work for very large n
        n_sqrtf = np.sqrt(n)
        n_sqrt = int(np.ceil(n_sqrtf))

        if n_sqrtf == n_sqrt:
            # Perfect square, we're done
            x, y = n_sqrt, n_sqrt
            
        elif n <= n_sqrt * (n_sqrt - 1):
            # An n_sqrt x n_sqrt - 1 grid is close enough to look pretty
            # square, so if n is less than that value, will use that rather
            # than jumping all the way to a square grid.
            x, y = n_sqrt, n_sqrt - 1

        elif not (n_sqrt % 2) and n % 2:
            # If the square root is even and the number of axes is odd, in
            # order to keep the arrangement horizontally symmetrical, using a
            # grid of size (n_sqrt + 1 x n_sqrt - 1) looks best and guarantees
            # symmetry.
            x, y = (n_sqrt + 1, n_sqrt - 1)

        else:
            # It's not a perfect square, but a square grid is best
            x, y = n_sqrt, n_sqrt

        # If exactly one of these is odd, make it the rows
        if (x % 2) != (y % 2) and (x % 2):
            x, y = y, x

        return x,y
