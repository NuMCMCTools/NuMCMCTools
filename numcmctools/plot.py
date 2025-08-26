import numpy as np
import matplotlib.pyplot as plt
import logging
from .jacobiangraph import JacobianGraph
from .empirical_priors import EmpiricalPrior
from typing import Union, List

logger = logging.getLogger(__name__)

class Plot:
    def __init__(self, variables, jacobians, empirical_priors, bins, axrange=None, mo_option=False):
        """
        Initialise a Plot instance.

        :variables: Array of strings indicating the variables to be plotted.
                    Only 1D and 2D are currently supported. Custom variables can be declared
                    when through the mcmcsamples class
        :jacobians: A dictionary of jacobian transform functions (can be None for
                 each parameter if no transformations are needed)
        :empirical_priors: A dictionary of empirical_priors to be applied to the plot
                      in the form of EmpiricalPriors objects
        :bins: number of bins or bin edges, formatted for either 1 or 2D as in the numpy
               documentation for histogram (https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)
               or histogram2D (https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html)
        :axrange: See as for bins
        :mo_option: When creating intervals, calculate intervals jointly over the mass hierarchies (True) or 
                    marginalized over hierarchies (False). Default is False; to be implemented.
        """
        self.variables = variables
        self.jacobian_funcs = jacobians
        self.empirical_priors = empirical_priors
        self.bins = bins
        self.axrange = axrange
        self.mo_option = mo_option
        self.nvar = len(self.variables)
        self.edges = []
        self.finalized = False

        if(self.nvar==1):
            self.hist, edges = np.histogram([], self.bins, self.axrange, weights=[])
            self.edges.append(edges)
            if self.mo_option:
                self.hist_no = np.zeros(np.shape(self.hist))
                self.hist_io = np.zeros(np.shape(self.hist))
        elif(self.nvar==2):
            self.hist, edgesx, edgesy = np.histogram2d([], [], self.bins, self.axrange, weights=[])
            self.edges.append(edgesx)
            self.edges.append(edgesy)
            if self.mo_option:
                self.hist_no = np.zeros(np.shape(self.hist))
                self.hist_io = np.zeros(np.shape(self.hist))
        else:
            logger.error("too many or too few variables!")
            return
                
    def fill_plot(self, data, weights=None):
        """
        Fill the plot. If the plot has been finalized, no more filling is allowed.

        :data: An array of the data
        :weights: Data weights
        """

        if weights==None:
            weights = np.ones(np.shape(data[self.variables[0]]))

        # Apply the jacobian transformation functions
        for var in self.jacobian_funcs:
            if not self.jacobian_funcs[var]:
                continue
            # Apply the weight
            weights *= self.jacobian_funcs[var](data[var])
        
        # Apply the empirical priors
        for empirical_prior in self.empirical_priors:
            if not isinstance(self.empirical_priors[empirical_prior], EmpiricalPrior):
                raise TypeError(f"Empirical prior {empirical_prior} is not an instance of EmpiricalPrior")
            
            if not self.empirical_priors[empirical_prior].is_inverted and not self.empirical_priors[empirical_prior].is_normal:
                raise ValueError(f"Empirical prior {empirical_prior} must be either inverted and/or normal")
            
            weights_tmp = self.empirical_priors[empirical_prior](data)

            if self.empirical_priors[empirical_prior].is_inverted and self.empirical_priors[empirical_prior].is_normal:
                weights *= weights_tmp
            elif self.empirical_priors[empirical_prior].is_inverted:
                weights *= np.where(np.less_equal(data["Deltam2_32"],0.0), weights_tmp, 1.0)
            elif self.empirical_priors[empirical_prior].is_normal:
                weights *= np.where(np.greater_equal(data["Deltam2_32"],0.0), weights_tmp, 1.0)
        
        if not self.finalized:
            if(self.nvar==1):
                if self.mo_option:
                    hist, edges = np.histogram(data[self.variables[0]], self.bins, self.axrange, weights = weights*np.greater_equal(data["Deltam2_32"],0))
                    self.hist_no += hist
                    hist, edges = np.histogram(data[self.variables[0]], self.bins, self.axrange, weights = weights*np.less_equal(data["Deltam2_32"],0))
                    self.hist_io += hist
                else:
                    hist, edges = np.histogram(data[self.variables[0]], self.bins, self.axrange, weights = weights)
                    self.hist += hist
            elif(self.nvar==2):
                if self.mo_option:
                    hist, edgesx, edgesy = np.histogram2d(data[self.variables[0]], data[self.variables[1]], self.bins, self.axrange,
                                                              weights = weights*np.greater_equal(data["Deltam2_32"],0))
                    self.hist_no += hist
                    hist, edgesx, edgesy = np.histogram2d(data[self.variables[0]], data[self.variables[1]], self.bins, self.axrange,
                                                              weights = weights*np.less_equal(data["Deltam2_32"],0))
                    self.hist_io += hist
                else:
                    hist, edgesx, edgesy = np.histogram2d(data[self.variables[0]], data[self.variables[1]], self.bins, self.axrange, weights = weights)
                    self.hist += hist
        else:
            logger.warn("histogram was finalized already! No filling allowed!")
            return

        def throw_error():
            """
            Throws an error if the filled histograms are empty
            """
            for varidx, var in enumerate(self.variables):
                min = np.min(data[var])
                max = np.max(data[var])

                ax = self.axrange if self.nvar == 1 else self.axrange[varidx]

                if not (ax[0] <= min <= ax[1]) or not (ax[0] <= max <= ax[1]):
                    raise ValueError(f"Variable {var} has posterior max and max values of {min} and {max}, "
                                     f"outside of the plotting range {ax} for plot with variables {self.variables}.")
        
        if self.mo_option:
            if np.sum(self.hist_no) == 0 or np.sum(self.hist_io) == 0:
                throw_error()
        else:
            if np.sum(self.hist) == 0:
                throw_error()

    def finalize_histogram(self):
        """
        Finalize the plot. The plot is finalized to make a probability density function.
        """
        if(not self.finalized):
            if(self.nvar==1):
                self.areas = np.diff(self.edges[0])
            if(self.nvar==2):
                self.areas = np.outer(np.diff(self.edges[0]),np.diff(self.edges[1]))

            if self.mo_option:
                self.hist = np.concatenate((self.hist_no,self.hist_io))
                self.areas = np.concatenate((self.areas,self.areas))
            total = np.sum(self.hist)
            self.hist = self.hist/self.areas/total
            if self.mo_option:
                if(self.nvar==1):
                    sh = np.shape(self.hist_no)
                    self.hist_no = self.hist[:sh[0]]
                    self.hist_io = self.hist[sh[0]:]
                if(self.nvar==2):
                    sh = np.shape(self.hist_no)
                    self.hist_no = self.hist[:sh[0],:]
                    self.hist_io = self.hist[sh[0]:,:]
            self.finalized = True
            
    def draw_plot(self, saxis: Union[plt.Axes, List[plt.Axes]]):
        """
        Draw the plot. 
        :saxis: matplotlib axis, or a list of axes to draw the plot on

        returns 1 or 2 matplotlib subplots for further manipulation
        """

        # Make sure the input either plt.Axes or array of plt.Axes
        is_axis_array = False
        if isinstance(saxis, plt.Axes):
            if self.nvar > 1 and self.mo_option:
                raise ValueError(f"Cannot plot 2 2D heatmaps on top of each other, need to provide 2 plt.Axes in a list")
        elif isinstance(saxis, (list, np.ndarray)) and all(isinstance(ax, plt.Axes) for ax in saxis):
            if len(saxis) > 2:
                raise ValueError(f"Maximum length of the axis list is 2, one per mass ordering of Deltam2_32")
            is_axis_array = True
        else:
            raise ValueError("Wrong instance of axes, should be either plt.Axes or List[plt.axes]!")

        if(self.nvar==1):
            if self.mo_option:
                if is_axis_array:
                    saxis[0].stairs(self.hist_no,self.edges[0])
                    saxis[0].set_xlabel(self.variables[0]+" NO")
                    saxis[1].stairs(self.hist_io,self.edges[0])
                    saxis[1].set_xlabel(self.variables[0]+" IO")
                else:
                    saxis.stairs(self.hist_no,self.edges[0])
                    saxis.stairs(self.hist_io,self.edges[0])
                    saxis.set_xlabel(self.variables[0])

            else:
                saxis.stairs(self.hist,self.edges[0])
                saxis.set_xlabel(self.variables[0])

        if(self.nvar==2):
            if self.mo_option:
                if is_axis_array:
                    cm = saxis[0].pcolormesh(self.edges[0], self.edges[1], self.hist_no.T)
                    saxis[0].set_xlabel(self.variables[0]+" NO")
                    saxis[0].set_ylabel(self.variables[1])
                    cm = saxis[1].pcolormesh(self.edges[0], self.edges[1], self.hist_io.T)
                    saxis[1].set_xlabel(self.variables[0]+" IO")
                else:
                    logger.error("2D color orverlap looks terrible, not drawn!")
                    #cm = saxis.pcolormesh(self.edges[0], self.edges[1], self.hist_no.T)
                    #saxis[0].set_xlabel(self.variables[0]+" NO")
                    #saxis[0].set_ylabel(self.variables[1])
                    #cm = saxis[1].pcolormesh(self.edges[0], self.edges[1], self.hist_io.T)
                    #saxis[1].set_xlabel(self.variables[0]+" IO")
            else:
                cm = saxis.pcolormesh(self.edges[0], self.edges[1], self.hist.T)
                saxis.set_xlabel(self.variables[0])
                saxis.set_ylabel(self.variables[1])

    def draw_interval(self, saxis: Union[plt.Axes, List[plt.Axes]]):
        """
        Draw the intervals.
        :saxes: matplotlib axes, or list of axes, to draw the plot on

        returns 1 or 2 matplotlib subplots for further manipulation
        """

        is_axis_array = False
        if isinstance(saxis, plt.Axes):
            is_axis_array = False
        elif isinstance(saxis, (list, np.ndarray)) and all(isinstance(ax, plt.Axes) for ax in saxis):
            if len(saxis) > 2:
                raise ValueError(f"Maximum length of the axis list is 2, one per mass ordering of Deltam2_32")
            is_axis_array = True
        else:
            raise ValueError("Wrong instance of axes, should be either plt.Axes or List[plt.axes]!")
        
        #need to put in a check if the intervals have been calculated
        if(self.nvar==1):
            if self.mo_option:

                if is_axis_array:
                    saxis[0].stairs(self.hist_no,self.edges[0], color='black')
                    for lev in self.prob_levels:
                        saxis[0].stairs(self.hist_no*np.greater_equal(self.hist_no,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
                    saxis[0].set_xlabel(self.variables[0]+" NO")

                    saxis[1].stairs(self.hist_io,self.edges[0], color='black')
                    for lev in self.prob_levels:
                        saxis[1].stairs(self.hist_io*np.greater_equal(self.hist_io,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
                    saxis[1].set_xlabel(self.variables[0]+" IO")
                else:
                    saxis.stairs(self.hist_no,self.edges[0], color='cornflowerblue')
                    for lev in self.prob_levels:
                        saxis.stairs(self.hist_no*np.greater_equal(self.hist_no,lev),self.edges[0], fill=True, color='cornflowerblue', alpha=0.3)
                    saxis.set_xlabel(self.variables[0])

                    saxis.stairs(self.hist_io,self.edges[0], color='lightcoral')
                    for lev in self.prob_levels:
                        saxis.stairs(self.hist_io*np.greater_equal(self.hist_io,lev),self.edges[0], fill=True, color='lightcoral', alpha=0.3)
            else:
                saxis.stairs(self.hist,self.edges[0], color='black')
                for lev in self.prob_levels:
                    saxis.stairs(self.hist*np.greater_equal(self.hist,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
                saxis.set_xlabel(self.variables[0])
        
        if(self.nvar==2):

            linestyles_base = ['solid','dashed','dashdot','dotted']
            linestyles=['']*len(self.prob_levels)
            for i in range(len(self.prob_levels)):
                linestyles[i] = linestyles_base[i%len(linestyles_base)]
            linestyles = list(reversed(linestyles))

            if self.mo_option:
                if is_axis_array:
                    saxis[0].contour(0.5*(self.edges[0][:-1]+self.edges[0][1:]), 0.5*(self.edges[1][:-1]+self.edges[1][1:]),self.hist_no.T, np.sort(self.prob_levels), linestyles=linestyles, colors='lightgrey')
                    saxis[0].set_xlabel(self.variables[0]+" NO")
                    saxis[0].set_ylabel(self.variables[1])
                    saxis[1].contour(0.5*(self.edges[0][:-1]+self.edges[0][1:]), 0.5*(self.edges[1][:-1]+self.edges[1][1:]),self.hist_io.T, np.sort(self.prob_levels), linestyles=linestyles, colors='lightgrey')
                    saxis[1].set_xlabel(self.variables[0]+" IO")
                else:
                    saxis.contour(0.5*(self.edges[0][:-1]+self.edges[0][1:]), 0.5*(self.edges[1][:-1]+self.edges[1][1:]),self.hist_no.T, np.sort(self.prob_levels), linestyles=linestyles, colors='cornflowerblue')
                    saxis.set_xlabel(self.variables[0]+" NO")
                    saxis.set_ylabel(self.variables[1])
                    saxis.contour(0.5*(self.edges[0][:-1]+self.edges[0][1:]), 0.5*(self.edges[1][:-1]+self.edges[1][1:]),self.hist_io.T, np.sort(self.prob_levels), linestyles=linestyles, colors='lightcoral')
                    saxis.set_xlabel(self.variables[0]+" IO")
            else:
                saxis.contour(0.5*(self.edges[0][:-1]+self.edges[0][1:]), 0.5*(self.edges[1][:-1]+self.edges[1][1:]),self.hist.T, np.sort(self.prob_levels), linestyles=linestyles, colors='lightgrey')
                saxis.set_xlabel(self.variables[0])
                saxis.set_ylabel(self.variables[1])
        
    def make_intervals(self,levels):
        """
        Create intervals. 
        :levels: an array of numbers between 0 and 1 representing the credible interval
                 levels required. 
        """
        if(not self.finalized):
            self.finalize_histogram()
        
        self.levels = -1*np.sort(-1*np.array(levels))
        nlev = len(self.levels)-1

        index_sort = np.argsort(-1*self.hist*self.areas, axis=None) #what a stupid hack--order high to low
        index_sort_unrav = np.unravel_index(index_sort,self.hist.shape)

        self.prob_levels = np.zeros(self.levels.shape)

        total = np.sum(self.hist*self.areas)
        index = 0
        process_sum = 0
        if(self.nvar==1):
            while nlev >= 0:
                process_sum+=self.hist[index_sort_unrav[0][index]]*self.areas[index_sort_unrav[0][index]]
                index+=1
                if(process_sum/total > self.levels[nlev]):
                    self.prob_levels[nlev]=self.hist[index_sort_unrav[0][index-1]]
                    nlev-=1
        if(self.nvar==2):       
            while nlev >= 0:
                process_sum+=self.hist[index_sort_unrav[0][index],index_sort_unrav[1][index]]*self.areas[index_sort_unrav[0][index],index_sort_unrav[1][index]]
                index+=1
                if(process_sum/total > self.levels[nlev]):
                    self.prob_levels[nlev]=self.hist[index_sort_unrav[0][index-1],index_sort_unrav[1][index-1]]
                    nlev-=1            