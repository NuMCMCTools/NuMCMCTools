import numpy as np
import matplotlib.pyplot as plt

class Plot:
        
    def __init__(self, variables, priors, bins, axrange=None, mo_option=False):
        """
        Initialise a Plot instance.

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
        self.variables = variables
        self.priors = priors
        self.bins = bins
        self.axrange = axrange
        self.mo_option = mo_option
        self.nvar = len(self.variables)
        self.edges = []
        self.finalized = False
                
        if(self.nvar==1):
            self.hist, edges = np.histogram([], self.bins, self.axrange)
            self.edges.append(edges)
        elif(self.nvar==2):
            self.hist, edgesx, edgesy = np.histogram2d([], [],self.bins, self.axrange)
            self.edges.append(edgesx)
            self.edges.append(edgesy)
        else:
            print("too many or too few variables!")
            return
        
        self.hist_interval = []
        
    def fill_plot(self, data, weights=None):
        """
        Fill the plot. If the plot has been finalized, no more filling is allowed.

        :data: An array of the data
        :weights: Data weights
        """
        if not self.finalized:
            if(self.nvar==1):
                hist, edges = np.histogram(data[self.variables[0]], self.bins, self.axrange, weights = weights)
                self.hist += hist
            elif(self.nvar==2):
                hist, edgesx, edgesy = np.histogram2d(data[self.variables[0]], data[self.variables[1]], self.bins, self.axrange, weights = weights)
                self.hist += hist
        else:
            print("histogram was finalized already! No filling allowed!")

    def finalize_histogram(self):
        """
        Finalize the plot. The plot is finalized to make a probability density function.
        """
        if(not self.finalized):
            if(self.nvar ==1):
                self.areas = np.diff(self.edges[0])
            if(self.nvar==2):
                self.areas = np.outer(np.diff(self.edges[0]),np.diff(self.edges[1]))
            total = np.sum(self.hist)
            self.hist = self.hist/self.areas/total

            
    def draw_plot(self, ax):
        """
        Draw the plot. 
        :ax: matplotlib axes to draw the plot on
        """
        if(self.nvar==1):
            ax.stairs(self.hist,self.edges[0])
            ax.set_xlabel(self.variables[0])

        if(self.nvar==2):
            ax.pcolormesh(self.edges[0], self.edges[1], self.hist.T)
            ax.set_xlabel(self.variables[0])
            ax.set_ylabel(self.variables[1])
            

    def draw_interval(self, ax):
        """
        Draw the intervals. To be improved
        :ax: matplotlib axes to draw the plot on
        """
        if(self.nvar==1):
            ax.stairs(self.hist,self.edges[0], color='black')
            for lev in self.levels:
                ax.stairs(self.hist*np.less_equal(self.intervals,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
            ax.set_xlabel(self.variables[0])
        
        if(self.nvar==2):
            ax.pcolormesh(self.edges[0], self.edges[1], self.intervals.T)
            ax.contour(self.intervals.T, np.sort(self.levels))
            ax.set_xlabel(self.variables[0])
            ax.set_ylabel(self.variables[1])
        
            
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

        self.intervals = np.ones(self.hist.shape)

        total = np.sum(self.hist*self.areas)
        index = 0
        process_sum = 0
        if(self.nvar==1):
            self.intervals[index_sort_unrav[0][index]] = self.levels[nlev]
            while nlev >= 0:
                process_sum+=self.hist[index_sort_unrav[0][index]]*self.areas[index_sort_unrav[0][index]]
                index+=1
                if(process_sum/total < self.levels[nlev]):
                    self.intervals[index_sort_unrav[0][index]] = self.levels[nlev]
                else:
                    nlev-=1
                    if(nlev>=0):
                        self.intervals[index_sort_unrav[0][index]] = self.levels[nlev]
        if(self.nvar==2):       
            self.intervals[index_sort_unrav[0][index],index_sort_unrav[1][index]] = self.levels[nlev]
            while nlev >= 0:
                process_sum+=self.hist[index_sort_unrav[0][index],index_sort_unrav[1][index]]*self.areas[index_sort_unrav[0][index],index_sort_unrav[1][index]]
                index+=1
                if(process_sum/total < self.levels[nlev]):
                    self.intervals[index_sort_unrav[0][index],index_sort_unrav[1][index]] = self.levels[nlev]
                else:
                    nlev-=1
                    if(nlev>=0):
                        self.intervals[index_sort_unrav[0][index],index_sort_unrav[1][index]] = self.levels[nlev]
                        
        
        

