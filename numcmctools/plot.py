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
            self.hist, edges = np.histogram([], self.bins, self.axrange, weights=[])
            self.edges.append(edges)
            if(self.mo_option):
                self.hist_no = np.zeros(np.shape(self.hist))
                self.hist_io = np.zeros(np.shape(self.hist))
        elif(self.nvar==2):
            self.hist, edgesx, edgesy = np.histogram2d([], [],self.bins, self.axrange, weights=[])
            self.edges.append(edgesx)
            self.edges.append(edgesy)
        else:
            print("too many or too few variables!")
            return
                
    def fill_plot(self, data, weights=None):
        """
        Fill the plot. If the plot has been finalized, no more filling is allowed.

        :data: An array of the data
        :weights: Data weights
        """

        if weights==None:
            weights = np.ones(np.shape(data[self.variables[0]]))
        
        if not self.finalized:
            if(self.nvar==1):
                if(self.mo_option):
                    hist, edges = np.histogram(data[self.variables[0]], self.bins, self.axrange, weights = weights*np.greater_equal(data["Deltam2_32"],0)*np.ones(np.shape(data["Deltam2_32"])))
                    self.hist_no += hist
                    hist, edges = np.histogram(data[self.variables[0]], self.bins, self.axrange, weights = weights*np.less_equal(data["Deltam2_32"],0)*np.ones(np.shape(data["Deltam2_32"])))
                    self.hist_io += hist
                else:
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

            if(self.nvar==1):
                self.areas = np.diff(self.edges[0])
                if(self.mo_option):
                    self.hist = np.concatenate((self.hist_no,self.hist_io))
                    self.areas = np.concatenate((self.areas,self.areas))
            if(self.nvar==2):
                self.areas = np.outer(np.diff(self.edges[0]),np.diff(self.edges[1]))
            total = np.sum(self.hist)
            self.hist = self.hist/self.areas/total
            if self.mo_option:
                self.hist_no = self.hist[:len(self.hist_no)]
                self.hist_io = self.hist[len(self.hist_no):]
            self.finalized = True
            
    def draw_plot(self, sfig):
        """
        Draw the plot. 
        :ax: matplotlib axes to draw the plot on
        """
        if self.mo_option:
            ax = sfig.subplots(1,2, sharey=True)
        else:
            ax = sfig.subplots(1,1)
        
        if(self.nvar==1):
            if self.mo_option:
                ax[0].stairs(self.hist_no,self.edges[0])
                ax[0].set_xlabel(self.variables[0]+" NO")
                ax[1].stairs(self.hist_io,self.edges[0])
                ax[1].set_xlabel(self.variables[0]+" IO")
                                
            else:
                ax.stairs(self.hist,self.edges[0])
                ax.set_xlabel(self.variables[0])

        if(self.nvar==2):
            cm = ax.pcolormesh(self.edges[0], self.edges[1], self.hist.T)
            ax.set_xlabel(self.variables[0])
            ax.set_ylabel(self.variables[1])

        sfig.subplots_adjust(wspace=0, hspace=0)

    def draw_interval(self, sfig):
        """
        Draw the intervals. To be improved
        :ax: matplotlib axes to draw the plot on
        """

        if self.mo_option:
            ax = sfig.subplots(1,2, sharey=True)
        else:
            ax = sfig.subplots(1,1)
        
        #need to put in a check if the intervals have been calculated
        if(self.nvar==1):
            if self.mo_option:
                ax[0].stairs(self.hist_no,self.edges[0], color='black')
                for lev in self.prob_levels:
                    ax[0].stairs(self.hist_no*np.greater_equal(self.hist_no,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
                ax[0].set_xlabel(self.variables[0]+" NO")
                
                ax[1].stairs(self.hist_io,self.edges[0], color='black')
                for lev in self.prob_levels:
                    ax[1].stairs(self.hist_io*np.greater_equal(self.hist_io,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
                ax[1].set_xlabel(self.variables[0]+" IO")
            else:
                ax.stairs(self.hist,self.edges[0], color='black')
                for lev in self.prob_levels:
                    ax.stairs(self.hist*np.greater_equal(self.hist,lev),self.edges[0], fill=True, color='grey', alpha=0.3)
                ax.set_xlabel(self.variables[0])
        
        if(self.nvar==2):

            linestyles_base = ['solid','dashed','dashdot','dotted']
            linestyles=['']*len(self.prob_levels)
            for i in range(len(self.prob_levels)):
                linestyles[i] = linestyles_base[i%len(linestyles_base)]
            linestyles = list(reversed(linestyles))

            cm = ax.pcolormesh(self.edges[0], self.edges[1], self.hist.T)
            ax.contour(0.5*(self.edges[0][:-1]+self.edges[0][1:]), 0.5*(self.edges[1][:-1]+self.edges[1][1:]),self.hist.T, np.sort(self.prob_levels), linestyles=linestyles, colors='lightgrey')
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
                        
        
        

