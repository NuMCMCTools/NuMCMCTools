#!/usr/bin/env python3
import uproot
import os
import numpy as np
import matplotlib.pyplot as plt

from numcmctools import PlotStack, MCMCSamples

# Get the path of the test MCMC chain file
script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "testchaindata.root")

# Create the MCMCSamples object
samples = MCMCSamples(root_file_path, "mcmc")

def AbsSinDcp(DeltaCP) ->np.ndarray:
    return np.abs(np.sin(DeltaCP))

# Define custom variable: ssth23
def ssth23(Theta23) ->np.ndarray:
    return np.power(np.sin(Theta23), 2)

# Define custom variable: ss2th13
def ss2th13(Theta13) ->np.ndarray:
    return np.power(np.sin(2*Theta13), 2)

# Define custom variable: Ue3
def AbsUe3(Theta13, DeltaCP):
    return np.abs(np.sin(Theta13) * np.exp(-1j * DeltaCP))

def AbsDm2(Deltam2_32):
    return np.abs(Deltam2_32)

# Add the new variables to the MCMCSamples object
samples.add_variable("AbsSinDeltaCP", AbsSinDcp)
samples.add_variable("SinSqTheta23", ssth23)
samples.add_variable("SinSq2Theta13", ss2th13)
samples.add_variable("AbsUe3", AbsUe3)
samples.add_variable("AbsDm2_32",AbsDm2)

# Create stack of plots
stack = PlotStack(samples)

priors= ["Gaussian(0.55, 0.01):sin^2(Theta23)"]

priors_sindcp = ["Uniform:sin(DeltaCP)",
                 "Gaussian(0.55, 0.01):sin^2(Theta23)"]

priors_expdcp = ["Uniform:sin(DeltaCP)",
                 "Gaussian(0.55, 0.01):sin^2(Theta23)"]

# 1D plots
stack.add_plot(["AbsUe3"], priors_expdcp, 50, [0.1, 0.25], True)
stack.add_plot(["JarlskogInvariant"], priors_sindcp, 50, [-0.05, 0.05], True)
stack.add_plot(["SinSqTheta23"], priors, 50, [0.35, 0.65])
stack.add_plot(["SinSqTheta23"], priors, 50, [0.35, 0.65],True)
stack.add_plot(["AbsDm2_32"],priors,50,[2.2E-3,2.8E-3],True)

# 2D plots
stack.add_plot(["JarlskogInvariant", "AbsUe3"], priors_expdcp, [50, 50], [[-0.05, 0.05], [0.1, 0.25]],True)
stack.add_plot(["SinSqTheta23", "SinSq2Theta13"], priors, [50, 50], [[0.35, 0.65], [0.04, 0.15]])
stack.add_plot(["JarlskogInvariant", "AbsSinDeltaCP"],priors_sindcp,[50, 50], [[-0.05, 0.05], [0, 1]])
stack.add_plot(["SinSqTheta23", "AbsDm2_32"],priors,[50, 50], [[0.35, 0.65], [2.2E-3,2.8E-3]], True)

# Fill the plots
stack.fill_plots()
stack.make_intervals([0.68,0.95])
stack.draw_plots()
fig, sfig, axs = stack.draw_intervals()


#all the things on the plot can be called with get_children
#this changes the color of the 1D histograms in a MO split plot
plot0_no = axs[0][0].get_children()
plot0_no[0].set_edgecolor('red')
plot0_no[1].set_facecolor('red')
plot0_no[2].set_facecolor('red')

plot0_io = axs[0][1].get_children()
plot0_io[0].set_edgecolor('blue')
plot0_io[1].set_facecolor('blue')
plot0_io[2].set_facecolor('blue')

print(axs[6])

#For a 2D plot, you can remove the color histogram and just leave
#the contours and also change their color
plot6 = axs[6][0].get_children()
plot6[0].remove()
plot6[1].set_edgecolor('magenta')
#why not improve the axis labels while we're at it
axs[6][0].set_xlabel(r'$\sin^2\theta_{23}$')
axs[6][0].set_ylabel(r'$\sin^2 2\theta_{13}$')

plt.show()
