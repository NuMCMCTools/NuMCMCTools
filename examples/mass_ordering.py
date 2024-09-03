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

# Define custom variable: ssth23
def ssth23(Theta23) ->np.ndarray:
    return np.power(np.sin(Theta23), 2)

# Define custom variable: ss2th13
def ss2th13(Theta13) ->np.ndarray:
    return np.power(np.sin(2*Theta13), 2)

# Define custom variable: Ue3
def AbsUe3(Theta13, DeltaCP):
    return np.abs(np.sin(Theta13) * np.exp(-1j * DeltaCP))

# Add the new variables to the MCMCSamples object
samples.add_variable("SinSqTheta23", ssth23)
samples.add_variable("SinSq2Theta13", ss2th13)
samples.add_variable("AbsUe3", AbsUe3)

# Create stack of plots
stack = PlotStack(samples)

# 1D plots
#stack.add_plot(["AbsUe3"], [], 50, [0.1, 0.25], True)
stack.add_plot(["JarlskogInvariant"], [], 50, [-0.05, 0.05], True)
stack.add_plot(["SinSqTheta23"], [], 50, [0.35, 0.65])
stack.add_plot(["SinSqTheta23"], [], 50, [0.35, 0.65],True)

# 2D plots
stack.add_plot(["JarlskogInvariant", "AbsUe3"], [], [50, 50], [[-0.05, 0.05], [0.1, 0.25]],True)
#stack.add_plot(["SinSqTheta23", "SinSq2Theta13"], [], [50, 50], [[0.35, 0.65], [0.075, 0.095]])
#stack.add_plot(["JarlskogInvariant", "DeltaCP"],[],[50, 50], [[-0.05, 0.05], [-np.pi, np.pi]])

# Fill the plots
stack.fill_plots()

# Plot with pyplot
fig = plt.figure()
axs = fig.subfigures(2, 2)
stack.plots[0].draw_plot(axs[0,0])
stack.plots[1].draw_plot(axs[0,1])
stack.plots[2].draw_plot(axs[1,0])
stack.plots[3].draw_plot(axs[1,1])
#stack.plots[4].draw_plot(axs[1,1])
#stack.plots[5].draw_plot(axs[1,2])
#plt.show()

stack.make_intervals([0.68,0.95])

# Plot with pyplot
fig2 = plt.figure()
axs2 = fig2.subfigures(2, 2)
stack.plots[0].draw_interval(axs2[0,0])
stack.plots[1].draw_interval(axs2[0,1])
stack.plots[2].draw_interval(axs2[1,0])
stack.plots[3].draw_interval(axs2[1,1])
#stack.plots[4].draw_interval(axs2[1,1])
#stack.plots[5].draw_interval(axs2[1,2])
plt.show()
