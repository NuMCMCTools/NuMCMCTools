#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

from numcmctools import PlotStack, MCMCSamples

# Define custom variable: ssth23
def ssth23(Theta23) ->np.ndarray:
    return  np.power(np.sin(Theta23), 2)

# Define custom variable: ss2th13
def ss2th13(Theta13) ->np.ndarray:
    return np.power(np.sin(2*Theta13), 2)

# Define custom variable: ss2th12
def ss2th12(Theta12) ->np.ndarray:
    return np.power(np.sin(2*Theta12), 2)

# Define custom variable: abs(dm32) (scaled)
def absdm32(Deltam2_32) ->np.ndarray:
    return np.abs(Deltam2_32 * 1e3)

# Define custom variable: abs(dm21) (scaled)
def scaleddm21(Deltam2_21) ->np.ndarray:
    return np.abs(Deltam2_21 * 1e5)

def main(file: str, chain_name: str):
  # Create the MCMCSamples object
  samples = MCMCSamples(file, chain_name)

  variable_names = ["DeltaCP_02pi", "SinSq2Theta13", "SinSqTheta23", "SinSq2Theta12", "AbsDm32", "Dm21"]
  labels = [r"$\delta_{CP}$", r"$\sin^{2}(2\theta_{13})$", r"$\sin^{2}(\theta_{23})$", r"$\sin^{2}(2\theta_{12})$", r"$|\Delta m^{2}_{32}|$ ($\times10^{-3}$eV)", r"$|\Delta m^{2}_{21}|$ ($\times10^{-5}$eV)"]
  bins = [50, 50, 50, 50, 50, 50]
  ranges = [[0, 2*np.pi], [0.03, 0.16], [0.35, 0.65], [0.78, 0.91], [2.15, 2.65], [6, 9]]

  priors = ["Uniform:DeltaCP",
            "Gaussian(0.085,0.003):sin^2(2Theta13)",
            "Uniform:sin^2(2Theta12)",
            "Uniform:sin^2(Theta23)"]

  # Add the new variables to the MCMCSamples object
  samples.add_variable("SinSqTheta23", ssth23)
  samples.add_variable("SinSq2Theta13", ss2th13)
  samples.add_variable("SinSq2Theta12", ss2th12)
  samples.add_variable("AbsDm32", absdm32)
  samples.add_variable("Dm21", scaleddm21)

  # Add all the 1D and 2D plots we want into the stack
  stack = PlotStack(samples)
  for i, pl_i in enumerate(variable_names):
    for j, pl_j in enumerate(variable_names):
      if i < j:
        continue
      if i == j:
        stack.add_plot([variable_names[i]], priors, bins[i], ranges[i], True)
      else:
        stack.add_plot([variable_names[j], variable_names[i]], priors, [bins[j], bins[i]], [ranges[j], ranges[i]], True)

  # Fill the plots
  stack.fill_plots()

  # Plot with pyplot
  n = len(variable_names)
  fig = plt.figure()

  # Create the NxN plot grid for the triangle plot
  gs = fig.add_gridspec(n, n, hspace=0, wspace=0)
  axes = gs.subplots(sharex='col')

  enum = 0
  for i, pl_i in enumerate(variable_names):
    for j, pl_j in enumerate(variable_names):
      # Don't show anything for the upper triangle half
      if i < j:
        axes[i, j].set_axis_off()
        continue
      
      # Get the axis and stack plot for easier use later
      ax = axes[i,j]
      plot = stack.plots[enum]

      # Make the 1 and 2 sigma credible intervals
      plot.make_intervals([0.6827,0.9545])
      plot.draw_interval(ax)

      # Set the labels for all the 1D and 2D plots
      if i == j:
           ax.set_xlabel(labels[i])
           ax.set_ylabel("Posterior probability")
      else:
           ax.set_xlabel(labels[j])
           ax.set_ylabel(labels[i])

      # Force the dcp axis to stay within 0--2pi
      if i == 0:
          ax.set_xlim(0, 2*np.pi)

      # Ticks stylistics
      ax.tick_params(bottom=True, top=True, left=True, right=True, direction='in')
      ax.minorticks_on()
      ax.tick_params(which='minor', bottom=True, top=True, left=True, right=True, direction='in')
      enum += 1

  # Only show lables for the outer plots
  for ax in fig.get_axes():
    ax.label_outer()
  
  # Save the plot
  fig.set_size_inches(10, 10)
  plt.savefig("triangle_noio.png", dpi=300)

def run():
  import argparse
  
  parser = argparse.ArgumentParser()

  parser.add_argument('-f', '--file', dest="file", required=True, help="Input root file with the MCMC samples")
  parser.add_argument('-c', '--chain', dest="chain", required=True, help="Chain TTree name in the MCMC file")

  args = parser.parse_args()
  main(str(args.file), str(args.chain))

if __name__ == "__main__":
    run()
