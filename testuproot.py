#!/usr/bin/env python3

import uproot
import PlotStack
import numpy as np
import matplotlib.pyplot as plt


tree = uproot.open("testchaindata.root:mcmc")
stack = PlotStack.PlotStack(tree)

stack.add_plot(["Deltam2_32", "Theta23"],[],[50,50],[[2.2E-3,2.7E-3],[0.7,0.9]])
#stack.add_plot(["DeltaCP"],[],100,[-3.14159,3.14159])

stack.add_plot(["DeltaCP"],[],[-np.pi, -0.9*np.pi, -0.8*np.pi, -0.7*np.pi, -0.6*np.pi, -0.5*np.pi,
                                -0.4*np.pi, -0.3*np.pi, -0.2*np.pi, -0.1*np.pi, 0, 0.2*np.pi, 0.4*np.pi, 0.6*np.pi, 0.8*np.pi, np.pi])
stack.fill_plots()

stack.plots[0].draw_plot()
stack.plots[1].draw_plot()

stack.plots[0].make_intervals([0.68,0.95])
stack.plots[0].draw_interval()

stack.plots[1].make_intervals([0.68,0.95])
stack.plots[1].draw_interval()
