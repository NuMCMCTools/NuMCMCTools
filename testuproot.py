#!/usr/bin/env python3

import uproot
import PlotStack
import numpy as np
import matplotlib.pyplot as plt


tree = uproot.open("testchaindata.root:mcmc")
stack = PlotStack.PlotStack(tree)

stack.add_plot(["Deltam2_32", "Theta23"],[],[50,50],[[2.2E-3,2.7E-3],[0.7,0.9]])
stack.add_plot(["DeltaCP"],[],100,[-3.14159,3.14159])

stack.fill_plots()

stack.plots[0].draw_plot()
stack.plots[1].draw_plot()

#stack.plots[0].make_intervals([0.68,0.95])
#stack.plots[0].draw_interval()

#stack.plots[1].make_intervals([0.68,0.95])
#stack.plots[1].draw_interval()
