# NuMCMCTools

This package is intended as a lightweight tool to assist in analysis using 
public releases of MCMC chains from neutrino oscilliation analyses. The 
following experiments have released data compatible with this software:

* Experiment 1
* Experiment 2


## File Format

The expected input format for this code is a ROOT file containg at
least one TTree. The tree can have any name. Inside the tree, there
will be at least six branches called:

* DeltaCP, in radians, across any 2pi interval
* Theta13, in radians, [0, pi/2]
* Theta23,  in radians, [0, pi/2]
* Theta12,  in radians, [0, pi/2]
* Deltam2_32, in eV^2
* Deltam2_21, in eV^2

which correspond to the standard parameters of neutrino mixing; the
first four are the PMNS angles, and the last two are the mass
splittings. For more information on the parameterization see the [PDG
review on neutrino mixing](https://pdg.lbl.gov/2024/web/viewer.html?file=../reviews/rpp2024-rev-neutrino-mixing.pdf). 

Each parameter has some prior set by the original analyzers. The
format of this information is TBD. 

There may be additional information contained in the file; please see
the data release from the particular analysis for more detail.

In this code, the class `mcmcsamples` is responsible for loading
in the required information from the file. The constructor takes the
filepath for the input file and the name of the oscillation parameters
tree, and an optional argument for additional variables to be loaded
(TO BE IMPLEMENTED). 

An example tree (`examples/testchaindata.root`) with the required
features is provided for testing the functionality of this software. 

A very basic script for loading the example tree is
`examples/load_mcmc.py` and can be run as follows:

 ```
 python -m examples.load_mcmc
```



## Posterior Density Functions and Contours

After the tree is loaded, an instance of `plotstack` can be created,
which is initialized with an `mcmcsamples` instance. Individual
`plot`s can be added to the plot stack. One and two dimensional plots
are currently supported. 

A plot is added to the plot stack with the `add_plot` function, which
takes as arguments an array of the name(s) of the variables to be
plotted as strings, any change in priors to be applied (TO BE
IMPLEMENTED), the bins and axis ranges for the plots, as defined in
[numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)
or
[numpy.histogram2d](https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html),
and finally whether to treat the two mass orderings separately
(`True`) or together (`False`, default) (TO BE IMPLEMENTED).

Any number of plots can be added to the stack. When all plots are
created, they can be filled by calling `fill_plots`, which has two
optional arguments. The first causes only the first
`n_steps` to be filled from the tree; this is useful in testing the
plot stack for very large input trees; the second restricts the number
of steps read from the tree to be `batchsize` at any given time, which
is helpful for memory management. Above a batch size of approximately
100000, there is no performance impact on the time it takes to read
large trees and fill plots; this value is therefore set at the
default. 

After the histograms are filled, they are normalized to create
probability _density_ functions, such that individual bins are
normalized to the width/area of the bin and the integral is normalized
to 1. 

Credible intervals can be created by calling `make_intervals`
on either the plot stack or individual plots within the stack. The
function takes an array of credible interval levels. 

PDFs and intervals can be drawn with `draw_plot` and `draw_interval`
respectively, and can take an `Axes` argument from matplotlib. (TO BE
IMPLEMENTED: automatic figure/axes)

An example of the plotting features is in `examples/testuproot.py` and
can be run as follows: 

Example that loads example chain and makes example plots
 (`examples/testuproot.py`) and can be run as follows:
 ```
 python -m examples.testuproot
```
The output of this example should look as follows

![testuproot](https://github.com/user-attachments/assets/dff6d98e-4600-4b43-9425-3d2aad3d0bbc)


## Derived Variables

Variables derived from the six required parameters—for example,
$\sin^2\theta_{23}$-can be defined and used for plotting. Each new
variable must be given a unique name and defined as a function which
uses the required parameters as input variables, spelled exactly as
they are defined above. 

The new variable is then attached to an instance of the `mcmcsamples`
class using the `add_variable` function, which takes as arguments the
name of the new variable and the function defining it. 

An example of this feature is in `examples/custom_variables.py` and
can be run as follows:

 ```
 python -m examples.custom_variables
```
The output of this example should look as follows

![custom_variables](https://github.com/user-attachments/assets/38ca3621-eba9-400c-a6ac-b4497a7e56e8)

## Changing Priors

Gotta do some work, gotta change some priors
