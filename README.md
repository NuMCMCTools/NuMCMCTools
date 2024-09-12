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

* DeltaCP, in radians, across any $2\pi$ interval
* Theta13, in radians, [0, $\pi$/2]
* Theta23,  in radians, [0, $\pi$/2]
* Theta12,  in radians, [0, $\pi$/2]
* Deltam2_32, in $\text{eV}^2$
* Deltam2_21, in $\text{eV}^2$

which correspond to the standard parameters of neutrino mixing; the
first four are the PMNS angles, and the last two are the mass
splittings. For more information on the parameterization see the [PDG
review on neutrino mixing](https://pdg.lbl.gov/2024/web/viewer.html?file=../reviews/rpp2024-rev-neutrino-mixing.pdf). 

Each parameter has some prior set by the original analyzers. The
format of this information is a TList containing a TNamed for each branch, 
which specifies the name of the branch and its prior. 

The priors are specified as either `Uniform` or `Gaussian(mean, sigma)`.
There is then a further specification as to which variable the functional form
applies to. For example, a specification of 

``` Theta23  Uniform:sin^2(Theta23)```

indicates that the prior for Theta23 is uniform in $\sin^2\theta_{23}$.

There may be additional information contained in the file; please see
the data release from the particular analysis for more detail.

In this code, the class `mcmcsamples` is responsible for loading
in the required information from the file. The constructor takes the
filepath for the input file and the name of the oscillation parameters
tree, and an optional argument for additional variables to be loaded. 

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
plotted as strings, any change in priors to be applied (see section below for details),
the bins and axis ranges for the plots, as defined in
[numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html)
or
[numpy.histogram2d](https://numpy.org/doc/stable/reference/generated/numpy.histogram2d.html),
and finally whether to treat the two mass orderings separately
(`True`) or together (`False`, default).

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
respectively, and can take an `Figure` argument from matplotlib. When
dividing a figure into subfigures, use the `subfigures` command, rather
than `subplots`, as `subplots` is used internally. 

The plot stack can also be used to automatically draw all of the plots
in the stack, choosing the optimal dimensions for the subfigure array.

An example of the plotting features is in `examples/testuproot.py` and
can be run as follows: 

Example that loads example chain and makes example plots
 (`examples/testuproot.py`) and can be run as follows:
 ```
 python -m examples.testuproot
```
The output of this example should look as follows

<img width="2869" alt="testuproot" src="https://github.com/user-attachments/assets/db68f411-2a06-4646-83c1-db8da26b5b41">


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

<img width="2904" alt="custom_variables" src="https://github.com/user-attachments/assets/3cfb6af4-e188-49c7-8b68-2ab15f2586ac">

## Changing Priors

Priors for parameters can be changed automatically when making plots.
A list of new priors is passed when adding a new `plot` to the `plotstack`.
This feature currently only works in 1D and with Uniform or Gaussian priors
in a limited set of transformations of variables. The list of currently
available transformations is:

```
x
sin(x)
sin^2(x)
cos(x)
cos^2(x)
2x
sin(2x)
sin^2(2x)
cos(2x)
cos^2(2x)
exp(-ix)
```
For example, an input chain that has a prior set to be uniform in DeltaCP can
have the following passed to the input chain to make plots with a prior uniform
in $\sin\delta_{CP}$:

```
["Uniform:sin(DeltaCP)"]
```
Examples of this in practice are in the `custom_variables` example. 
