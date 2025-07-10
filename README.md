# NuMCMCTools

This package is intended as a lightweight tool to assist in analysis using 
public releases of MCMC chains from neutrino oscillation analyses. The
following experiments have released data compatible with this software:

* Experiment 1
* Experiment 2

## Installation

To install `numcmctools` with PyPI, run the following command:

```bash
pip install numcmctools
```

To install directly from source, first clone numcmctools with git:

```bash
git clone https://github.com/NuMCMCTools/NuMCMCTools.git
```

enter the new NuMCMCTool directory and simply run:

```bash
pip install .
```

The pacakge does not need to be installed to use it, you can simply `git clone`
the package and add your own macros that include the `numcmcmtools` folder.

## File Format

### MCMC samples
The expected input format for this code is a ROOT file containing at
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

### Priors
Each parameter has some prior set by the original analyzers. The format of this
information is a `TList` object named `priors`, containing a `TNamed` for each
branch, which specifies the name of the branch and its prior.

The priors are specified as:
1. `Uniform` 
2. `Gaussian(mean, sigma)`
4. `BimodalGaussian(mean1, sigma1, mean2, sigma2, bias)` where bias is in % and optional, default 50%.
3. `Step(bias, boundary)` where bias is in % and boundary is optional, default 0

There is then a further specification as to which variable the functional form
applies to. For example, a specification of

``` Theta23  Uniform:sin^2(Theta23)```

indicates that the prior for Theta23 is uniform in $\sin^2\theta_{23}$.

### Constraints

Optionally, file can contain 1D or 2D constraints to be applied by default to the plots when the `constraints` object input to `add_plot` function is set to `None`. For this, there needs to be a `TDirectoryFile` object named `constraints`. Each constraint object inside that directory can be either `TGraph2D` or `THnD`. The object must have a name in a format:

1. unique name
2. parameter names separated by colon, e.g. `sin^2(2Theta13):abs(Deltam2_32)`
3. Optional string "NO" for constraint only applied to Normal Mass Ordering, IO for the Inverted mass ordering, or none for a constraint applied across both Mass Orderings.
4. 1 if the constraint should be applied automatically, or 0 if not.

Example:
1. `DayaBay2D_2024_IO:sin^2(2Theta13):abs(Deltam2_32):IO:1` for a 2D constraint automatically applied in the Inverted Mass Ordering
1. `some_constraint_from_theory:sin^2(2Theta13):sin^2(Theta23):1` for a 2D constraint automatically applied across both mass orderings.

### Other


The file additionally contains a citation to the original analysis that
produced the chain. The citation must be inside of a `TObjString` object named
`citation`.

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

```bash
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
default. The number of plots to be filled does have an impact on the run
time of the code.

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

An example of the plotting features is in `examples/simpleplots.py` and
can be run as follows:

```bash
 python -m examples.simpleplots
```

The output of this example should look as follows

![simpleplots](https://github.com/user-attachments/assets/e79060c5-42f4-43b5-90b9-6e89aa0b7f0a)

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

```bash
 python -m examples.custom_variables
```

The output of this example should have two canvases which look as follows

![custom_variable_1](https://github.com/user-attachments/assets/a6243f95-7cc4-4e21-bd27-27777a404b6e)
![custom_variable_2](https://github.com/user-attachments/assets/f377d5fa-6154-41d7-b3aa-1a3822b2a95b)

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
cos^4(x)
2x
sin(2x)
sin^2(2x)
cos(2x)
cos^2(2x)
cos^4(2x)
exp(-ix)
exp(ix)
abs(ix)
```

For example, an input chain that has a prior set to be uniform in DeltaCP can
have the following passed to the input chain to make plots with a prior uniform
in $\sin\delta_{CP}$:

```
["Uniform:sin(DeltaCP)"]
```

Examples of this in practice are in the `custom_variables` example.

## Plotting

Plots can be drawn directly from the `plotstack` (shown in the `custom_variables` 
example), where the arrangement of pdfs or intervals is done algorithmically. 
For more control, the figure and its subdivisions can be arranged by the user
and each plot drawn individually. An example of the latter is shown in both the
`simpleplots` example, and in the `plot_triangle` example, which creates a standard
triangle (corner) plot with the mass ordering separated. This example can be run as 
follows:

```bash
  python -m examples.plot_triangle -f examples/testchaindata.root -c mcmc    
```

or 

```bash
  plot_triangle -f examples/testchaindata.root -c mcmc    
```
if installed throgh pip, which produces the output below.

![triangle_noio](https://github.com/user-attachments/assets/3f5940dc-cf99-4244-b82e-608c4112cd84)
