# NuMCMCTools
MCMC tools for neutrino oscillations

## File Format

The expected input format for this code is a root file containg at
least one tree. The tree can have any name. Inside the tree, there
will be at least six branches called:

* DeltaCP, in radians, across any 2pi interval
* Theta13, in radians, [0, pi/2]
* Theta23,  in radians, [0, pi/2]
* Theta12,  in radians, [0, pi/2]
* Deltam2_32, in eV^2
* Deltam2_21, in eV^2

which correspond to the standard parameters of neutrino mixing; the
first four are the PMNS angles, and the last two are the mass
splittings. For more information on the parameterization see the PDG
review on neutrino mixing. 

Each parameter has some prior set by the original analyzers. The
format of this information is TBD. 

There may be additional information contained in the file; please see
the data release from the particular analysis for more detail.

In this code, the class `mcmcsamples` is responsible for loading
in the required information from the file. The constructor takes the
filepath for the input file and the name of the oscillation parameters
tree, and an optional argument for additional variables to be loaded
(TO BE IMPLEMENTED). 

A very basic script for loading the example tree is
`examples/load_mcmc.py` and can be run as follows:

 ```
 python -m examples.load_mcmc
```



## Posterior Density Functions




## Contours


## Derived Variables



## Changing Priors

Gotta do some work, gotta change some priors


## Examples:

Loads example chain into MCMCSamples object (`examples/load_mcmc.py`):
 ```
 python -m examples.load_mcmc
```

Example that loads example chain and makes example plots (`examples/testuproot.py`):
 ```
 python -m examples.testuproot
```

Example that loads example chain, creates custom Variables and makes example plots with them (`examples/custom_variables.py`):
 ```
 python -m examples.custom_variables
```
