import uproot
import numpy as np
from typing import Callable, Dict
from .variable import Variable

class MCMCSamples:
    compulsory_variables = [
            "DeltaCP",
            "Theta13",
            "Theta23",
            "Theta12",
            "Deltam2_32",
            "Deltam2_21",
            ]

    def __init__(self, _filepath, _treename):
        """
        Initialise the MCMCSamples object.

        :param _filepath: Path to the ROOT file.
        :param _treename: Name of the TTree inside of the ROOT file.
        """

        self.filepath = _filepath
        self.treename = _treename
        self.tree = uproot.open(f"{self.filepath}:{self.treename}")
        # Get self.variables
        self.variables = {}
        self.__check_and_extract_variables()
        # Get self.priors
        self.__check_and_extract_priors()

        print(f"Successfully initialised MCMCSamples object with ROOT file '{self.filepath}' and '{self.treename}' TTree inside")

    def __check_and_extract_variables(self):
        """
        Check if the compulsory variables exist, fill the "variables" map and
        create attributes for them.
        """

        keys = self.tree.keys()
        for var in self.compulsory_variables:
            if var in keys:
                # Create new variable
                self.variables[var] = Variable(var, lambda **kwargs: kwargs[var])
            else:
                raise ValueError(f"Compulsory variable '{var}' not found in the TTree.")

        self.__expand_variables()


    def __expand_variables(self):
        """
        Add extra variables that should be available to the user, that are not
        already in the input MCMC chain, e.g. Jarlskog-Invariant.
        """

        # Adding JarlskogInvariant
        def jarlskog_invariant(Theta12, Theta13, Theta23, DeltaCP):
            jarlskog = np.cos(Theta12) * np.pow(np.cos(Theta13), 2) * np.cos(Theta23) * np.sin(Theta12) * np.sin(Theta13) * np.sin(Theta23) * np.sin(DeltaCP)
            return jarlskog
        self.variables["JarlskogInvariant"] = Variable("JarlskogInvariant", jarlskog_invariant)

    def add_variable(self, name: str, function: Callable[..., np.ndarray]):
        """
        Add a custom variable to the MCMCSamples object. Required if want to
        plot a custom variable.

        :param name: Unique name of the custom variable.
        :param function: Function to compute the value fo the custom variable.
                         The function can accept keyqord arguments where keys match the names of
                         already defined variables.
        """
        if name in self.variables:
           raise ValueError(f"Variable '{name}' is already defined!")

        self.variables[name] = Variable(name, function)

    def __check_and_extract_priors(self):
        """
        Check if the default priors exist for each of the compulsory variables,
        and fill the self.priors map. Fill the transform map here too?
        """
        pass

    def __repr__(self):
       """
       Provide information about the class so it can be printed
       """
       return f"MCMCSamples(filepath='{self.filepath}', treename='{self.treename}', tree='{self.tree}', variables='{self.variables}')"
