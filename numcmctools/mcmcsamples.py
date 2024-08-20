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
        self._check_and_extract_variables()
        # Get self.priors
        self._check_and_extract_priors()

        print(f"Successfully initialised MCMCSamples object with ROOT file '{self.filepath}' and '{self.treename}' TTree inside")

    def _check_and_extract_variables(self):
        """
        Check if the compulsory variables exist, fill the "variables" map and
        create attributes for them.
        """

        keys = self.tree.keys()
        for var in self.compulsory_variables:
            if var in keys:
                # Create new variable
                self.variables[var] = Variable(var, lambda data, var=var: data[var])
            else:
                raise ValueError(f"Compulsory variable '{var}' not found in the TTree.")

        # Adding JarlskogInvariant
        def jarlskog_invariant(data: Dict[str, np.ndarray]) ->np.ndarray:
            return np.sin(2 * data["Theta12"]) * np.sin(2 * data["Theta13"]), * np.sin(2 * data["Theta23"]) * sin(data["DeltaCP"])
        self.variables["JarlskogInvariant"] = Variable("JarlskogInvariant", jarlskog_invariant)

    def add_variable(self, _name: str, _function: Callable[[Dict[str, np.ndarray]], np.ndarray]):
        self.variables[_name] = Variable(_name, _function)

    def _check_and_extract_priors(self):
        """
        Check if the default priors exist for each of the compulsory variables,
        and fill the self.priors map. Fill the transform map here too?
        """

    def __repr__(self):
       """
       Provide information about the class so it can be printed
       """
       return f"MCMCSamples(filepath='{self.filepath}', treename='{self.treename}', tree='{self.tree}', variables='{self.variables}')"
