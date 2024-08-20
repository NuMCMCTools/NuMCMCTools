import uproot

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
                print(f"Variable {var} exists in the chain")
                # TODO: Create the variable
                # self.variables[var] = Variable(XYZ)
                # setattr(self, var, self.variables[var])
            else:
                raise ValueError(f"Compulsory variable '{var}' not found in the TTree.")

    def _check_and_extract_priors(self):
        """
        Check if the default priors exist for each of the compulsory variables,
        and fill the self.priors map. Fill the transform map here too?
        """
