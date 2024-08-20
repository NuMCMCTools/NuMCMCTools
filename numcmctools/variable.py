import numpy as np
from typing import Callable, Dict

class Variable:
    def __init__(self, name: str, function: Callable[..., np.ndarray]):
        """
        Initialise a Variable instance.

        :param name: Unique name of the variable. It is used for
                     identification and should match the name used in the plot parameter
                     names.
        :param function: Function that transforms the input values data to
                         produce some transformed variable's values. Should accept keyword
                         arguments with names matching the keys in the input MCMC chain.
        """
        self.name = name
        self.function = function

    def evaluate(self, data: Dict[str, np.ndarray]) -> np.ndarray:
        """
        Evaluate the variable using provided values.

        :param data: A dictionary where keys are variable names and values are
                     numpy arrays of the corresponding data.
        :return: A numpy array of the computed values vased on the
                 self.function
        """
        # Get the argunent names from the transform function
        func_code = self.function.__code__
        arg_names = func_code.co_varnames[:func_code.co_argcount]

        # Get keyword args for the transform function
        #kwargs = {}
        #for arg in arg_names:
        #    if arg in data:
        #        kwargs[arg] = data[arg]
        #    else:
        #        raise KeyError(f"Missing required argument '{arg}' in the input data. Use MCMCSamples.add_variable!")

        kwargs = {arg : data[arg] for arg in arg_names if arg in data}

        # Call the transform function!
        try:
            return self.function(**kwargs)
        except Exception as e:
            raise RuntimeError(f"Error evaluating variable '{self.name}': {e}")

    def __repr__(self):
       return f"Variable(name='{self.name}', function='{self.function}')"
