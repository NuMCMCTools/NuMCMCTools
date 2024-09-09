import sympy as sp
import re
from typing import Callable, Optional, Dict
import numpy as np

class JacobianGraph:

    # Define symbolic variable
    x = sp.Symbol('x')

    # Pre-defined priors
    variables = {
            'x' : x,
            'sin(x)' : sp.sin(x),
            'sin^2(x)': sp.sin(x)**2,
            'cos(x)': sp.cos(x),
            'cos^2(x)': sp.cos(x)**2,
            '2x' : 2*x,
            'sin(2x)' : sp.sin(2*x),
            'sin^2(2x)': sp.sin(2*x)**2,
            'cos(2x)': sp.cos(2*x),
            'cos^2(2x)': sp.cos(2*x)**2,
            'exp(-ix)': sp.exp(-sp.I * x)
            }

    # Possible distribution functions (not including Uniform)
    distribution_functions: Dict[str, Callable[..., Callable[[np.ndarray], np.ndarray]]] = {
        'Gaussian': lambda mean, std: lambda x: (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
    }

    def __init__(self):
        self.transform_graph: Dict[str, Callable[[np.ndarray], np.ndarray]] = {}
        self.jacobian_graph = self.__build_graph()

    def gaussian_prior(self, mean: float, std: float) -> Callable[[np.ndarray], np.ndarray]:
        """
        Returns a Gaussian prior function with the specified mean and standard deviation.

        :param mean: Mean (mu) of the Gaussian distribution.
        :param std: Standard deviation (sigma) of the Gaussian distribution.
        :returns: A Numpy-compatible function that represents the Gaussian prior.
        """
        def gaussian(x: np.ndarray) -> np.ndarray:
            return (1 / (std * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std) ** 2)
        return gaussian

    def __jacobian_chain_rule(self, source: sp.Expr, target: sp.Expr):
        """
        Attempts a derivative from the source variable to the target variable,
        creates the derivative's modulus numpified function as the Jacobian
        determinant from the parameter transformation. Only works in 1D.
        TODO: Implement 2D, or ND version of this

        :param source: Source variable for prior transformation
        :param target: Target variable for prior transformation
        :returns: Numpy-compatible function that computes the Jacobian determinant
        """
        try:
            # Differentiate the source function w.r.t. x
            dsource_dx = sp.diff(source, self.x)

            # Differentiate the target function w.r.t. x
            dtarget_dx = sp.diff(target, self.x)

            # Apply chain rule to get absolute value of target function w.r.t. source function
            jacobian = sp.Abs(dtarget_dx / dsource_dx)

            if jacobian == 0:
                print(f"Jacobian transformation from {source} to {target} not possible, got 0")
                return None

            # Return numpified Jacobian transform function
            return sp.lambdify(self.x, target, modules=['numpy']), sp.lambdify(self.x, jacobian, modules=['numpy'])

        except Exception as e:
            print(f"Jacobian transformation from {source} to {target} error: {e}")
            return None

    def __build_graph(self):
        """
        Builds the Jacobian transformation prior dictionary from uniform in one
        variable to uniform in another, with pre-defined priors.

        :returns: Jacobian transformation dictionary
        """
        jacobian_graph = {}
        for from_prior, expr_from in self.variables.items():
            jacobian_graph[from_prior] = {}
            for to_prior, expr_to in self.variables.items():
                # Calculate the transformation function and Jacobian if possible
                transform, jacobian = self.__jacobian_chain_rule(expr_from, expr_to)
                if to_prior not in self.transform_graph:
                    self.transform_graph[to_prior] = transform

                # Store in the jacobian graph
                jacobian_graph[from_prior][to_prior] = jacobian
        return jacobian_graph

    def __parse_prior_string(self, prior: str):
        """
        Parses the prior string and extracts the prior type, it's parameter
        values (if any) and the target distribution expression.

        :param prior: String representing the prior, e.g. Gaussian(0.085, 0.003):sin^2(2x)
        :returns: Distribution type as string, distribution parameters & target distribution expression
        """

        match = re.match(r'(\w+)\((.*?)\):(.+)', prior)
        if match:
            dist_type = match.group(1)
            params_str = match.group(2)
            expr = match.group(3).strip()

            # Extract parameters as floats
            params = [float(p) for p in params_str.split(',') if p]

            return dist_type, params, expr
        else:
            return 'Uniform', [], prior.split(':')[-1].strip()

    def get_jacobian_func(self, from_prior: str, to_prior: str):
        """
        Returns the Jacobian prior transformation function from prior uniform
        in one variable to prior uniform in another. The function includes
        reweights to a different target distribution, e.g. Gaussian.

        :param from_prior: variable to the prior transform from
        :param to_prior: variable to the prior transform to
        :returns: Jacobian transformation function
        """

        # Get the distribution types, dist. parameters & expressions
        from_dist_type, from_params, from_expr = self.__parse_prior_string(from_prior)
        to_dist_type, to_params, to_expr = self.__parse_prior_string(to_prior)

        # Error checking
        if from_dist_type is not 'Uniform':
            print(f"Warning: Transforming from a gaussian distribution not supported, here be dragons!")
        if from_expr not in self.jacobian_graph or to_expr not in self.jacobian_graph[from_expr]:
            raise ValueError(f"No direct transformation available for {from_prior} to {to_prior}")

        jacobian = self.jacobian_graph[from_expr][to_expr]

        # If transforming to uniform distribution, simply return the jacobian
        if to_dist_type == 'Uniform':
            return jacobian

        # If transforming to non-uniform (e.g. Gaussian), create new function for it
        def transformed_jacobian(values: np.ndarray) -> np.ndarray:
            # Make sure the distribution we target is available
            if to_dist_type not in self.distribution_functions:
                raise TypeError(f"Distribution {to_dist_type} not available!")

            # Get the jacobian transform
            jacobian_values = jacobian(values)

            # Get the target distribution function
            dist_func = self.distribution_functions[to_dist_type](*to_params)

            # Get the values from the target distribution
            # TODO: This is likely already being calculated in MCMCSamples class! Need to integrate this better
            dist_values = dist_func(self.transform_graph[to_expr](values))

            # Return the product of jacobian & target distribution weights
            return dist_values * jacobian_values

        # Return the jacobian distribution multiplied by non-Uniform target distribution
        return transformed_jacobian
