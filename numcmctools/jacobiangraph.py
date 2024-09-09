import sympy as sp

class JacobianGraph:

    # Define symbolic variable
    x = sp.Symbol('x')

    # Pre-defined priors
    variables = {
            'Uniform:x' : x,
            'Uniform:sin(x)' : sp.sin(x),
            'Uniform:sin^2(x)': sp.sin(x)**2,
            'Uniform:cos(x)': sp.cos(x),
            'Uniform:cos^2(x)': sp.cos(x)**2,
            'Uniform:2x' : 2*x,
            'Uniform:sin(2x)' : sp.sin(2*x),
            'Uniform:sin^2(2x)': sp.sin(2*x)**2,
            'Uniform:cos(2x)': sp.cos(2*x),
            'Uniform:cos^2(2x)': sp.cos(2*x)**2,
            'Uniform:exp(-ix)': sp.exp(-sp.I * x)
            }

    def __init__(self):
        self.graph = self.__build_graph()

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
            return sp.lambdify(self.x, jacobian, modules=['numpy'])

        except Exception as e:
            print(f"Jacobian transformation from {source} to {target} error: {e}")
            return None

    def __build_graph(self):
        """
        Builds the Jacobian transformation prior dictionary from uniform in one
        variable to uniform in another, with pre-defined priors.

        :returns: Jacobian transformation dictionary
        """
        graph = {}
        for from_prior, expr_from in self.variables.items():
            graph[from_prior] = {}
            for to_prior, expr_to in self.variables.items():
                # Calculate the transformation function and Jacobian if possible
                jacobian = self.__jacobian_chain_rule(expr_from, expr_to)

                # Store in the graph
                graph[from_prior][to_prior] = jacobian
        return graph

    def get_jacobian_func(self, from_prior: str, to_prior: str):
        """
        Returns the Jacobian prior transformation function from prior uniform
        in one variable to prior uniform in another.

        :param from_prior: variable to the prior transform from
        :param to_prior: variable to the prior transform to
        :returns: Jacobian transformation function
        """
        if from_prior in self.graph and to_prior in self.graph[from_prior]:
            return self.graph[from_prior][to_prior]
        else:
            raise ValueError(f"No direct transformation available for {from_prior} to {to_prior}")
