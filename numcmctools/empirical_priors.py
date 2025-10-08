import uproot
import numpy as np
import logging
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator
from types import FunctionType

logger = logging.getLogger(__name__)

class EmpiricalPrior:
    _INTERPOLATORS = ['regular', 'linear']

    def __init__(self, root_obj, variables, interpolator_type: str ='linear'):
        type_name = type(root_obj).__name__
        logger.debug(f"Initializing EmpiricalPrior with input root object type: {type_name}")
        if "THnT" in type_name: 
            self._init_thnd(root_obj, interpolator_type)
        elif "TGraph2D" in type_name: 
            self._init_tgraph2d(root_obj, interpolator_type)
        else:
            raise ValueError(f"The provided object is not a valid THnT or TGraph2D object. Got {type_name}.")

        if interpolator_type not in self._INTERPOLATORS:
            raise ValueError(f"Invalid interpolator type. Choose from {self._INTERPOLATORS}. Got {interpolator_type}.") 
        
        self._init_interpolator(interpolator_type)
        logger.debug(f"EmpiricalPrior object initialized with {self.dimensions} dimensions "\
                    f"and {self.interpolate.__class__.__name__} interpolation.")

        self.variables = variables
    
    def _init_thnd(self, thnd, interpolator_type):
        # Throw if the number of dimensions is greater than 2
        self.dimensions: int = thnd.member("fNdimensions")
        if self.dimensions > 2:
            raise ValueError("External empirical prior with more than 2 dimensions are not supported yet.")
        logger.debug(f"Initializing THnD-empirical prior with {self.dimensions} dimensions.")

        # Get the axes and their properties
        axes =  thnd.member(f"fAxes")
        shape  = tuple(int(ax.member("fNbins")) for ax in axes)
        mins = tuple(ax.member("fXmin") for ax in axes)
        maxs = tuple(ax.member("fXmax") for ax in axes)

        # Expand the mins and maxes to take the underflow and overflow bins into the account,
        # and move the bins to the center
        step =  tuple((mx - m) / (2*s) for m, mx, s in zip(mins, maxs, shape))
        mins = tuple(m - step[i] for i, m in enumerate(mins))
        maxs = tuple(mx + step[i] for i, mx in enumerate(maxs))
        shape = tuple(s + 2 for s in shape)

        bin_centers = [np.linspace(mn, mx, s) for mn, mx, s in zip(mins, maxs, shape)]

        # Load the data from the THnT object
        data = np.array(thnd.member("fArray").member("fData"), dtype=float).reshape(shape)

        if interpolator_type == 'linear':
            # Cartesian product of bin centers
            mesh = np.meshgrid(*bin_centers, indexing="ij")
            data = data.flatten()                                   # (npoints,)
            points = np.stack([m.flatten() for m in mesh], axis=-1) # (npoints, ndim)
            points_shape = points.shape
        elif interpolator_type == 'regular':
            points = tuple(bin_centers) # ndim arrays
            points_shape = [len(p) for p in points]
        else:
            raise ValueError(f"Unknown interpolator_type: {interpolator_type}")

        self.data = data
        self.points = points
        logger.debug(f"Points shape: {points_shape}, z data shape: {self.data.shape}")

    def _init_tgraph2d(self, tgraph2d, interpolator_type):
        ''' Initialize from a TGraph2D object with scattered data. 
            Only 'linear' interpolation is supported.'''
        
        self.dimensions = 2
    
        x_array = tgraph2d.member('fX')
        y_array = tgraph2d.member('fY')
        
        if interpolator_type == 'linear':
            self.data = np.array(tgraph2d.member("fZ"))
            self.points = np.column_stack((x_array, y_array))
        elif interpolator_type == 'regular':
            raise NotImplementedError("'regular' interpolation is not supported for TGraph2D-style scattered data.")
        else:
            raise ValueError(f"Unknown interpolator_type: {interpolator_type}")
        logger.debug(f"Points shape: {self.points.shape}, z data shape: {self.data.shape}")
    
    def _init_interpolator(self, interpolator_type: str):

        logger.debug(f"Initialising the interpolator of type: {interpolator_type}")
        if interpolator_type == 'linear':
                self.interpolate = LinearNDInterpolator(self.points, self.data, fill_value=0.0, rescale=True)
        elif interpolator_type == 'regular':
                self.interpolate = RegularGridInterpolator(self.points, self.data, bounds_error=False, fill_value=0.0)

        logger.debug(f"Interpolator initialized: {self.interpolate.__class__.__name__}")
    
    def __call__(self, data):
        for var in self.variables:
            if var not in data:
                raise KeyError(f"Variable '{var}' not found in the provided data. Available: {list(data.keys())}")
        
        # Extract the values for the variables
        values = np.asarray([data[var] for var in self.variables])
        return self.interpolate(values.T)
        

    def __repr__(self):
        return f"EmpiricalPrior(interp_type={self.interpolate.__class__.__name__}, " \
               f"dims={self.dimensions})"
