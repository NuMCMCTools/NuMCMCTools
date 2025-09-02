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
            self._init_thnd(root_obj)
        elif "TGraph2D" in type_name:
            self._init_tgraph2d(root_obj)
        else:
            raise ValueError(f"The provided object is not a valid THnT or TGraph2D object. Got {type_name}.")

        if interpolator_type not in self._INTERPOLATORS:
            raise ValueError(f"Invalid interpolator type. Choose from {self._INTERPOLATORS}. Got {interpolator_type}.") 
        
        self._init_interpolator(interpolator_type)
        logger.debug(f"EmpiricalPrior object initialized with {self.dimensions} dimensions "\
                    f"and {self.interpolate.__class__.__name__} interpolation.")

        self.variables = variables
    
    def _init_thnd(self, thnd):
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

        # Expand the mins and maxes to take the underflow and overflow bins intothe account, 
        # and move the bins to the center
        step =  tuple((mx - m) / (2*s) for m, mx, s in zip(mins, maxs, shape))
        mins = tuple(m - step[i] for i, m in enumerate(mins))
        maxs = tuple(mx + step[i] for i, mx in enumerate(maxs))
        shape = tuple(s + 2 for s in shape)

        # Load the data from the THnT object
        self.data: np.ndarray = np.array(thnd.member("fArray").member("fData"), dtype=float)
        self.data = self.data.reshape(shape)

        # Create the bins
        self.points = np.asarray(tuple(np.linspace(m, mx, s) for m, mx, s in zip(mins, maxs, shape)))
        logger.debug(f"Points shape: {self.points.shape}, z data shape: {self.data.shape}")
    
    def _init_tgraph2d(self, tgraph2d):
        self.dimensions = 2

        x_increasing = not np.all(np.diff(tgraph2d.member('fX')) >= 0)
        y_increasing = not np.all(np.diff(tgraph2d.member('fY')) >= 0)
        
        if x_increasing and y_increasing:
            raise ValueError("The X and Y axis of the TGraph2D are not increasing.")

        points = []
        for axis, is_increasing in zip(['fX', 'fY'], [x_increasing, y_increasing]):
            if is_increasing:
                npoints = np.where(np.diff(tgraph2d.member(axis)) < 0)[0] + 1
                points.append(tgraph2d.member(axis)[0:npoints[0]])
            else:
                npoints = np.where(np.diff(tgraph2d.member(axis)) > 0)[0] + 1
                points.append(tgraph2d.member(axis)[::npoints[0]])

        self.data = np.array(tgraph2d.member("fZ"))
        self.data = self.data.reshape(int(np.sqrt(len(self.data))), int(np.sqrt(len(self.data))))
        self.points = np.asarray(points)

        logger.debug(f"Points shape: {self.points.shape}, z data shape: {self.data.shape}")
    
    def _init_interpolator(self, interpolator_type: str):

        logger.debug(f"Initialising the interpolator of type: {interpolator_type}")
        if interpolator_type == 'linear':
            # Linear interpolation requires a meshgrid of points for the Delaunay triangulation
            self.points = np.column_stack([p.ravel() for p in np.meshgrid(*self.points)])
            self.interpolate = LinearNDInterpolator(self.points,
                                                    np.ravel(self.data.T),
                                                    fill_value=0.0,
                                                    rescale=True)

        elif interpolator_type == 'regular':
            self.interpolate = RegularGridInterpolator(self.points,
                                                       self.data,
                                                       method='linear',
                                                       bounds_error=False,
                                                       fill_value=0.0)
        
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