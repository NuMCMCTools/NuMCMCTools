import numpy as np
from typing import Callable, Dict

class Variable:
    def __init__(self, _name: str, _function: Callable[[Dict[str, np.ndarray]], np.ndarray]):
        self.name = _name
        self.function = _function

    def evaluate(self, data: Dict[str, np.ndarray]) -> np.ndarray:
        return self.function(data)

    def __repr__(self):
       return f"Variable(name='{self.name}', function='{self.function}')"

