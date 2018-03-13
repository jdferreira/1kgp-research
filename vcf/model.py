"""
TODO
"""

import numpy as np

from two_way_dict import TwoWayDict

class Model:
    """
    Represents a model that can be used to classify an individual based on its
    distance to individuals in the various groups of a population
    """
    
    def __init__(self):
        self.within = {}
        self.between = TwoWayDict()
    
    
    def add_within_values(self, group, values):
        """
        Insert the distance values between individuals of the same group in
        this model
        """
        
        self.within[group] = {
            'min':  min(values),
            'max':  max(values),
            'mean': np.mean(values),
            'std':  np.std(values),
        }
    
    
    def add_between_values(self, group1, group2, values):
        """
        Insert the distance values between individuals of two distinct groups
        in this model
        """
        
        self.between[group1, group2] = {
            'min':  min(values),
            'max':  max(values),
            'mean': np.mean(values),
            'std':  np.std(values),
        }
