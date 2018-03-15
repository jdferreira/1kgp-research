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
        
        self.within[group] = Model.compute_stats(values)
    
    
    def add_between_values(self, group1, group2, values):
        """
        Insert the distance values between individuals of two distinct groups
        in this model
        """
        
        self.between[group1, group2] = Model.compute_stats(values)
    
    
    @staticmethod
    def compute_stats(values):
        """
        Return a dictionary of statistics that describe the given list of values.
        Namely, this dictionary contains:
        
            * the minimum with the key 'min'
            * the maximum with the key 'max'
            * the mean with the key 'mean'
            * the standard deviation with the key 'stdev'
        """
        
        return {
            'min':  min(values),
            'max':  max(values),
            'mean': np.mean(values),
            'stdev':  np.std(values),
        }
