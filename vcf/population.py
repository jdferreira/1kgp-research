# class Individual:
    
#     def __init__(self, identifier, group):
#         self.identifier = identifier
#         self.group = group

from collections import defaultdict

class Population:
    """
    This class represents a population, which is a list of individuals along
    with the group they belong to
    """
    
    def __init__(self):
        self.individual_to_group = {}
        self.group_to_individuals = defaultdict(list)
    
    def add_individual(self, identifier, group):
        """
        Add a new individual associated with the given group. Individual
        identifiers must be unique and can only belong to one group.
        """
        
        self.individual_to_group[identifier] = group
        self.group_to_individuals[group].append(identifier)
    
    def individuals(self):
        """
        Return a list of the individuals in this population
        """
        
        return list(self.individual_to_group)
    
    def groups(self):
        """
        Return a list of the groups in this population
        """
        
        return list(self.group_to_individuals)
    
    def items(self):
        """
        Return an item-like iterator that yields the various groups in this
        population. Each item is of the form
            (group_name, list_of_individuals)
        """
        
        return self.group_to_individuals.items()
    
    def has_individual(self, individual):
        """
        Determines whether this population contains a specific individual.
        """
        
        return individual in self.individual_to_group
