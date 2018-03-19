"""
TODO
"""

import string
import random
import itertools
from collections import defaultdict

import vcf.handler

__all__ = [
    'DefaultComparer',
    'RandomComparer',
    'ALL_COMPARERS',
]

def extract_genotype(variant):
    """
    Return exactly the GT data for a particular VCF field
    """
    
    # Get the genotype of this individual
    genotype = variant.split(':')[0]
    
    # Split for diploid (or higher ploidy) chromosomes
    if '/' in genotype:
        alleles = genotype.split('/')
    elif '|' in genotype:
        alleles = genotype.split('|')
    else:
        alleles = [genotype]
    
    # Remove extra spaces surrounding the alleles
    return [i.strip() for i in alleles]

class Comparer(vcf.handler.Handler):
    """
    Bare bone implementation of a `Handler` that returns a comparison of two
    individuals based on the information in one VCF file
    """
    
    def __init__(self):
        super().__init__()
        self.input_individuals = None
        self.individual_indeces = {}
    
    def set_individuals(self, input_individuals):
        """
        Inform this comparer that the individuals to compare are the ones in
        the provided list. This method must run *before* the `self.run` method
        is executed.
        """
        
        self.input_individuals = input_individuals
        self.individual_indeces = {}
    
    def process_individuals(self, individuals):
        for idx, individual in enumerate(individuals):
            if individual in self.input_individuals:
                self.individual_indeces[individual] = idx
        
        if len(self.input_individuals) > len(self.individual_indeces):
            # Find the first that is not present
            for individual in self.input_individuals:
                if individual not in individuals:
                    raise Exception(
                        f'Individual {individual} not present in the VCF file'
                    )
    
    def process_data(self, fixed_fields, format_field, data_fields):
        for id1, id2 in itertools.combinations(self.individual_indeces, 2):
            variant1 = data_fields[self.individual_indeces[id1]]
            variant2 = data_fields[self.individual_indeces[id2]]
            
            self.update_comparison(fixed_fields, format_field, id1, id2, variant1, variant2)
    
    def update_comparison(self, fixed_fields, format_field, id1, id2, variant1, variant2):
        """
        Update the state of this comparer based on this piece of VCF information
        where `fixed_fields` is the list of metadata about the particular VCF
        being compared (the first 8 columns of the line), `format_field` is
        either the FORMAT field for this line or an empty string, if the column
        does not exist, `id1` and `id2` are identifier strings of the two
        individuals being compared, and `variant1` and `variant2` are the data
        for this VCF for these two individuals.
        """
        
        raise NotImplementedError
    
    
    def compare(self, id1, id2):
        """
        Uses the data read from the VCF file and returns the comparison value
        between the two individuals
        """
        
        raise NotImplementedError


class DefaultComparer(Comparer):
    """
    This comparer assigns a difference `n` to a pair of individuals if the
    number of VCF's that are diffferent between the two is `n`.
    """
    
    def __init__(self):
        super().__init__()
        self.result = defaultdict(int)
    
    def process_meta(self, line):
        pass
    
    def update_comparison(self, fixed_fields, format_field, id1, id2, variant1, variant2):
        variant1 = extract_genotype(variant1)
        variant2 = extract_genotype(variant2)
        
        # For each individual, remove the non-mutated alleles, sort and remove
        # duplicates. Equality is defined here as whether the resulting sets
        # are equal
        variant1 = {i for i in variant1 if i != '0'}
        variant2 = {i for i in variant2 if i != '0'}
        
        if variant1 != variant2:
            value = 1
        else:
            value = 0
        
        self.result[id1, id2] += value
    
    def compare(self, id1, id2):
        if (id1, id2) in self.result:
            return self.result[id1, id2]
        elif (id2, id1) in self.result:
            return self.result[id2, id1]
        else:
            raise ValueError(
                f'Cannot find individuals {id1!r} or {id2!r} in the result'
            )


class RandomComparer(Comparer):
    """
    This comparer assigns a random comparison value to a pair of individuals.
    The actual mechanism is just 
    """
    
    def __init__(self):
        super().__init__()
        self.salt = ''.join(random.choices(string.ascii_uppercase, k=10))
    
    def process_meta(self, line):
        pass
    
    def update_comparison(self, fixed_fields, format_field, id1, id2, variant1, variant2):
        pass
    
    def compare(self, id1, id2):
        hash1 = hash(self.salt + id1)
        hash2 = hash(self.salt + id2)
        
        result = 0
        for char1, char2 in zip(hash1, hash2):
            result += abs(ord(char1) - ord(char2))
        
        return result


ALL_COMPARERS = {
    'default': DefaultComparer,
    'random': RandomComparer,
}

# def main():
#     """
#     Test the functionality of this module
#     """
    
#     import gzip
    
#     # handler = vcf.Printer()
#     handler = DefaultComparer()
#     handler.set_individuals(['HG00114', 'HG00115', 'HG00116', 'HG00145'])
    
#     filename = r'data\ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz'
#     with gzip.open(filename, 'rt') as stream:
#         handler.run(stream)
    
#     from pprint import pprint
#     pprint(handler.result)

# if __name__ == '__main__':
#     main()
