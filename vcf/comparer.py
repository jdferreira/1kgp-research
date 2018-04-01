"""
TODO
"""

import string
import random
import itertools

import vcf.handler

__all__ = [
    'DefaultComparer',
    'RandomComparer',
    'ALL_COMPARERS',
]

class Comparer(vcf.handler.Handler):
    """
    Bare bone implementation of a `Handler` that returns a comparison of two
    individuals based on the information in one VCF file
    """
    
    def __init__(self, init_factory=int):
        """
        Initialize a comparer by configuring the initial value that all pairs of
        individuals will start with. This is given not as an actual value but as
        a 0-argument function that produces a new value each time, just like the
        default value in `collections.defaultdict`.
        """
        
        super().__init__()
        self.input_individuals = None
        self.individual_indeces = {}
        self.pairs = []
        self.idx_pairs = []
        self.values = []
        self.init_factory = init_factory
    
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
        
        # Each line read from now on will be used to update the comparison
        # of all pairs of individuals of interest. To speed things up, we need
        # to be able to quickly find all interesting pairs and compare them
        self.pairs = list(itertools.combinations(self.individual_indeces.keys(), 2))
        self.idx_pairs = list(itertools.combinations(self.individual_indeces.values(), 2))
        self.values = [self.init_factory() for _ in self.pairs]
    
    def process_variant(self, fixed_fields, format_field, variant_field):
        """
        For a given variant value, transform a variant field value into
        something that can be directly compared (to be used by the
        `update_comparison` method method. This transformation is useful if that
        method needs to perform a transformation, but ensures that the
        transformation is only performed once per field, instead of every time
        the value is needed to update the comparison values.
        
        `fixed_fields` is the list of metadata fields in this line (the first 8
        columns), `format_field` is either the FORMAT field for this line or an
        empty string, if the column does not exist, and `variant_field` is the
        actual data to transform.
        """
        
        #pylint: disable=W0613,R0201
        
        # By default, do not do any transformation
        return variant_field
    
    def process_data(self, fixed_fields, format_field, data_fields):
        # Start by transforming the data where it matters
        for idx in self.individual_indeces.values():
            data_fields[idx] = self.process_variant(
                fixed_fields,
                format_field,
                data_fields[idx]
            )
        
        # For each pair to be compared, update the comparison value
        for i, idx_pair in enumerate(self.idx_pairs):
            idx1, idx2 = idx_pair
            
            self.update_comparison(i, data_fields[idx1], data_fields[idx2])
    
    def update_comparison(self, pair_idx, variant1, variant2):
        """
        Update the state of this comparer based on this piece of VCF
        information. Each line of the VCF file is used to update the _current_
        version of the comparison value for all pairs of interest. The pairs of
        interest are actually computed with the `set_individuals` method. In
        `update_comparison`, a pair of variants `(variant1, variant2)` is used
        to update the comparison value that currently lives in
        `self.value[pair_idx]`. As such, any implementation of this method must
        use the values of the arguments `variant1` and `variant2` and then
        somehow write to `self.value[pair_idx]`. The arguments passed on to this
        method are the result of applying `process_variant` to the actual
        variant fields in the VCF file.
        """
        
        raise NotImplementedError
    
    
    def terminate(self):
        """
        Produce a dictionary of all pairs of individuals and place there the
        values computed from the VCF file, and then release some resources.
        """
        
        new_values = {}
        for i, pair in enumerate(self.pairs):
            id1, id2 = pair
            new_values[id1, id2] = new_values[id2, id1] = self.values[i]
        self.values = new_values
        
        # Release resources
        self.pairs = None
        self.idx_pairs = None
    
    def compare(self, id1, id2):
        """
        Uses the data read from the VCF file and returns the comparison value
        between the two individuals
        """
        
        return self.values[id1, id2]


class DefaultComparer(Comparer):
    """
    This comparer assigns a difference `n` to a pair of individuals if the
    number of VCF's that are diffferent between the two is `n`.
    """
    
    def process_meta(self, line):
        pass
    
    def process_variant(self, fixed_fields, format_field, variant_field):
        """
        Return exactly the GT data for a particular VCF field
        """
        
        # Get the genotype of this individual
        genotype = variant_field.split(':')[0]
        
        # Split for diploid (or higher ploidy) chromosomes
        if '/' in genotype:
            alleles = genotype.split('/')
        elif '|' in genotype:
            alleles = genotype.split('|')
        else:
            alleles = [genotype]
        
        # Remove extra spaces surrounding the alleles
        return [i.strip() for i in alleles]
    
    def update_comparison(self, pair_idx, variant1, variant2):
        # For each individual, remove the non-mutated alleles, sort and remove
        # duplicates. Equality is defined here as whether the resulting sets
        # are equal
        variant1 = {i for i in variant1 if i != '0'}
        variant2 = {i for i in variant2 if i != '0'}
        
        if variant1 != variant2:
            value = 1
        else:
            value = 0
        
        self.values[pair_idx] += value


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
    
    def update_comparison(self, pair_idx, variant1, variant2):
        pass
    
    def compare(self, id1, id2):
        import sys
        mask = (1 << sys.hash_info.width) - 1
        
        hash1 = hex(hash(self.salt + id1) & mask)
        hash2 = hex(hash(self.salt + id2) & mask)
        
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
