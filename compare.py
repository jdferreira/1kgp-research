#!/usr/bin/python3

"""
Given a pickle file containing the polymorphism data of genomes, compare
two individuals based on their polymorphisms. The result is the number
of said polymorphisms where the genome of the two individuals differ.

More than two individuals can be provided, in which case the tool will output
the comparison between any pair of individuals.

Several metrics can be used (but for now only the one explained above!)
TODO: Add new metrics
"""

import itertools

from vcf.comparer import ALL_COMPARERS as all_comparers
from vcf.population import Population
from vcf.model import Model

import vcf.comparer

def parse_population(file) -> Population:
    """
    Read a population file, where each non-empty line has two space-separated
    fields: an individual identifier and the group they belong to.
    """
    
    result = Population()
    
    for line in file:
        line = line.rstrip('\n')
        
        # Remove possible comments
        comment_start = line.find('#')
        if comment_start >= 0:
            line = line[:comment_start]
        line = line.strip()
        
        if not line:
            # Ignore empty lines
            continue
        
        fields = line.split()
        if len(fields) != 2:
            raise Exception(
                f'Line {line!r} is invalid: needs 2 fields, found {len(fields)}'
            )
        
        identifier, group = fields
        result.add_individual(identifier, group)
    
    return result
        

def make_model(population, comparer):
    """
    Create a model that can be used to classify individuals in groups based on
    the distances between them and the individuals in the groups.
    """
    
    model = Model()
    
    for group, identifiers in population.items():
        values = [
            comparer.compare(id1, id2)
            for id1, id2 in itertools.combinations(identifiers, 2)
        ]
        
        model.add_within_values(group, values)
    
    for group1, group2 in itertools.combinations(population.groups(), 2):
        group1_ids = population.group_to_individuals[group1]
        group2_ids = population.group_to_individuals[group2]
        
        values = [
            comparer.compare(id1, id2)
            for id1, id2 in itertools.product(group1_ids, group2_ids)
        ]
        
        model.add_between_values(group1, group2, values)
    
    return model


def main():
    """
    Processes the command line arguments and performs the comparisons
    """
    
    import argparse
    import gzip
    
    parser = argparse.ArgumentParser(
        description='Output the comparison values of the genomes of several '
                    'individuals'
    )
    parser.add_argument(
        'vcf_file', metavar='VCF_FILE', type=argparse.FileType('rb'),
        help='The VCF file where the polymorphism data is stored'
    )
    parser.add_argument(
        'training_set', metavar='TRAIN_FILE', type=argparse.FileType('rt'),
        help='The population file containing the individuals to train the model'
    )
    parser.add_argument(
        'testing_set', metavar='TEST_FILE', type=argparse.FileType('rt'),
        help='The population file containing the individuals to test the model'
    )
    parser.add_argument(
        '-c', '--comparer', choices=all_comparers, default='default',
        help='The method used to compare two individuals. Right now the only '
             'valid options are "default" and "random".'
    )
    
    args = parser.parse_args()
    
    # We use argparse's machinery to detect that the input VCF file is valid
    # but we actaully want its name, not the stream, as this is what we need
    # to give to the gzip library
    args.vcf_file = args.vcf_file.name
    
    # Read the populations
    training_set = parse_population(args.training_set)
    testing_set = parse_population(args.testing_set)
    
    # Instantiate the comparer with the correct set of individuals
    comparer = all_comparers[args.comparer]()
    assert isinstance(comparer, vcf.comparer.Comparer)
    
    all_individual_identifiers = (
        list(training_set.individual_to_group) +
        list(testing_set.individual_to_group)
    )
    comparer.set_individuals(all_individual_identifiers)
    with gzip.open(args.vcf_file, 'rt') as stream:
        comparer.run(stream)
    
    # Compute distances within groups and between different groups
    # and create a model for classification
    model = make_model(training_set, comparer)
    
    # # Compare the test individuals with all the known
    # test_distances = get_test_distances(vcf, individuals, compare)
    
    from pprint import pprint
    pprint(model.within)
    pprint(model.between.dict)


if __name__ == '__main__':
    main()
