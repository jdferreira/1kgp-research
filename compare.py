#!/usr/bin/python3

"""
Given a VCF file (or a gzipped VSF file) containing the polymorphism data of
genomes, compare two individuals based on their polymorphisms. The result is the
number of said polymorphisms where the genome of the two individuals differ.

More than two individuals can be provided, in which case the tool will output
the comparison between any pair of individuals.

Several metrics can be used (but for now only the one explained above!)
TODO: Add new metrics
"""

import argparse
import gzip
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
        

def get_distances_by_group(population: Population, comparer):
    """
    Return a dictionary that associates a pair of groups in the provided
    population with the list of comparison values calculated between the
    individuals in those two groups.
    """
    
    result = {}
    
    pairs = itertools.combinations_with_replacement(population.groups(), 2)
    for group1, group2 in pairs:
        group1_ids = population.group_to_individuals[group1]
        group2_ids = population.group_to_individuals[group2]
        
        values = [
            comparer.compare(id1, id2)
            for id1, id2 in itertools.product(group1_ids, group2_ids)
            if group1 != group2 or id1 != id2
            # We need to ensure tht, in case the two groups are the same,
            # then the individual is not compared with itself
        ]
        
        result[group1, group2] = values
    
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


def compare_test_individuals(test_population, train_population, comparer):
    """
    Compare test individuals with the train individuals, and return a dictionary
    where each key is a pair of (test individual, group) and the values are the
    statistics of the comparison values 
    """
    
    result = {}
    
    for individual in test_population.individuals():
        result[individual] = {}
        for group, group_individuals in train_population.items():
            values = [comparer.compare(individual, i) for i in group_individuals]
            stats = Model.compute_stats(values)
            result[individual][group] = stats
    
    return result


def print_table(table):
    """
    Print the dictionary that is created with the `get_distances_by_group`
    function
    """
    
    for (group1, group2), values in table.items():
        print(f'{group1}\t{group2}\t{values!r}')


def print_model(model):
    """
    Print the details of a model to the stardard output following a
    non-specified reference
    """
    
    print('# Model')
    print('\t'.join(('Group pair', 'Min', 'Max', 'Mean', 'Stdev')))
    
    print('# Distances within groups')
    for group, stats in model.within.items():
        print(group + '/' + group + '\t' + format_stats(stats))
    
    print('# Distances between groups')
    for group_pair, stats in model.between.items():
        group1, group2 = group_pair
        print(group1 + '/' + group2 + '\t' + format_stats(stats))


def print_test_distances(test_distances: dict, test_population: Population):
    """
    Print the details of the comparison between the test and the train
    individuals following a non-specified reference
    """
    
    print('# Test individuals')
    print('\t'.join(('Individual/group', 'Min', 'Max', 'Mean', 'Stdev')))
    for individual, group_stats in test_distances.items():
        expected = test_population.individual_to_group[individual]
        print(f'# {individual} -- expected group {expected}')
        for group, stats in group_stats.items():
            print(individual + '/' + group + '\t' + format_stats(stats))


def format_stats(stats):
    """
    Return a tab-separated line with the given statistics 
    """
    
    fields = (
        f'{stats["min"]}',
        f'{stats["max"]}',
        f'{stats["mean"]:.2f}',
        f'{stats["stdev"]:.2f}',
    )
    return '\t'.join(fields)


def get_arguments():
    """
    Creates an argument parser that can handle command line arguments and
    returns the parsed arguments
    """
    
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
        '-c', '--comparer', choices=all_comparers, default='default',
        help='The method used to compare two individuals. Right now the only '
             'valid options are "default" and "random".'
    )
    
    subparsers = parser.add_subparsers(
        title='Subcommands',
        dest='action',
        metavar='ACTION'
    )
    subparsers.required = True
    
    predict_parser = subparsers.add_parser(
        'predict',
        help='Predict the population of a set of individuals'
    )
    predict_parser.add_argument(
        'testing_set', metavar='TEST_FILE', type=argparse.FileType('rt'),
        help='The population file containing the individuals to test the model'
    )
    
    subparsers.add_parser(
        'table',
        help='Output a table with the comparison values of all possible pairs '
             'of individuals in the training dataset'
    )
    
    args = parser.parse_args()
    
    # We use argparse's machinery to detect that the input VCF file is valid
    # but we actaully want its name, not the stream, as this is what we need
    # to give to the gzip library
    args.vcf_file.close()
    args.vcf_file = args.vcf_file.name
    
    return args


def main():
    """
    Processes the command line arguments and performs the comparisons
    """
    
    args = get_arguments()
    
    if args.action == 'predict':
        make_prediction(args)
    elif args.action == 'table':
        make_table(args)


def make_prediction(args):
    """
    Perform the comparisons in the training set and report on the trained model
    and the possible classification of the provided test individuals.
    """
    
    # Read the populations
    training_set = parse_population(args.training_set)
    testing_set = parse_population(args.testing_set)
    
    # Instantiate the comparer with the correct set of individuals
    comparer = all_comparers[args.comparer]()
    assert isinstance(comparer, vcf.comparer.Comparer) # pylint hint!
    
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
    
    # Compare the test individuals with all the trainig individuals
    test_distances = compare_test_individuals(testing_set, training_set, comparer)
    
    print_model(model)
    print()
    print_test_distances(test_distances, testing_set)


def make_table(args):
    """
    Perform the comparisons between all possible pairs of individuals and
    print those values
    """
    
    # Read the populations
    population = parse_population(args.training_set)
    
    # Instantiate the comparer with the correct set of individuals
    comparer = all_comparers[args.comparer]()
    assert isinstance(comparer, vcf.comparer.Comparer) # pylint hint!
    
    comparer.set_individuals(list(population.individual_to_group))
    with gzip.open(args.vcf_file, 'rt') as stream:
        comparer.run(stream)
    
    # Compute the distances between every pair of individuals and report
    # the actual distance values within each group and between groups
    distances = get_distances_by_group(population, comparer)
    
    print_table(distances)


if __name__ == '__main__':
    main()
