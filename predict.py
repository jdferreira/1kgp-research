#!/usr/bin/python3

"""
Given a VCF file (or a gzipped VSF file) containing the polymorphism data of
genomes, predict the population of a set of individuals based on the
similarities of their genomes with the genomes of other individuals with known
population
"""

import argparse
import gzip
import random

from vcf.predictor import ALL_PREDICTORS as all_predictors
from vcf.population import Population, parse_population

import vcf.predictor

def get_arguments():
    """
    Creates an argument parser that can handle command line arguments and
    returns the parsed arguments
    """
    
    parser = argparse.ArgumentParser(
        description='Predict the population of individuals'
    )
    parser.add_argument(
        'vcf_file', metavar='VCF_FILE', type=argparse.FileType('rb'),
        help='The VCF file where the polymorphism data is stored'
    )
    parser.add_argument(
        'population', metavar='POPULATION', type=argparse.FileType('rt'),
        help='The population file containing the individuals to train the model'
    )
    parser.add_argument(
        '-i', '--individuals', metavar='INDIVIDUALS',
        type=argparse.FileType('rt'),
        help='The file containing the individuals to test the model. '
             'If not provided, 10% of the individuals in the population are '
             'used to test the model. See also -t.'
    )
    parser.add_argument(
        '-p', '--polymorphisms', metavar='POLYMORPHISMS',
        type=argparse.FileType('rt'),
        help='The file containing the polymorphism identifiers to base the '
             'prediciton on. Each line must be a single identifier'.
    )
    parser.add_argument(
        '-t', '--test-fraction', metavar='FRACTION', type=float, default=0.1,
        help='If no list of individuals is given, a fraction of the population '
             'is randomly selected to serve as testing data. This flag '
             'controls the fraction to use. Defaults to 0.1.'
    )
    parser.add_argument(
        '-p', '--predictor', choices=all_predictors, default='PredictorJoao',
        help='The method used to compare two individuals.'
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
    make_prediction(args)


def read_polymorphisms(file):
    result = []
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
        
        result.append(line)
    
    return result


def make_prediction(args):
    """
    Classify based on the given arguments.
    """
    
    # Read the populations
    population = parse_population(args.population)
    
    if args.individuals is None:
        individuals = Population()
        
        size = int(len(population.individuals()) * args.test_fraction)
        for identifier in random.sample(population.individuals(), size):
            group = population.individual_to_group[identifier]
            individuals.add_individual(identifier, group)
    else:
        individuals = parse_population(args.individuals)
    
    for i in individuals.individuals():
        population.remove_individual(i)
    
    # Instantiate the comparer with the correct set of individuals
    predictor: vcf.predictor.Predictor = all_predictors[args.predictor](0, 0)
    
    predictor.set_populations(population, individuals.individuals())
    
    if args.polymorphisms:
        polymorphisms = read_polymorphisms(args.polymorphisms)
        predictor.set_polymorphisms(polymorphisms)
    
    with gzip.open(args.vcf_file, 'rt') as stream:
        predictor.run(stream)
    
    print('\n'.join(f'{i} : {j} : {individuals.individual_to_group[i]}' for i, j in predictor.labels.items()))

if __name__ == '__main__':
    main()
