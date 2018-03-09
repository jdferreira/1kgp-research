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

import sys

import argparse
import itertools
import pickle
import numpy as np
import comparers
import parse_vcf

from collections import defaultdict
from two_way_dict import TwoWayDict


def parse_individuals(file):
    identifiers = {}
    test_set = {}
    superpopulations = defaultdict(list)
    
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
        
        identifier, superpopulation = fields
        if superpopulation[0] == '?':
            # This individual has an unknown superpopulation
            # (or is to be used as a test individual)
            superpopulation = superpopulation[1:]
            if not superpopulation:
                superpopulation = None
            test_set[identifier] = superpopulation
        
        else:
            identifiers[identifier] = {
                'superpopulation': superpopulation
            }
            superpopulations[superpopulation].append(identifier)
    
    if not test_set:
        # No unknown individuals. Is this an error for this study?
        print(
            'WARNING: The individuals file does not contain any individuals '
            'with an unknown superpopulation',
            file=sys.stderr
        )
    
    return {
        'identifiers': identifiers,
        'test_set': test_set,
        'superpopulations': superpopulations,
    }

def get_all_distances(vcf, individuals, compare):
    """
    Compute and return the distance values between all possible pairs of
    individuals, given a comparison function 
    """    
    
    pairs = itertools.combinations(individuals['identifiers'], 2)
    return TwoWayDict(
        {(id1, id2): compare(vcf, id1, id2) for id1, id2 in pairs}
    )


def get_test_distances(vcf, individuals, compare):
    """
    Compute and return the distance values between the individuals with
    an unknown superpopulation and the other ones
    """
    
    pairs = itertools.product(individuals['test_set'], individuals['identifiers'])
    return {
        (id1, id2): compare(vcf, id1, id2) for id1, id2 in pairs
    }


def make_model(distances, individuals):
    """
    For each superpopulation, find the internal distances and compute
    mean and standard deviation.
    """
    
    model = TwoWayDict()
    
    superpopulations = individuals['superpopulations']
    for superpopulation, identifiers in superpopulations.items():
        values = []
        for id1, id2 in itertools.combinations(identifiers, 2):
            values.append(distances[id1, id2])
        
        model[superpopulation, superpopulation] = {
            'min':  min(values),
            'max':  max(values),
            'mean': np.mean(values),
            'std':  np.std(values),
        }
    
    for sup1, sup2 in itertools.combinations(superpopulations, 2):
        sup1_ids = superpopulations[sup1]
        sup2_ids = superpopulations[sup2]
        values = []
        for id1, id2 in itertools.product(sup1_ids, sup2_ids):
            values.append(distances[id1, id2])
        
        model[sup1, sup2] = {
            'min':  min(values),
            'max':  max(values),
            'mean': np.mean(values),
            'std':  np.std(values),
        }
    
    return model


def main():
    """
    Processes the command line arguments and performs the comparisons
    """
    
    parser = argparse.ArgumentParser(
        description='Output the comparison values of the genomes of several '
                    'individuals'
    )
    parser.add_argument(
        'file', metavar='FILE', type=argparse.FileType('rb'),
        help='The pickle file where the polymorphism data is stored'
    )
    parser.add_argument(
        'individuals', metavar='IDENTIFIERS_FILE', type=argparse.FileType('rt'),
        help='The file containing the identifiers of the individuals to compare'
    )
    parser.add_argument(
        '-c', '--comparer', choices=comparers.ALL, default='default',
        help='The method used to compare two individuals. Right now the only '
             'valid option is "default".'
    )
    
    args = parser.parse_args()
    
    # Load the pickle file containing the dictionary where each key is an
    # individual and each value is the list of their polymorphisms. Two
    # polymorphisms are equal if their numbers are equal.
    vcf = parse_vcf.parse(args.file)
    
    # Parse the individuals file
    individuals = parse_individuals(args.individuals)
    
    # Grab the actual comparer function
    compare = comparers.ALL[args.comparer]
    
    # Compute distances (all vs. all)
    distances = get_all_distances(vcf, individuals, compare)
    
    # Compute distances within superpopulations and between different
    # superpopulations in order to create a model for prediction
    model = make_model(distances, individuals)
    
    # Compare the test individuals with all the known
    test_distances = get_test_distances(vcf, individuals, compare)
    
    from pprint import pprint
    pprint(model.dict)
    pprint(test_distances)
    


if __name__ == '__main__':
    main()
