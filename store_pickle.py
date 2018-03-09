#!/usr/bin/python3

"""
Parses a VCF file from in the data/ directory and constructs a dictionary in the
format

    data = {identifier: [value]}

where `identifier` is a string identifying each individual in the original file
and `[value]` is a list of integers where each integer corresponds to a
polymorphism of this individual. The list is sorted in the same order as the
original file, thus comparison of lists is meaningful.

The dictionary is saved in a 'pickle' file.
"""

import argparse
import pickle
import parse_vcf


def main():
    """
    Parse the file(s) specified with the command line arguments and save them
    in a more available format (using python's pickle format)
    """
    
    parser = argparse.ArgumentParser(
        description='Parse the specified file(s) from the 1000 genome '
                    'project and save the data in a pickle file'
    )
    parser.add_argument(
        'filenames', metavar='FILE', nargs='+',
        help='The name of the file to parse'
    )
    
    args = parser.parse_args()
    
    for filename in args.filenames:
        obj = parse_vcf.parse(filename)
        with open(filename + '.pickle', 'wb') as target:
            pickle.dump(obj, target)


if __name__ == '__main__':
    main()
