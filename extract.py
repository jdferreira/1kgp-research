#!/usr/bin/python3

"""
Extracts a file from the FTP server containing the 1000 genome project files
"""

import argparse
import ftplib
import os
import pathlib

HOST = 'ftp.1000genomes.ebi.ac.uk'
REMOTE_DIR = '/vol1/ftp/release/20130502/'


def extract(filenames):
    """
    Connect to the FTP server and retrieve the file(s) specified
    with the command line arguments
    """
    
    ftp = ftplib.FTP(HOST)
    ftp.login()
    ftp.cwd(REMOTE_DIR)
    
    for filename in filenames:
        destination = os.path.join('data', filename)
        with open(destination, 'wb') as f:
            ftp.retrbinary('RETR ' + filename, f.write)


def create_data_dir():
    """
    Creates the 'data/' directory if it does not yet exists
    """

    pathlib.Path('data').mkdir(exist_ok=True)


def main():
    """
    Connect to the FTP server and retrieve the file(s) specified
    with the command line arguments
    """
    
    parser = argparse.ArgumentParser(
        description='Extract the specified file(s) from the 1000 genome '
                    'project'
    )
    parser.add_argument(
        'filenames', metavar='FILE', nargs='+',
        help='The name of the file to download'
    )
    
    args = parser.parse_args()
    
    if any('/' in name for name in args.filenames):
        first_bad = [name for name in args.filenames if '/' in name][0]
        msg = f'File names cannot contain slashes: {first_bad!r}'
        parser.error(msg)
    
    create_data_dir()
    extract(args.filenames)


if __name__ == '__main__':
    main()
