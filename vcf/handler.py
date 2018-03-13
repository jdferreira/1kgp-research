"""
This module contains classes to deal with VCF files. The `Handler` is the base
class which does not perform any computation. There is a generic `Comparer`
class that returns a simple count of the number of VCFs that are different
between two individuals. More sophisticated methdos may be implemented based on
that example.
"""

__all__ = [
    'Handler',
    'Printer',
]

FIXED_FIELDS_LENGTH = 8

class Handler(object):
    """
    A class that defines how to handle a line of a VCF file. This class itself
    does not do anything other than pass through the lines with no effect.
    Subclasses can be created to properly perform any computation that can fit
    the VCF data model
    """
    
    def __init__(self):
        """
        Initializes this `Handler` with the provided stream, corresponding to
        a VCF file
        """
        
        self.has_genotype_field = None # Still undetermined
    
    
    def run(self, stream):
        """
        Run this handler, which will read the contents of the stream one line at
        a time, and will process metadata lines, column names and actual data
        lines in sequence, with specific methods being executed for each one
        """
        
        # Run any preparations that need to exist before the handler runs
        self.prepare()
        
        for line in stream:
            line = line.rstrip('\n')
            
            if line.startswith('##'):
                self.process_meta(line)
            else:
                fields = line.split('\t')
                
                # Do we have a genotype field?
                if line.startswith('#'):
                    self.has_genotype_field = 'FORMAT' in fields
                    if self.has_genotype_field:
                        individuals = fields[FIXED_FIELDS_LENGTH+1:]
                    else:
                        individuals = fields[FIXED_FIELDS_LENGTH:]
                    self.process_individuals(individuals)
                
                else:
                    fixed_fields = fields[:FIXED_FIELDS_LENGTH]
                    if self.has_genotype_field:
                        format_field = fields[FIXED_FIELDS_LENGTH]
                        data_fields = fields[FIXED_FIELDS_LENGTH+1:]
                    else:
                        format_field = ''
                        data_fields = fields[FIXED_FIELDS_LENGTH:]
                    
                    self.process_data(fixed_fields, format_field, data_fields)
        
        self.terminate()
    
    
    def prepare(self):
        """
        Perform any computations before starting to read the file
        """
        
        pass
    
    def process_meta(self, line):
        """
        Receives a line of meta information from the VCF file
        
        Arguments
        - line: str containing the exact line
        """
        
        raise NotImplementedError
    
    def process_individuals(self, individuals):
        """
        Receives a list of individuals in the VCF file.
        
        Arguments:
            individuals: list of str containing the identifiers in the file
        """
        
        raise NotImplementedError
    
    def process_data(self, fixed_fields, format_field, data_fields):
        """
        Receives a line of data from the VCF file
        
        Arguments:
            fixed_fields:
                The data for the fixed columns of the VCF file format
                (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
            format_field:
                The data for the format field, if it exists, or an
                empty string 
            data_fields:
                A list of str containing the actual data for the
                individuals in the file, in the order that was
                presented in the `self.process_individuals` method
        """
        
        raise NotImplementedError
    
    
    def terminate(self):
        """
        Perform any computations after reading and processing the last line in
        the file
        """
        
        pass


class Printer(Handler):
    """
    A `Handler` that prints the VCF file back to an output stream
    """
    
    def __init__(self, output):
        super().__init__()
        self.output = output
    
    def process_meta(self, line):
        print(line, file=self.output)
    
    def process_individuals(self, individuals):
        headers = [
            '#CHROM',
            'POS',
            'ID',
            'REF',
            'ALT',
            'QUAL',
            'FILTER',
            'INFO',
        ]
        if self.has_genotype_field:
            headers.append('FORMAT')
        
        print(*headers, sep='\t', end='\t', file=self.output)
        print(*individuals, sep='\t', file=self.output)
    
    def process_data(self, fixed_fields, format_field, data_fields):
        print(*fixed_fields, sep='\t', end='\t', file=self.output)
        if self.has_genotype_field:
            print(format_field, end='\t', file=self.output)
        print(*data_fields, sep='\t', file=self.output)


# def main():
#     """
#     Test the functionality of this module
#     """
    
#     import gzip
#     import sys
    
#     handler = Printer(sys.stdout)
    
#     filename = r'data\ALL.chrY.phase3_integrated_v2a.20130503.genotypes.vcf.gz'
#     with gzip.open(filename, 'rt') as stream:
#         handler.run(stream)

# if __name__ == '__main__':
#     main()
