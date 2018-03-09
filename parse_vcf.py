
import gzip

def parse(filename):
    """
    Parse a stream of bytes containing a gzipped vcf-file, and return a
    dictionary contianing the interesting bits
    """
    
    with gzip.open(filename, 'rt') as source:
        for line in source:
            if line.startswith('##'):
                continue
            
            fields = line.rstrip('\n').split('\t')
            
            if line.startswith('#'):
                # This is the column header; extract column names
                
                # The format column is optional. If it exists, the number of
                # meta keys is 9; otherwise, it is 8
                if 'FORMAT' in fields:
                    first_individual = 9
                else:
                    first_individual = 8
                
                meta_keys = fields[:first_individual]
                identifiers = fields[first_individual:]
                
                meta_values = []
                data = {identifier: [] for identifier in identifiers}
            
            else:
                # This is the actual data, with polymorphism information
                # on the first columns and individuals' data on the rest
                meta_values.append(fields[:first_individual])
                
                for identifier, value in zip(identifiers, fields[first_individual:]):
                    data[identifier].append(value)
    
    return {
        'meta_keys': meta_keys,
        'meta_values': meta_values,
        'data': data
    }
