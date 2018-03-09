__all__ = [
    'default_compare',
    'random_compare',
    'ALL',
]

def extract_gt(variant):
    return variant.split(':')[0]


def default_compare(vcf, id1, id2):
    data = vcf['data']
    
    if id1 not in data:
        raise Exception(f'{id1!r} is not a known individual')
    if id2 not in data:
        raise Exception(f'{id2!r} is not a known individual')
    
    counter = 0
    for variant1, variant2 in zip(data[id1], data[id2]):
        # We simply compare the alternative value associated with the GT format
        # which is always the first one
        if extract_gt(variant1) != extract_gt(variant2):
            counter += 1
    
    return counter


def random_compare(_vcf, id1, id2):
    return int(hash(id1) < hash(id2))


ALL = {
    'default': default_compare,
    'random': random_compare,
}
