__all__ = [
    'TwoWayDict',
]

class TwoWayDict(object):
    
    def __init__(self, items=None):
        if items is None:
            self.dict = {}
            return
        
        if isinstance(items, dict):
            items = items.items()
        
        self.dict = {normalize_key(key): value for key, value in items}
    
    def __setitem__(self, key, value):
        self.dict[normalize_key(key)] = value
    
    def __getitem__(self, key):
        return self.dict[normalize_key(key)]
    
def normalize_key(key):
    key1, key2 = key
    if key1 < key2:
        return (key1, key2)
    else:
        return (key2, key1)
