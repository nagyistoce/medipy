import re

import config

def configure_file(source, destination, **kwargs) :
    data = open(source).read()
    
    for key, value in kwargs.items() :
        data = data.replace("@{0}@".format(key), value)
    
    open(destination, "w").write(data)

def filter_dimensions(dimensions_constraint):
    dimensions = []
    
    string_constraint = (isinstance(dimensions_constraint, basestring) and 
                         re.match(r"^[0-9]+\+$", dimensions_constraint))
    
    if string_constraint :
        min_dim = int(re.sub(r"^([0-9]+)\+$", r"\1", dimensions_constraint))
        max_disallowed = min_dim-1
        var_name = ""
        for d in config.dimensions :
            if d >= min_dim :
                dimensions.append(d)
    else :
        # The condition is just a list of dims. Return the intersection of these
        # dims with the selected ones.
        s1 = set(config.dimensions)
        s2 = set(dimensions_constraint)
        dimensions = list(s1.intersection(s2))
    
    return dimensions