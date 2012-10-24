import os
import re

def get_instantiated_types(wrapitk_root):
    """ Return the types WrapITK was compiled for.
    """
    
    result = []
    
    config_file = open(os.path.join(wrapitk_root, "WrapITKConfig.cmake"))
    for line in config_file.readlines() :
        match = re.match(r"SET\(WRAP_[a-z_]+ ON CACHE BOOL "
                         "\"Wrap ([a-z ]+) type\"\)", line)
        if match :
            result.append(match.group(1))
    return result
    

def get_unsigned_int_types(wrapitk_root):
    """ Return the unsigned integer types for which WrapITK was compiled.
    """
    
    result = []
    for type in get_instantiated_types(wrapitk_root) :
        if type.startswith("unsigned ") :
           result.append(type)
    return result

def get_signed_int_types(wrapitk_root):
    """ Return the signed integer types for which WrapITK was compiled.
    """
    
    result = []
    for type in get_instantiated_types(wrapitk_root) :
        if type.startswith("signed ") :
           result.append(type)
    return result

def get_real_types(wrapitk_root):
    """ Return the real types for which WrapITK was compiled.
    """
    
    result = []
    for type in get_instantiated_types(wrapitk_root) :
       if type in ["float", "double"]:
           result.append(type)
    return result

def get_int_types(wrapitk_root):
    """ Return the integer types for which WrapITK was compiled.
    """
    
    return get_unsigned_int_types(wrapitk_root)+get_signed_int_types(wrapitk_root)

def get_scalar_types(wrapitk_root):
    """ Return the scalar types for which WrapITK was compiled.
    """
    
    return get_int_types(wrapitk_root)+get_real_types(wrapitk_root)

def get_dimensions(wrapitk_root):
    """ Return the dimensions for which WrapITK was compiled.
    """
    
    result = []
    
    config_file = open(os.path.join(wrapitk_root, "WrapITKConfig.cmake"))
    for line in config_file.readlines() :
        match = re.match(r"SET\(WRAP_ITK_DIMS \"([0-9;]+)\" CACHE STRING "
                         "\"dimensions available separated by semicolons "
                         "\(;\)\"\)", line)
        if match :
            value = match.group(1)
            dimensions = [int(x) for x in value.split(";")]
            result.extend(dimensions)
    return result
