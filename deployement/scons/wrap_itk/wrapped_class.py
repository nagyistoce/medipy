import config
from type import Type
from utils import filter_dimensions

class WrappedClass(object):
    def __init__(self, name, wrap_method="", instantiations=None):
        self.name = name
        self.wrap_method = wrap_method
        self.instantiations = instantiations or []
    
    def _get_wrap_method(self):
        return self._wrap_method
    
    def _set_wrap_method(self, wrap_method):
        wrap_methods = ["", "POINTER", "POINTER_WITH_SUPERCLASS", 
            "POINTER_WITH_2_SUPERCLASSES", "EXPLICIT_SPECIALIZATION",
            "POINTER_WITH_EXPLICIT_SPECIALIZATION", "ENUM", "AUTOPOINTER"]
        if wrap_method not in wrap_methods :
            raise Exception("Invalid wrap method '{0}'. "
                            "Possible values are {1}.".format(
                            wrap_method,", ".join("'{0}'".format(x) for x in wrap_methods)))
        self._wrap_method = wrap_method
    
    wrap_method = property(_get_wrap_method, _set_wrap_method)

def ImageFilter(types, parameters_count, dimensions_constraint=""):
    if dimensions_constraint :
        dimensions = filter_dimensions(dimensions_constraint)
    else :
        dimensions = config.dimensions
    
    instantiations = []
    
    for t in types :
        for d in dimensions :
            instantiation = []
            for i in range(parameters_count) :
                instantiation.append(Type("itk::Image", [t, d]))
            instantiations.append(instantiation)
    
    return instantiations