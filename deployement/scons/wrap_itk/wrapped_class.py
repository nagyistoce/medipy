import config
from type import Type
from utils import filter_dimensions

class WrappedClass(object):
    def __init__(self, name, wrap_method="", instantiations=None, includes=None):
        self.name = name
        self.wrap_method = wrap_method
        self.instantiations = instantiations or []
        self.includes = includes or ["{0}.h".format(Type(name).mangled_name)]
    
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
    
    def _get_typedefs(self):
        typedefs = []
        for instantiation in self.instantiations :
            typedefs.append(Type(self.name, instantiation))

            if "2_SUPERCLASSES" in self.wrap_method :
                typedefs.append(Type(self.name, instantiation, "Superclass::Superclass"))
                typedefs.append(Type(self.name, instantiation, "Superclass::Superclass::Pointer"))
            
            if "SUPERCLASS" in self.wrap_method :
                typedefs.append(Type(self.name, instantiation, "Superclass"))
                typedefs.append(Type(self.name, instantiation, "Superclass::Pointer"))
            
            if "POINTER" in self.wrap_method :
                if self.wrap_method == "AUTOPOINTER" :
                    typedefs.append(Type(self.name, instantiation, "SelfAutoPointer"))
                else :
                    typedefs.append(Type(self.name, instantiation, "Pointer"))
        
        return typedefs
    
    wrap_method = property(_get_wrap_method, _set_wrap_method)
    typedefs = property(_get_typedefs)

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