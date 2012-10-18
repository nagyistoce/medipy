import itkTypes

class Instantiation(object) :
    """ Instantiation of a template type.
    
        >>> rgb_pixel = Instantiation("itk::RGBPixel", "float")
        >>> image = Instantiation("itk::Image", rgb_pixel, 3)
        >>> itk_filter = Instantiation("itk::BETImageFilter", image, image)
    """
    
    def __init__(self, class_name, *parameters):
        self.class_name = class_name
        self.parameters = parameters
    
    def _get_cpp_name(self):
        if self.parameters :
            parameters = []
            for parameter in self.parameters :
                if hasattr(parameter, "cpp_name") :
                    parameters.append(parameter.cpp_name)
                else :
                    parameters.append(str(parameter))
            return "{0}< {1} >".format(self.class_name, ",".join(parameters))
        else :
            return self.class_name
    
    def _get_mangled_name(self):
        if self.parameters :
            parameters = []
            for parameter in self.parameters :
                if hasattr(parameter, "mangled_name") :
                    parameters.append(parameter.mangled_name)
                else :
                    mangled_parameter = mangled_type_map.get(parameter, str(parameter))
                    parameters.append(mangled_parameter)
            mangled_class_name = mangled_type_map.get(
                self.class_name, get_swig_class_name(self.class_name))
            return "{0}{1}".format(mangled_class_name, "".join(parameters))
        else : 
            return mangled_type_map.get(self.class_name, self.class_name)
    
    cpp_name = property(_get_cpp_name)
    mangled_name = property(_get_mangled_name)

def get_includes(instantiation):
    """ Return a list of include files for the given instantiation.
    """
    
    includes = set()
    
    if isinstance(instantiation, Instantiation) :
        # Try to guess the include of the class name
        if "::" in instantiation.class_name :
            namespace, class_name = instantiation.class_name.split("::", 1)
            if namespace == "itk" :
                include = "itk{0}".format(instantiation.class_name.split("::")[-1])
            elif namespace == "std" :
                include = class_name
            else :
                raise Exception("No rule for namespace {0}".format(namespace))
        includes.add(include)
        # Recursively add the parameters' includes
        for parameter in instantiation.parameters :
            includes.update(get_includes(parameter))
    # Otherwise do nothing
    
    return includes

def get_typedefs(instantiation, pointer):
    """ Return a list of pair for typedef declarations.
    
        pointer can be None, "pointer" or "pointer_with_superclass", following
        WrapITK's semantics (Documentation/Guide.txt) :
            If POINTER is passed, then the class and the typedef'd class::Pointer
            type is wrapped. (Class::Pointer had better be a SmartPointer
            instantiation, or things won't work. This is always the case for
            ITK-style code.) If POINTER_WITH_SUPERCLASS is provided, then
            class::Pointer, class::Superclass and class::Superclass::Pointer are
            all wrapped. (Again, this only works for ITK-style code where the
            class has a typedef'd Superclass, and the superclass has Self and
            Pointer typedefs.)
    """
    
    typedefs = []
    
    cpp_name, mangled_name = instantiation.cpp_name, instantiation.mangled_name 
    
    typedefs.append((cpp_name, mangled_name))
        
    if pointer is not None :
        if "pointer" in pointer :
            typedefs.append(("{0}::Pointer::SmartPointer".format(cpp_name),
                            "{0}_Pointer".format(mangled_name)))
        if "superclass" in pointer :
            typedefs.append(("{0}::Superclass::Self".format(cpp_name), 
                            "{0}_Superclass".format(mangled_name)))
        if pointer == "pointer_with_superclass" :
            typedefs.append(("{0}::Superclass::Pointer::SmartPointer".format(cpp_name), 
                            "{0}_Superclass_Pointer".format(mangled_name)))
    
    return typedefs

def get_swig_class_name(class_name) :
    """ Return the Swig name of a C++ class by deleting the scope separator ("::")
        operators from the name
        
        >>> get_swig_class_name("itk::Numerics::MyMarvelousOptimizer")
        'itkNumericsMyMarvelousOptimizer'
    """
    # Split the class name according to its namespaces
    class_name_elements = class_name.split("::")
    swig_class_name = "".join(class_name_elements)
    
    return swig_class_name

# Dictionary to map from C++ name to mangled WrapITK name, 
# e.g. "unsigned short" -> "US"
mangled_type_map = {}
    
# Build the mangled type map from the names in the itkTypes module. Objects
# which are type should have a name and a short_name. 
for name in dir(itkTypes) :
    value = getattr(itkTypes, name)
    if hasattr(value, "name") and hasattr(value, "short_name") :
        mangled_type_map[value.name] = value.short_name

# From WrapITKTypes.cmake
mangled_type_map["itk::Offset"] = "O"
mangled_type_map["itk::Vector"] = "V"
mangled_type_map["itk::CovariantVector"] = "CV"
mangled_type_map["itk::ContinuousIndex"] = "CI"
mangled_type_map["itk::Array"] = "A"
mangled_type_map["itk::FixedArray"] = "FA"
mangled_type_map["itk::RGBPixel"] = "RGB"
mangled_type_map["itk::RGBAPixel"] = "RGBA"
mangled_type_map["std::complex"] = "C"
mangled_type_map["itk::Image"] = "I"
mangled_type_map["itk::VectorImage"] = "VI"
mangled_type_map["itk::VariableLenghtVector"] = "VLV"
mangled_type_map["itk::Point"] = "P"
mangled_type_map["itk::LevelSetNode"] = "LSN"
mangled_type_map["itk::FlatStructuringElement"] = "SE"
mangled_type_map["itk::SpatialObject"] = "SO"
mangled_type_map["itk::Statistics::Histogram"] = "H"
mangled_type_map["itk::Matrix"] = "M"