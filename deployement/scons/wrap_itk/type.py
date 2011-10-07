import copy

import config

class Type(object):
    def __init__(self, name, template_parameters=None, nested_type=None):
        self.name = name
        self.template_parameters = template_parameters or []
        self.nested_type = nested_type
    
    def __eq__(self, other):
        return self.full_name == other.full_name
    
    def _get_full_name(self):
        if not self.template_parameters :
            return self.name
        
        parameters = []
        for p in self.template_parameters :
            if isinstance(p, basestring) :
                parameters.append(p)
            elif isinstance(p, int) :
                parameters.append(str(p))
            else :
                parameters.append(p.full_name)
        
        parameters = ", ".join(parameters)
        
        if self.nested_type :
            if isinstance(self.nested_type, basestring) :
                nested_type = "::{0}".format(self.nested_type)
            else :
                nested_type = "::{0}".format(self.nested_type.full_name)
        else :
            nested_type = ""
        
        return "{0}<{1}{2}>{3}".format(
            self.name,
            parameters, " " if parameters.endswith(">") else "", 
            nested_type)
                
    def _get_mangled_name(self):
        
        name = Type._mangled_names.get(self.name, self.name)
        name = name.replace("::", "")
        
        if not self.template_parameters :
            return name
        
        parameters = []
        for p in self.template_parameters :
            if isinstance(p, basestring) :
                parameters.append(Type._mangled_names[p])
            elif isinstance(p, int) :
                parameters.append(str(p))
            else :
                parameters.append(p.mangled_name)
        
        parameters = "".join(parameters)
        
        if self.nested_type :
            if isinstance(self.nested_type, basestring) :
                nested_type = "_{0}".format(self.nested_type)
            else :
                nested_type = "_{0}".format(self.nested_type.mangled_name)
        else :
            nested_type = ""
            
        return "{0}{1}{2}".format(name, parameters, nested_type)
    
    full_name = property(_get_full_name)
    mangled_name = property(_get_mangled_name)
    
    _mangled_names = {
        "unsigned char" : "UC", "unsigned short" : "US", "unsigned int" : "UI", 
        "unsigned long" : "UL", "signed char" : "SC", "signed short" : "SS", 
        "signed int" : "SI", "signed long" : "SL", "float" : "F", 
        "double" : "D", "long double" : "LD", "bool" : "B",
        
        "itk::Offset" : "O", "itk::Vector" : "V", "itk::CovariantVector" : "CV",
        "itk::ContinuousIndex" : "CI", "itk::Array" : "A", 
        "itk::FixedArray" : "FA", "itk::RGBPixel" : "RGB",
        "itk::RGBAPixel" : "RGBA", "std::complex" : "C",
        "itk::SymmetricSecondRankTensor" : "SSRT", "itk::Image" : "I",
        "itk::VectorImage" : "VI", "itk::VariableLengthVector" : "VLV",
        "itk::Point" : "P", "itk::LevelSetNode" : "LSN",
        "itk::FlatStructuringElement" : "SE", "itk::SpatialObject" : "SO",
        "itk::Statistics::Histogram" : "H", "itk::LabelMap" : "LM",
    }
    