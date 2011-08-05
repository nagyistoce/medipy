import logging
import re

import config
from gccxml import GccXML

class Wrapper(object):
    
    def __init__(self, env):
        self.env = env
        
        self.default_include = ["itkCommand.h",
                                # itkStatisticsLabelObject.h does not exist in 
                                # ITK 3.14
                                #"itkStatisticsLabelObject.h"
                                ]+config.type_wrapper.includes
        
        self.languages = ["GccXML", "SwigInterface", "Python"]
        allowed_languages = ["GccXML", "SwigInterface", "Doc", 
                             "Python", 
                             #"TCL", "Ruby", "Java",
                             #"Explicit"
                            ]
        
        self.library_name = None
        self.module_name = None
        self.auto_include_headers = None
        self.include_files = None
        self.wrap_method = None
        self.class_ = None
        self.swig_name = None
        
        self.language_wrappers = {}
        for language in allowed_languages :
            if language in globals() :
                self.language_wrappers[language] = globals()[language](env)
            else :
                logging.warning("No builder for language \"{0}\"".format(language))
    
    def wrap_libraries(self) :
        pass
    
    def end_wrap_libraries(self):
        pass
    
#    AUTO_INCLUDE_MODULES
#    INCLUDE_MODULE
    
    def wrap_library(self, library, languages=None):
        self.library_name = library
        
        if languages is not None :
            self.languages = languages
        
        allowed_languages = ["GccXML", "SwigInterface", "Doc", 
                             "Python", "TCL", "Ruby", "Java",
                             "Explicit"]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_library(library)
    
    def end_wrap_library(self):
        allowed_languages = ["GccXML", "SwigInterface", "Doc", 
                             "Python", "TCL", "Ruby", "Java",
                             "Explicit"]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].end_wrap_library()
    
    def wrap_module(self, module):
        upper_module = module.upper()
        upper_lib = self.library_name.upper()
        if upper_module == upper_lib :
            raise Exception("The module {0} can't have the same name than its "
                            "library. Note that the names are not case "
                            "sensitive.".format(module))
        
        self.module_name = module
        
        allowed_languages = ["GccXML", "SwigInterface", "Doc", 
                             "Python", "TCL", "Ruby", "Java",
                             "Explicit"]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_module(module)
        
        self.include_files = []
        for inc in self.default_include :
            self.wrap_include(inc)
        self.auto_include_headers = True
    
    def end_wrap_module(self, module):
        allowed_languages = ["GccXML", "SwigInterface", "Doc", 
                             "Python", "TCL", "Ruby", "Java",
                             "Explicit"]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].end_wrap_module(self.module_name)
    
    def wrap_class(self, class_, wrap_method=""):
        if "::" in class_ :
            elements = class_.split("::")
            top_namespace, base_name = elements[0], elements[-1]
            swig_name = top_namespace+base_name
        else :
            swig_name = class_
        
        self.wrap_named_class(class_, swig_name, wrap_method)
        
        if self.auto_include_headers :
            self.wrap_include("{0}.h".format(swig_name))

        allowed_languages = ["Doc", 
                             "Python", "TCL", "Ruby",
                            ]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].end_wrap_module(self.module_name)
    
    def end_wrap_class(self):
        pass
    
    def wrap_named_class(self, class_, swig_name, wrap_method=""):
        wrap_methods =["", "POINTER", "POINTER_WITH_SUPERCLASS", 
                       "POINTER_WITH_2_SUPERCLASSES", "EXPLICIT_SPECIALIZATION",
                       "POINTER_WITH_EXPLICIT_SPECIALIZATION", "ENUM", 
                       "AUTOPOINTER"]
        if wrap_method not in wrap_methods :
            raise Exception("WRAP_CLASS: Invalid option '{0}'. Possible values "
                "are POINTER, POINTER_WITH_SUPERCLASS, "
                "POINTER_WITH_2_SUPERCLASSES, EXPLICIT_SPECIALIZATION, "
                "POINTER_WITH_EXPLICIT_SPECIALIZATION, ENUM and AUTOPOINTER".format(wrap_method))
        self.wrap_method = wrap_method
        
        self.class_ = class_
        self.swig_name = swig_name
        # TODO
#        set(WRAPPER_WARN_ABOUT_NO_TEMPLATE ON)

        allowed_languages = ["Doc", 
                             "Python", "TCL", "Ruby",
                            ]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_named_class(class_, swig_name)
        
    
    def wrap_non_template_class(self, class_, wrap_method=""):
        self.wrap_class(class_, wrap_method)
        # TODO
#        # to avoid useless warning: no template can be defined in
#        set(WRAPPER_WARN_ABOUT_NO_TEMPLATE OFF)
        self.add_one_typedef(self.class_, self.swig_name, wrap_method)
        self.end_wrap_class()
        
        allowed_languages = []
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_named_class(class_, swig_name)
    
    def wrap_named_non_template_class(self, class_, swig_name, wrap_method=""):
        self.wrap_named_class(class_, swig_name, wrap_method)
        # TODO
#        # to avoid useless warning: no template can be defined in
#        set(WRAPPER_WARN_ABOUT_NO_TEMPLATE OFF)
        self.add_one_typedef(self.class_, self.swig_name, wrap_method)
        self.end_wrap_class()
        
        allowed_languages = []
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_named_non_template_class(class_, swig_name)
    
    def wrap_include(self, include_file):
        if include_file not in self.include_files :
            self.include_files.append(include_file)
        
        allowed_languages = ["GccXML", "SwigInterface",
                             "Explicit"]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_include(include_file)
    
    def add_simple_typedef(self, wrap_class, swig_name):
        allowed_languages = ["GccXML",
                             "Python", "TCL", "Ruby", "Java"
                            ]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].add_simple_typedef(wrap_class, swig_name)
    
    def add_one_typedef(self, wrap_class, swig_name, wrap_method="", template_parameters=None):
        template_parameters = template_parameters or []
        
        elements = wrap_class.split("::")
        base_name = elements[-1]
        
        wrap_pointer = 0
        if template_parameters :
            full_class_name = "{0}< {1} >".format(wrap_class, ", ".join(template_parameters))
        else :
            full_class_name = wrap_class
        
        allowed_languages = ["GccXML", "SwigInterface", "Doc"
                             "Python", "TCL", "Ruby", "Java",
                             "Explicit"
                            ]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].add_one_typedef(wrap_class, swig_name, wrap_method, template_parameters)
    
        if "2_SUPERCLASSES" in wrap_method :
            self.add_simple_typedef("{0}::Superclass::Superclass".format(full_class_name),
                                    "{0}_Superclass_Superclass".format(swig_name))
            self.add_simple_typedef("{0}::Superclass::Superclass::Pointer".format(full_class_name),
                                    "{0}_Superclass_Superclass_Pointer".format(swig_name))
        
        if "SUPERCLASS" in wrap_method :
            self.add_simple_typedef("{0}::Superclass".format(full_class_name),
                                    "{0}_Superclass".format(swig_name))
            self.add_simple_typedef("{0}::Superclass::Pointer".format(full_class_name),
                                    "{0}_Superclass_Pointer".format(swig_name))
        
        self.add_simple_typedef(full_class_name, swig_name)
        
        if "POINTER" in wrap_method :
            if wrap_method == "AUTOPOINTER" :
                self.add_simple_typedef("{0}::SelfAutoPointer".format(full_class_name),
                                        "{0}_AutoPointer".format(swig_name))
            else :
                self.add_simple_typedef("{0}::Pointer".format(full_class_name),
                                        "{0}_Pointer".format(swig_name))
    
    def wrap_template(self, name, types):
        # TODO
#        set(WRAPPER_WARN_ABOUT_NO_TEMPLATE OFF)
        self.add_one_typedef(self.class_, self.swig_name+name, self.wrap_method, types)
        
        allowed_languages = ["Python", "TCL", "Ruby"
                            ]
        for language in self.languages :
            if language in allowed_languages and language in self.language_wrappers :
                self.language_wrappers[language].wrap_template(name, types)
    
    def wrap_image_filter(self, param_types, param_count, dim_cond=None):
        for param_type in param_types :
            param_list = []
            for i in range(param_count) :
                param_list.append(param_type)
            
            self.wrap_image_filter_types(param_list, dim_cond)
    
    def wrap_image_fitler_combinations(self, *args):
        pass
    
    def wrap_image_filter_types(self, param_list, dim_cond=None):
        if dim_cond :
            dims = self.filter_dims(dim_cond)
        else :
            dims = env["WRAP_ITK_DIMS"]
        
        for d in dims :
            template_params = []
            mangled_name = []
            
            for type in param_list :
                if type in config.WRAP_ITK_VECTOR :
                    type = "{0}{1}".format(type, d)
                image_type = config.type_wrapper.ITKT["I{0}{1}".format(type,d)]
                mangle_type = config.type_wrapper.ITKM["I{0}{1}".format(type,d)]
                template_params.append(image_type)
                mangled_name.append(mangle_type)
            self.wrap_template("".join(mangled_name), template_params)
    
    def filter_dims(self, dimension_condition):
        dims = []
        
        if isinstance(dimension_condition, basestring) and re.match(r"^[0-9]+\+$", dimension_condition) :
            min_dim = int(re.sub(r"^([0-9]+)\+$", r"\1", dimension_condition))
            max_disallowed = min_dim-1
            var_name = ""
            for d in env["WRAP_ITK_DIMS"] :
                if d >= min_dim :
                    dims.append(d)
        else :
            # The condition is just a list of dims. Return the intersection of these
            # dims with the selected ones.
            s1 = set(env["WRAP_ITK_DIMS"])
            s2 = set(dimension_condition)
            dims = list(s1.intersection(s2))
        
        return dims
    
if __name__ == "SCons.Script" :
    env = Environment()
    env.AppendUnique(CPPPATH=[
        "/usr/include/InsightToolkit/",
        "/usr/include/InsightToolkit/BasicFilters/", 
        "/usr/include/InsightToolkit/Common/",
        "/usr/include/InsightToolkit/Numerics/Statistics/",
        "/usr/include/InsightToolkit/Review/",
        "/usr/include/InsightToolkit/SpatialObject/",
        "/usr/include/InsightToolkit/Utilities/vxl/core/",
        "/usr/include/InsightToolkit/Utilities/vxl/vcl/",
    ])
    
    env["WRAPPER_LIBRARY_OUTPUT_DIR"] = "."
    env["WRAP_ITK_DIMS"] = [2,3]
    env["WRAP_ITK_USE_CCACHE"] = False
    env["WRAP_ITK_GCCXML_SOURCE_DIR"] = "/usr/lib/InsightToolkit/WrapITK/Configuration/Languages/GccXML"
    
    wrapper = Wrapper(env)
    wrapper.wrap_library("PixelMath")
    wrapper.wrap_module("itkAndImageFilter")
    wrapper.wrap_class("itk::AndImageFilter", "POINTER_WITH_SUPERCLASS")
    wrapper.wrap_image_filter(config.WRAP_ITK_USIGN_INT, 3)
    wrapper.wrap_image_filter(config.WRAP_ITK_SIGN_INT, 3)
    wrapper.end_wrap_class()
    wrapper.end_wrap_module("itkAndImageFilter")
    wrapper.end_wrap_library()