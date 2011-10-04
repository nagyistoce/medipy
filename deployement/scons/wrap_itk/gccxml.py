import os
import re

from SCons.Builder import Builder

from utils import configure_file

class GccXML(object) :
    
    def __init__(self, env) :
        self.env = env
        
        self.library_depends = None
        
        self.includes = ""
        self.typedefs = ""
        self.force_instantiate = ""
    
    def wrap_library(self, library) :
        gccxml_inc_file = os.path.join(self.env["WRAPPER_LIBRARY_OUTPUT_DIR"], "gcc_xml.inc")
        
        CONFIG_GCCXML_INC_CONTENTS = ""
        
        if self.env["WRAP_ITK_USE_CCACHE"] :
            pass
            # TODO
        else :
            for dir in self.env["CPPPATH"] :
                CONFIG_GCCXML_INC_CONTENTS += "-I{0}\n".format(dir)
                def gcc_xml_inc_action(target, source, env) :
                    configure_file(os.path.join(self.env["WRAP_ITK_GCCXML_SOURCE_DIR"], "gcc_xml.inc.in"),
                                   target[0].path,
                                   CONFIG_GCCXML_INC_CONTENTS=CONFIG_GCCXML_INC_CONTENTS)
                Builder(action=gcc_xml_inc_action)(self.env, gccxml_inc_file, "")
    
    def end_wrap_library(self) :
        pass
    
    def wrap_module(self, module) :
        self.typedefs = ""
        self.includes = ""
        self.force_instantiate = ""
    
    def end_wrap_module(self, module) :
        file_name = "{0}.cxx".format(module)
        cxx_file = os.path.join(self.env["WRAPPER_LIBRARY_OUTPUT_DIR"], file_name)
        
        def cxx_file_action(target, source, env) :
            configure_file(os.path.join(self.env["WRAP_ITK_GCCXML_SOURCE_DIR"], "wrap_.cxx.in"),
                target[0].path,
                CONFIG_WRAPPER_INCLUDES=self.includes,
                CONFIG_WRAPPER_MODULE_NAME=module,
                CONFIG_WRAPPER_TYPEDEFS=self.typedefs,
                CONFIG_WRAPPER_FORCE_INSTANTIATE=self.force_instantiate,
                CONFIG_WRAPPER_EXTRAS="")
        
        Builder(action=cxx_file_action)(self.env, cxx_file, "")
        
        gccxml_inc_file = os.path.join(self.env["WRAPPER_LIBRARY_OUTPUT_DIR"], "gcc_xml.inc")
        xml_file = os.path.join(self.env["WRAPPER_LIBRARY_OUTPUT_DIR"], "{0}.xml".format(module))
        
        if self.env["WRAP_ITK_USE_CCACHE"] :
            pass
            # TODO
        else :
            xml_file_action = (
                "gccxml -fxml-start=_cable_ -fxml=$TARGET "
                "--gccxml-gcc-options {gccxml_inc_file} -DCSWIG "
                "-DCABLE_CONFIGURATION -DITK_MANUAL_INSTANTIATION "
                "$SOURCE".format(gccxml_inc_file=gccxml_inc_file)
            )
        
        self.env.Depends(xml_file, gccxml_inc_file)
        Builder(action=xml_file_action)(self.env, xml_file, cxx_file)
    
    def wrap_include(self, include_file) :
        if re.match(r"<.*>", include_file) :
            self.includes += "#include {0}\n".format(include_file)
        else :
            self.includes += "#include \"{0}\"\n".format(include_file)
    
    def add_one_typedef(self, wrap_class, swig_name, wrap_method="", template_params=None) :
        self.typedefs += "\n"
        self.force_instantiate += "\n"
    
    def add_simple_typedef(self, wrap_class, swig_name) :
        self.typedefs += "      typedef {0} {1};\n".format(wrap_class, swig_name)
        self.force_instantiate += "  sizeof({0});\n".format(swig_name)
