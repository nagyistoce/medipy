import os
import re

from SCons.Builder import Builder

import config
from utils import configure_file

class Python(object) :
    def __init__(self, env) :
        self.env = env
        
        self.library_depends = None
        
        self.WRAP_ITK_PYTHON_CONFIGURATION_TEMPLATES = None
        self.WRAP_ITK_PYTHON_LIBRARY_CALLS = None
        self.WRAP_ITK_PYTHON_LIBRARY_DECLS = None
        self.WRAP_ITK_PYTHON_LIBRARY_IMPORTS = None
        self.WRAP_ITK_PYTHON_SWIG_EXT = None
        
        self._cxx_files = None
        self._current_class = None
        self._current_swig_name = None
        self._deps = None
        
    def wrap_library(self, library) :
        self.WRAP_ITK_PYTHON_CONFIGURATION_TEMPLATES = ""
        self.WRAP_ITK_PYTHON_LIBRARY_IMPORTS = ""
        self._deps = []
        self.WRAP_ITK_PYTHON_LIBRARY_DECLS = ""
        self.WRAP_ITK_PYTHON_LIBRARY_CALLS = ""
        self._cxx_files = ""
    
    def end_wrap_library(self) :
        pass
    
    def wrap_module(self, module) :
        self.WRAP_ITK_PYTHON_SWIG_EXT = ""
        self.WRAP_ITK_PYTHON_LIBRARY_IMPORTS += "from {0}Python import *\n".format(module)
    
    def end_wrap_module(self, group_name) :
        base_name = group_name
        
        if group_name == "itkVectorContainer" :
            text = ""
            for d in config.WRAP_ITK_DIMS :
                text += "%template(vectoritkVectorD{0}) std::vector< itkVectorD{0} >;\n".format(d)
                self.add_config_template("vector", "std::vector", "vectoritkVectorD{0}".format(d), config.type_wrapper.ITKT["VD{0}".format(d)])
                text += "%template(vectoritkPointD{0}) std::vector< itkPointD{0} >;\n".format(d)
                self.add_config_template("vector", "std::vector", "vectoritkPointD{0}".format(d), config.type_wrapper.ITKT["PD{0}".format(d)])
                text += "%template(vectoritkOffset{0}) std::vector< itkOffset{0} >;\n".format(d)
                self.add_config_template("vector", "std::vector", "vectoritkOffset{0}".format(d), config.type_wrapper.ITKT["O{0}".format(d)])
                text += "%template(vectoritkContinuousIndexD{0}) std::vector< itkContinuousIndexD{0} >;\n".format(d)
                self.add_config_template("vector", "std::vector", "vectoritkContinuousIndexD{0}".format(d), config.type_wrapper.ITKT["CID{0}".format(d)])
            text += "%template(vectoritkVectorUC1) std::vector< itkVectorUC1 >;\n"
            self.add_config_template("vector", "std::vector", "vectoritkVectorUC1", config.type_wrapper.ITKT["VUC1"])
            self.WRAP_ITK_PYTHON_SWIG_EXT += text
        
        if group_name == "itkMapContainer" :
            text = ""
            for d in config.WRAP_ITK_DIMS :
                text += "%template(mapULitkVectorD{0}) std::map< unsigned long, itkVectorD{0}, std::less< unsigned long > >;\n".format(d)
                self.add_config_template("map", "std::map", "mapULitkVectorD{0}".format(d), "unsigned long, {0}".format(config.type_wrapper.ITKT["VD{0}".format(d)]))

                text += "%template(mapULitkPointD{0}) std::map< unsigned long, itkPointD{0}, std::less< unsigned long > >;\n".format(d)
                self.add_config_template("map", "std::map", "mapULitkPointD{0}".format(d), "unsigned long, {0}".format(config.type_wrapper.ITKT["PD{0}".format(d)]))
            self.WRAP_ITK_PYTHON_SWIG_EXT += text
        
        if group_name == "ITKQuadEdgeMeshBase" :
            text = ""
            for d in config.WRAP_ITK_DIMS :
                text += "%template(mapULitkQuadEdgeMeshPointF{0}) std::map< unsigned long, itkQuadEdgeMeshPointF{0}, std::less< unsigned long > >;\n".format(d)
                self.add_config_template("map", "std::map", "mapULitkQuadEdgeMeshPointF{0}".format(d), "unsigned long, itk::QuadEdgeMeshPoint< float, {0} >".format(d))

                text += "%traits_swigtype(itkCellInterfaceDQEMCTI{0});\n".format(d)
                text += "%fragment(SWIG_Traits_frag(itkCellInterfaceDQEMCTI{0}));\n".format(d)
                text += "%template(mapULitkCellInterfaceDQEMCTI{0}) std::map< unsigned long, itkCellInterfaceDQEMCTI{0} *, std::less< unsigned long > >;\n".format(d)
                self.add_config_template("map", "std::map", "mapULitkCellInterfaceDQEMCTI{0}".format(d), "unsigned long, itk::CellInterface< double, itk::QuadEdgeMeshCellTraitsInfo< {0} > >*".format(d))
            self.WRAP_ITK_PYTHON_SWIG_EXT += text
        
        if group_name == "ITKFastMarchingBase" :
            text = ""
            for d in config.WRAP_ITK_DIMS :
                for t in config.WRAP_ITK_SCALARS :
                    text += "%template(vectoritkLevelSetNode{0}{1}) std::vector< itkLevelSetNode{0}{1} >;\n".format(t, d)
                    self.add_config_template("vector", "std::vector", "vectoritkLevelSetNode{0}{1}".format(t, d), config.type_wrapper.ITKT["LSN{0}{1}".format(t,d)])
            self.WRAP_ITK_PYTHON_SWIG_EXT += text
        
        # TODO
#        if(WRAP_ITK_DOC AND NOT WRAP_ITK_PYTHON_PROCCESS_SWIG_INPUTS)
#            # yes. Include the docstring file
#            set(doc_file "${WRAPPER_MASTER_INDEX_OUTPUT_DIR}/${WRAPPER_MODULE_NAME}_doc.i")
#            set(WRAP_ITK_PYTHON_SWIG_EXT "%include ${WRAPPER_MODULE_NAME}_doc.i\n\n${WRAP_ITK_PYTHON_SWIG_EXT}")
#        else(WRAP_ITK_DOC AND NOT WRAP_ITK_PYTHON_PROCCESS_SWIG_INPUTS)
#            # no. Clear the doc_file var
#            set(doc_file "")
#        endif(WRAP_ITK_DOC AND NOT WRAP_ITK_PYTHON_PROCCESS_SWIG_INPUTS)
        doc_file = ""
        self.WRAP_ITK_PYTHON_SWIG_EXT = "%import pyBase.i\n\n" + self.WRAP_ITK_PYTHON_SWIG_EXT
        
        lib = "{0}Python".format(group_name)
        python_file = "{0}.py".format(lib)
        cpp_file = "{0}Python.cpp".format(base_name)
        
        def module_ext_i_action(target, source, env):
            configure_file(os.path.join(self.env["WRAP_ITK_PYTHON_SOURCE_DIR"], "module_ext.i.in"),
                           target[0].path,
                           WRAP_ITK_PYTHON_SWIG_EXT=self.WRAP_ITK_PYTHON_SWIG_EXT)
        Builder(action=module_ext_i_action)(self.env, "{0}_ext.i".format(group_name), "")
        
        self._deps = []
        for dep in self.library_depends :
            self._deps.append(dep)
    
    def wrap_named_class(self, class_, swig_name) :
        self._current_class = class_
        self._current_swig_name = swig_name
    
    def wrap_template(self, name, types) :
        if self._current_class == "itk::ImageSource" :
            image_source = "{0}{1}".format(self._current_swig_name, name)
            image = config.type_wrapper.ITKN[name]
            
            text = "\n\n"
            text += "%typemap(in) {0} * {\n".format(image)
            text += "  {0} * imgsrc;\n".format(image_source)
            text += "  {0} * img;\n".format(image)
            text += "  if( $input != Py_None && SWIG_ConvertPtr($input,(void **)(&imgsrc),\$descriptor({0} *), 0) == 0 )\n".format(image_source)
            text += "    {\n"
            text += "    \$1 = imgsrc->GetOutput(0);\n"
            text += "    }\n"
            text += "  else if( SWIG_ConvertPtr($input,(void **)(&img),\$descriptor({0} *), 0) == 0 )\n".format(image)
            text += "    {\n"
            text += "    \$1 = img;\n"
            text += "    }\n"
            text += "  else\n"
            text += "    {\n"
            text += "    PyErr_SetString(PyExc_TypeError, \"Expecting argument of type {0} or {1}.\");\n".format(image, image_source)
            text += "    SWIG_fail;\n"
            text += "    }\n"
            text += "}\n"
            text += "\n"
            text += "\n"
            text += "%typemap(typecheck) {0} * {\n".format(image)
            text += "  {0} * imgsrc;\n".format(image_source)
            text += "  {0} * img;\n".format(image)
            text += "  if( $input != Py_None && SWIG_ConvertPtr($input,(void **)(&imgsrc),\$descriptor({0} *), 0) == 0 )\n".format(image_source)
            text += "    {\n"
            text += "    \$1 = 1;\n"
            text += "    }\n"
            text += "  else if( SWIG_ConvertPtr($input,(void **)(&img),\$descriptor({0} *), 0) == 0 )\n".format(image)
            text += "    {\n"
            text += "    \$1 = 1;\n"
            text += "    }\n"
            text += "  else\n"
            text += "    {\n"
            text += "    PyErr_Clear();\n"
            text += "    \$1 = 0;\n"
            text += "    }\n"
            text += "}\n"
            self.WRAP_ITK_PYTHON_SWIG_EXT += text
    
    def add_one_typedef(self, wrap_class, swig_name, wrap_method="", template_params=None) :
        base_name = wrap_class.split("::")[-1]
        if wrap_class != "MetaEvent" and "ENUM" not in wrap_method :
            self.add_config_template(base_name, wrap_class, swig_name, template_params)
        
        if wrap_class == "vcl_complex" :
            self.add_config_template("complex", "std::complex", swig_name, template_params)
    
    def add_simple_typedef(self, wrap_class, swig_name) :
        if re.match(r".*<.*>", wrap_class) :
            cpp_name = re.sub(r"^([^<]+)< *(.+) *>([^>]*)$", r"\1", wrap_class)
            template_params = re.sub(r"^([^<]+)< *(.+) *>([^>]*)$", r"\2", wrap_class)
            ext_def = re.sub(r"^([^<]+)< *(.+) *>([^>]*)$", r"\3",  wrap_class)
        else :
            cpp_name = wrap_class
            template_params = None
            ext_def = ""
        
        simple_name = re.sub(".*::", "", cpp_name)
        
        if swig_name.endswith("_Pointer") :
            smart_pointed = re.sub("_Pointer$", "", swig_name)
            self.add_pointer_typemap(smart_pointed)
        
        if cpp_name == "vcl_complex" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_VCL_COMPLEX_CLASS({0})\n".format(swig_name)
        if swig_name == "itkLightObject" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(listitkLightObject) std::list< itkLightObject_Pointer >;\n\n"
            self.add_config_template("list", "std::list", "listitkLightObject", "itk::LightObject")
        if swig_name == "itkObject" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_OBJECT_CLASS({0})\n".format(swig_name)
        if swig_name == "itkProcessObject" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_PROCESSOBJECT_CLASS({0})\n".format(swig_name)
        if swig_name == "itkDataObject" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(vectoritkDataObject) std::vector< itkDataObject_Pointer >;\n"
            self.add_config_template("vector", "std::vector", "vectoritkDataObject", "itk::DataObject")
        if swig_name == "itkObjectFactoryBase" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(listitkObjectFactoryBase) std::list< itkObjectFactoryBase * >;\n"
            self.add_config_template("list", "std::list", "listitkObjectFactoryBase", "itk::ObjectFactoryBase")
        if swig_name == "itkMetaDataDictionary" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(vectoritkMetaDataDictionary) std::vector< itkMetaDataDictionary * >;\n"
            self.add_config_template("vector", "std::vector", "vectoritkMetaDataDictionary", "itk::MetaDataDictionary")
        if swig_name == "itkCommand" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%feature(\"director\") itkCommand;\n"
        if cpp_name == "itk::ImageBase" and not swig_name.endswith("Pointer") :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_IMAGEBASE_CLASS({0}, {1})\n".format(swig_name, ", ".join(template_params))
            self.WRAP_ITK_PYTHON_SWIG_EXT += "${WRAP_ITK_PYTHON_SWIG_EXT}%inline %{\n"
            self.WRAP_ITK_PYTHON_SWIG_EXT += "#include \"itkContinuousIndexSwigInterface.h\"\n"
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%}\n"
        if cpp_name == "itk::StatisticsLabelObject" and not swig_name.endswith("Pointer") :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(map{0}) std::map< unsigned long, {1}_Pointer, std::less< unsigned long > >;\n".format(swig_name, swig_name)
            self.add_config_template("map", "std::map", "map{0}".format(swig_name), "unsigned long, {0}< {1} >".format(cpp_name, ", ".join(template_params)))
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(vector{swig_name}) std::vector< {swig_name}_Pointer >;\n".format(swig_name, swig_name)
            self.add_config_template("vector", "std::vector", "vector{0}".format(swig_name), "{0}< {1} >".format(cpp_name, ", ".join(template_params)))
        if cpp_name == "itk::LabelMap" and not swig_name.endswith("Pointer") :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_LABELMAP_CLASS(${swig_name})\n"
        if cpp_name == "itk::ComponentTreeNode" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(list{0}) std::list< {0}* >;\n".format(swig_name)
            self.add_config_template("list", "std::list", "list{0}".format(swig_name), "{0}< {1} > *".format(cpp_name, ", ".join(template_params)))
        if cpp_name == "itk::ImageRegion" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_IMAGEREGION_CLASS(${swig_name})\n"
        if cpp_name == "itk::Index" :
            self.add_seq_typemap(swig_name, template_params)
        if cpp_name == "itk::Size" :
            self.add_seq_typemap(swig_name, template_params)
        if cpp_name == "itk::RGBPixel" :
            self.add_variable_length_seq_typemap(swig_name, template_params)
        if cpp_name == "itk::RGBAPixel" :
            self.add_variable_length_seq_typemap(swig_name, template_params)
        if cpp_name == "itk::Offset" :
            self.add_seq_typemap(swig_name, template_params)
        if cpp_name == "itk::FixedArray" :
            self.add_vec_typemap(swig_name, template_params)
        if cpp_name == "itk::Vector" :
            self.add_vec_typemap(swig_name, template_params)
        if cpp_name == "itk::CovariantVector" :
            self.add_vec_typemap(swig_name, template_params)
        if cpp_name == "itk::Point" :
            self.add_vec_typemap(swig_name, template_params)
        if cpp_name == "itk::ContinuousIndex" :
            self.add_vec_typemap(swig_name, template_params)
        if cpp_name == "itk::Array" :
            self.add_variable_length_seq_typemap(swig_name, template_params)
        if swig_name == "itkTransformBase" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(list{0}_Pointer) std::list< {1}_Pointer >;\n".format(swig_name, swig_name) 
            self.add_config_template("list", "std::list", "list{0}_Pointer".format(swig_name), "itk::TransformBase")
        if cpp_name == "itk::SpatialObjectPoint" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_SPATIALOBJECTPPOINT_CLASS({0})%template(vector{0}) std::vector< {0} >;\n".format(swig_name)
            self.add_config_template("vector", "std::vector", "vector{0}".format(swig_name), "{0}< {1} >".format(cpp_name, template_params))
        if cpp_name == "itk::ContourSpatialObjectPoint" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(vector{0}) std::vector< {0} >;\n".format(swig_name)
            self.add_config_template("vector", "std::vector", "vector{0}".format(swig_name), "{0}< {1} >".format(cpp_name, template_params))
        if cpp_name == "itk::LineSpatialObjectPoint" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(vector{0}) std::vector< {0} >;\n".format(swig_name)
            self.add_config_template("vector", "std::vector", "vector{0}".format(swig_name), "{0}< {1} >".format(cpp_name, template_params))
        if cpp_name == "itk::SurfaceSpatialObjectPoint" :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(vector{0}) std::vector< {0} >;\n".format(swig_name)
            self.add_config_template("vector", "std::vector", "vector{0}".format(swig_name), "{0}< {1} >".format(cpp_name, template_params))
        if cpp_name == "itk::SpatialObject" and Pointer not in ext_def :
            self.WRAP_ITK_PYTHON_SWIG_EXT += "%template(list{0}_Pointer) std::list< {0}_Pointer >;\n".format(swig_name)
            self.add_config_template("list", "std::list", "list{0}_Pointer".format(swig_name), "{0}< {1} >".format(cpp_name, template_params))
        
    def add_config_template(self, base_name, wrap_class, swig_name, template_params):
        if not template_params :
            self.WRAP_ITK_PYTHON_CONFIGURATION_TEMPLATES += "  ('{0}', '{1}', '{2}'),\n".format(
                base_name, wrap_class, swig_name)
        else :
            self.WRAP_ITK_PYTHON_CONFIGURATION_TEMPLATES += "  ('{0}', '{1}', '{2}', '{3}'),\n".format(
                base_name, wrap_class, swig_name, ", ".join(template_params))
    
    def add_pointer_typemap(self, template_params) :
        text = "DECLARE_REF_COUNT_CLASS({0})\n".format(template_params)
        self.WRAP_ITK_PYTHON_SWIG_EXT = text+self.WRAP_ITK_PYTHON_SWIG_EXT
    
    def add_seq_typemap(self, swig_name, dim) :
        self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_SEQ_TYPEMAP({0}, {1})\n".format(swig_name, dim)
        
    def add_seq_typemap(self, swig_name, value_type) :
        self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_VARLEN_SEQ_TYPEMAP({0}, {1})\n".format(type, value_type)
    
    def add_vec_typemap(self, swig_name, template_params) :
        type, dim = template_params
        self.WRAP_ITK_PYTHON_SWIG_EXT += "DECL_PYTHON_VEC_TYPEMAP({0}, {1}, {2})\n".format(swig_name, type, dim)
