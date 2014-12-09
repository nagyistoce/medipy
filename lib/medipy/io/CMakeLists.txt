add_subdirectory("dicom")
add_subdirectory("schemes")

file(GLOB files "*.py")
install_py_files(${files})

itk_wrap_module("PyArrayIO")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon ITKIOImageBase)
    
    set(modules "PyArrayFileReader;PyArrayFileWriter")
    
    find_swig_library_files("PyArrayFileReader" "PyArrayFileWriter")
    
    foreach(module ${modules})
        itk_wrap_submodule("itk${module}")
            itk_wrap_class("itk::${module}" POINTER)
                set(dimensions "2;3;4")
                foreach(type ${WRAP_ITK_SCALAR})
                    foreach(dimension ${dimensions})
                        set(parameters "${ITKT_${type}},${dimension}")
                        set(mangled_parameters "${type}${dimension}")
                        itk_wrap_template("${mangled_parameters}" "${parameters}")
                    endforeach()
                endforeach()
            itk_end_wrap_class()
        itk_end_wrap_submodule()
    endforeach()
itk_end_wrap_module()
