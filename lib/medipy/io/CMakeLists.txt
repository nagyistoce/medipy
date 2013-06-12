add_subdirectory("dicom")
add_subdirectory("schemes")

find_package(WrapITK REQUIRED)

file(GLOB files "*.py")
install_py_files(${files})

prepare_wrap_library()
wrap_library("PyArrayIO" python)
    set(WRAPPER_LIBRARY_DEPENDS Base IO)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    set(modules "PyArrayFileReader;PyArrayFileWriter")
    
    find_swig_library_files("PyArrayFileReader" "PyArrayFileWriter")
    
    foreach(module ${modules})
        find_swig_library_files(${module})
        
        wrap_module("itk${module}")
            wrap_class("itk::${module}" POINTER)
                set(dimensions "2;3;4")
                foreach(type ${WRAP_ITK_SCALAR})
                    foreach(dimension ${dimensions})
                        set(parameters "${ITKT_${type}},${dimension}")
                        set(mangled_parameters "${type}${dimension}")
                        wrap_template("${mangled_parameters}" "${parameters}")
                    endforeach()
                endforeach()
            end_wrap_class()
        end_wrap_module()
    endforeach()
end_wrap_library()
wrap_ikt_post_install("PyArrayIO")
