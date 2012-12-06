add_subdirectory("dicom")
add_subdirectory("schemes")

file(GLOB files "*.py")
install_py_files(${files})

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

wrap_library("PyArrayIO" python)
    set(WRAPPER_LIBRARY_DEPENDS Base IO)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    set(modules "PyArrayFileReader;PyArrayFileWriter")
    foreach(module ${modules})
        set(WRAPPER_SWIG_LIBRARY_FILES 
            ${WRAPPER_SWIG_LIBRARY_FILES} "${CMAKE_CURRENT_SOURCE_DIR}/itk${module}.swg")
        
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
