file(GLOB files "*.py")
install_py_files(${files})

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

wrap_library("MediPyBridge" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    set(WRAPPER_SWIG_LIBRARY_FILES ${WRAPPER_SWIG_LIBRARY_FILES} 
        "${CMAKE_CURRENT_SOURCE_DIR}/itkNumpyBridge.swg"
        "${CMAKE_CURRENT_SOURCE_DIR}/itkMatrixBridge.swg")
    
    wrap_module("itkNumpyBridge")
        wrap_class("itk::NumpyBridge")
            foreach(image "I" "VI")
                foreach(type ${WRAP_ITK_SCALAR})
                    foreach(dimension ${WRAP_ITK_DIMS})
                        set(instantiation ${image}${type}${dimension})
                        wrap_template(${ITKM_${instantiation}} ${ITKT_${instantiation}})
                    endforeach()
                endforeach()
            endforeach()
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkMatrixBridge")
        wrap_class("itk::MatrixBridge")
            foreach(type "F" "D")
                foreach(dimension ${WRAP_ITK_DIMS})
                    set(instantiation "M${type}${dimension}${dimension}")
                    wrap_template("${ITKM_${instantiation}}" "${ITKT_${instantiation}}")
                endforeach()
            endforeach()
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()

wrap_ikt_post_install("MediPyBridge")
