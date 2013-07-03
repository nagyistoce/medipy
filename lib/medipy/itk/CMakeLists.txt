find_package(WrapITK REQUIRED)

file(GLOB files "*.py")
install_py_files(${files})

prepare_wrap_library()
wrap_library("MediPyBridge" python)

    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    find_swig_library_files("NumpyBridge" "MatrixBridge")
    
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
                    wrap_template(${ITKM_${instantiation}} ${ITKT_${instantiation}})
                endforeach()
            endforeach()
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPyBridge")
