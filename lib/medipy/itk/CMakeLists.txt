# WrapITK is already included: do NOT re-include it, otherwise the patched
# functions will disappear

file(GLOB files "*.py")
install_py_files(${files})

itk_wrap_module("MediPyBridge")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon)
    
    find_swig_library_files("NumpyBridge" "MatrixBridge")
    
    itk_wrap_submodule("itkNumpyBridge")
        itk_wrap_class("itk::NumpyBridge")
            foreach(image "I" "VI")
                foreach(type ${WRAP_ITK_SCALAR})
                    foreach(dimension ${ITK_WRAP_DIMS})
                        set(instantiation ${image}${type}${dimension})
                        itk_wrap_template(${ITKM_${instantiation}} ${ITKT_${instantiation}})
                    endforeach()
                endforeach()
            endforeach()
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    itk_wrap_submodule("itkMatrixBridge")
        itk_wrap_class("itk::MatrixBridge")
            foreach(type "F" "D")
                foreach(dimension ${ITK_WRAP_DIMS})
                    set(instantiation "M${type}${dimension}${dimension}")
                    itk_wrap_template(${ITKM_${instantiation}} ${ITKT_${instantiation}})
                endforeach()
            endforeach()
        itk_end_wrap_class()
    itk_end_wrap_submodule()
itk_end_wrap_module()
