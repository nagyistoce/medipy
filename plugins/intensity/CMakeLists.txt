file(GLOB files "*.py")
install_py_files(${files})

prepare_wrap_library()
wrap_library("MediPyIntensity" python)
    
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon ITKStatistics)
    
    wrap_module("itkJointHistogramNormalizationFilter")
        wrap_class("itk::JointHistogramNormalizationFilter" "POINTER")
            wrap_image_filter("${WRAP_ITK_SCALAR}" 3)
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPyIntensity")
