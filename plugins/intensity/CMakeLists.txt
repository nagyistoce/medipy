find_package(VTK REQUIRED)
include(${USE_VTK_FILE})

file(GLOB files "*.py")
install_py_files(${files})

prepare_wrap_library()
wrap_library("MediPyIntensity" python)
    
    set(WRAPPER_LIBRARY_DEPENDS Base Numerics)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon ITKNumerics ITKStatistics)
    
    find_swig_library_files("JointHistogramTransferFunctionCalculator")
    
    wrap_module("itkJointHistogramCalculator")
        wrap_class("itk::JointHistogramCalculator" "POINTER")
            wrap_image_filter("${WRAP_ITK_SCALAR}" 2)
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkJointHistogramTransferFunctionCalculator")
        wrap_class("itk::JointHistogramTransferFunctionCalculator" "POINTER")
            wrap_template("HD2" ${ITKT_HF2})
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkTransformIntensityImageFilter")
        wrap_class("itk::TransformIntensityImageFilter" "POINTER")
            wrap_image_filter("${WRAP_ITK_SCALAR}" 2)
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPyIntensity")
