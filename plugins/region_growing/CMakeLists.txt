file(GLOB_RECURSE files "*.py")
install_py_files(${files})

itk_wrap_module("MediPyRegionGrowing")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon)
    
    itk_wrap_submodule("itkConnectedThresholdWithRadiusImageFilter")
        itk_wrap_class("itk::ConnectedThresholdWithRadiusImageFilter" "POINTER")
            itk_wrap_image_filter("${WRAP_ITK_SCALAR}" 2)
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
itk_end_wrap_module()
