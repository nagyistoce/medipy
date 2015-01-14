file(GLOB_RECURSE files "*.py")
install_py_files(${files})

itk_wrap_module("MediPyLogic")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon)
    
    itk_wrap_submodule("itkMaskWithValueImageFilter")
        itk_wrap_class("itk::MaskWithValueImageFilter" "POINTER")
            itk_wrap_image_filter("${WRAP_ITK_SCALAR}" 3)
        itk_end_wrap_class()
    itk_end_wrap_submodule()

itk_end_wrap_module()
