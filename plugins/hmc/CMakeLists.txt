file(GLOB_RECURSE files "*.py")
install_py_files(${files})

add_library(medipy_hmc SHARED EM.cpp HMCInitializer.cpp)

itk_wrap_module("MediPyHMC")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon)
    
    itk_wrap_submodule("itkHMCSegmentationImageFilter")
        itk_wrap_class("itk::HMCSegmentationImageFilter" POINTER)
            itk_wrap_template("IF3IF3IF3" "${ITKT_IF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        itk_end_wrap_class()
    itk_end_wrap_submodule()

itk_end_wrap_module()

add_dependencies(${lib} medipy_hmc)
    
