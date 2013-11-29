find_package(WrapITK REQUIRED)

file(GLOB_RECURSE files "*.py")
install_py_files(${files})

add_library(medipy_hmc SHARED EM.cpp HMCInitializer.cpp)

prepare_wrap_library()
wrap_library("MediPyHMC" python)
    set(WRAPPER_LIBRARY_DEPENDS Base Numerics)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon medipy_hmc)
    
    wrap_module("itkHMCSegmentationImageFilter")
        wrap_class("itk::HMCSegmentationImageFilter" POINTER)
            wrap_template("IF3IF3IF3" "${ITKT_IF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        end_wrap_class()
    end_wrap_module()

end_wrap_library()
wrap_ikt_post_install("MediPyHMC")
