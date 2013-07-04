find_package(WrapITK REQUIRED)

file(GLOB_RECURSE files "*.py")
install_py_files(${files})

prepare_wrap_library()
wrap_library("MediPyHMC" python)
    set(WRAPPER_LIBRARY_DEPENDS Base Numerics)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    wrap_module("itkHiddenMarkovChainFilter")
        wrap_class("itk::HiddenMarkovChainFilter" POINTER)
            wrap_template("VIF3IF3" "${ITKT_VIF3}, ${ITKT_IF3}")
        end_wrap_class()
    end_wrap_module()

end_wrap_library()
wrap_ikt_post_install("MediPyHMC")
