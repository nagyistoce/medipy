file(GLOB_RECURSE files "*.py")
install_py_files(${files})

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

wrap_library("MediPyRegionGrowing" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    wrap_module("itkConnectedThresholdWithRadiusImageFilter")
        wrap_class("itk::ConnectedThresholdWithRadiusImageFilter" "POINTER")
            wrap_image_filter("${WRAP_ITK_SCALAR}" 2)
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPyRegionGrowing")
