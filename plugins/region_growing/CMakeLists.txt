file(GLOB_RECURSE files "*.py")
install_py_files(${files})

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

wrap_library("MediPyRegionGrowing" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    set(WRAPPER_SWIG_LIBRARY_FILES ${WRAPPER_SWIG_LIBRARY_FILES}
        "${CMAKE_CURRENT_SOURCE_DIR}/itkConnectedThresholdWithRadiusImageFilter.swg"
    )
    auto_include_modules()
end_wrap_library()
wrap_ikt_post_install("MediPyRegionGrowing")
