file(GLOB files "*.py")
install_py_files(${files})

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)
#set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")

wrap_library("MediPyITK" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    set(WRAPPER_SWIG_LIBRARY_FILES ${WRAPPER_SWIG_LIBRARY_FILES} 
        "${CMAKE_CURRENT_SOURCE_DIR}/itkNumpyBridge.swg"
        "${CMAKE_CURRENT_SOURCE_DIR}/itkMatrixBridge.swg")
    auto_include_modules()
end_wrap_library()

wrap_ikt_post_install("MediPyITK")
