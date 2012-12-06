file(GLOB files "*.py")
install_py_files(${files})

# Build "parse" Python module
find_package(ITK)
if(${ITK_FOUND})
    if(${ITK_USE_SYSTEM_GDCM})
        find_package(GDCM REQUIRED)
    else(${ITK_USE_SYSTEM_GDCM})
        include_directories(${ITK_INCLUDE_DIRS})
        link_directories(${ITK_LIBRARY_DIRS})
    endif(${ITK_USE_SYSTEM_GDCM})
else(${ITK_FOUND})
    find_package(GDCM REQUIRED)
endif(${ITK_FOUND})

set_source_files_properties(parse.i PROPERTIES CPLUSPLUS ON)
swig_python_module(parse parse.i parse.cpp)
swig_link_libraries(parse "itkgdcm")
install_swig_python_module(parse)

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

wrap_library("AssembleTilesImageFilter" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    set(WRAPPER_SWIG_LIBRARY_FILES ${WRAPPER_SWIG_LIBRARY_FILES} 
        "${CMAKE_CURRENT_SOURCE_DIR}/itkAssembleTilesImageFilter.swg")
    auto_include_modules()
end_wrap_library()
wrap_ikt_post_install("AssembleTilesImageFilter")
