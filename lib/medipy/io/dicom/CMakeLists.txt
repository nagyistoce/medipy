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

# TODO : build "AssembleTilesImageFilter" WrapITK filters
