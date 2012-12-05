function(install_py_files)
    # Install the given .py files to a directory mirroring their location.
    foreach(python_file ${ARGV})
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${python_file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${python_file} DESTINATION ${destination})
    endforeach(python_file ${ARGV})
endfunction(install_py_files)

macro(swig_python_module name)
    # Build a SWIG Python module. This macro is a wrapper around swig_add_module.
    find_package(SWIG REQUIRED)
    include(${SWIG_USE_FILE})

    find_package(PythonLibs)
    include_directories(${PYTHON_INCLUDE_PATH})
    
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    
    swig_add_module(${name} python ${ARGV})
endmacro(swig_python_module name)

macro(install_swig_python_module name)
    # Install the files of a SWIG Python module
    
    # Get path to the generated library
    get_target_property(lib_location ${SWIG_MODULE_${name}_REAL_NAME} LOCATION)
    # Install the library and its companion .py file
    set(files ${CMAKE_CURRENT_BINARY_DIR}/${name}.py ${lib_location})
    
    # Destination : mirror current location
    file(RELATIVE_PATH destination ${CMAKE_BINARY_DIR} ${lib_location})
    get_filename_component(destination ${destination} PATH)
    
    install(FILES ${files} DESTINATION ${destination})
endmacro(install_swig_python_module name)
