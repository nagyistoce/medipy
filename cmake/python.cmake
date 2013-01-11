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

    set(cplusplus "OFF")

    foreach(it ${ARGN})
        if(${it} MATCHES ".*\\.i$")
            set(swig_dot_i_sources ${swig_dot_i_sources} "${it}")
        elseif(${it} MATCHES ".*\\.cpp$")
            set(cplusplus "ON")
        elseif(${it} MATCHES ".*\\.cxx$")
            set(cplusplus "ON")
        endif()
    endforeach()
    
    foreach(it ${swig_dot_i_sources})
        get_source_file_property(value ${it} CPLUSPLUS)
        if(${value} STREQUAL "NOTFOUND")
            set_source_files_properties(${it} PROPERTIES CPLUSPLUS ${cplusplus})
        endif()
    endforeach()
    
    find_package(PythonLibs)
    include_directories(${PYTHON_INCLUDE_PATH})
    
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    
    swig_add_module(${name} python ${ARGN})
    
    swig_link_libraries(${name} ${PYTHON_LIBRARIES})
endmacro(swig_python_module name)

macro(install_swig_python_module name)
    # Install the files of a SWIG Python module
    
    # Get path to the generated library
    get_target_property(lib_location ${SWIG_MODULE_${name}_REAL_NAME} LOCATION)
    # Install the library and its companion .py file
    set(files ${CMAKE_CURRENT_BINARY_DIR}/${name}.py ${lib_location})
    
    # Destination : mirror current location
    file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
    get_filename_component(destination ${destination} PATH)
    
    install(FILES ${files} DESTINATION ${destination})
endmacro(install_swig_python_module name)
