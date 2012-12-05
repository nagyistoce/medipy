cmake_minimum_required(VERSION 2.8)
project("MediPy")

macro(install_py_files)
    # Install the given .py files to a directory mirroring their location.
    foreach(_file ${ARGV})
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${_file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${_file} DESTINATION ${destination})
    endforeach(_file ${ARGV})
endmacro(install_py_files)

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

macro(vtk_python_module name sources)
    # TODO : Use hint file ?
    
    find_package(VTK REQUIRED)
    include(${USE_VTK_FILE})
    
    get_filename_component(vtk_cmake_path ${VTK_USE_FILE} PATH)
    include(${vtk_cmake_path}/CMake/vtkWrapPython.cmake)

    find_package(PythonLibs)
    include_directories(${PYTHON_INCLUDE_PATH})

    include_directories(${CMAKE_CURRENT_SOURCE_DIR})

    set(wrappers )
    VTK_WRAP_PYTHON3(${name} wrappers "${sources}")
    
    # Fix the module name in the Init file
    add_custom_command(
        OUTPUT ${name}InitModified.cxx
        DEPENDS ${name}Init.cxx
        COMMAND python -c 
            "import sys; open(sys.argv[2], 'w').write(open(sys.argv[1]).read().replace('lib${name}', '${name}'))" 
            ${name}Init.cxx ${name}InitModified.cxx
        VERBATIM
    )

    # Build the module using the modified Init file
    python_add_module(${name} ${sources} ${wrappers} ${name}InitModified.cxx)
    target_link_libraries(${name} vtkHybrid vtkHybridPythonD)
    set_target_properties(${name} PROPERTIES PREFIX "")
endmacro(vtk_python_module name sources)

macro(install_vtk_python_module name)
    # Get path to the generated library
    get_target_property(lib_location ${name} LOCATION)
    # Destination : mirror current location
    file(RELATIVE_PATH destination ${CMAKE_BINARY_DIR} ${lib_location})
    get_filename_component(destination ${destination} PATH)
    
    install(TARGETS ${name} DESTINATION ${destination})
endmacro(install_vtk_python_module name sources)

# Call this /before/ overloading WRAP_ITK_INSTALL
find_package(WrapITK REQUIRED)

macro(WRAP_ITK_INSTALL path)
    # Overload the WrapITK macro so that we install the files in the same
    # directory as the regular Python modules
    foreach(_file ${ARGN})
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
        get_filename_component(destination ${destination} PATH)
        message("New WRAP_ITK_INSTALL ${ARGN} ${destination}")
        #message("New WRAP_ITK_INSTALL ${CMAKE_CURRENT_SOURCE_DIR} ${_file} ${destination}")
        #message("New WRAP_ITK_INSTALL ${path} ${_file} ${destination}")
        install(FILES ${_file} DESTINATION ${destination})
    endforeach(_file ${ARGN})
endmacro(WRAP_ITK_INSTALL)

macro(wrap_ikt_post_install library)
    # Post-install to make sure that everything is installed in the same 
    # directory as the regular Python modules.
    file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
    get_filename_component(destination ${destination} PATH)
    
    get_target_property(lib_location "${library}Python" LOCATION)
    install(CODE "
        file(COPY 
            ${lib_location} DESTINATION ${CMAKE_INSTALL_PREFIX}/${destination} 
            USE_SOURCE_PERMISSIONS)
        file(GLOB_RECURSE foo \"${CMAKE_INSTALL_PREFIX}/${destination}/*.py\")
        message(\"${CMAKE_INSTALL_PREFIX}/${destination}/*.py\")
        message(\"To remove : ${foo}\")#file(REMOVE ${to_remove})
        file(REMOVE_RECURSE ${CMAKE_INSTALL_PREFIX}/${WRAP_ITK_INSTALL_PREFIX})
    ")
endmacro(wrap_ikt_post_install)

add_subdirectory("lib")
add_subdirectory("plugins")
