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
