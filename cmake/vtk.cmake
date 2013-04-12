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
    
    if("${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" STRLESS "5.5")
        # Fix the module name in the Init file
        add_custom_command(
            OUTPUT ${name}InitModified.cxx
            DEPENDS ${name}Init.cxx
            COMMAND python -c 
                "import sys; open(sys.argv[2], 'w').write(open(sys.argv[1]).read().replace('lib${name}', '${name}'))" 
                ${name}Init.cxx ${name}InitModified.cxx
            VERBATIM
        )
        set(init_file ${name}InitModified.cxx)
        
        # Old versions of VTK (Ubuntu 10.04, Debian 5.0) do not define VTK_LIBRARIES
        set(VTK_LIBRARIES vtkCommon vtkFiltering vtkGenericFiltering vtkGraphics
                          vtkHybrid vtkImaging vtkInfovis vtkIO vtkParallel
                          vtkRendering vtkViews vtkVolumeRendering vtkWidgets)
    else()
        set(init_file ${name}Init.cxx)
    endif()

    # Build the module using the modified Init file
    python_add_module(${name} ${sources} ${wrappers} ${init_file})
    
    set(vtk_python_libraries)
    foreach(lib ${VTK_LIBRARIES})
        set(vtk_python_libraries ${vtk_python_libraries} ${lib}PythonD)
    endforeach()
    
    target_link_libraries(${name} ${VTK_LIBRARIES} ${vtk_python_libraries})
    set_target_properties(${name} PROPERTIES PREFIX "")
endmacro(vtk_python_module)

macro(install_vtk_python_module name)
    # Get path to the generated library
    get_target_property(lib_location ${name} LOCATION)
    # Destination : mirror current location
    file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
    get_filename_component(destination ${destination} PATH)
    
    install(TARGETS ${name} DESTINATION ${destination})
endmacro(install_vtk_python_module)
