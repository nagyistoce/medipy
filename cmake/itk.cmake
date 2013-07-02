# Call this /before/ overloading WRAP_ITK_INSTALL
find_package(WrapITK REQUIRED)

macro(prepare_wrap_library)
    # Overload the WrapITK macro to tweak some settings
    
    # Disable generation of documentation
    set(WRAP_ITK_DOC OFF) 
    
    # Old WrapITK does not work with g++-4.6
    set(CMAKE_CXX_COMPILER g++-4.4)
endmacro()

macro(wrap_ikt_post_install library)

    # This is not really post-install, but used to make sure that the module
    # config file is in the right place.
    configure_file("${WRAP_ITK_PYTHON_SOURCE_DIR}/ModuleConfig.py.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${WRAPPER_LIBRARY_NAME}Config.py"
        @ONLY IMMEDIATE)

    # Post-install to make sure that everything is installed in the same 
    # directory as the regular Python modules.
    file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
    get_filename_component(destination ${destination} PATH)
    
    get_target_property(lib_location "${library}Python" LOCATION)
    install(CODE "
        file(COPY 
            ${lib_location} DESTINATION ${CMAKE_INSTALL_PREFIX}/${destination} 
            USE_SOURCE_PERMISSIONS)
        file(REMOVE_RECURSE ${CMAKE_INSTALL_PREFIX}/${WRAP_ITK_INSTALL_PREFIX})
    ")
endmacro()

macro(find_swig_library_files)
    foreach(module ${ARGN})
        set(WRAPPER_SWIG_LIBRARY_FILES 
            ${WRAPPER_SWIG_LIBRARY_FILES} "${CMAKE_CURRENT_SOURCE_DIR}/itk${module}.swg")
    endforeach()
endmacro()

wrap_type("itk::Matrix" "M")
    unique(types "${WRAP_ITK_SCALAR};D")
    foreach(d ${WRAP_ITK_DIMS})
        foreach(type ${types})
            add_template("${ITKM_${type}}${d}${d}"  "${ITKT_${type}},${d},${d}")
        endforeach()
    endforeach()
end_wrap_type()
set(itk_Wrap_Matrix ${WRAPPER_TEMPLATES})

# Reset the INCLUDE_DIRECTORIES to avoid propagation to sub-directories
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "")
