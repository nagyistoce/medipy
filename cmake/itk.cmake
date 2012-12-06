# Call this /before/ overloading WRAP_ITK_INSTALL
find_package(WrapITK REQUIRED)

macro(WRAP_ITK_INSTALL path)
    # Overload the WrapITK macro so that we install the files in the same
    # directory as the regular Python modules
    foreach(_file ${ARGN})
        # Install only Python-related files
        if(NOT (("${_file}" MATCHES ".*\\.i$") OR
                ("${_file}" MATCHES ".*\\.idx$") OR
                ("${_file}" MATCHES ".*\\.includes$") OR
                ("${_file}" MATCHES ".*\\.mdx$") OR
                ("${_file}" MATCHES ".*\\.swg$")
          )) 
            file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
            get_filename_component(destination ${destination} PATH)
            install(FILES ${_file} DESTINATION ${destination})
        endif()
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
        file(REMOVE_RECURSE ${CMAKE_INSTALL_PREFIX}/${WRAP_ITK_INSTALL_PREFIX})
    ")
endmacro(wrap_ikt_post_install)

wrap_type("itk::Matrix" "M")
    unique(types "${WRAP_ITK_SCALAR};D")
    foreach(d ${WRAP_ITK_DIMS})
        foreach(type ${types})
            add_template("${ITKM_${type}}${d}${d}"  "${ITKT_${type}},${d},${d}")
        endforeach()
    endforeach()
end_wrap_type()
set(itk_Wrap_Matrix ${WRAPPER_TEMPLATES})
