# Call this /before/ overloading WRAP_ITK_INSTALL
find_package(WrapITK REQUIRED)

macro(WRAP_LIBRARY library_name)
    # Overload the WrapITK macro to tweak some settings
    
    # Disable generation of documentation
    set(WRAP_ITK_DOC OFF) 
    
    # Old WrapITK does not work with g++-4.6
    set(CMAKE_CXX_COMPILER g++-4.4)
    
    # Original WrapITK code
    SET(WRAPPER_LIBRARY_NAME "${library_name}")

    MESSAGE(STATUS "${WRAPPER_LIBRARY_NAME}: Creating library.")

    # Mark the current source dir for inclusion because it may contain header files.
    INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")
    INCLUDE_DIRECTORIES(BEFORE ${WRAPPER_LIBRARY_INCLUDE_DIRECTORIES})

    # WRAPPER_LIBRARY_INCLUDE_DIRECTORIES. List of other include directories that
    # contain the desired header files.
    #SET(WRAPPER_LIBRARY_INCLUDE_DIRECTORIES )

    # WRAPPER_LIBRARY_SOURCE_DIR. Directory to be scanned for wrap_*.cmake files. 
    SET(WRAPPER_LIBRARY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

    # WRAPPER_LIBRARY_OUTPUT_DIR. Directory in which generated cxx, xml, and idx
    # files will be placed. 
    SET(WRAPPER_LIBRARY_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")

    # WRAPPER_LIBRARY_DEPENDS. List of names of other wrapper libraries that
    # define symbols used by this wrapper library.
    SET(WRAPPER_LIBRARY_DEPENDS )

    # WRAPPER_LIBRARY_LINK_LIBRARIES. List of other libraries that should
    # be linked to the wrapper library.
    SET(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)

    # WRAPPER_LIBRARY_GROUPS. List of wrap_*.cmake groups in the source dir
    # that should be included/wrapped before the rest. Just the group name is needed,
    # not the full path or file name. 
    SET(WRAPPER_LIBRARY_GROUPS )

    # WRAPPER_LIBRARY_CABLESWIG_INPUTS. List of C++ source files to be used
    # as input for CableSwig. This list is then appended to by
    # WRAPPER_LIBRARY_AUTO_INCLUDE_WRAP_FILES. A full path to each input is required.
    SET(WRAPPER_LIBRARY_CABLESWIG_INPUTS )

    # WRAPPER_SWIG_LIBRARY_FILES. List of swig .swg files to pass to cswig to control
    # type handling and so forth. A full path to each include is required.
    # The itk.swg file and the library file for the current library are implicitly added.
    SET(WRAPPER_SWIG_LIBRARY_FILES )

    # WRAPPER_LIBRARY_SWIG_INPUTS. SWIG input files to be fed to swig (not
    # CableSwig). A full path to each input is required.
    SET(WRAPPER_LIBRARY_SWIG_INPUTS ) 

    # WRAPPER_LIBRARY_CXX_SOURCES. C++ sources to be compiled and linked in
    # to the wrapper library (with no prior processing by swig, etc.)
    # A full path to each input is required.
    SET(WRAPPER_LIBRARY_CXX_SOURCES ) 

    IF("${ARGC}" EQUAL 2)
        FOREACH(lang ${WRAP_ITK_LANGUAGES})
            STRING(TOUPPER ${lang} LANG)
            SET(WRAPPER_LIBRARY_${LANG} OFF)
        ENDFOREACH(lang)
        FOREACH(lang ${ARGV1})
            STRING(TOUPPER ${lang} LANG)
            SET(WRAPPER_LIBRARY_${LANG} ON)
        ENDFOREACH(lang)
    ELSE("${ARGC}" EQUAL 2)
        FOREACH(lang ${WRAP_ITK_LANGUAGES})
            STRING(TOUPPER ${lang} LANG)
            SET(WRAPPER_LIBRARY_${LANG} ON)
        ENDFOREACH(lang)
    ENDIF("${ARGC}" EQUAL 2)

    IF("${WRAPPER_LIBRARY_WRAP_LIBRARIES_STATUS}" STREQUAL "NOT_EXECUTED")
        WRAP_LIBRARIES()
        # change the status of WRAPPER_LIBRARY_WRAP_LIBRARIES_STATUS, so we can call END_WRAP_LIBRARIES when
        # END_WRAP_LIBRARY will be called
        SET(WRAPPER_LIBRARY_WRAP_LIBRARIES_STATUS "EXECUTED_IN_WRAP_LIBRARY" CACHE INTERNAL "status var used to avoid the use of WRAP_LIBRARIES in simple contributions.")
    ENDIF("${WRAPPER_LIBRARY_WRAP_LIBRARIES_STATUS}" STREQUAL "NOT_EXECUTED")

    # Call the language support initialization function
    WRAP_LIBRARY_ALL_LANGUAGES("${library_name}")
endmacro()

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
