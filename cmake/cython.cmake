find_package(Cython)
if(CYTHON_FOUND)
    include(UseCython)
    
    function(medipy_cython_add_module _name )
        cython_add_module(${_name} ${ARGN})
        # Install the module in the same directory as its source file, and make
        # sure no "lib" prefix is added to the target file.
        set_target_properties(${_name} PROPERTIES 
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            PREFIX ""
        )
    endfunction()

endif()
