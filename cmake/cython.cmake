find_package(Cython)
if(CYTHON_FOUND)
    include(UseCython)
    
    function(medipy_cython_add_module _name )
        cython_add_module(${_name} ${ARGN})

        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/dummy)
        get_filename_component(destination ${destination} PATH)
        set_target_properties(${_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${destination})
    endfunction()

endif()
