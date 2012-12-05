file(GLOB_RECURSE files ${CMAKE_CURRENT_SOURCE_DIR} "*")

foreach(file ${files})
    
    get_filename_component(leaf ${file} NAME)
    if(NOT ((${leaf} STREQUAL "CMakeLists.txt") OR (${leaf} STREQUAL "SConstruct")))
        #message("${file} -> ${leaf}")
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${file} DESTINATION ${destination})
    endif()
endforeach(file)