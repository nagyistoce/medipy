cmake_minimum_required(VERSION 2.8)
project("MediPy")

macro(install_py_files)
    foreach(_file ${ARGV})
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${_file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${_file} DESTINATION ${destination})
    endforeach(_file ${ARGV})
endmacro(install_py_files)

add_subdirectory("lib")
#add_subdirectory("plugins")
