file(GLOB files "*.py")
install_py_files(${files})

file(GLOB plugins "*/CMakeLists.txt")
foreach(plugin ${plugins})
    get_filename_component(directory ${plugin} PATH)
    add_subdirectory(${directory})
endforeach(plugin ${plugins})
