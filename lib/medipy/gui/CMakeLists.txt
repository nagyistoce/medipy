set(subdirectories annotations base brushes control dicom function_gui_builder
                   image io shapes xrc_wrapper)
foreach(subdirectory ${subdirectories})
    add_subdirectory("${subdirectory}")
endforeach()

file(GLOB files "*.py")
install_py_files(${files})
