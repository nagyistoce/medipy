set(subdirectories base gui io itk network resources vtk)
foreach(sudirectory ${subdirectories})
    add_subdirectory(${sudirectory})
endforeach()

file(GLOB files "*.py")
install_py_files(${files})
