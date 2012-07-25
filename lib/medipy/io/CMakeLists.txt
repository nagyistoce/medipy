add_subdirectory("dicom")
add_subdirectory("schemes")

file(GLOB files "*.py")
install_py_files(${files})
