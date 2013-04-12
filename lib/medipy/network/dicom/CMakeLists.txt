add_subdirectory("scu")

file(GLOB files "*.py")
install_py_files(${files})
