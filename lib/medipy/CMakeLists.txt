add_subdirectory("base")
add_subdirectory("gui")
add_subdirectory("io")
add_subdirectory("itk")
add_subdirectory("resources")
add_subdirectory("vtk")

file(GLOB files "*.py")
install_py_files(${files})
