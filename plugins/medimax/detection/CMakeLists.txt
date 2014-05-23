file(GLOB files "*.py")
install_py_files(${files})

set_source_files_properties(detection.i PROPERTIES CPLUSPLUS ON)
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/..")
swig_python_module(detection detection.i)
swig_link_libraries(detection medimax)
install_swig_python_module(detection)
