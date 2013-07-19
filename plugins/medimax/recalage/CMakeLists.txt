file(GLOB files "*.py")
install_py_files(${files})

set_source_files_properties(recalage.i PROPERTIES CPLUSPLUS ON)
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/..")
swig_python_module(recalage recalage.i)
swig_link_libraries(recalage medimax)
install_swig_python_module(recalage)
