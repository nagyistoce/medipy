file(GLOB files "*.py")
install_py_files(${files})

set(sources vtkNormalsByPlaneFitting.cpp)
vtk_python_module(normals "${sources}")

find_package("ITK" REQUIRED)
target_link_libraries(normals itkvnl_algo)

install_vtk_python_module(normals)
