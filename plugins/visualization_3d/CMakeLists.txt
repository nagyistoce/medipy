file(GLOB files "*.py")
install_py_files(${files})

set(sources vtkNormalsByPlaneFitting.cpp)
vtk_python_module(normals "${sources}")
install_vtk_python_module(normals)
