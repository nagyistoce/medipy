find_package("ITK" REQUIRED)
find_package("VTK" REQUIRED)

file(GLOB files "*.py")
install_py_files(${files})

include(${ITK_USE_FILE})

set(sources vtkNormalsByPlaneFitting.cpp)
vtk_python_module(normals "${sources}")



target_link_libraries(normals itkvnl_algo)

install_vtk_python_module(normals)
