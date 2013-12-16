file(GLOB files "*.py")
install_py_files(${files})

set(sources vtkColorTransferFunctionWithAlpha.cxx vtkEnhancedLookupTable.cpp
    vtkMisc.cpp vtkOrientationAnnotation.cxx)
vtk_python_module(vtk_addons "${sources}")
install_vtk_python_module(vtk_addons)
