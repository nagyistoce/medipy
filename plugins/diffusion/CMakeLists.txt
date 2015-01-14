file(GLOB files "*.py" "gui/*.py")
install_py_files(${files})

# The underlying filters are VectorImage->VectorImage, and 
# itk::ImageToImageFilter is not wrapped for this case. Simple Swig modules are
# created to solve this.
set(swig_modules itkutils spectral_analysis)
foreach(swig_module ${swig_modules})
    swig_python_module(${swig_module} ${swig_module}.i ${swig_module}.cpp)
    swig_link_libraries(${swig_module} ${ITK_LIBRARIES})
    install_swig_python_module(${swig_module})
endforeach()

find_package(Cython)
if(CYTHON_FOUND)
    include(UseCython)
    medipy_cython_add_module(fiber_clustering fiber_clustering.pyx)
endif()

set(instantiations )
foreach(type ${WRAP_ITK_REAL})
    foreach(dimension ${ITK_WRAP_DIMS})
        set(instantiations ${instantiations} "${type}${dimension}")
    endforeach()
endforeach()

itk_wrap_module("MediPyDiffusion")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon ITKIOImageBase)
    
    find_swig_library_files("TractographyAlgorithm")
    
    # IO : (no input) -> VectorImage
    foreach(filter "Tensor2ImageFileReader" "Tensor2ImageFileWriter")
        itk_wrap_submodule("itk${filter}")
            itk_wrap_class("itk::${filter}" POINTER)
                foreach(instantiation ${instantiations})
                    itk_wrap_template("VI${instantiation}" ${ITKT_VI${instantiation}})
                endforeach()
            itk_end_wrap_class()
        itk_end_wrap_submodule()
    endforeach()

    # Estimation : Image -> (VectorImage, Image)
    foreach(filter "TensorReconstructionImageFilter" "WeightedLeastSquaresImageFilter")
        itk_wrap_submodule("itk${filter}")
            itk_wrap_class("itk::${filter}" POINTER)
                foreach(instantiation ${instantiations})
                    itk_wrap_template(
                        "I${instantiation}VI${instantiation}I${instantiation}" 
                        "${ITKT_I${instantiation}}, ${ITKT_VI${instantiation}}, ${ITKT_I${instantiation}}")
                endforeach()
            itk_end_wrap_class()
        itk_end_wrap_submodule()
    endforeach()
    
    # WeightedMeanImageFilter : Image -> VectorImage (float,3) only
    itk_wrap_submodule("itkWeightedMeanImageFilter")
        itk_wrap_class("itk::WeightedMeanImageFilter" POINTER)
            itk_wrap_template("IF3VIF3VIF3" "${ITKT_IF3}, ${ITKT_VIF3}, ${ITKT_VIF3}")
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    # Scalars : VectorImage -> Image
    foreach(scalar "AxialDiffusivity" "FractionalAnisotropy" "MeanDiffusivity" "RadialDiffusivity")
        set(filter "${scalar}ImageFilter")
        itk_wrap_submodule("itk${filter}")
            itk_wrap_class("itk::${filter}" POINTER_WITH_SUPERCLASS)
                foreach(instantiation ${instantiations})
                    itk_wrap_template(
                        "VI${instantiation}I${instantiation}" 
                        "${ITKT_VI${instantiation}}, ${ITKT_I${instantiation}}")
                endforeach()
            itk_end_wrap_class()
        itk_end_wrap_submodule()
    endforeach()
    
    # Tractography : VectorImage -> Image, (float,3) only
    foreach(filter "TractographyAlgorithm" "StreamlineTractographyAlgorithm")
        itk_wrap_submodule("itk${filter}")
            itk_wrap_class("itk::${filter}" POINTER)
                itk_wrap_template("VIF3IF3" "${ITKT_VIF3}, ${ITKT_IF3}")
            itk_end_wrap_class()
        itk_end_wrap_submodule()
    endforeach()
    
    # Decomposition Spectrale : VectorImage -> VectorImage, (float,3) only
    itk_wrap_submodule("itkComposeSpectralImageFilter")
        itk_wrap_class("itk::ComposeSpectralImageFilter" POINTER)
            itk_wrap_template("VIF3VIF3" "${ITKT_VIF3}, ${ITKT_VIF3}")
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    
    # DWIStatistics : Vector -> (VectorImage, Image) (float,3) for Spatial
    # and Image -> ( VectorImage, Image)(float,3) only for Bootstrap
    itk_wrap_submodule("itkDWIStatisticsImageFilter")
        itk_wrap_class("itk::DWIStatisticsImageFilter" POINTER)
            itk_wrap_template("VIF3VIF3IF3IF3" "${ITKT_VIF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
            itk_wrap_template("IF3VIF3IF3IF3" "${ITKT_IF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    # SpatialDWIStatistics : Vector -> (VectorImage, Image) (float,3) only
    itk_wrap_submodule("itkSpatialDWIStatisticsImageFilter")
        itk_wrap_class("itk::SpatialDWIStatisticsImageFilter" POINTER)
            itk_wrap_template("VIF3VIF3IF3IF3" "${ITKT_VIF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    # BootstrapDWIStatistics : Image -> (VectorImage, Image) (float,3) only
    itk_wrap_submodule("itkBootstrapDWIStatisticsImageFilter")
        itk_wrap_class("itk::BootstrapDWIStatisticsImageFilter" POINTER)
            itk_wrap_template("IF3VIF3IF3IF3" "${ITKT_IF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
itk_end_wrap_module()
