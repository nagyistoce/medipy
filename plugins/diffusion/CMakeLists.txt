file(GLOB files "*.py" "gui/*.py")
install_py_files(${files})

# The underlying filters are VectorImage->VectorImage, and 
# itk::ImageToImageFilter is not wrapped for this case. Simple Swig modules are
# created to solve this.
set(swig_modules itkutils spectral_analysis)
foreach(swig_module ${swig_modules})
    swig_python_module(${swig_module} ${swig_module}.i ${swig_module}.cpp)
    swig_link_libraries(${swig_module} ${ITK_LIBS})
    install_swig_python_module(${swig_module})
endforeach()

find_package(Cython)
if(CYTHON_FOUND)
    include(UseCython)
    medipy_cython_add_module(fiber_clustering fiber_clustering.pyx)
endif()

set(instantiations )
foreach(type ${WRAP_ITK_REAL})
    foreach(dimension ${WRAP_ITK_DIMS})
        set(instantiations ${instantiations} "${type}${dimension}")
    endforeach()
endforeach()

wrap_library("MediPyDiffusion" python)
    set(WRAPPER_LIBRARY_DEPENDS Base IO)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    find_swig_library_files("TractographyAlgorithm")
    
    # IO : (no input) -> VectorImage
    foreach(filter "Tensor2ImageFileReader" "Tensor2ImageFileWriter")
        wrap_module("itk${filter}")
            wrap_class("itk::${filter}" POINTER)
                foreach(instantiation ${instantiations})
                    wrap_template("VI${instantiation}" ${ITKT_VI${instantiation}})
                endforeach()
            end_wrap_class()
        end_wrap_module()
    endforeach()
    
    # Estimation : Image -> VectorImage
    foreach(filter "SecondOrderSymmetricTensorReconstructionFilter")
        wrap_module("itk${filter}")
            wrap_class("itk::${filter}" POINTER)
                foreach(instantiation ${instantiations})
                    wrap_template(
                        "I${instantiation}VI${instantiation}" 
                        "${ITKT_I${instantiation}}, ${ITKT_VI${instantiation}}")
                endforeach()
            end_wrap_class()
        end_wrap_module()
    endforeach()
    
    # Estimation : Image -> VectorImage
    foreach(filter "WeightedLeastSquaresImageFilter")
        wrap_module("itk${filter}")
            wrap_class("itk::${filter}" POINTER)
                foreach(instantiation ${instantiations})
                    wrap_template(
                        "I${instantiation}VI${instantiation}" 
                        "${ITKT_I${instantiation}}, ${ITKT_VI${instantiation}}")
                endforeach()
            end_wrap_class()
        end_wrap_module()
    endforeach()
    
    # Scalars : VectorImage -> Image
    foreach(filter "FractionalAnisotropyImageFilter" "MeanDiffusivityImageFilter")
        wrap_module("itk${filter}")
            wrap_class("itk::${filter}" POINTER)
                foreach(instantiation ${instantiations})
                    wrap_template(
                        "VI${instantiation}I${instantiation}" 
                        "${ITKT_VI${instantiation}}, ${ITKT_I${instantiation}}")
                endforeach()
            end_wrap_class()
        end_wrap_module()
    endforeach()
    
    # Tractography : VectorImage -> Image, (float,3) only
    foreach(filter "TractographyAlgorithm" "StreamlineTractographyAlgorithm")
        wrap_module("itk${filter}")
            wrap_class("itk::${filter}" POINTER)
                wrap_template("VIF3IF3" "${ITKT_VIF3}, ${ITKT_IF3}")
            end_wrap_class()
        end_wrap_module()
    endforeach()
    
    # DWIStatistics : Vector -> (VectorImage, Image) (float,3) for Spatial
    # and Image -> ( VectorImage, Image)(float,3) only for Bootstrap
    wrap_module("itkDWIStatisticsImageFilter")
        wrap_class("itk::DWIStatisticsImageFilter" POINTER)
            wrap_template("VIF3VIF3IF3IF3" "${ITKT_VIF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
            wrap_template("IF3VIF3IF3IF3" "${ITKT_IF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        end_wrap_class()
    end_wrap_module()
    
    # SpatialDWIStatistics : Vector -> (VectorImage, Image) (float,3) only
    wrap_module("itkSpatialDWIStatisticsImageFilter")
        wrap_class("itk::SpatialDWIStatisticsImageFilter" POINTER)
            wrap_template("VIF3VIF3IF3IF3" "${ITKT_VIF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        end_wrap_class()
    end_wrap_module()
    
    # BootstrapDWIStatistics : Image -> (VectorImage, Image) (float,3) only
    wrap_module("itkBootstrapDWIStatisticsImageFilter")
        wrap_class("itk::BootstrapDWIStatisticsImageFilter" POINTER)
            wrap_template("IF3VIF3IF3IF3" "${ITKT_IF3}, ${ITKT_VIF3}, ${ITKT_IF3}, ${ITKT_IF3}")
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPyDiffusion")
