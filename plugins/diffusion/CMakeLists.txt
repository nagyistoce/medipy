file(GLOB files "*.py" "gui/*.py")
install_py_files(${files})

swig_python_module(spectral_analysis spectral_analysis.i spectral_analysis.cpp)
swig_link_libraries(spectral_analysis ${ITK_LIBS})
install_swig_python_module(spectral_analysis)

set(instantiations )
foreach(type ${WRAP_ITK_REAL})
    foreach(dimension ${WRAP_ITK_DIMS})
        set(instantiations ${instantiations} "${type}${dimension}")
    endforeach()
endforeach()

wrap_library("MediPyDiffusion" python)
    set(WRAPPER_LIBRARY_DEPENDS Base IO)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
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
    
end_wrap_library()
wrap_ikt_post_install("MediPyDiffusion")
