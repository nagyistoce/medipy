add_subdirectory("dataset_io")
add_subdirectory("routing")

find_package(WrapITK REQUIRED)

file(GLOB files "*.py")
install_py_files(${files})

prepare_wrap_library()
wrap_library("AssembleTilesImageFilter" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    wrap_module("itkAssembleTilesImageFilter")
        wrap_class("itk::AssembleTilesImageFilter" POINTER)
            foreach(type ${WRAP_ITK_SCALAR})
                set(dimensions "2")
                foreach(dimension ${dimensions})
                    math(EXPR next_dimension "${dimension}+1")
                    
                    set(input_image_type ${ITKT_I${type}${dimension}})
                    set(output_image_type ${ITKT_I${type}${next_dimension}})
                    
                    set(template_params "${input_image_type}, ${output_image_type}")
                    
                    set(input_mangled_type ${ITKM_I${type}${dimension}})
                    set(output_mangled_type ${ITKM_I${type}${next_dimension}})
                    set(mangled_type "${input_mangled_type}${output_mangled_type}")
                    
                    wrap_template("${mangled_type}" "${template_params}")
                endforeach()
            endforeach(type)
        end_wrap_class()
    end_wrap_module()
end_wrap_library()
wrap_ikt_post_install("AssembleTilesImageFilter")
