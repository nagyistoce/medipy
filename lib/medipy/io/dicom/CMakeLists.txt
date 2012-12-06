file(GLOB files "*.py")
install_py_files(${files})

# Build "parse" Python module
find_package(ITK)
if(${ITK_FOUND})
    if(${ITK_USE_SYSTEM_GDCM})
        find_package(GDCM REQUIRED)
    else(${ITK_USE_SYSTEM_GDCM})
        include_directories(${ITK_INCLUDE_DIRS})
        link_directories(${ITK_LIBRARY_DIRS})
    endif(${ITK_USE_SYSTEM_GDCM})
else(${ITK_FOUND})
    find_package(GDCM REQUIRED)
endif(${ITK_FOUND})

set_source_files_properties(dataset_io.i PROPERTIES CPLUSPLUS ON)
swig_python_module(dataset_io dataset_io.i dataset_io.cpp DataSetBridge.cpp)
swig_link_libraries(dataset_io "itkgdcm")
install_swig_python_module(dataset_io)

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

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
