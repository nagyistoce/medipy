add_subdirectory("dataset_io")
add_subdirectory("routing")

find_package(WrapITK REQUIRED)

file(GLOB files "*.py")
install_py_files(${files})

add_definitions(-DHAVE_CONFIG_H)

find_path(OFLOG "dcmtk/oflog/oflog.h")

if(NOT ${OFLOG} STREQUAL "OFLOG-NOTFOUND")
    swig_python_module(generate_dicomdir 
        generate_dicomdir.i generate_dicomdir.cpp DicomDirGenerator.cpp)
endif()

get_directory_property(include_directories INCLUDE_DIRECTORIES)
set_directory_properties(PROPERTIES
    INCLUDE_DIRECTORIES "${include_directories};${DCMTK_INCLUDE_DIRS};${PYTHON_INCLUDE_DIRS}"
)

if("${DCMTK_LIBRARIES}" STREQUAL "")
    set(DCMTK_LIBRARIES dcmdata dcmimage dcmimgle dcmjpeg dcmnet dcmpstat 
        dcmqrdb dcmsr dcmtls ijg12 ijg16 ijg8 ofstd)
endif()
find_library(oflog oflog)
if(NOT ${oflog} STREQUAL "oflog-NOTFOUND")
    set(DCMTK_LIBRARIES ${DCMTK_LIBRARIES} oflog)
endif()

if(NOT ${OFLOG} STREQUAL "OFLOG-NOTFOUND")
    swig_link_libraries(generate_dicomdir 
        dcmdsig ${DCMTK_LIBRARIES} gdcmMSFF wrap pthread)
    install_swig_python_module(generate_dicomdir)
endif()

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
