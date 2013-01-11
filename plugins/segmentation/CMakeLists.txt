file(GLOB_RECURSE files "*.py")
install_py_files(${files})

file(GLOB_RECURSE files "${CMAKE_CURRENT_SOURCE_DIR}/resources/*")
foreach(file ${files})
    get_filename_component(leaf ${file} NAME)
    if(NOT ((${leaf} STREQUAL "CMakeLists.txt") OR (${leaf} STREQUAL "SConstruct")))
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${file} DESTINATION ${destination})
    endif()
endforeach()

find_package("VTK" REQUIRED)
get_directory_property(include_directories INCLUDE_DIRECTORIES)
set(include_directories ${include_directories} ${VTK_INCLUDE_DIRS})
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${include_directories}")

wrap_library("MediPySegmentation" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    wrap_module("itkBETImageFilter")
        wrap_class("itk::BETImageFilter" "POINTER")
            # Wrap only for 3D images
            wrap_image_filter("${WRAP_ITK_SCALAR}" 2 3)
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkChangeDetectionClusteringImageFilter")
        wrap_class("itk::ChangeDetectionClusteringImageFilter" "POINTER")
            wrap_image_filter("${WRAP_ITK_REAL}" 3)
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkClustersToAnnotationsCalculator")
        wrap_class("itk::ClustersToAnnotationsCalculator" "POINTER")
            wrap_image_filter("${WRAP_ITK_REAL}" 1)
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkMonoModalStatisticalChangeDetectionImageFilter")
        wrap_class("itk::MonoModalStatisticalChangeDetectionImageFilter" "POINTER")
            wrap_image_filter("${WRAP_ITK_REAL}" 3)
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPySegmentation")

target_link_libraries("MediPySegmentationPython" vtkHybrid)
