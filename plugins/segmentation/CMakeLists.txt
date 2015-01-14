file(GLOB_RECURSE files "*.py")
install_py_files(${files})

file(GLOB_RECURSE files "${CMAKE_CURRENT_SOURCE_DIR}/resources/*")
foreach(file ${files})
    get_filename_component(leaf ${file} NAME)
    if(NOT (${leaf} STREQUAL "CMakeLists.txt"))
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${file} DESTINATION ${destination})
    endif()
endforeach()

find_package("VTK" REQUIRED)
get_directory_property(include_directories INCLUDE_DIRECTORIES)
set(include_directories ${include_directories} ${VTK_INCLUDE_DIRS})
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${include_directories}")

itk_wrap_module("MediPySegmentation")
    set(WRAPPER_LIBRARY_DEPENDS ITKCommon)
    
    itk_wrap_submodule("itkBETImageFilter")
        itk_wrap_class("itk::BETImageFilter" "POINTER")
            # Wrap only for 3D images
            itk_wrap_image_filter("${WRAP_ITK_SCALAR}" 2 3)
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    itk_wrap_submodule("itkChangeDetectionClusteringImageFilter")
        itk_wrap_class("itk::ChangeDetectionClusteringImageFilter" "POINTER")
            itk_wrap_image_filter("${WRAP_ITK_REAL}" 3)
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    itk_wrap_submodule("itkClustersToAnnotationsCalculator")
        itk_wrap_class("itk::ClustersToAnnotationsCalculator" "POINTER")
            itk_wrap_image_filter("${WRAP_ITK_REAL}" 1)
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
    itk_wrap_submodule("itkMonoModalStatisticalChangeDetectionImageFilter")
        itk_wrap_class("itk::MonoModalStatisticalChangeDetectionImageFilter" "POINTER")
            itk_wrap_image_filter("${WRAP_ITK_REAL}" 3)
        itk_end_wrap_class()
    itk_end_wrap_submodule()
    
itk_end_wrap_module()

target_link_libraries("MediPySegmentationPython" vtkHybrid)
