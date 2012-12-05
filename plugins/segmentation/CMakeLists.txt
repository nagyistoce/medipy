file(GLOB_RECURSE files "*.py")
install_py_files(${files})

file(GLOB_RECURSE files ${CMAKE_CURRENT_SOURCE_DIR} "*")
foreach(file ${files})
    
    get_filename_component(leaf ${file} NAME)
    if(NOT ((${leaf} STREQUAL "CMakeLists.txt") OR (${leaf} STREQUAL "SConstruct")))
        #message("${file} -> ${leaf}")
        file(RELATIVE_PATH destination ${CMAKE_SOURCE_DIR} ${file})
        get_filename_component(destination ${destination} PATH)
        install(FILES ${file} DESTINATION ${destination})
    endif()
endforeach(file)

# Old WrapITK does not work with g++-4.6
set(CMAKE_CXX_COMPILER g++-4.4)

set(WRAP_ITK_DOC OFF)

find_package("VTK" REQUIRED)
get_directory_property(include_directories INCLUDE_DIRECTORIES)
set(include_directories ${include_directories} ${VTK_INCLUDE_DIRS})
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${include_directories}")

wrap_library("MediPySegmentation" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    set(WRAPPER_SWIG_LIBRARY_FILES ${WRAPPER_SWIG_LIBRARY_FILES}
        "${CMAKE_CURRENT_SOURCE_DIR}/itkBETImageFilter.swg" 
        "${CMAKE_CURRENT_SOURCE_DIR}/itkChangeDetectionClusteringImageFilter.swg"
        "${CMAKE_CURRENT_SOURCE_DIR}/itkClustersToAnnotationsCalculator.swg"
        "${CMAKE_CURRENT_SOURCE_DIR}/itkMonoModalStatisticalChangeDetectionImageFilter.swg"
    )
    auto_include_modules()
end_wrap_library()
wrap_ikt_post_install("MediPySegmentation")

target_link_libraries("MediPySegmentationPython" vtkHybrid)