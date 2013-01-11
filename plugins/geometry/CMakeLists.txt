file(GLOB files "*.py")
install_py_files(${files})

wrap_library("MediPyGeometry" python)
    set(WRAPPER_LIBRARY_DEPENDS Base)
    set(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon)
    
    wrap_module("itkPlane")
        wrap_class("itk::Plane")
            foreach(type "F" "D")
                wrap_template(${type} ${ITKT_${type}})
            endforeach()
        end_wrap_class()
    end_wrap_module()
    
    wrap_module("itkSymmetryPlaneCalculator")
        wrap_class("itk::SymmetryPlaneCalculator")
            wrap_image_filter("${WRAP_ITK_SCALAR}" 1 3)
        end_wrap_class()
    end_wrap_module()
    
end_wrap_library()
wrap_ikt_post_install("MediPyGeometry")
