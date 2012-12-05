WRAP_CLASS("itk::AssembleTilesImageFilter" POINTER)
    foreach(type ${WRAP_ITK_SCALAR})
        set(input_image_type ${ITKT_I${type}2})
        set(output_image_type ${ITKT_I${type}3})
        
        set(template_params "${input_image_type}, ${output_image_type}")
        
        set(input_mangled_type ${ITKM_I${type}2})
        set(output_mangled_type ${ITKM_I${type}3})
        set(mangled_type "${input_mangled_type}${output_mangled_type}")
        
        WRAP_TEMPLATE("${mangled_type}" "${template_params}")
    endforeach(type)
END_WRAP_CLASS()

