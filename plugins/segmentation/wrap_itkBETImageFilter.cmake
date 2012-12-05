WRAP_CLASS("itk::BETImageFilter" "POINTER")
    # Wrap only for 3D images
    WRAP_IMAGE_FILTER("${WRAP_ITK_SCALAR}" 2 3)
END_WRAP_CLASS()
