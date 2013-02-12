##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from numpy_bridge import (array_to_itk_matrix, itk_matrix_to_array,
                          array_to_itk_image, array_to_itk_vector_image,
                          itk_image_to_array, itk_vector_image_to_array, 
                          medipy_image_to_itk_image, itk_image_to_medipy_image, 
                          itk_to_dtype, dtype_to_itk, itk_image_type)
from loader import load_wrapitk_module
import types

import os.path
load_wrapitk_module(os.path.dirname(__file__), "MediPyBridge")

__all__ = ["array_to_itk_matrix", "itk_matrix_to_array",
           "array_to_itk_image", "array_to_itk_vector_image",
           "itk_image_to_array", "itk_vector_image_to_array", 
           "medipy_image_to_itk_image", "itk_image_to_medipy_image", 
           "itk_to_dtype", "dtype_to_itk", "itk_image_type",
           "load_wrapitk_module", "types"]