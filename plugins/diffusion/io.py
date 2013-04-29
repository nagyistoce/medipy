##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import medipy.io.schemes
import itk_io
import vtk_io

# Add the IO classes to the known classes
medipy.io.schemes.file.io_classes.insert(0, itk_io.ITK)
medipy.io.schemes.file.io_classes.insert(0, vtk_io.Vtk)
