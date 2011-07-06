##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from vtk import *

def quadric_decimation(input, reduction, output):
    """ Use vtkQuadricDecimation to decimate the input object. Only triangle
        cells of the input object are processed.
    """
    
    decimation_filter = vtkQuadricDecimation()
    decimation_filter.SetInput(input.dataset)
    decimation_filter.SetTargetReduction(reduction)
    
    decimation_filter.Update()
    output.dataset = decimation_filter.GetOutput()
