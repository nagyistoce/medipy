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
