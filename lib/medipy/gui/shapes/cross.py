##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from vtk import vtkCellArray, vtkPoints, vtkPolyDataMapper, vtkPolyData

from medipy.gui.shapes import Shape

class Cross(Shape):
    """ 2D X-cross
    """
    
    def __init__(self, position=None, color=None, size=None, filled=None):
        """ Size is the half-length of the branches
        """
        Shape.__init__(self, position, color, size, filled)
        
        points = vtkPoints()
        points.InsertNextPoint(-0.707, -0.707, 0)
        points.InsertNextPoint(+0.707, +0.707, 0)
        points.InsertNextPoint(-0.707, +0.707, 0)
        points.InsertNextPoint(+0.707, -0.707, 0)
        
        lines = vtkCellArray()
        lines.InsertNextCell(2)
        lines.InsertCellPoint(0)
        lines.InsertCellPoint(1)
        lines.InsertNextCell(2)
        lines.InsertCellPoint(2)
        lines.InsertCellPoint(3)
        
        polydata = vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)
        
        mapper = vtkPolyDataMapper()
        mapper.SetInput(polydata)
        self._actor.SetMapper(mapper)
        # Use flat lightning
        self._actor.GetProperty().SetAmbient(1)
        
        if size is not None :
            self._set_size(size)
    
    ##############
    # Properties #
    ##############
    
    def _set_size(self, size):
        Shape._set_size(self, size)
        self._actor.SetScale(2.*size)
