##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from vtk import vtkPolyDataMapper, vtkRegularPolygonSource 

from medipy.gui.shapes import Shape

class Disk(Shape):
    """ 2D disk
    """
    
    def __init__(self, position=None, color=None, size=None, filled=None):
        """ Size is the radius of the disk
        """
        Shape.__init__(self, position, color, size, filled)
        
        source = vtkRegularPolygonSource()
        source.SetNumberOfSides(30)
        source.SetCenter(0,0,0)
        
        mapper = vtkPolyDataMapper()
        mapper.SetInput(source.GetOutput())
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