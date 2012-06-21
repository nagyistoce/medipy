##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

def set_origin_as_index(image, index):
    """ Set the origin of the image to the specified index.
        
        <gui>
            <item name="image" type="Image" label="Image" />
            <item name="index" type="Coordinates" 
                  initializer="image=${image}, display_coordinates='index'"
                  label="Index" />
        </gui>
    """
    
    image.origin = -numpy.dot(image.direction, index*image.spacing)

def set_origin_as_point(image, point):
    """ Set the origin of the image to the specified point.
        
        <gui>
            <item name="image" type="Image" label="Image" />
            <item name="point" type="Coordinates" 
                  initializer="image=${image}, display_coordinates='physical'" 
                  label="Point" />
        </gui>
    """
    
    pass
    # TODO : transform point to index, call set_origin_as_index(image, index)