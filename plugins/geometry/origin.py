##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

def set_origin_as_index(image, index):
    """ Set the origin of the image to the specified index.
        
        <gui>
            <item name="image" type="Image" label="Image" />
            <item name="index" type="Coordinates" 
                  initializer="image=${image}, display_coordinates='index'"
                  label="Index" />
        </gui>
    """
    
    image.origin = -(index*image.spacing)
