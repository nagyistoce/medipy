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
