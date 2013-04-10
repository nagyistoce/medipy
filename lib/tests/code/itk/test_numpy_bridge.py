import unittest

import itk
import numpy

import medipy.base
import medipy.itk

class TestNumpyBridge(unittest.TestCase):
    
    def setUp(self):
        self.shape = [32,31,21]
        self.shape_4d = [16,15,7]
        self.dtype = numpy.single
        
        self.origin = (1.1, 1.2, 1.3)
        self.spacing = (0.7, 0.3, 1.5)
        self.direction = medipy.base.coordinate_system.RAS
        
        self.PixelType = medipy.itk.dtype_to_itk[self.dtype]
        self.Dimension = len(self.shape)
        self.ImageType = itk.Image[self.PixelType, self.Dimension]
        self.VectorImageType = itk.VectorImage[self.PixelType, self.Dimension]
    
    def array_to_itk_image(self, transferOwnership):
        array = self._create_random_array()
        
        itk_image = medipy.itk.array_to_itk_image(array, transferOwnership)
        
        # Check result type
        self._test_type(array, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if transferOwnership :
            self.assertFalse(array.flags.owndata)
        else :
            self.assertTrue(array.flags.owndata)
        
        # Check content
        self._test_scalar_content(array, itk_image)
    
    def test_array_to_itk_image_without_transfer(self):
        self.array_to_itk_image(False)
    
    def test_array_to_itk_image_with_transfer(self):
        self.array_to_itk_image(True)
    
    def itk_image_to_array(self, transferOwnership):
        itk_image = self._create_random_itk_image()
        
        array = medipy.itk.itk_image_to_array(itk_image, transferOwnership)
        
        # Check result type
        self._test_type(array, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if transferOwnership :
            self.assertTrue(array.flags.owndata)
        else :
            self.assertFalse(array.flags.owndata)
        
        # Check content
        self._test_scalar_content(array, itk_image)
    
    def test_array_to_itk_image_without_transfer(self):
        self.itk_image_to_array(False)
    
    def test_array_to_itk_image_with_transfer(self):
        self.itk_image_to_array(True)
    
    def array_to_itk_vector_image(self, transferOwnership):
        array = self._create_random_array_4d()
        
        itk_image = medipy.itk.array_to_itk_vector_image(array, transferOwnership)
        
        # Check result type
        self._test_type(array, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if transferOwnership :
            self.assertFalse(array.flags.owndata)
        else :
            self.assertTrue(array.flags.owndata)
        
        # Check content
        self._test_vector_content(array, itk_image)
    
    def test_array_to_itk_vector_image_without_transfer(self):
        self.array_to_itk_vector_image(False)
    
    def test_array_to_itk_vector_image_with_transfer(self):
        self.array_to_itk_vector_image(True)

    def itk_vector_image_to_array(self, transferOwnership):
        itk_image = self._create_random_itk_vector_image()
        
        array = medipy.itk.itk_vector_image_to_array(itk_image, transferOwnership)
        
        # Check result type
        self._test_type(array, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if transferOwnership :
            self.assertTrue(array.flags.owndata)
        else :
            self.assertFalse(array.flags.owndata)
        
        # Check content
        self._test_vector_content(array, itk_image)

    def test_itk_vector_image_to_array_without_transfer(self):
        self.itk_vector_image_to_array(False)
    
    def test_itk_vector_image_to_array_with_transfer(self):
        self.itk_vector_image_to_array(True)

    def medipy_image_to_itk_image(self, transferOwnership):
        medipy_image = self._create_random_medipy_image()
        
        itk_image = medipy.itk.medipy_image_to_itk_image(medipy_image, transferOwnership)
        
        # Check result type
        self._test_type(medipy_image, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if transferOwnership :
            self.assertFalse(medipy_image.data.flags.owndata)
        else :
            self.assertTrue(medipy_image.data.flags.owndata)
        
        # Check origin, spacing and direction
        self._test_grid(medipy_image, itk_image)
        
        # Check content
        self._test_scalar_content(medipy_image, itk_image)
    
    def test_medipy_image_to_itk_image_without_transfer(self):
        self.medipy_image_to_itk_image(False)
    
    def test_medipy_image_to_itk_image_with_transfer(self):
        self.medipy_image_to_itk_image(True)
    
    def itk_image_to_medipy_image(self, provide_image, transferOwnership):
        itk_image = self._create_random_itk_image()
        
        if provide_image :
            medipy_image = medipy.base.Image()
            medipy.itk.itk_image_to_medipy_image(
                itk_image, medipy_image, transferOwnership)
        else :
            medipy_image = medipy.itk.itk_image_to_medipy_image(
                itk_image, None, transferOwnership)
        
        # Check result type
        self._test_type(medipy_image, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if not provide_image :
            self.assertTrue(medipy_image.data.flags.owndata)
        else :
            if transferOwnership :
                self.assertTrue(medipy_image.data.flags.owndata)
            else :
                self.assertFalse(medipy_image.data.flags.owndata)
        
        # Check origin, spacing and direction
        self._test_grid(medipy_image, itk_image)
        
        # Check content
        self._test_scalar_content(medipy_image, itk_image)
    
    def test_itk_image_to_medipy_image_with_image_without_transfer(self):
        self.itk_image_to_medipy_image(True, False)
    
    def test_itk_image_to_medipy_image_without_image_with_transfer(self):
        self.itk_image_to_medipy_image(False, True)
    
    def test_itk_image_to_medipy_image_with_image_with_transfer(self):
        self.itk_image_to_medipy_image(True, True)
    
    def medipy_image_to_itk_vector_image(self, transferOwnership):
        medipy_image = self._create_random_medipy_image_4d()
        
        itk_image = medipy.itk.medipy_image_to_itk_image(medipy_image, transferOwnership)
        
        # Check result type
        self._test_type(medipy_image, itk_image)
        
        # Check ownership. TODO : check ownership of ITK image
        if transferOwnership :
            self.assertFalse(medipy_image.data.flags.owndata)
        else :
            self.assertTrue(medipy_image.data.flags.owndata)
        
        # Check origin, spacing and direction
        self._test_grid(medipy_image, itk_image)

        # Check content
        self._test_vector_content(medipy_image, itk_image)
    
    def test_medipy_image_to_itk_vector_image_without_transfer(self):
        self.medipy_image_to_itk_vector_image(False)
#    
#    def test_medipy_image_to_itk_vector_image_with_transfer(self):
#        self.medipy_image_to_itk_vector_image(True)
#    
#    def itk_vector_image_to_medipy_image(self, provide_image, transferOwnership):
#        itk_image = self._create_random_itk_vector_image()
#        
#        if provide_image :
#            medipy_image = medipy.base.Image()
#            medipy.itk.itk_image_to_medipy_image(
#                itk_image, medipy_image, transferOwnership)
#        else :
#            medipy_image = medipy.itk.itk_image_to_medipy_image(
#                itk_image, None, transferOwnership)
#        
#        # Check result type
#        self._test_type(medipy_image, itk_image)
#        
#        # Check ownership. TODO : check ownership of ITK image
#        if not provide_image :
#            self.assertTrue(medipy_image.data.flags.owndata)
#        else :
#            if transferOwnership :
#                self.assertTrue(medipy_image.data.flags.owndata)
#            else :
#                self.assertFalse(medipy_image.data.flags.owndata)
#        
#        # Check origin, spacing and direction
#        self._test_grid(medipy_image, itk_image)
#        
#        # Check content
#        self._test_vector_content(medipy_image, itk_image)
#    
#    def test_itk_vector_image_to_medipy_image_with_image_without_transfer(self):
#        self.itk_vector_image_to_medipy_image(True, False)
#    
#    def test_itk_vector_image_to_medipy_image_without_image_with_transfer(self):
#        self.itk_vector_image_to_medipy_image(False, True)
#    
#    def test_itk_vector_image_to_medipy_image_with_image_with_transfer(self):
#        self.itk_vector_image_to_medipy_image(True, True)
    
    ####################
    # Internal helpers #
    ####################
    
    def _create_random_array(self):
        array = (1024*numpy.random.random(self.shape)).astype(self.dtype)
        self.assertTrue(array.flags.owndata)
        return array
    
    def _create_random_medipy_image(self):
        image = medipy.base.Image(data=self._create_random_array(),
            origin=self.origin, spacing=self.spacing, direction=self.direction)
        return image
    
    def _create_random_array_4d(self):
        array = (1024*numpy.random.random(self.shape_4d+[6])).astype(self.dtype)
        self.assertTrue(array.flags.owndata)
        return array
    
    def _create_random_medipy_image_4d(self):
        image = medipy.base.Image(data=self._create_random_array_4d(),
            origin=self.origin, spacing=self.spacing, direction=self.direction,
            data_type="vector"
        )
        return image
    
    def _create_random_itk_image(self):
        itk_image = self.ImageType.New()
        region = itk.ImageRegion[self.Dimension]((0,0,0), list(reversed(self.shape)))
        itk_image.SetRegions(region)
        itk_image.Allocate()
        for index in numpy.ndindex(*self.shape) :
            itk_index = list(reversed(index))
            itk_image.SetPixel(itk_index, 1024*numpy.random.random())
        
        itk_image.SetOrigin(list(reversed(self.origin)))
        itk_image.SetSpacing(list(reversed(self.spacing)))
        itk_image.SetDirection(medipy.itk.array_to_itk_matrix(self.direction, True))

        return itk_image
    
    def _create_random_itk_vector_image(self):
        itk_image = self.VectorImageType.New()
        region = itk.ImageRegion[self.Dimension]((0,0,0), list(reversed(self.shape_4d)))
        itk_image.SetRegions(region)
        itk_image.SetNumberOfComponentsPerPixel(6)
        itk_image.Allocate()
        
        for index in numpy.ndindex(*self.shape_4d) :
            itk_index = list(reversed(index))
            itk_pixel = itk_image.GetPixel(itk_index)
            
            for i in range(itk_image.GetNumberOfComponentsPerPixel()) :
                itk_pixel.SetElement(i, 1024*numpy.random.random())
        
        itk_image.SetOrigin(list(reversed(self.origin)))
        itk_image.SetSpacing(list(reversed(self.spacing)))
        itk_image.SetDirection(medipy.itk.array_to_itk_matrix(self.direction, True))

        return itk_image
    
    ###########################
    # Internal test functions #
    ###########################
    
    def _test_type(self, array, itk_image):
        """ array can be any array-like Python object (including medipy.base.Image)
        """
        
        itk_dtype = medipy.itk.dtype_to_itk[array.dtype.type]
        
        if itk_image.GetNameOfClass() == "Image" :
            self.assertEqual(itk_image.__class__, 
                itk.Image[itk_dtype, array.ndim])
            self.assertEqual(
                list(array.shape),
                list(reversed(itk_image.GetRequestedRegion().GetSize()))
            )
        elif itk_image.GetNameOfClass() == "VectorImage" :
            array_ndim = None
            if isinstance(array, medipy.base.Image) :
                array_ndim = array.ndim
            else :
                array_ndim = array.ndim-1
            
            if isinstance(array, medipy.base.Image) :
                self.assertEqual(array.data_type, "vector")
            self.assertEqual(itk_image.__class__, 
                itk.VectorImage[itk_dtype, array_ndim])
            if isinstance(array, medipy.base.Image) :
                self.assertEqual(
                    list(array.shape),
                    list(reversed(itk_image.GetRequestedRegion().GetSize()))
                )
                self.assertEqual(
                     array.number_of_components, itk_image.GetNumberOfComponentsPerPixel())
            else :
                self.assertEqual(
                    list(array.shape[:-1]),
                    list(reversed(itk_image.GetRequestedRegion().GetSize()))
                )
                self.assertEqual(
                     array.shape[-1], itk_image.GetNumberOfComponentsPerPixel())
        else :
            raise Exception("Unknown ITK image type: {0}".format(itk_image.GetNameOfClass()))
    
    def _test_scalar_content(self, array, itk_image):
        """ array can be any array-like Python object (including medipy.base.Image)
        """
        
        for index in numpy.ndindex(*array.shape) :
            itk_index = list(reversed(index))
            self.assertEqual(array[index], itk_image.GetPixel(itk_index))
    
    def _test_grid(self, medipy_image, itk_image):
        # Check origin, spacing and direction
        numpy.testing.assert_array_almost_equal(
            medipy_image.origin, list(reversed(itk_image.GetOrigin())))
        numpy.testing.assert_array_almost_equal(
            medipy_image.spacing, list(reversed(itk_image.GetSpacing())))
        numpy.testing.assert_array_almost_equal(
            medipy_image.direction, 
            medipy.itk.itk_matrix_to_array(itk_image.GetDirection(), True))
    
    def _test_vector_content(self, array, itk_image):
        """ array can be any array-like Python object (including medipy.base.Image)
        """
        
        if isinstance(array, medipy.base.Image) :
            shape = array.shape
        else :
            shape = array.shape[:-1]
        
        for index in numpy.ndindex(*shape) :
            itk_value = itk_image.GetPixel(list(reversed(index)))
            itk_value = [itk_value.GetElement(i) 
                         for i in range(itk_image.GetNumberOfComponentsPerPixel())]
            
            numpy.testing.assert_array_almost_equal(list(array[index]), itk_value)

if __name__ == "__main__" :
    unittest.main()