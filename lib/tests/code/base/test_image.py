import unittest

import numpy

import medipy.base

class TestImage(unittest.TestCase):
    
    def setUp(self):
        self.un_initialized_image = medipy.base.Image((4,5,6), numpy.uint16)
        self.constant_image = medipy.base.Image((4,5,6), numpy.uint16, value=1)
        self.initialized_image = medipy.base.Image(
            data=numpy.arange(0, 4*5*6, dtype=numpy.uint16).reshape(4,5,6))
    
    def test_uninitialized_constructor(self):
        self.assertEqual(type(self.un_initialized_image.data), numpy.ndarray)
        self.assertEqual(self.un_initialized_image.data.shape, (4,5,6))
        self.assertEqual(self.un_initialized_image.data.dtype, numpy.uint16)
    
    def test_constant_constructor(self):
        self.assertEqual(type(self.constant_image.data), numpy.ndarray)
        self.assertEqual(self.constant_image.data.shape, (4,5,6))
        self.assertEqual(self.constant_image.data.dtype, numpy.uint16)
        numpy.testing.assert_equal(
           self.constant_image.data, 
           numpy.ones(self.constant_image.data.shape, self.constant_image.data.dtype))
    
    def test_initialized_constructor(self):
        self.assertEqual(type(self.initialized_image.data), numpy.ndarray)
        self.assertEqual(self.initialized_image.data.shape, (4,5,6))
        self.assertEqual(self.initialized_image.data.dtype, numpy.uint16)
        numpy.testing.assert_equal(
           self.initialized_image.data, 
           numpy.arange(0, 4*5*6, dtype=numpy.uint16).reshape(4,5,6))
    
    def test_as_array(self):
        array = numpy.asarray(self.initialized_image)
        self.assertEqual(type(array), numpy.ndarray)
        self.assertEqual(array.shape, (4,5,6))
        self.assertEqual(array.dtype, numpy.uint16)
        numpy.testing.assert_equal(
           array, numpy.arange(0, 4*5*6, dtype=numpy.uint16).reshape(4,5,6))
    
    def test_element_read(self):
        self.assertEqual(self.initialized_image[3,2,4], 4+6*(2+5*3))
    
    def test_element_write(self):
        self.initialized_image[3,2,4] = 0
        self.assertEqual(self.initialized_image[3,2,4], 0)
    
    def test_slice_read(self):
        numpy.testing.assert_equal(
           self.initialized_image[2:4,3:5,4:6], 
           [
            [[82, 83], 
             [88, 89]],
            [[112, 113], 
             [118, 119]]])
    
    def test_slice_write(self):
        self.initialized_image[2:4,3:5,4:6] = numpy.arange(
           0, 8, dtype=numpy.uint16).reshape(2,2,2)
        numpy.testing.assert_equal(
           self.initialized_image[2:4,3:5,4:6], 
           [
            [[0, 1], 
             [2, 3]],
            [[4, 5], 
             [6, 7]]])
    
    def test_inside(self):
        self.assertTrue(self.un_initialized_image.is_inside((1,2,3)))
        self.assertFalse(self.un_initialized_image.is_inside((4,2,3)))
    
    def test_index_to_physical(self):
        image = medipy.base.Image((64,64,64), numpy.uint8, 
                                  spacing=(1,2,3), origin=(4,5,6))
        index = (8,9,10)
        physical = image.index_to_physical(index)
        numpy.testing.assert_array_almost_equal(physical, (12, 23, 36))
    
    def test_physical_to_index(self):
        image = medipy.base.Image((64,64,64), numpy.uint8, 
                                  spacing=(1,2,3), origin=(4,5,6))
        physical = (12.,23.,36.)
        index = image.physical_to_index(physical)
        numpy.testing.assert_array_almost_equal(index, (8, 9, 10))
    
if __name__ == '__main__':
    unittest.main()