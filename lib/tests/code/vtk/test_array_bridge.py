import unittest

import numpy
from vtk import vtkImageData
import vtk.util.numpy_support

import medipy.vtk

class TestArrayBridge(unittest.TestCase) :
    
    def array_to_vtk_image_scalar(self, copy_data) :
        array = numpy.ndarray((91, 109, 60), dtype=numpy.uint16)
        array[10,20,30] = 42
        
        vtk_image = medipy.vtk.array_to_vtk_image(array, copy_data)
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            array.shape, [x for x in reversed(vtk_shape)])
        
        self.assertEqual(
            array.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(vtk_image.GetScalarComponentAsFloat(30,20,10,0), 42)
    
    def test_array_to_vtk_image_scalar_with_copy(self) :
        self.array_to_vtk_image_scalar(True)
    
    def test_array_to_vtk_image_scalar_without_copy(self) :
        self.array_to_vtk_image_scalar(False)
    
    def array_to_vtk_image_complex(self, copy_data) :
        array = numpy.ndarray((91, 109, 60), dtype=numpy.complex64)
        array[10,20,30] = 5+42j
        
        vtk_image = medipy.vtk.array_to_vtk_image(array, copy_data)
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            array.shape, [x for x in reversed(vtk_shape)])
        numpy.testing.assert_array_equal(
            vtk_image.GetNumberOfScalarComponents(), 2)
        
        self.assertEqual(vtk_image.GetScalarComponentAsFloat(30,20,10,0), 5)
        self.assertEqual(vtk_image.GetScalarComponentAsFloat(30,20,10,1), 42)
    
    def test_array_to_vtk_image_complex_with_copy(self) :
        self.array_to_vtk_image_complex(True)
    
    def test_array_to_vtk_image_complex_without_copy(self) :
        self.array_to_vtk_image_complex(False)
    
    def array_to_vtk_image_vector(self, copy_data) :
        array = numpy.ndarray((91, 109, 60, 6), dtype=numpy.uint16)
        array[10,20,30,3] = 42
        
        vtk_image = medipy.vtk.array_to_vtk_image(array, copy_data, "vector")
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            array.shape[:-1], [x for x in reversed(vtk_shape)])
        numpy.testing.assert_array_equal(
            array.shape[-1], vtk_image.GetNumberOfScalarComponents())
        
        self.assertEqual(
            array.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(vtk_image.GetScalarComponentAsFloat(30,20,10,3), 42)
    
    def test_array_to_vtk_image_vector_with_copy(self) :
        self.array_to_vtk_image_vector(True)
    
    def test_array_to_vtk_image_vector_without_copy(self) :
        self.array_to_vtk_image_vector(False)
     
    def array_to_vtk_image_complex_vector(self, copy_data) :
        array = numpy.ndarray((91, 109, 60, 6), dtype=numpy.complex64)
        array[10,20,30,3] = 5+42j
        
        vtk_image = medipy.vtk.array_to_vtk_image(array, copy_data, "vector")
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            array.shape[:-1], [x for x in reversed(vtk_shape)])
        numpy.testing.assert_array_equal(
            2*array.shape[-1], vtk_image.GetNumberOfScalarComponents())
        
        self.assertEqual(
            array[0,0,0,0].real.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(vtk_image.GetScalarComponentAsFloat(30,20,10,6), 5)
        self.assertEqual(vtk_image.GetScalarComponentAsFloat(30,20,10,7), 42)
        
    def test_array_to_vtk_image_complex_vector_with_copy(self) :
        self.array_to_vtk_image_complex_vector(True)
    
    def test_array_to_vtk_image_complex_vector_without_copy(self) :
        self.array_to_vtk_image_complex_vector(False)

    def get_vtk_image(self, number_of_scalar_components) :
        vtk_image = vtkImageData()
        vtk_image.SetScalarTypeToUnsignedShort()
        vtk_image.SetNumberOfScalarComponents(number_of_scalar_components)
        
        extent = (0, 59, 0, 108, 0, 90)
        vtk_image.SetWholeExtent(extent)
        vtk_image.SetUpdateExtentToWholeExtent()
        vtk_image.SetExtent(extent)
        
        vtk_image.AllocateScalars()
        
        vtk_image.SetOrigin(1,2,3)
        vtk_image.SetSpacing(4,5,6)
        
        return vtk_image
        
    def test_vtk_image_to_array_scalar(self) :
        vtk_image = self.get_vtk_image(1)
        vtk_image.SetScalarComponentFromFloat(30,20,10,0, 42)
        
        array = medipy.vtk.vtk_image_to_array(vtk_image)
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            array.shape, [x for x in reversed(vtk_shape)])
        
        self.assertEqual(
            array.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(array[10,20,30], 42)
    
    def test_vtk_image_to_array_vector(self) :
        vtk_image = self.get_vtk_image(6)
        vtk_image.SetScalarComponentFromFloat(30,20,10,3, 42)
        
        array = medipy.vtk.vtk_image_to_array(vtk_image)
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            array.shape[:-1], [x for x in reversed(vtk_shape)])
        numpy.testing.assert_array_equal(
            array.shape[-1], vtk_image.GetNumberOfScalarComponents())
        
        self.assertEqual(
            array.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(array[10,20,30][3], 42)
    
    def test_vtk_image_to_medipy_image_scalar(self) :
        vtk_image = self.get_vtk_image(1)
        vtk_image.SetScalarComponentFromFloat(30,20,10,0, 42)
        
        image = medipy.vtk.vtk_image_to_medipy_image(vtk_image, None)
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            image.shape, [x for x in reversed(vtk_shape)])
        
        self.assertEqual(
            image.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(image[10,20,30], 42)
        
        self.assertEqual(image.data_type, "scalar")
        
        numpy.testing.assert_array_equal(
            image.origin, 
            [x for x in reversed(vtk_image.GetOrigin())]
        )
        
        numpy.testing.assert_array_equal(
            image.spacing, 
            [x for x in reversed(vtk_image.GetSpacing())]
        )
    
    def test_vtk_image_to_medipy_image_vector(self) :
        vtk_image = self.get_vtk_image(6)
        vtk_image.SetScalarComponentFromFloat(30,20,10,3, 42)
        
        image = medipy.vtk.vtk_image_to_medipy_image(vtk_image, None)
        
        vtk_shape = [1+x for x in vtk_image.GetExtent()[1::2]]
        numpy.testing.assert_array_equal(
            image.shape, [x for x in reversed(vtk_shape)])
        
        self.assertEqual(
            image.dtype, 
            vtk.util.numpy_support.get_numpy_array_type(vtk_image.GetScalarType()))
        
        self.assertEqual(image[10,20,30][3], 42)
        
        self.assertEqual(image.data_type, "vector")
        
        numpy.testing.assert_array_equal(
            image.origin, 
            [x for x in reversed(vtk_image.GetOrigin())]
        )
        
        numpy.testing.assert_array_equal(
            image.spacing, 
            [x for x in reversed(vtk_image.GetSpacing())]
        )

if __name__ == "__main__" :
    unittest.main()
