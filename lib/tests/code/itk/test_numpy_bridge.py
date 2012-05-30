import unittest

import itk
import numpy

import medipy.itk

class TestNumpyBridge(unittest.TestCase):
    
    def setUp(self):
        shape = [4,5,6]
        dtype = numpy.uint16
        
        PixelType = medipy.itk.dtype_to_itk[dtype]
        Dimension = len(shape)
        ImageType = itk.Image[PixelType, Dimension]

        self.numpy_array = numpy.ndarray(shape, dtype)        
        
        self.itk_image = ImageType.New()
        region = itk.ImageRegion[Dimension]((0,0,0), list(reversed(shape)))
        self.itk_image.SetRegions(region)
        self.itk_image.Allocate()
        
        value = 0
        for index in numpy.ndindex(*shape) :
            self.numpy_array[index] = value
            self.itk_image.SetPixel(list(reversed(index)), value)
            value += 1
    
    def test_numpy_to_itk(self):
        result = medipy.itk.array_to_itk_image(self.numpy_array, False)
        for index in numpy.ndindex(*self.numpy_array.shape) :
            self.assertEqual(self.numpy_array[index],
                             result.GetPixel(list(reversed(index))))
    
    def test_itk_to_numpy(self):
        result = medipy.itk.itk_image_to_array(self.itk_image, False)
        for index in numpy.ndindex(*self.itk_image.GetBufferedRegion().GetSize()) :
            self.assertEqual(self.itk_image.GetPixel(index),
                             result[tuple(reversed(index))])

if __name__ == "__main__" :
    unittest.main()