import unittest

import itk
import numpy

import medipy.itk

class TestNumpyBridge(unittest.TestCase):
    
    def test_array_to_itk_matrix_no_flip(self):
        array = (1024*numpy.random.random((3,3))).astype(numpy.single)
        itk_matrix = medipy.itk.array_to_itk_matrix(array, False)
        for index in numpy.ndindex(*array.shape) :
            self.assertEqual(array[index], itk_matrix(*index))
    
    def test_array_to_itk_matrix_flip(self):
        array = (1024*numpy.random.random((3,3))).astype(numpy.single)
        itk_matrix = medipy.itk.array_to_itk_matrix(array, True)
        for index in numpy.ndindex(*array.shape) :
            self.assertEqual(numpy.fliplr(numpy.flipud(array))[index], 
                             itk_matrix(*index))
    
    def test_itk_matrix_to_array_no_flip(self):
        source_array = (1024*numpy.random.random((3,3))).astype(numpy.single)
        itk_matrix = medipy.itk.array_to_itk_matrix(source_array, False)
                
        array = medipy.itk.itk_matrix_to_array(itk_matrix, False)
        for index in numpy.ndindex(*array.shape) :
            self.assertEqual(array[index], itk_matrix(*index))
    
    def test_itk_matrix_to_array_flip(self):
        source_array = (1024*numpy.random.random((3,3))).astype(numpy.single)
        itk_matrix = medipy.itk.array_to_itk_matrix(source_array, False)
                
        array = medipy.itk.itk_matrix_to_array(itk_matrix, True)
        
        for index in numpy.ndindex(*array.shape) :
            self.assertEqual(numpy.fliplr(numpy.flipud(array))[index], 
                             itk_matrix(*index))

if __name__ == "__main__" :
    unittest.main()
