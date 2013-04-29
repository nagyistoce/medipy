#import os
#import shutil
#import tempfile
import unittest

import numpy

#import medipy.base
#import medipy.io
import medipy.diffusion

class TestIO(unittest.TestCase):
    
    def test_dti6to33(self):
        upper_diagonal = numpy.ndarray([1,1,1,6], dtype=numpy.single)
        upper_diagonal[0,0,0] = (1,2,3,4,5,6)
        
        matrix = medipy.diffusion.utils.dti6to33(upper_diagonal)
        
        self.assertEqual(upper_diagonal.shape[:-1], matrix.shape[:-2])
        self.assertEqual(matrix.shape[-2:], (3,3))
        self.assertEqual(upper_diagonal.dtype, matrix.dtype)
        numpy.testing.assert_array_equal(matrix[0,0,0], [[1,2,3],
                                                         [2,4,5],
                                                         [3,5,6]])
    
    def test_dti33to66(self):
        matrix = numpy.ndarray([1,1,1,3,3], dtype=numpy.single)
        matrix[0,0,0] = [[1,2,3],
                         [2,4,5],
                         [3,5,6]]
        
        upper_diagonal = medipy.diffusion.utils.dti33to6(matrix)
        
        self.assertEqual(matrix.shape[:-2], upper_diagonal.shape[:-1])
        self.assertEqual(upper_diagonal.shape[-1:], (6,))
        self.assertEqual(matrix.dtype, upper_diagonal.dtype)
        numpy.testing.assert_array_equal(upper_diagonal[0,0,0], [1,2,3,4,5,6])

if __name__ == "__main__" :
    unittest.main()
