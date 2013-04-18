import unittest

import numpy

import medipy.base
import medipy.intensity

class TestReplace(unittest.TestCase):
    def test_replace_labels(self):
        image = medipy.base.Image((2,2), dtype=numpy.int8, value=0)
        image[0,0] = 1
        image[1,1] = 42
        modified_indices = medipy.intensity.replace_labels(image, {1: 2, 42: 3}, True)
        
        # Transform the modified indices from numpy-friendly to user-friendly
        modified_indices = dict([[key, numpy.transpose(value).tolist()] 
                                  for key, value in modified_indices.items()])
        
        self.assertEqual(modified_indices, {1: [[0, 0]], 42: [[1, 1]]})
        self.assertEqual(image[0,0], 2)
        self.assertEqual(image[1,1], 3)
        self.assertEqual(image[0,1], 0)
        self.assertEqual(image[1,0], 0)
    
if __name__ == "__main__" :
    unittest.main()
