import unittest

import numpy

import medipy.base

class TestCoordinateSystem(unittest.TestCase):
    
    def test_best_fitting_axes_aligned_matrix(self):
        matrix = numpy.asarray([[1.1,  0.2, -0.1],
                                [0.1, -0.2, -0.99],
                                [0.0,  1.2, -0.05]])
        axis_aligned = medipy.base.coordinate_system.best_fitting_axes_aligned_matrix(matrix)
        numpy.testing.assert_array_equal(axis_aligned,
                                         numpy.asarray([[1, 0, 0],
                                                        [0, 0, -1],
                                                        [0, 1, 0]]))

if __name__ == "__main__" :
    unittest.main()
