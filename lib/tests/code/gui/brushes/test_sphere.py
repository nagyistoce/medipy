import unittest

import numpy
import numpy.testing

import medipy.base
import medipy.gui.brushes

class TestSphere(unittest.TestCase):
    def test_indices(self):
        radius = 3
        
        image = medipy.base.Image(3*(256,))
        brush = medipy.gui.brushes.Sphere(1, image, radius)
        
        diameter = 2*radius+1
        bounding_box = numpy.transpose(
            numpy.indices((diameter, diameter, diameter))
        ).reshape(diameter*diameter*diameter, 3)-(radius, radius, radius)
        
        expected_indices = []
        for index in bounding_box :
            if numpy.linalg.norm(index) <= radius :
                expected_indices.append(index)
        expected_indices = numpy.sort(expected_indices)
        
        numpy.testing.assert_almost_equal(numpy.sort(brush.indices()), 
                                          expected_indices)

if __name__ == '__main__':
    unittest.main()
