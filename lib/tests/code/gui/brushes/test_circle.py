import unittest

import numpy
import numpy.testing

import medipy.base
import medipy.gui.brushes

class TestCircle(unittest.TestCase):
    def test_indices(self):
        radius = 3
        
        image = medipy.base.Image(3*(256,))
        brush = medipy.gui.brushes.Circle(1, image, radius, [1,0,0])
        
        diameter = 2*radius+1
        bounding_box = numpy.transpose(
            numpy.indices((diameter, diameter))
        ).reshape(diameter*diameter, 2)-(radius,radius)
        
        expected_indices = []
        for index in bounding_box :
            if numpy.linalg.norm(index) <= radius :
                expected_indices.append(numpy.hstack([[0], index]))
        expected_indices = numpy.sort(expected_indices)
        
        numpy.testing.assert_almost_equal(numpy.sort(brush.indices()), 
                                          expected_indices)

if __name__ == '__main__':
    unittest.main()
