import unittest

import itk
import numpy

import medipy.geometry

class TestItkPlane(unittest.TestCase) :
    def test_origin_and_normal_constructor(self) :
        origin, normal = (1,0,0), (0,1,0)
        plane = itk.Plane[itk.F](origin, normal)
        
        self.assertEqual(plane.GetOrigin(), origin)
        self.assertEqual(plane.GetNormal(), normal)
        
        self.assertEqual(plane.GetDistance(plane.GetP1()), 0)
        self.assertEqual(plane.GetDistance(plane.GetP2()), 0)
        self.assertEqual(plane.GetDistance(plane.GetP3()), 0)
    
    def test_three_points_constructor(self) :
        p1, p2, p3 = (0,0,0), (1,0,0), (0,1,0)
        plane = itk.Plane[itk.F](p1, p2, p3)
        
        self.assertEqual(plane.GetP1(), p1)
        self.assertEqual(plane.GetP2(), p2)
        self.assertEqual(plane.GetP3(), p3)
        
        v1 = numpy.subtract(p2, p1)
        v2 = numpy.subtract(p3, p1)
        
        self.assertEqual(numpy.dot(plane.GetNormal(), v1), 0)
        self.assertEqual(numpy.dot(plane.GetNormal(), v2), 0)
    
    def test_distance(self) :
        p1, p2, p3 = (1, 0, 0), (0,1,0), (0,0,1)
        plane = itk.Plane[itk.F](p1, p2, p3)
        
        numpy.testing.assert_almost_equal(plane.GetDistance((0,0,0)), -(3**0.5)/3)
    
    def test_reflect(self) :
        p1, p2, p3 = (1, 0, 0), (0,1,0), (0,0,1)
        plane = itk.Plane[itk.F](p1, p2, p3)
        
        numpy.testing.assert_almost_equal(plane.Reflect((0,0,0)), 3*[2./3.])

if __name__ == "__main__" :
    unittest.main()
