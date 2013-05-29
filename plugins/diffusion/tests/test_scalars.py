import unittest
import numpy

import medipy.base
import medipy.diffusion

class TestScalars(unittest.TestCase):
    
    def setUp(self) :
        tensor = [ 6.06370799e-04, -4.86828321e-05, -5.83471810e-05,
                   4.17223084e-04,  1.12409536e-04,  7.96857639e-04]
        self.image = medipy.base.Image(shape=(1,1,1), dti="tensor_2")
        self.image[0,0,0] = tensor
    
    def test_fractional_anisotropy(self) :
        fa = medipy.diffusion.scalars.fractional_anisotropy(self.image)
        self.assertEqual(fa.shape, self.image.shape)
        self.assertEqual(fa.dtype, self.image.dtype)
        numpy.testing.assert_almost_equal(fa[0,0,0], 0.36685354)
    
    def test_mean_diffusivity(self) :
        md = medipy.diffusion.scalars.mean_diffusivity(self.image)
        self.assertEqual(md.shape, self.image.shape)
        self.assertEqual(md.dtype, self.image.dtype)
        numpy.testing.assert_almost_equal(md[0,0,0], 0.00060681719)

if __name__ == "__main__" :
    unittest.main()
