import unittest
import numpy

import medipy.base
import medipy.diffusion

class TestUtils(unittest.TestCase):
    
    def setUp(self) :
        tensor = [ 6.06370799e-04, -4.86828321e-05, -5.83471810e-05,
                   4.17223084e-04,  1.12409536e-04,  7.96857639e-04]
        self.image = medipy.base.Image(shape=(1,1,1), dti="tensor_2")
        self.image[0,0,0] = tensor
    
    def test_log_transformation(self) :
        result = medipy.diffusion.utils.log_transformation(self.image)
    
        self.assertEqual(result.shape, self.image.shape)
        self.assertEqual(result.image_type, "tensor_2")
        self.assertEqual(result.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(result[0,0,0], 
            [-7.41512775, -0.08842814, -0.07699567, 
             -7.80961609, 0.19082165, -7.15261126])
    
    def test_exp_transformation(self) :
        result = medipy.diffusion.utils.exp_transformation(self.image)
    
        self.assertEqual(result.shape, self.image.shape)
        self.assertEqual(result.image_type, "tensor_2")
        self.assertEqual(result.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(result[0,0,0], 
            [ 1.00060642e+00, -4.86820936e-05, -5.84162772e-05,
              1.00041723e+00,  1.12473965e-04,  1.00079715e+00])
    
    def test_spectral_decomposition(self) :
    
        eigenvalues, eigenvectors = medipy.diffusion.utils.spectral_decomposition(
            self.image)
        
        self.assertEqual(eigenvalues.shape, self.image.shape)
        self.assertEqual(len(eigenvalues[0,0,0]), 3)
        self.assertEqual(eigenvalues.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(eigenvalues[0,0,0], 
            [0.00038178, 0.00059103, 0.00084763])
        
        self.assertEqual(eigenvectors.shape, self.image.shape)
        self.assertEqual(len(eigenvectors[0,0,0]), 9)
        self.assertEqual(eigenvectors.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(eigenvectors[0,0,0], 
            [-0.1458625 , -0.95988911,  0.23945151, -0.94952697, 0.06789147,
             -0.3062503 , -0.27770963,  0.27203611,  0.92134345])
    
    def test_dti6to33(self) :
        result = medipy.diffusion.utils.dti6to33(self.image.data)
        
        self.assertEqual(result.shape, self.image.shape+(3,3))
        self.assertEqual(result.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(result[0,0,0], 
            [[  6.06370799e-04, -4.86828321e-05, -5.83471810e-05],
             [ -4.86828321e-05,  4.17223084e-04,  1.12409536e-04],
             [ -5.83471810e-05,  1.12409536e-04,  7.96857639e-04]])
    
    def test_dti33to6(self) :
        result = medipy.diffusion.utils.dti33to6(
            medipy.diffusion.utils.dti6to33(self.image.data))
        
        self.assertEqual(result.shape, self.image.shape+(6,))
        self.assertEqual(result.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(result[0,0,0], self.image[0,0,0])
    
    # TODO : rotation
    
    def test_compose_spectral(self) :
        eigenvalues, eigenvectors = medipy.diffusion.utils.spectral_decomposition(
            self.image)
        
        # FIXME : function should be modified to match spectral_decomposition
        result = medipy.diffusion.utils.compose_spectral(
            eigenvectors.data.reshape(self.image.shape+(3,3)), eigenvalues)
        
        self.assertEqual(result.shape, self.image.shape+(6,))
        self.assertEqual(result.dtype, self.image.dtype)
        
        numpy.testing.assert_almost_equal(result[0,0,0], self.image[0,0,0])

if __name__ == "__main__" :
    unittest.main()
