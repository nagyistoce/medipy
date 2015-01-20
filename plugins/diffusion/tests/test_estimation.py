import unittest
import numpy

import medipy.base
import medipy.io.dicom

import medipy.diffusion

class TestEstimation(unittest.TestCase):
    
    def test_least_squares(self) :
        signal = [1035.0, 555.0, 558.0, 597.0, 661.0, 503.0, 503.0, 469.0, 
                  690.0, 522.0, 648.0, 580.0, 430.0, 697.0, 557.0, 560.0, 721.0, 
                  624.0, 534.0, 654.0, 502.0, 703.0, 637.0, 491.0, 478.0, 637.0,
                  557.0, 521.0, 433.0, 469.0, 666.0, 554.0, 482.0, 585.0]

        b_values = [
            0, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 
            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
        directions = [
            [0.0, 0.0, 0.0], [0.875, -0.48, 0.062], [-0.281, -0.814, -0.508], 
            [0.83, -0.325, 0.454], [-0.601, 0.678, -0.423], 
            [0.683, 0.374, 0.627], [-0.614, -0.362, 0.701], 
            [-0.424, -0.094, -0.901], [0.286, -0.905, 0.315], 
            [0.053, -0.218, -0.975], [0.262, -0.607, 0.75], [1.0, 0.0, 0.0], 
            [-0.075, 0.638, 0.767], [0.094, -0.995, -0.035], 
            [0.904, 0.401, 0.148], [0.854, 0.401, -0.331], [0.629, 0.71, 0.317],
            [-0.887, -0.066, -0.456], [0.923, -0.022, -0.384], 
            [-0.48, -0.733, 0.481], [-0.823, 0.45, 0.348], 
            [-0.677, -0.728, 0.111], [-0.198, -0.552, 0.81], 
            [-0.094, 0.22, -0.971], [-0.565, 0.271, -0.779], 
            [0.313, 0.95, 0.0], [0.519, -0.687, -0.509], 
            [0.144, -0.895, -0.422], [-0.348, -0.102, 0.932], 
            [0.443, -0.393, -0.806], [0.09, 0.873, -0.479], 
            [0.314, 0.493, 0.811], [0.699, -0.113, -0.706], 
            [-0.566, 0.822, 0.063]]
        
        self.assertTrue(len(signal) == len(b_values) == len(directions))
        
        images = []
        for index in xrange(len(signal)) :
            image = medipy.base.Image(shape = (1,1,1), dtype=numpy.single)
            image[0,0,0] = signal[index]
            
            diffusion_dataset = medipy.io.dicom.DataSet(
                diffusion_bvalue=b_values[index],
                diffusion_gradient_direction_sequence = [medipy.io.dicom.DataSet(
                    diffusion_gradient_orientation = directions[index])])
            
            image.metadata["mr_diffusion_sequence"] = [diffusion_dataset]
            
            images.append(image)
        
        tensors = medipy.diffusion.estimation.least_squares(images)
        
        self.assertEqual(tensors.shape, images[0].shape)
        self.assertEqual(tensors.dtype, images[0].dtype)
        numpy.testing.assert_almost_equal(tensors[0,0,0],
            [6.03082415e-04, -4.93411935e-05, -5.63119102e-05,
             4.15756309e-04,  1.05372004e-04,  7.93865882e-04])

if __name__ == "__main__" :
    unittest.main()
