import numpy
import os
import unittest

import medipy.io.dicom
import medipy.io.dicom.normalize

class TestDiffusion(unittest.TestCase):
    
    def setUp(self):
        self.data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data"))
    
    def test_siemens_0(self):
        """ Normalization of a DWI Siemens image with a b-value of 0.
        """
        
        dataset = medipy.io.dicom.read(
            os.path.join(self.data_directory, "input", "siemens_dwi_0.dcm"))
        
        normalized = medipy.io.dicom.normalize.normalize(dataset)
        
        for dataset in normalized :
            self.assertTrue("mr_diffusion_sequence" in dataset)
            diffusion = dataset.mr_diffusion_sequence.value[0]
            self.assertEqual(diffusion.diffusion_bvalue.value, 0)
            self.assertEqual(diffusion.diffusion_directionality.value, "DIRECTIONAL")
            self.assertTrue("diffusion_gradient_direction_sequence" in diffusion)
            direction = diffusion.diffusion_gradient_direction_sequence.value[0].\
                diffusion_gradient_orientation
            numpy.testing.assert_array_equal(direction.value, [0,0,0])
    
    def test_siemens_1000(self):
        """ Normalization of a DWI Siemens image with a b-value of 1000.
        """
        
        dataset = medipy.io.dicom.read(
            os.path.join(self.data_directory, "input", "siemens_dwi_1000.dcm"))
        
        normalized = medipy.io.dicom.normalize.normalize(dataset)
        
        for dataset in normalized :
            self.assertTrue("mr_diffusion_sequence" in dataset)
            diffusion = dataset.mr_diffusion_sequence.value[0]
            self.assertEqual(diffusion.diffusion_bvalue.value, 1000)
            self.assertEqual(diffusion.diffusion_directionality.value, "DIRECTIONAL")
            self.assertTrue("diffusion_gradient_direction_sequence" in diffusion)
            direction = diffusion.diffusion_gradient_direction_sequence.value[0].\
                diffusion_gradient_orientation
            self.assertAlmostEqual(numpy.linalg.norm(direction.value), 1)

if __name__ == "__main__" :
    unittest.main()
