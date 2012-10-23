import os
import shutil
import tarfile
import tempfile
import unittest

import numpy

import medipy.io

class TestSerie(unittest.TestCase):
    
    def extract(self, filename):
        data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "data","input"))
        
        tempdir = tempfile.mkdtemp()
        
        archive = tarfile.open(os.path.join(data_directory, filename), "r")
        archive.extractall(tempdir)
        
        return tempdir
    
    def cleanup(self, tempdir):
        shutil.rmtree(tempdir)
    
    def setUp(self):
        self.data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "data"))
        
    def test_load_dicom(self):
        tempdir = self.extract("ep2d_DTI_20_dir_RL_pitch_yaw_22.tar.bz2")
        
        series_instance_uid = "1.3.12.2.1107.5.2.32.35389.2010072309200364977723827.0.0.0"
        url = "dicom:{0}#series_instance_uid={1}".format(
            tempdir, series_instance_uid)
        images = medipy.io.load_serie(url, numpy.uint32)
        
        self.cleanup(tempdir)
        
        self._test_dti_serie(images)
    
    def test_load_nifti(self):
        filename = os.path.join(self.data_directory, "input",
                                "ep2d_DTI_20_dir_RL_pitch_yaw_22.nii.gz")
        images = medipy.io.load_serie(filename)
        self._test_dti_serie(images)
    
    #####################
    # Private interface #
    #####################    
    
    def _test_dti_serie(self, images):
        self.assertEqual(len(images), 42)
        
        b_0_images = []
        b_non_0_images = []
        for image in images :
            self.assertTrue("mr_diffusion_sequence" in image.metadata)
            diffusion = image.metadata["mr_diffusion_sequence"][0]
            if diffusion.diffusion_bvalue == 0 :
                b_0_images.append(image)
            else :
                b_non_0_images.append(image)
        
        self.assertEqual(len(b_0_images), 2)
        self.assertEqual(len(b_non_0_images), 40)
        
        for image in b_0_images :
            diffusion = image.metadata["mr_diffusion_sequence"][0]
            gradient = diffusion.diffusion_gradient_direction_sequence[0].\
                diffusion_gradient_orientation
            self.assertAlmostEqual(numpy.linalg.norm(gradient), 0, 3)
            
            self.assertEqual(image.shape, (60, 102, 102))
        
        for image in b_non_0_images :
            diffusion = image.metadata["mr_diffusion_sequence"][0]
            self.assertEqual(diffusion.diffusion_bvalue, 1000)
            
            gradient = diffusion.diffusion_gradient_direction_sequence[0].\
                diffusion_gradient_orientation
            self.assertAlmostEqual(numpy.linalg.norm(gradient), 1, 3)
            
            self.assertEqual(image.shape, (60, 102, 102))
        

if __name__ == "__main__" :
    unittest.main()