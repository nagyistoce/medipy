import os
import shutil
import tarfile
import tempfile
import unittest

import numpy.testing

import medipy.base
import medipy.io

class TestSchemes(unittest.TestCase):
    
    def setUp(self):
        self.data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "data","input"))
    
    def extract_dicom(self):
        destination = tempfile.mkdtemp()
        archive = tarfile.open(
            os.path.join(self.data_directory, "brainix.tgz"), "r")
        archive.extractall(destination)
        return destination
    
    def clean_dicom(self, location):
        shutil.rmtree(location)
    
    def test_default(self) :
        path = os.path.join(self.data_directory, "avg152T1_LR_nifti.nii.gz")
        image = medipy.io.load(path)
        self.assertTrue(isinstance(image, medipy.base.Image))
        self.assertEqual(image.metadata["loader"]["url"], path)
    
    def test_file(self) :
        url = "file:{0}".format(
            os.path.join(self.data_directory, "avg152T1_LR_nifti.nii.gz"))
        image = medipy.io.load(url)
        self.assertTrue(isinstance(image, medipy.base.Image))
        self.assertEqual(image.metadata["loader"]["url"], url) 
    
    def test_dicomdir(self):
        location = self.extract_dicom()
        
        series_instance_uid = "1.3.46.670589.11.0.0.11.4.2.0.8743.5.5396.2006120114285654497"
        
        url = "dicomdir:{0}#series_instance_uid={1}".format(
            os.path.join(location, "BRAINIX", "DICOMDIR"), series_instance_uid)
        image1 = medipy.io.load(url)
        self.assertTrue(isinstance(image1, medipy.base.Image))
        self.assertEqual(image1.metadata["loader"]["url"], url)
        
        tags = ["0x0020000e", "0020000e",
                "(0020,000e)", "(0x0020,0x000e)"]
        for tag in tags :
            url = "dicomdir:{0}#{1}={2}".format(
                os.path.join(location, "BRAINIX", "DICOMDIR"), 
                tag, series_instance_uid)
            image2 = medipy.io.load(url)
            numpy.testing.assert_array_equal(image1, image2)
        
        self.clean_dicom(location)
    
    def test_dicom(self):
        location = self.extract_dicom()
        
        url = "dicom:{0}".format(
            os.path.join(location, "BRAINIX", "2182114", "401"))
        image1 = medipy.io.load(url)
        self.assertTrue(isinstance(image1, medipy.base.Image))
        self.assertEqual(image1.metadata["loader"]["url"], url)
        
        tags = ["series_instance_uid", 
                "0x0020000e", "0020000e",
                "(0020,000e)", "(0x0020,0x000e)"]
        series_instance_uid = "1.3.46.670589.11.0.0.11.4.2.0.8743.5.5396.2006120114285654497"
        for tag in tags :
            url = "dicom:{0}#{1}={2}".format(
                os.path.join(location, "BRAINIX", "2182114", "401"), 
                tag, series_instance_uid)
            image2 = medipy.io.load(url)
            numpy.testing.assert_array_equal(image1, image2)
        
        self.clean_dicom(location)
        
if __name__ == '__main__':
    unittest.main()
