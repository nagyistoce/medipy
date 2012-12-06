import numpy
import os
import shutil
import tarfile
import tempfile
import unittest

import medipy.base
import medipy.io
import medipy.io.dicom
import medipy.io.dicom.normalize

class TestReconstruction(unittest.TestCase):
    
    def extract(self, filename):
        data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data","input"))
        
        tempdir = tempfile.mkdtemp()
        
        archive = tarfile.open(os.path.join(data_directory, filename), "r")
        archive.extractall(tempdir)
        
        return tempdir
    
    def cleanup(self, tempdir):
        shutil.rmtree(tempdir)
    
    def setUp(self):
        self.data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data","input"))
    
    def test_brainix(self):
        tempdir = self.extract("brainix.tgz")
        
        series_instance_uid = "1.3.46.670589.11.0.0.11.4.2.0.8743.5.5396.2006120114285654497"
        url = "dicomdir:{0}#series_instance_uid={1}".format(
            os.path.join(tempdir, "BRAINIX", "DICOMDIR"), series_instance_uid)
        image = medipy.io.load(url)
        
        self.cleanup(tempdir)
        
        self.assertEqual(image.shape, (22, 288, 288))
        numpy.testing.assert_array_almost_equal(image.spacing, 
                                                (6., 0.7986111, 0.7986111))
    
    def test_nema_mf_mr_mammo(self):
        tempdir = self.extract("nemamfmr.imagesDG.tar.bz2")
        dataset = medipy.io.dicom.read(os.path.join(tempdir, "DISCIMG/IMAGES/DYNMAMMO"))
        self.cleanup(tempdir)
        
        self.assertEqual(dataset.number_of_frames, 168)
        
        normalized = medipy.io.dicom.normalize.normalize(dataset)
        self.assertEqual(len(normalized), 168)
        self.assertRaises(medipy.base.Exception, medipy.io.dicom.image, normalized)
        
        stacks = medipy.io.dicom.stacks(normalized)
        self.assertEqual(len(stacks), 6)
        
        for stack in stacks :
            image = medipy.io.dicom.image(stack)
            self.assertEqual(image.shape, (28, 512, 512))
            numpy.testing.assert_almost_equal(image.spacing,
                                              (2.5, 0.703125, 0.703125))
    
    def test_nm(self):
        tempdir = self.extract("nm.tar.gz")
        image = medipy.io.load(os.path.join(tempdir, "nm.dcm"))
        self.cleanup(tempdir)
        
        self.assertEqual(image.shape, (43, 128, 128))
        numpy.testing.assert_array_almost_equal(image.spacing, 
                                                (3.935541, 3.895369, 3.895369))
        
    
if __name__ == "__main__" :
    unittest.main()