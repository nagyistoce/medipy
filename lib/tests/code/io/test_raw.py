import gzip
import os
import shutil
import tempfile
import unittest

import numpy

import medipy.base
import medipy.io

class TestSchemes(unittest.TestCase):
    def setUp(self):
        self.data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "data","input"))
        self.filename = "ep2d_DTI_20_dir_RL_pitch_yaw_22.nii"
        
        self.shape = (42, 60, 102, 102)
        self.dtype = numpy.int16
        self.offset = 352
    
    def extract(self):
        destination = tempfile.mkdtemp()
        
        source = gzip.open(
            os.path.join(self.data_directory, self.filename+".gz"), "rb")
        target = open(os.path.join(destination, self.filename), "wb")
        
        target.write(source.read())
        
        return destination
    
    def clean(self, location):
        shutil.rmtree(location)
    
    def test_can_load(self):
        directory = self.extract()
        
        can_load = medipy.io.raw.can_load(
            os.path.join(directory, self.filename), 
            self.shape, self.dtype, offset=self.offset)
        
        self.assertTrue(can_load)
        
        self.clean(directory)
    
    def test_load(self):
        directory = self.extract()
        
        image = medipy.io.raw.load(
            os.path.join(directory, self.filename), 
            self.shape, self.dtype, offset=self.offset)
        
        self.assertTrue(isinstance(image, medipy.base.Image))
        self.assertEqual(image.shape, self.shape)
        self.assertEqual(image.dtype, self.dtype)
        
        numpy.testing.assert_array_equal(image.origin, image.ndim*(0,))
        numpy.testing.assert_array_equal(image.spacing, image.ndim*(1,))
        numpy.testing.assert_array_equal(image.direction, numpy.identity(image.ndim))
        
        self.assertTrue("loader" in image.metadata)
        self.assertTrue("url" in image.metadata["loader"])
        self.assertEqual(image.metadata["loader"]["url"], 
            os.path.join(directory, self.filename))
        
        self.clean(directory)

if __name__ == "__main__" :
    unittest.main()
