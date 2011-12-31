import unittest

import medipy.fsl

class TestSienax(unittest.TestCase):
    
    def setUp(self):
        self.sienax = medipy.fsl.Sienax()
        self.sienax.input = "foo.nii"
    
    def test_basic(self):
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii"])
    
    def test_output_directory(self):
        self.sienax.output_directory = "bar"
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-o", "bar"])
    
    def test_two_class_segmentation(self):
        self.sienax.two_class_segmentation = True
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-2"])
    
    def test_t2_weighted_input(self):
        self.sienax.t2_weighted_input = True
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-t2"])
    
    def test_ignore_upwards(self):
        self.sienax.ignore_upwards = 10
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-t", "10"])
    
    def test_ignore_downwards(self):
        self.sienax.ignore_downwards = 10
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-b", "10"])
    
    def test_regional(self):
        self.sienax.regional = True
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-r"])
    
    def test_lesion_mask(self):
        self.sienax.lesion_mask = "bar.nii"
        self.assertEqual(self.sienax.command, ["sienax", "foo.nii", "-lm", "bar.nii"])
    
if __name__ == "__main__" :
    unittest.main()