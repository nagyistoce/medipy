import os
import unittest

import medipy.fsl

class TestSiena(unittest.TestCase):
    
    def setUp(self):
        self.siena = medipy.fsl.Siena()
        self.siena.input1 = "foo1.nii"
        self.siena.input2 = "foo2.nii"
    
    def test_basic(self):
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii"])
    
    def test_output_directory(self):
        self.siena.output_directory = "bar"
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii", "-o", "bar"])
    
    def test_two_class_segmentation(self):
        self.siena.two_class_segmentation = True
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii", "-2"])
    
    def test_t2_weighted_input(self):
        self.siena.t2_weighted_input = True
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii", "-t2"])
    
    def test_masking(self):
        self.siena.masking = True
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii", "-m"])
    
    def test_ignore_upwards(self):
        self.siena.ignore_upwards = 10
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii", "-t", "10"])
    
    def test_ignore_downwards(self):
        self.siena.ignore_downwards = 10
        self.assertEqual(self.siena.command, ["siena", "foo1.nii", "foo2.nii", "-b", "10"])
    
#    def test_parse_report(self):
#        report = self.sienax.parse_report(os.path.join(os.path.dirname(__file__), "report.sienax"))
#        
#        self.assertEqual(report.get("version", None), "2.6")
#        self.assertEqual(report.get("vscale", None), 1.2580182223)
#        
#        self.assertEqual(report.get("grey", {}).get("normalized", None), 724097.23) 
#        self.assertEqual(report.get("grey", {}).get("raw", None), 575585.64)
#        
#        self.assertEqual(report.get("white", {}).get("normalized", None), 639077.80) 
#        self.assertEqual(report.get("white", {}).get("raw", None), 508003.61)
#        
#        self.assertEqual(report.get("brain", {}).get("normalized", None), 1363175.03) 
#        self.assertEqual(report.get("brain", {}).get("raw", None), 1083589.25)
#        
#        self.assertEqual(report.get("pgrey", {}).get("normalized", None), 556299.24) 
#        self.assertEqual(report.get("pgrey", {}).get("raw", None), 442202.84)
#        
#        self.assertEqual(report.get("vcsf", {}).get("normalized", None), 75373.30) 
#        self.assertEqual(report.get("vcsf", {}).get("raw", None), 59914.32)
    
    def test_default_output_directory(self):
        self.siena.input1 = "/foo/bar1.nii"
        self.siena.input2 = "/foo/bar2.nii"
        self.assertEqual(self.siena.default_output_directory, "bar1_to_bar2_siena")
    
if __name__ == "__main__" :
    unittest.main()