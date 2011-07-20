import unittest

import numpy

import medipy.io.dicom

class TestDataSet(unittest.TestCase):
    
    def setUp(self):
        self.dictionary = { (0x0010,0x0010) : "Doe^John",
                            (0x0010,0x0020) : "DOE_Joh" }
        self.dataset =  medipy.io.dicom.DataSet.from_dict(self.dictionary)
    
    def test_numeric_tag_access(self):
        self.assertEqual(self.dataset[(0x0010,0x0010)], "Doe^John")
        self.assertEqual(self.dataset[0x00100020], "DOE_Joh")
        
        self.assertTrue((0x0010,0x0010) in self.dataset)
        self.assertTrue(0x00100020 in self.dataset)
        
        self.assertFalse(0x00100040 in self.dataset)
        self.assertFalse((0x0010,0x0040) in self.dataset)
        
        self.dataset[0x00100040] = "O"
        self.assertEqual(self.dataset[0x0010,0x0040], "O")
    
    def test_named_tag_access(self):
        self.assertEqual(self.dataset["patients_name"], "Doe^John")
        self.assertEqual(self.dataset.patient_id, "DOE_Joh")
        
        self.assertTrue("patients_name" in self.dataset)
        
        self.assertFalse("patients_sex" in self.dataset)
        
        self.dataset.patients_sex = "O"
        self.assertEqual(self.dataset["patients_sex"], "O")
    
    def test_dict_interface(self):
        self.assertEqual(set(self.dataset.keys()), set([0x00100010, 0x00100020]))
        self.assertEqual(self.dataset.get((0x0010,0x0010), ""), "Doe^John")
        self.assertEqual(self.dataset.get(0x00100010, ""), "Doe^John")
        self.assertEqual(self.dataset.get("patients_name", ""), "Doe^John")
        
        self.assertEqual(self.dataset.get((0x0010,0x0040), "O"), "O")
        self.assertEqual(self.dataset.get(0x00100040, "O"), "O")
        self.assertEqual(self.dataset.get("patients_sex", "O"), "O")
    
    def test_tags(self):
        self.assertEqual(set(self.dataset.tags()), 
                         set(["patients_name", "patient_id"]))
    
    def test_pixel_array(self):
        array = numpy.arange(12, dtype=numpy.uint16).reshape((4,3))
        
        self.dataset.pixel_data = array.flatten()
        self.dataset.rows = 4
        self.dataset.columns = 3
        self.dataset.bits_allocated = 16
        self.dataset.pixel_representation = 0
        self.dataset.transfer_syntax_uid = None
        
        numpy.testing.assert_equal(self.dataset.pixel_array, array)

if __name__ == '__main__':
    unittest.main()
