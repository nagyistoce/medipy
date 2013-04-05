import unittest

import numpy

import medipy.io.dicom

class TestDataSet(unittest.TestCase):
    
    def setUp(self):
        self.dataset = medipy.io.dicom.DataSet()
        dict.__setitem__(self.dataset, medipy.io.dicom.Tag(0x0010,0x0010), 
                         medipy.io.dicom.PN("Doe^John"))
        dict.__setitem__(self.dataset, medipy.io.dicom.Tag(0x0010,0x0020), 
                         medipy.io.dicom.LO("DOE_Joh"))
    
    def test_numeric_tag_access(self):
        self.assertEqual(self.dataset[(0x0010,0x0010)].value, "Doe^John")
        self.assertEqual(self.dataset[0x00100020].value, "DOE_Joh")
        
        self.assertTrue((0x0010,0x0010) in self.dataset)
        self.assertTrue(0x00100020 in self.dataset)
        
        self.assertFalse(0x00100040 in self.dataset)
        self.assertFalse((0x0010,0x0040) in self.dataset)
        self.assertFalse(0x00100040 in self.dataset)
        self.assertFalse((0x0010,0x0040) in self.dataset)
        
        self.dataset[0x00100040] = medipy.io.dicom.CS("O")
        self.assertEqual(self.dataset[0x0010,0x0040].value, "O")
        
        del self.dataset[0x00100020]
        self.assertFalse(0x00100020 in self.dataset)
    
    def test_named_tag_access(self):
        self.assertEqual(self.dataset["patients_name"].value, "Doe^John")
        self.assertEqual(self.dataset.patient_id.value, "DOE_Joh")
        
        self.assertTrue("patients_name" in self.dataset)
        
        self.assertFalse("patients_sex" in self.dataset)
        self.assertFalse("foobar" in self.dataset)
        
        self.dataset.patients_sex = medipy.io.dicom.CS("O")
        self.assertEqual(self.dataset["patients_sex"].value, "O")
        
        delattr(self.dataset, "patients_name")
        self.assertFalse("patients_name" in self.dataset)
    
    def test_dict_interface(self):
        self.assertEqual(set(self.dataset.keys()), set([0x00100010, 0x00100020]))
        self.assertEqual(
            self.dataset.get((0x0010,0x0010), medipy.io.dicom.PN("")).value, "Doe^John")
        self.assertEqual(
            self.dataset.get(0x00100010, medipy.io.dicom.PN("")).value, "Doe^John")
        self.assertEqual(
            self.dataset.get("patients_name", medipy.io.dicom.LO("")).value, "Doe^John")
        
        self.assertEqual(
            self.dataset.get((0x0010,0x0040), medipy.io.dicom.CS("O")).value, "O")
        self.assertEqual(
            self.dataset.get(0x00100040, medipy.io.dicom.CS("O")).value, "O")
        self.assertEqual(
            self.dataset.get("patients_sex", medipy.io.dicom.CS("O")).value, "O")
    
    def test_tags(self):
        self.assertEqual(set(self.dataset.tags()), 
                         set(["patients_name", "patient_id"]))
    
    def test_pixel_array(self):
        array = numpy.arange(12, dtype=numpy.uint16).reshape((4,3))
        
        self.dataset.pixel_data = medipy.io.dicom.OB(array.flatten())
        self.dataset.rows = medipy.io.dicom.US(4)
        self.dataset.columns = medipy.io.dicom.US(3)
        self.dataset.bits_allocated = medipy.io.dicom.US(16)
        self.dataset.pixel_representation = medipy.io.dicom.US(0)
        self.dataset.transfer_syntax_uid = medipy.io.dicom.UI(None)
        
        numpy.testing.assert_equal(self.dataset.pixel_array, array)
    
    def test_constructor_parameters(self):
        dataset = medipy.io.dicom.DataSet(patients_name="Doe^John")
        self.assertTrue("patients_name" in dataset)
        self.assertTrue(0x00100010 in dataset)
        self.assertTrue(isinstance(dataset.patients_name, medipy.io.dicom.PN))
        self.assertEqual(dataset.patients_name.value, "Doe^John")

if __name__ == '__main__':
    unittest.main()
