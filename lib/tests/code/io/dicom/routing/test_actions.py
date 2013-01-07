import unittest

import medipy.io.dicom

class TestActions(unittest.TestCase):
    def setUp(self):
        self.dataset = medipy.io.dicom.DataSet()
    
    def test_set_element_absent(self):
        action = medipy.io.dicom.routing.SetElement("patients_name", "Doe^John")
        action(self.dataset)
        
        self.assertTrue("patients_name" in self.dataset)
        self.assertEqual(self.dataset.patients_name, "Doe^John")
    
    def test_set_element_present(self):
        self.dataset.patients_name = "Brouchard^Georges"
        
        action = medipy.io.dicom.routing.SetElement("patients_name", "Doe^John")
        action(self.dataset)
        
        self.assertTrue("patients_name" in self.dataset)
        self.assertEqual(self.dataset.patients_name, "Doe^John")
    
    def test_delete_element_absent(self):
        action = medipy.io.dicom.routing.DeleteElement("patients_name")
        action(self.dataset)
        
        self.assertFalse("patients_name" in self.dataset)
    
    def test_delete_element_present(self):
        self.dataset.patients_name = "Brouchard^Georges"
        
        action = medipy.io.dicom.routing.DeleteElement("patients_name")
        action(self.dataset)
        
        self.assertFalse("patients_name" in self.dataset)
    
    def test_empty_element(self):
        action = medipy.io.dicom.routing.EmptyElement("patients_name")
        action(self.dataset)
        
        self.assertTrue("patients_name" in self.dataset)
        self.assertEqual(self.dataset.patients_name, "")
        
        action = medipy.io.dicom.routing.EmptyElement("language_code_sequence")
        action(self.dataset)
        
        self.assertTrue("language_code_sequence" in self.dataset)
        self.assertEqual(self.dataset.language_code_sequence, [])

if __name__ == "__main__" :
    unittest.main()