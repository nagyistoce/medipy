import unittest

import medipy.io.dicom

class TestConditions(unittest.TestCase):
    def setUp(self):
        self.dataset = medipy.io.dicom.DataSet()
        
    def test_always_true(self):
        condition = medipy.io.dicom.routing.AlwaysTrue()
        self.assertTrue(condition(self.dataset))
    
    def test_always_false(self):
        condition = medipy.io.dicom.routing.AlwaysFalse()
        self.assertFalse(condition(self.dataset))
    
    def test_and(self):
        condition = medipy.io.dicom.routing.And(
            medipy.io.dicom.routing.AlwaysTrue(),
            medipy.io.dicom.routing.AlwaysTrue())
        
        self.assertTrue(condition(self.dataset))
        
        condition = medipy.io.dicom.routing.And(
            medipy.io.dicom.routing.AlwaysTrue(),
            medipy.io.dicom.routing.AlwaysFalse())
        
        self.assertFalse(condition(self.dataset))
    
    def test_or(self):
        condition = medipy.io.dicom.routing.Or(
            medipy.io.dicom.routing.AlwaysTrue(),
            medipy.io.dicom.routing.AlwaysFalse())
        
        self.assertTrue(condition(self.dataset))
        
        condition = medipy.io.dicom.routing.Or(
            medipy.io.dicom.routing.AlwaysFalse(),
            medipy.io.dicom.routing.AlwaysFalse())
        
        self.assertFalse(condition(self.dataset))
    
    def test_not(self):
        condition = medipy.io.dicom.routing.Not(
            medipy.io.dicom.routing.AlwaysTrue())
        self.assertFalse(condition(self.dataset))
    
    def test_element_match(self):
        condition = medipy.io.dicom.routing.ElementMatch(
            "patients_name", "Doe^John")
        self.assertFalse(condition(self.dataset))
        
        self.dataset.patients_name = "Foo"
        self.assertFalse(condition(self.dataset))
        
        self.dataset.patients_name = "Doe^John"
        self.assertTrue(condition(self.dataset))

if __name__ == "__main__" :
    unittest.main()