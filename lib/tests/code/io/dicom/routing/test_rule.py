import unittest

import medipy.io.dicom

class TestRule(unittest.TestCase):
    def setUp(self):
        self.dataset = medipy.io.dicom.DataSet()
        self.dataset.patients_name = "Doe^John"
    
    def test_rule_no_match(self):
        condition = medipy.io.dicom.routing.AlwaysFalse()
        action = medipy.io.dicom.routing.SetElement("patient_id", "FOO")
        
        rule = medipy.io.dicom.routing.Rule(condition, [action])
        rule(self.dataset)
        
        self.assertFalse("patient_id" in self.dataset)
    
    def test_rule_match(self):
        condition = medipy.io.dicom.routing.AlwaysTrue()
        action = medipy.io.dicom.routing.SetElement("patient_id", "FOO")
        
        rule = medipy.io.dicom.routing.Rule(condition, [action])
        rule(self.dataset)
        
        self.assertTrue("patient_id" in self.dataset)
        self.assertEqual(self.dataset.patient_id, "FOO")

if __name__ == "__main__" :
    unittest.main()