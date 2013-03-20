import unittest

import medipy.io.dicom

class TestActions(unittest.TestCase):
    def setUp(self):
        self.dataset = medipy.io.dicom.DataSet()
    
    def test_set_element_absent(self):
        action = medipy.io.dicom.routing.SetElement("patients_name", "Doe^John")
        action(self.dataset)
        
        self.assertTrue("patients_name" in self.dataset)
        self.assertEqual(self.dataset.patients_name.value, "Doe^John")
    
    def test_set_element_present(self):
        self.dataset.patients_name = "Brouchard^Georges"
        
        action = medipy.io.dicom.routing.SetElement("patients_name", "Doe^John")
        action(self.dataset)
        
        self.assertTrue("patients_name" in self.dataset)
        self.assertEqual(self.dataset.patients_name.value, "Doe^John")
    
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
        self.assertEqual(self.dataset.patients_name.value, "")
        
        action = medipy.io.dicom.routing.EmptyElement("language_code_sequence")
        action(self.dataset)
        
        self.assertTrue("language_code_sequence" in self.dataset)
        self.assertEqual(self.dataset.language_code_sequence.value, [])
    
    def test_modify_dataset_dict(self):
        dataset = medipy.io.dicom.DataSet()
        dataset.patients_name = "Doe^John"
        dataset.series_description = "3D T1 SPGR 180 slices 1mm"
        action = medipy.io.dicom.routing.ModifyDataSet(
            {"series_description": "T1_3D"}, "source", "system", "COERCE")
        action(dataset)
        self.assertEqual(dataset.series_description.value, "T1_3D")
    
    def test_modify_dataset_callable(self):
        def modify_series_description(dataset) :
            series_description = dataset.get(
                "series_description", medipy.io.dicom.LO(None)).value
            if series_description == "3D T1 SPGR 180 slices 1mm" :
                return { "series_description" : "T1_3D" }
            else :
                return {}
        
        dataset = medipy.io.dicom.DataSet()
        dataset.patients_name = "Doe^John"
        dataset.series_description = "3D T1 SPGR 180 slices 1mm"
        action = medipy.io.dicom.routing.ModifyDataSet(
            modify_series_description, "source", "system", "COERCE")
        action(dataset)
        
        self.assertEqual(dataset.series_description.value, "T1_3D")
    
    def test_restore_dataset(self):
        dataset = medipy.io.dicom.DataSet()
        dataset.patients_name = "Doe^John"
        dataset.series_description = "3D T1 SPGR 180 slices 1mm"
        modify = medipy.io.dicom.routing.ModifyDataSet(
            {"series_description": "T1_3D"}, "source", "system", "COERCE")
        modify(dataset)
        
        restore = medipy.io.dicom.routing.RestoreDataSet("source", "system", "COERCE")
        restore(dataset)
        
        self.assertEqual(dataset.series_description.value, "3D T1 SPGR 180 slices 1mm")
        self.assertFalse("original_attributes_sequence" in dataset)
    
if __name__ == "__main__" :
    unittest.main()