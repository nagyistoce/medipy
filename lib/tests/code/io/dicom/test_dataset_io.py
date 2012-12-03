import os
import unittest

import medipy.io.dicom

class testParse(unittest.TestCase):
            
    def setUp(self):
        data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data"))
        
        self.dataset = medipy.io.dicom.read(
            os.path.join(data_directory, "input", "siemens_mosaic.dcm"))
        
    def test_type(self) :
        self.assertTrue(isinstance(self.dataset, medipy.io.dicom.DataSet))
    
    def test_single_value(self):
        self.assertEqual(self.dataset.header.media_storage_sop_instance_uid, 
                         "1.3.12.2.1107.5.2.32.35389.2010072309084257478998732")
    
    def test_multiple_value(self):
        self.assertEqual(self.dataset.pixel_spacing, [2., 2.])
    
    def test_sequence(self):
        self.assertTrue(isinstance(self.dataset.referenced_image_sequence, list))
        self.assertEqual(len(self.dataset.referenced_image_sequence), 3)
        self.assertTrue(isinstance(self.dataset.referenced_image_sequence[0],
                                   medipy.io.dicom.DataSet))
        self.assertEqual(
             self.dataset.referenced_image_sequence[0].referenced_sop_instance_uid,
             "1.3.12.2.1107.5.2.32.35389.201007230834425475364627")

if __name__ == '__main__':
    unittest.main()