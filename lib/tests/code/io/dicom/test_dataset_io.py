import os
import shutil
import tempfile
import unittest

import medipy.io.dicom

class testDataSetIO(unittest.TestCase):
            
    def setUp(self):
        data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data"))
        
        self.dataset = medipy.io.dicom.read(
            os.path.join(data_directory, "input", "siemens_mosaic.dcm"))
        
    def test_type(self) :
        self.assertTrue(isinstance(self.dataset, medipy.io.dicom.DataSet))
    
    def test_single_value(self):
        self.assertEqual(self.dataset.header.media_storage_sop_instance_uid.value, 
                         "1.3.12.2.1107.5.2.32.35389.2010072309084257478998732")
    
    def test_multiple_value(self):
        self.assertEqual(self.dataset.pixel_spacing.value, [2., 2.])
    
    def test_sequence(self):
        self.assertTrue(isinstance(self.dataset.referenced_image_sequence.value, list))
        self.assertEqual(len(self.dataset.referenced_image_sequence.value), 3)
        self.assertTrue(isinstance(self.dataset.referenced_image_sequence.value[0],
                                   medipy.io.dicom.DataSet))
        self.assertEqual(
             self.dataset.referenced_image_sequence.value[0].referenced_sop_instance_uid.value,
             "1.3.12.2.1107.5.2.32.35389.201007230834425475364627")
    
    def test_write(self):
        
        tempdir = tempfile.mkdtemp()
        
        medipy.io.dicom.write(self.dataset, os.path.join(tempdir, "foo.dcm"))
        other_dataset = medipy.io.dicom.read(os.path.join(tempdir, "foo.dcm"))
        
        shutil.rmtree(tempdir)
        
        # Check that other_dataset is included in self.dataset
        for tag in sorted(self.dataset) :
            if tag == (0x0008,0x0005) :
                # Specific Character Set is not tested, since file might be re-encoded
                continue
            elif tag.private :
                # TODO : private tags
                continue
            
            self.assertTrue(tag in other_dataset)
            self.assertEqual(self.dataset[tag], other_dataset[tag])
        # Check that dataset is included in other_dataset
        for tag in sorted(other_dataset) :
            if tag.private :
                # TODO : private tags
                continue
            self.assertTrue(tag in self.dataset)
            # No need to check for equality, we just did it.

if __name__ == '__main__':
    unittest.main()