import os
import shutil
import tempfile
import unittest

import numpy

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
        self.assertTrue(isinstance(self.dataset.referenced_image_sequence, medipy.io.dicom.SQ))
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
        
        self._compare_datasets(self.dataset, other_dataset)
    
    def _compare_datasets(self, dataset1, dataset2):
        for tag in sorted(dataset1) :
            if tag == (0x0008,0x0005) :
                # Specific Character Set is not tested, since file might be re-encoded
                continue
            
            self.assertTrue(tag in dataset2)
            if isinstance(dataset1[tag], medipy.io.dicom.SQ) :
                self.assertEqual(len(dataset1[tag].value), len(dataset2[tag].value))
                for index, sub_dataset1 in enumerate(dataset1[tag].value) :
                    sub_dataset2 = dataset2[tag].value[index]
                    self._compare_datasets(sub_dataset1, sub_dataset2)
            else :
                self.assertEqual(dataset1[tag].__class__, dataset2[tag].__class__)
                
                value1 = dataset1[tag].value
                value2 = dataset2[tag].value
                
                float_vrs = (
                    medipy.io.dicom.DS, medipy.io.dicom.FD, medipy.io.dicom.FL)
                binary_vrs = (
                    medipy.io.dicom.OB, medipy.io.dicom.OF, medipy.io.dicom.OW,
                    medipy.io.dicom.UN)
                
                if isinstance(dataset1[tag], float_vrs) :
                    numpy.testing.assert_array_almost_equal(
                        value1, value2,
                        err_msg="Failure for tag {0}".format(tag))
                elif isinstance(dataset1[tag], binary_vrs) :
                    numpy.testing.assert_array_equal(value1, value2)
                else :
                    self.assertEqual(value1, value2)
        # Check that dataset is included in other_dataset
        for tag in sorted(dataset2) :
            self.assertTrue(tag in dataset1)
            # No need to check for equality, we just did it.

if __name__ == '__main__':
    unittest.main()
