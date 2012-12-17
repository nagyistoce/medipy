import datetime
import unittest

import medipy.io.dicom
import medipy.io.dicom.misc

class TestMisc(unittest.TestCase):
    
    def test_parse_da(self):
        parsed = medipy.io.dicom.misc.parse_da("20110721")
        expected = datetime.datetime(2011, 07, 21)
        self.assertEqual(parsed, expected)
    
    def test_parse_tm(self):
        
        parsed = medipy.io.dicom.misc.parse_tm("17")
        expected = datetime.time(17, 00, 00)
        self.assertEqual(parsed, expected)
        
        parsed = medipy.io.dicom.misc.parse_tm("1733")
        expected = datetime.time(17, 33, 00)
        self.assertEqual(parsed, expected)
        
        parsed = medipy.io.dicom.misc.parse_tm("173325")
        expected = datetime.time(17, 33, 25)
        self.assertEqual(parsed, expected)
        
        parsed = medipy.io.dicom.misc.parse_tm("173325.1234")
        expected = datetime.time(17, 33, 25, 123400)
        self.assertEqual(parsed, expected)
    
    def test_generate_uid(self):
        
        uid_vr = medipy.io.dicom.generate_uid()
        self.assertTrue(isinstance(uid_vr, medipy.io.dicom.UI))
        self.assertTrue(isinstance(uid_vr.value, str))
        
        uid_str = medipy.io.dicom.generate_uid(False)
        self.assertTrue(isinstance(uid_str, str))
        
        uid1 = medipy.io.dicom.generate_uid()
        uid2 = medipy.io.dicom.generate_uid()
        self.assertNotEqual(uid1.value, uid2.value)

if __name__ == '__main__':
    unittest.main()