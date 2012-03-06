import datetime
import unittest

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

if __name__ == '__main__':
    unittest.main()