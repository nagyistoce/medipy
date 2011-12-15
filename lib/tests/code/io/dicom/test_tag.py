import unittest

import medipy.io.dicom

class TestTag(unittest.TestCase):
    
    def test_from_int(self):
        tag = medipy.io.dicom.Tag(0x00100020)
        self.assertEqual(tag.group, 0x0010)
        self.assertEqual(tag.element, 0x0020)
        self.assertFalse(tag.private)
    
    def test_from_tuple(self):
        tag = medipy.io.dicom.Tag(0x0010,0x0020)
        self.assertEqual(tag.group, 0x0010)
        self.assertEqual(tag.element, 0x0020)
        self.assertFalse(tag.private)
    
    def test_private(self):
        public_tag = medipy.io.dicom.Tag(0x0010,0x0020)
        self.assertFalse(public_tag.private)
        
        private_tag = medipy.io.dicom.Tag(0x0029,0x1010)
        self.assertTrue(private_tag.private)

if __name__ == '__main__':
    unittest.main()