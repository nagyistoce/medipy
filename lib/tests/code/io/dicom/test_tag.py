import unittest

import medipy.io.dicom

class TestTag(unittest.TestCase):
    
    def test_from_int(self):
        tag = medipy.io.dicom.Tag(0x00100020)
        self.assertEqual(tag.group, 0x0010)
        self.assertEqual(tag.element, 0x0020)
        self.assertFalse(tag.private)
    
    def test_high_from_int(self):
        # Test tag > (7fff,ffff) (i.e. > 2^31) for 32 bits systems
        tag = medipy.io.dicom.Tag(0xfffee00d)
        self.assertEqual(tag.group, 0xfffe)
        self.assertEqual(tag.element, 0xe00d)
        self.assertFalse(tag.private)
    
    def test_from_tuple(self):
        tag = medipy.io.dicom.Tag(0x0010,0x0020)
        self.assertEqual(tag.group, 0x0010)
        self.assertEqual(tag.element, 0x0020)
        self.assertFalse(tag.private)
    
    def test_high_from_tuple(self):
        # Test tag > (7fff,ffff) (i.e. > 2^31) for 32 bits systems
        tag = medipy.io.dicom.Tag(0xfffe, 0xe00d)
        self.assertEqual(tag.group, 0xfffe)
        self.assertEqual(tag.element, 0xe00d)
        self.assertFalse(tag.private)
    
    def test_private(self):
        public_tag = medipy.io.dicom.Tag(0x0010,0x0020)
        self.assertFalse(public_tag.private)
        
        private_tag = medipy.io.dicom.Tag(0x0029,0x1010)
        self.assertTrue(private_tag.private)
    
    def test_equality(self):
        tag1 = medipy.io.dicom.Tag(0x00100010)
        tag2 = medipy.io.dicom.Tag(0x00100010)
        tag3 = medipy.io.dicom.Tag(0x00100020)
        
        self.assertEqual(tag1, tag2)
        self.assertNotEqual(tag1, tag3)
        
        self.assertEqual(tag1, (0x0010,0x0010))
        self.assertEqual(tag1, 0x00100010)
        
        self.assertNotEqual(tag1, (0x0010,0x0020))
        self.assertNotEqual(tag1, 0x00100020)

if __name__ == '__main__':
    unittest.main()