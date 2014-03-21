import unittest
import medipy.base

class TestEnum(unittest.TestCase):
    
    def test_int_creation(self):
        Colors = medipy.base.enum("Colors", "red", "green", "blue")
        self.assertEqual(Colors.red, 0)
    
    def test_value_creation(self):
        Colors = medipy.base.enum("Colors", "red", green="green", blue=3.14)
        self.assertEqual(Colors.red, 0)
        self.assertEqual(Colors.green, "green")
        self.assertEqual(Colors.blue, 3.14)
    
    def test_name(self):
        Colors = medipy.base.enum("Colors", "red", green="green", blue=3.14)
        self.assertEqual(Colors[0], "red")
        self.assertEqual(Colors["green"], "green")
        self.assertEqual(Colors[3.14], "blue")
    
    def test_members(self):
        Colors = medipy.base.enum("Colors", "red", green="green", blue=3.14)
        self.assertEqual(Colors.members, ["red", "green", "blue"])
    
    def test_type(self):
        Colors = medipy.base.enum("Colors", "red", green="green", blue=3.14)
        self.assertTrue(isinstance(Colors, medipy.base.Enum))

if __name__ == "__main__" :
    unittest.main()
