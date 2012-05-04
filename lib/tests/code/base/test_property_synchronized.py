import unittest

import medipy.base

class TestPropertySynchronized(unittest.TestCase):
    
    class TestClass(medipy.base.PropertySynchronized):
        def __init__(self):
            medipy.base.PropertySynchronized.__init__(self, ["value"])
            self._value = None
        
        def _get_value(self):
            return self._value
        
        def _set_value(self, value):
            self._value = value
            self.notify_observers("value")
        
        value = property(_get_value, _set_value)
    
    def test_synchronization(self):
        object1 = TestPropertySynchronized.TestClass()
        object2 = TestPropertySynchronized.TestClass()
        
        object1.value = 1
        object2.value = 2
        
        self.assertEqual(object1.value, 1)
        self.assertEqual(object2.value, 2)
        
        # Changing object1 will change object2
        object1.append_child(object2)
        object1.value = 3
        self.assertEqual(object1.value, 3)
        self.assertEqual(object2.value, 3)
        
        # Changing object2 will not change object1
        object2.value = 4
        self.assertEqual(object1.value, 3)
        self.assertEqual(object2.value, 4)
    
    def test_two_way_synchronization(self):
        object1 = TestPropertySynchronized.TestClass()
        object2 = TestPropertySynchronized.TestClass()
        
        object1.value = 1
        object2.value = 2
        
        self.assertEqual(object1.value, 1)
        self.assertEqual(object2.value, 2)
        
        object1.append_child(object2)
        object2.append_child(object1)
        
        # Changing object1 will change object2
        object1.append_child(object2)
        object1.value = 3
        self.assertEqual(object1.value, 3)
        self.assertEqual(object2.value, 3)
        
        # Changing object2 will change object1
        object2.value = 4
        self.assertEqual(object1.value, 4)
        self.assertEqual(object2.value, 4)
    
    def test_children(self):
        object1 = TestPropertySynchronized.TestClass()
        object2 = TestPropertySynchronized.TestClass()
        
        self.assertFalse(object1.has_child(object2))
        self.assertFalse(object2.has_child(object1))

        object1.append_child(object2)
        
        self.assertTrue(object1.has_child(object2))
        self.assertEqual(object1.child_index(object2), 0)
        self.assertRaises(ValueError, lambda : object2.child_index(object1))
        self.assertFalse(object2.has_child(object1))

if __name__ == '__main__':
    unittest.main()