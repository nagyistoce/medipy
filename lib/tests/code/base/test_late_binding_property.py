import unittest

import medipy.base

class TestLateBindingProperty(unittest.TestCase):
    
    def test(self):
        class Base(object):
            
            def __init__(self):
                self._x = 0
            
            def get_x(self):
                self.counter += 1
                return self._x
            
            x = property(get_x)
        
        class Derived(Base):
            
            def __init__(self):
                super(Derived, self).__init__()
                self.counter = 0
                
            def get_x(self):
                self.counter += 1
                return super(Derived, self).get_x()
        
        d = Derived()
        d.x
        self.assertEqual(d.counter, 1)

if __name__ == '__main__':
    unittest.main()
