import unittest

import medipy.base

class TestObservable(unittest.TestCase):

    def test_allowed_events(self):
        observable = medipy.base.Observable(["foo"])
        self.assertEqual(set(observable.allowed_events), set(["foo", "any"]))
        self.assertTrue(observable.is_allowed_event("foo"))
        self.assertFalse(observable.is_allowed_event("bar"))
        observable.add_allowed_event("bar")
        self.assertEqual(set(observable.allowed_events), set(["foo", "bar", "any"]))
    
    def test_observers(self):
        
        class Observer(object):
            def f(self, event):
                pass
        
        def observer_function(event):
            pass
        
        observable = medipy.base.Observable(["foo"])
        
        observer_object = Observer()
        observable.add_observer("foo", observer_object.f)
        self.assertTrue(observable.has_observer("foo", observer_object.f))
        self.assertFalse(observable.add_observer("foo", observer_object.f))
        observable.remove_observer("foo", observer_object.f)
        self.assertFalse(observable.has_observer("foo", observer_object.f))
        self.assertRaises(Exception, observable.add_observer, ("bar", observer_object.f))
        
        observable.add_observer("foo", observer_function)
        self.assertTrue(observable.has_observer("foo", observer_function))
        self.assertFalse(observable.add_observer("foo", observer_function))
        observable.remove_observer("foo", observer_function)
        self.assertFalse(observable.has_observer("foo", observer_function))
        self.assertRaises(Exception, observable.add_observer, ("bar", observer_function))
    
    def test_events(self):
        
        class Observer(object):
            def __init__(self):
                self.counter_1 = 0
                self.counter_2 = 0
                self.counter_3 = 0
            
            def f1(self, event):
                self.counter_1 += 1
            
            def f2(self, event):
                self.counter_2 += 1
            
            def f3(self, event):
                self.counter_3 += 1
        
        counters = {1 : 0, 2 : 0, 3 : 0}
        def observer_function(event, counter):
            counters[counter] += 1
        
        observable = medipy.base.Observable(["foo", "bar"])
        
        observer_object = Observer()
        observable.add_observer("foo", observer_object.f1)
        observable.add_observer("bar", observer_object.f2)
        observable.add_observer("any", observer_object.f3)
        
        observable.notify_observers("foo")
        self.assertEqual(observer_object.counter_1, 1)
        self.assertEqual(observer_object.counter_2, 0)
        self.assertEqual(observer_object.counter_3, 1)
        
        observable.notify_observers("bar")
        self.assertEqual(observer_object.counter_1, 1)
        self.assertEqual(observer_object.counter_2, 1)
        self.assertEqual(observer_object.counter_3, 2)
        
        observable.add_observer("foo", lambda event : observer_function(event, 1))
        observable.add_observer("bar", lambda event : observer_function(event, 2))
        observable.add_observer("any", lambda event : observer_function(event, 3))
        
        observable.notify_observers("foo")
        self.assertEqual(counters, {1 : 1, 2 : 0, 3 : 1})
        
        observable.notify_observers("bar")
        self.assertEqual(counters, {1 : 1, 2 : 1, 3 : 2})
        
if __name__ == '__main__':
    unittest.main()