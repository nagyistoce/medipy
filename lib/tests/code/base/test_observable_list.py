import unittest

import medipy.base

class TestObservableList(unittest.TestCase):
    
    def setUp(self):
        self.events_fired = {}
        self.event = None
    
    def observer(self, event):
        self.events_fired.setdefault(event.event, 0)
        self.events_fired[event.event] += 1
        self.event = event
    
    def test_init_empty(self):
        observable_list = medipy.base.ObservableList()
        self.assertEqual(len(observable_list), 0)
    
    def test_init_sequence(self):
        sequence = [1,2,3]
        observable_list = medipy.base.ObservableList(sequence)
        self.assertEqual(observable_list, sequence)
    
    def test_setitem(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list[1] = 4
        
        # Check value
        self.assertEqual(observable_list, [1,4,3])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("set_item", 0), 1)
        # Check event content
        self.assertEqual(self.event.index, 1)
        self.assertEqual(self.event.old_value, 2)
    
    def test_delitem(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        del observable_list[1]
        
        # Check value
        self.assertEqual(observable_list, [1,3])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("delete_item", 0), 1)
        # Check event content
        self.assertEqual(self.event.index, 1)
        self.assertEqual(self.event.old_value, 2)
    
    def test_append(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list.append(4)
        
        # Check value
        self.assertEqual(observable_list, [1,2,3,4])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("append", 0), 1)
    
    def test_pop(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list.pop()
        
        # Check value
        self.assertEqual(observable_list, [1,2])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("pop", 0), 1)
        # Check event content
        self.assertEqual(self.event.index, 2)
        self.assertEqual(self.event.popped_value, 3)
    
    def test_extend(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        
        extension = [4,5,6]
        observable_list.extend(extension)
        
        # Check value
        self.assertEqual(observable_list, [1,2,3,4,5,6])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("extend", 0), 1)
        # Check event content
        self.assertEqual(self.event.value, extension)
    
    def test_insert(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list.insert(1, 4)
        
        # Check value
        self.assertEqual(observable_list, [1,4,2,3])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("insert", 0), 1)
        # Check event content
        self.assertEqual(self.event.index, 1)
        self.assertEqual(self.event.value, 4)
    
    def test_remove(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list.remove(1)
        
        # Check value
        self.assertEqual(observable_list, [2,3])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("remove", 0), 1)
        # Check event content
        self.assertEqual(self.event.index, 0)
        self.assertEqual(self.event.value, 1)
    
    def test_reverse(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list.reverse()
        
        # Check value
        self.assertEqual(observable_list, [3,2,1])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("reverse", 0), 1)
    
    def test_sort(self):
        observable_list = medipy.base.ObservableList([1,2,3])
        observable_list.add_observer("any", self.observer)
        observable_list.sort(reverse=True)
        
        # Check value
        self.assertEqual(observable_list, [3,2,1])
        # Check fired events
        self.assertEqual(len(self.events_fired), 1)
        self.assertEqual(self.events_fired.get("sort", 0), 1)
        # Check event content
        self.assertTrue(self.event.reverse)

if __name__ == '__main__':
    unittest.main()
