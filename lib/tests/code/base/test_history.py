import unittest

import medipy.base

class TestHistory(unittest.TestCase):
    class Command(medipy.base.UndoableCommand) :
        def __init__(self, label) :
            self.execution_count = 0
            self.undo_count = 0
            medipy.base.UndoableCommand.__init__(self, label)
        
        def execute(self) :
            self.execution_count += 1
        
        def undo(self) :
            self.undo_count += 1
            
    def setUp(self):
        self.old_old_command = TestHistory.Command("Old old")
        self.old_command = TestHistory.Command("Old")
        self.recent_command = TestHistory.Command("Recent")
        
        self.history = medipy.base.History()
    
    def test_empty(self) :
        # Test history state
        self.assertTrue(self.history.empty)
        self.assertFalse(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, [])
    
    def test_one_command(self) :
        self.history.add(self.recent_command)
        
        # Test command state
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 0)
        
        # Test history state
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent"])
        
        # Test undo
        self.history.undo()
        
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 1)
        
        self.assertFalse(self.history.empty)
        self.assertFalse(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent"])
        
        # Test redo
        self.history.redo()
        
        self.assertEqual(self.recent_command.execution_count, 2)
        self.assertEqual(self.recent_command.undo_count, 1)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent"])
    
    def test_multiple_commands(self) :
        self.history.add(self.old_old_command)
        self.history.add(self.old_command)
        self.history.add(self.recent_command)
        
        # Test commands state
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 0)
        
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 0)
        
        # Test history state
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        
        # Test undo
        self.history.undo()
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 1)
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 0)
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        
        # Test redo
        self.history.redo()
        self.assertEqual(self.recent_command.execution_count, 2)
        self.assertEqual(self.recent_command.undo_count, 1)
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 0)
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        
        # Test multiple undo
        self.history.undo(2)
        
        self.assertEqual(self.recent_command.execution_count, 2)
        self.assertEqual(self.recent_command.undo_count, 2)
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 1)
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        
        # Test multiple redo
        self.history.redo(2)
        
        self.assertEqual(self.recent_command.execution_count, 3)
        self.assertEqual(self.recent_command.undo_count, 2)
        self.assertEqual(self.old_command.execution_count, 2)
        self.assertEqual(self.old_command.undo_count, 1)
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        
if __name__ == '__main__':
    unittest.main()