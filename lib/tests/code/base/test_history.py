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
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        
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
    
    def test_limited_size(self):
        self.history = medipy.base.History(2)
        
        self.history.add(self.old_old_command)
        self.history.add(self.old_command)
        
        # Test commands state
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 0)
        
        # Test history state
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Old", "Old old"])
        
        # Add new command, discarding earliest one
        self.history.add(self.recent_command)
        
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
        self.assertEqual(self.history.labels, ["Recent", "Old"])
    
    def test_add_after_undo(self):
        self.history.add(self.old_old_command)
        self.history.add(self.old_command)
        
        # Test commands state
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 0)
        
        self.assertEqual(self.recent_command.execution_count, 0)
        self.assertEqual(self.recent_command.undo_count, 0)
        
        # Test history state
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Old", "Old old"])
        
        # Undo
        self.history.undo()
        
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 1)
        
        self.assertEqual(self.recent_command.execution_count, 0)
        self.assertEqual(self.recent_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Old", "Old old"])
        
        # Add new command
        self.history.add(self.recent_command)
        
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 1)
        
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old old"])
    
    def test_cursor(self):
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
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        self.assertEqual(self.history.cursor, 0)
        
        # Test undo
        self.history.cursor = 2
        
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 1)
        
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 1)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        self.assertEqual(self.history.cursor, 2)
        
        # Test redo
        self.history.cursor = 1
        
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 0)
        
        self.assertEqual(self.old_command.execution_count, 2)
        self.assertEqual(self.old_command.undo_count, 1)
        
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 1)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        self.assertEqual(self.history.cursor, 1)
    
    def test_cursor_limit(self):
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
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        self.assertEqual(self.history.cursor, 0)
        
        # Test undo everything
        self.history.cursor = self.history.steps_count
        
        self.assertEqual(self.old_old_command.execution_count, 1)
        self.assertEqual(self.old_old_command.undo_count, 1)
        
        self.assertEqual(self.old_command.execution_count, 1)
        self.assertEqual(self.old_command.undo_count, 1)
        
        self.assertEqual(self.recent_command.execution_count, 1)
        self.assertEqual(self.recent_command.undo_count, 1)
        
        self.assertFalse(self.history.empty)
        self.assertFalse(self.history.can_undo)
        self.assertTrue(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        self.assertEqual(self.history.cursor, 3)
        
        # Test redo everything
        self.history.cursor = 0
        
        self.assertEqual(self.old_old_command.execution_count, 2)
        self.assertEqual(self.old_old_command.undo_count, 1)
        
        self.assertEqual(self.old_command.execution_count, 2)
        self.assertEqual(self.old_command.undo_count, 1)
        
        self.assertEqual(self.recent_command.execution_count, 2)
        self.assertEqual(self.recent_command.undo_count, 1)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Recent", "Old", "Old old"])
        self.assertEqual(self.history.cursor, 0)
        
        # Test undo everything, then add a new step
        newest_command = TestHistory.Command("Only me")
        self.history.cursor = self.history.steps_count
        self.history.add(newest_command)
        
        self.assertEqual(self.old_old_command.execution_count, 2)
        self.assertEqual(self.old_old_command.undo_count, 2)
        
        self.assertEqual(self.old_command.execution_count, 2)
        self.assertEqual(self.old_command.undo_count, 2)
        
        self.assertEqual(self.recent_command.execution_count, 2)
        self.assertEqual(self.recent_command.undo_count, 2)
        
        self.assertEqual(newest_command.execution_count, 1)
        self.assertEqual(newest_command.undo_count, 0)
        
        self.assertFalse(self.history.empty)
        self.assertTrue(self.history.can_undo)
        self.assertFalse(self.history.can_redo)
        self.assertEqual(self.history.labels, ["Only me"])
        self.assertEqual(self.history.cursor, 0)
        
if __name__ == '__main__':
    unittest.main()