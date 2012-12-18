# encoding: utf-8

import os
import shutil
import sys
import tarfile
import tempfile
import unittest

import wx
import medipy.gui.dicom
import medipy.io.dicom

if sys.platform == "linux2" and "DISPLAY" in os.environ :

    class TestHierarchyTree(unittest.TestCase):
        
        def extract(self, filename):
            data_directory = os.path.abspath(os.path.join(
                os.path.dirname(__file__), "..", "..", "..", "data","input"))
            
            tempdir = tempfile.mkdtemp()
            
            archive = tarfile.open(os.path.join(data_directory, filename), "r")
            archive.extractall(tempdir)
            
            return tempdir
        
        def cleanup(self, tempdir):
            shutil.rmtree(tempdir)
        
        def setUp(self):
            self.app = wx.PySimpleApp()
            self.frame = wx.Frame(None)
            self.hierarchy_tree = medipy.gui.dicom.HierarchyTree(self.frame)
            
            sizer = wx.BoxSizer()
            sizer.Add(self.hierarchy_tree, 1, wx.EXPAND)
            self.frame.SetSizer(sizer)
            self.frame.Show()
            
        def tearDown(self):
            self.frame.Close()
            self.app.Exit()
        
        def test_content(self):
            tempdir = self.extract("brainix.tgz")
            dicomdir = medipy.io.dicom.read(os.path.join(tempdir, "BRAINIX", "DICOMDIR"))
            
            self.hierarchy_tree.set_datasets([dicomdir])
            
            self.assertEqual(self.hierarchy_tree.GetCount(), 9)
            
            labels = ["BRAINIX", 
                      u"IRM cérébrale, neuro-crâne", 
                      "sT2/TSE/T", "sT2W/FLAIR", "T2W/FE-EPI", 
                      "T1/SE/extrp", "SOUS", "T1/SE/extrp", 
                      "T1/3D/FFE/C"]
            
            for index, node in enumerate(traverse(self.hierarchy_tree)) :
                self.assertEqual(self.hierarchy_tree.GetItemText(node), labels[index])
                
                if labels[index] not in ["BRAINIX", u"IRM cérébrale, neuro-crâne"] : 
                    py_data = self.hierarchy_tree.GetItemPyData(node)
                    self.assertTrue(isinstance(py_data, list))
                    self.assertEqual(len(py_data), 1)
                    
                    record = py_data[0]
                    self.assertTrue(isinstance(record, medipy.io.dicom.DataSet))
                    self.assertEqual(
                        record.get("directory_record_type", medipy.io.dicom.CS("")).value,
                        "SERIES")
                    self.assertEqual(
                        record.get("modality", medipy.io.dicom.CS("")).value,
                        "MR")
            
            self.cleanup(tempdir)
        
def traverse(treectrl, root=None):
    if root is None:
        root = treectrl.GetRootItem()
    
    item, cookie = treectrl.GetFirstChild(root)
    while item.IsOk():
        yield item 
        for child in traverse(treectrl, item) :
            yield child
        item, cookie = treectrl.GetNextChild(root, cookie)

if __name__ == "__main__" :
    unittest.main()

