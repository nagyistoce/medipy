import os
import shutil
import tarfile
import tempfile
import unittest

import medipy.base
import medipy.io.dicom

class TestDicomSeries(unittest.TestCase):
    def extract(self, filename):
        data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data","input"))
        
        tempdir = tempfile.mkdtemp()
        
        archive = tarfile.open(os.path.join(data_directory, filename), "r")
        archive.extractall(tempdir)
        
        return tempdir
    
    def cleanup(self, tempdir):
        shutil.rmtree(tempdir)
    
    def test_root(self):
        path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), 
            "..", "..", "..", "data", "input", "dicom_series_brainix"))
        
        dicom_series = medipy.io.dicom.DicomSeries.read(path)
        self.assertEqual(dicom_series.root, "dicomdir:./BRAINIX/DICOMDIR")
        
    def test_uid(self):
        path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), 
            "..", "..", "..", "data", "input", "dicom_series_brainix"))
        
        tempdir = self.extract("brainix.tgz")
        old_dir = os.getcwd()
        os.chdir(tempdir)
        
        shutil.copy(path, os.path.join(tempdir, "dicom_series"))
        path = os.path.join(tempdir, "dicom_series")
        
        uid = "1.3.46.670589.11.0.0.11.4.2.0.8743.5.5396.2006120114285654497"
        
        dicom_series = medipy.io.dicom.DicomSeries.read(path)
        self.assertTrue(dicom_series.has_uid(uid))
        self.assertFalse(dicom_series.has_uid(""))
         
        image = medipy.io.load(dicom_series.url_from_uid(uid))
        self.assertEqual(image.metadata["series_description"], "sT2W/FLAIR")
        
        os.chdir(old_dir)
        self.cleanup(tempdir)
    
    def test_description(self):
        path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), 
            "..", "..", "..", "data", "input", "dicom_series_brainix"))
        
        tempdir = self.extract("brainix.tgz")
        old_dir = os.getcwd()
        os.chdir(tempdir)
        
        shutil.copy(path, os.path.join(tempdir, "dicom_series"))
        path = os.path.join(tempdir, "dicom_series")
        
        description = "sT2W/FLAIR"
        
        dicom_series = medipy.io.dicom.DicomSeries.read(path)
        self.assertTrue(dicom_series.has_description(description))
        self.assertFalse(dicom_series.has_description(""))
         
        image = medipy.io.load(dicom_series.url_from_description(description))
        self.assertEqual(image.metadata["series_instance_uid"], 
                         "1.3.46.670589.11.0.0.11.4.2.0.8743.5.5396.2006120114285654497")
        
        os.chdir(old_dir)
        self.cleanup(tempdir)
    
    def test_custom_name(self):
        path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), 
            "..", "..", "..", "data", "input", "dicom_series_brainix"))
        
        tempdir = self.extract("brainix.tgz")
        old_dir = os.getcwd()
        os.chdir(tempdir)
        
        shutil.copy(path, os.path.join(tempdir, "dicom_series"))
        path = os.path.join(tempdir, "dicom_series")
        
        custom_name = "T2 Flair"
        
        dicom_series = medipy.io.dicom.DicomSeries.read(path)
        self.assertTrue(dicom_series.has_custom_name(custom_name))
        self.assertFalse(dicom_series.has_custom_name("foo"))
         
        image = medipy.io.load(dicom_series.url_from_custom_name(custom_name))
        self.assertEqual(image.metadata["series_instance_uid"], 
                         "1.3.46.670589.11.0.0.11.4.2.0.8743.5.5396.2006120114285654497")
        
        os.chdir(old_dir)
        self.cleanup(tempdir)

if __name__ == "__main__" :
    unittest.main()
