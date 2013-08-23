import os
import shutil
import tarfile
import tempfile
import unittest

import medipy.io.dicom
import medipy.io.dicom.split

class TestSplit(unittest.TestCase) :
    def extract(self, filename):
        data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data","input"))
        
        tempdir = tempfile.mkdtemp()
        
        archive = tarfile.open(os.path.join(data_directory, filename), "r")
        archive.extractall(tempdir)
        
        return tempdir
    
    def cleanup(self, tempdir):
        shutil.rmtree(tempdir)
            
    def test_series(self) :
        tempdir = self.extract("brainix.tgz")
        dicomdir = medipy.io.dicom.read(os.path.join(tempdir, "BRAINIX", "DICOMDIR"))
        
        series = medipy.io.dicom.series([dicomdir])
        
        self.assertEqual(len(series), 7)
        
        self.cleanup(tempdir)

if __name__ == "__main__" :
    unittest.main()
