import os
import shutil
import tempfile
import unittest

import numpy

import medipy.base
import medipy.io
import medipy.diffusion

class TestIO(unittest.TestCase):
    def test_itk_io(self):
        tensors = medipy.base.Image(
            data=numpy.arange(10*20*30*6, dtype=numpy.single).reshape((10,20,30,6)),
            dti="tensor_2",
            origin=(1,2,3), spacing=(4,5,6), direction=medipy.base.coordinate_system.RAS)
        
        directory = tempfile.mkdtemp()
        medipy.io.save(tensors, os.path.join(directory, "tensors.nii"))
        
        other_tensors = medipy.io.load(
            os.path.join(directory, "tensors.nii"), None)
        
        self.assertTrue(isinstance(
            other_tensors.metadata["loader"]["loader"], 
            medipy.diffusion.itk_io.ITK))
        
        self.assertEqual(other_tensors.shape, (10,20,30))
        self.assertEqual(other_tensors.dtype, numpy.single)
        self.assertEqual(other_tensors.data_type, "vector")
        self.assertEqual(other_tensors.image_type, "tensor_2")
        
        numpy.testing.assert_array_almost_equal(other_tensors.origin, (1,2,3))
        numpy.testing.assert_array_almost_equal(other_tensors.spacing, (4,5,6))
        numpy.testing.assert_array_almost_equal(
            other_tensors.direction, medipy.base.coordinate_system.RAS)
        
        numpy.testing.assert_array_almost_equal(
            other_tensors.data, 
            numpy.arange(10*20*30*6, dtype=numpy.single).reshape((10,20,30,6)))

        shutil.rmtree(directory)
    
    def test_vtk_io(self):
        tensors = medipy.base.Image(
            data=numpy.arange(10*20*30*6, dtype=numpy.single).reshape((10,20,30,6)),
            dti="tensor_2",
            origin=(1,2,3), spacing=(4,5,6), direction=medipy.base.coordinate_system.LPS)
        
        directory = tempfile.mkdtemp()
        medipy.io.save(tensors, os.path.join(directory, "tensors.vtk"))
        
        other_tensors = medipy.io.load(
            os.path.join(directory, "tensors.vtk"), None)
        
        self.assertTrue(isinstance(
            other_tensors.metadata["loader"]["loader"], 
            medipy.diffusion.vtk_io.Vtk))
        
        self.assertEqual(other_tensors.shape, (10,20,30))
        self.assertEqual(other_tensors.dtype, numpy.single)
        self.assertEqual(other_tensors.data_type, "vector")
        self.assertEqual(other_tensors.image_type, "tensor_2")
        
        numpy.testing.assert_array_almost_equal(other_tensors.origin, (1,2,3))
        numpy.testing.assert_array_almost_equal(other_tensors.spacing, (4,5,6))
        numpy.testing.assert_array_almost_equal(
            other_tensors.direction, medipy.base.coordinate_system.LPS)
        
        numpy.testing.assert_array_almost_equal(
            other_tensors.data, 
            numpy.arange(10*20*30*6, dtype=numpy.single).reshape((10,20,30,6)))

        shutil.rmtree(directory)

if __name__ == "__main__" :
    unittest.main()
