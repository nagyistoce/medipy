import math
import os
import unittest

import itk
import numpy.testing

import medipy.io.dicom
import medipy.io.dicom.csa2
import medipy.itk

class testMosaic(unittest.TestCase):
    def setUp(self):
        self.data_directory = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "data"))
        
        dataset = medipy.io.dicom.parse(
            os.path.join(self.data_directory, "input", "siemens_mosaic.dcm"))
        self.assertTrue(0x00291010 in dataset)
        image_csa = medipy.io.dicom.csa2.parse_csa(dataset[0x0029,0x1010])
        
        self.number_of_tiles = image_csa["NumberOfImagesInMosaic"][0]
        
        mosaic_size = math.ceil(math.sqrt(self.number_of_tiles))
        self.tile_size = (int(dataset.columns/mosaic_size),
                          int(dataset.rows/mosaic_size))
        
    def testAssembleTilesImageFilter(self):
        
        PixelType = itk.US
        TileImageType = itk.Image[PixelType, 2]
        VolumeImageType = itk.Image[PixelType, 3]
        
        reader = itk.ImageFileReader[TileImageType].New(
            FileName = os.path.join(self.data_directory, "input", "siemens_mosaic.dcm"))
        assemble_tiles = itk.AssembleTilesImageFilter[TileImageType, VolumeImageType].New(
            Input=reader[0], TileSize=self.tile_size, 
            NumberOfTiles=self.number_of_tiles, EmptyTiles=1, AssemblyOrder=0)
        result = medipy.itk.itk_image_to_array(assemble_tiles()[0], False)
        
        nifti_reader = itk.ImageFileReader[VolumeImageType].New(
            FileName = os.path.join(self.data_directory, "baseline", "io", 
                                    "dicom", "siemens_mosaic.nii"))
        expected = medipy.itk.itk_image_to_array(nifti_reader()[0], False)
        