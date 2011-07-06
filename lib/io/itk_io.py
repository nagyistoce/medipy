##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import os

import itk
import numpy

from medipy.io.io_base import IOBase
from medipy.itk import itk_image_to_array, medipy_image_to_itk_image

class ITK(IOBase) :
    
    _io_classes = [itk.BMPImageIO, itk.GDCMImageIO, itk.JPEGImageIO, 
                   itk.MetaImageIO, itk.NiftiImageIO, itk.NrrdImageIO, 
                   itk.PNGImageIO, itk.TIFFImageIO, itk.VTKImageIO]
    
    _component_type_to_type = {
        itk.ImageIOBase.UCHAR : itk.UC,
        itk.ImageIOBase.CHAR : itk.SC,
        itk.ImageIOBase.USHORT : itk.US,
        itk.ImageIOBase.SHORT : itk.SS,
        itk.ImageIOBase.UINT : itk.UI,
        itk.ImageIOBase.INT : itk.SI,
        itk.ImageIOBase.ULONG : itk.UL,
        itk.ImageIOBase.LONG : itk.SL,
        itk.ImageIOBase.FLOAT : itk.F,
        itk.ImageIOBase.DOUBLE : itk.D,
    }
    
    # Merge all supported read extensions, add a "*" before each of them       
    filenames = ["*"+str(x) for x in 
                 reduce(lambda l1, l2:list(l1)+list(l2), 
                        [x.New().GetSupportedReadExtensions() 
                         for x in _io_classes], 
                        [])]
    
    def __init__(self, filename=None, report_progress=None):
        
        self._loader = None
        self._image_informations_read = False
        self._saver = None
        
        IOBase.__init__(self, filename, report_progress)
        
        if filename is not None :
            self._set_filename(filename)
        
    def can_load(self):
        return self._loader is not None
    
    def number_of_images(self) :
        return 1
    
    def load_data(self, index=0) :
        if not self._image_informations_read :
            self._loader.SetFileName(self._filename)
            self._loader.ReadImageInformation()
            self._image_informations_read = True
        
        # TODO : NumberOfComponents != 1
        PixelType = ITK._component_type_to_type[self._loader.GetComponentType()]
        Dimension = self._loader.GetNumberOfDimensions()
        ImageType = itk.Image[PixelType, Dimension]
        
        reader = itk.ImageFileReader[ImageType].New()
        reader.SetImageIO(self._loader)
        reader.SetFileName(self._filename)
        reader.Update()
        
        array = itk_image_to_array(reader.GetOutput(), True)
        return array
    
    def load_metadata(self, index=0) :
        
        metadata = {}
        
        if not self._image_informations_read :
            self._loader.SetFileName(self._filename)
            self._loader.ReadImageInformation()
            self._image_informations_read = True
    
        ndim = self._loader.GetNumberOfDimensions()
        # GetDirection returns columns of the direction matrix 
        # (cf. itkImageFileReader.txx), so we need to transpose the numpy array.
        # We then need to change the axes order from ITK to numpy
        metadata["direction"] = numpy.asarray([self._loader.GetDirection(i) 
                                               for i in range(ndim)]).T
        metadata["direction"] = numpy.fliplr(numpy.flipud(metadata["direction"]))
        
        metadata["origin"] = numpy.asarray([self._loader.GetOrigin(ndim-i-1) 
                                            for i in range(ndim)])
        metadata["spacing"] = numpy.asarray([self._loader.GetSpacing(ndim-i-1) 
                                             for i in range(ndim)])
        
        # TODO : other metadata from dictionary
        
        return metadata
    
    def can_save(self):
        return self._saver is not None
    
    def save(self, image) :
        itk_image = medipy_image_to_itk_image(image, False)
        writer = itk.ImageFileWriter[itk_image].New(
            ImageIO = self._saver, FileName = self._filename, Input = itk_image)
        writer.Update()
    
    ##############
    # Properties #
    ##############
    
    def _get_filename(self) :
        return self._filename
    
    def _set_filename(self, filename):
        self._filename = str(filename)
        if os.path.isfile(self._filename) :
            self._loader = self._find_loader()
        self._saver = self._find_saver()
    
    filename = property(_get_filename, _set_filename)
    
    #####################
    # Private interface #
    #####################
    
    def _find_loader(self) :
        """ Return an instance of a subclass of itk.ImageIOBase that can read
            the current filename.
        """
        
        loader = None
        for load_class in self._io_classes :
            l = load_class.New()
            if l.CanReadFile(self._filename) :
                l.SetFileName(self._filename)
                l.ReadImageInformation()
                dimensions_ok = (l.GetNumberOfDimensions() in 
                                 [x[1] for x in itk.Image.__template__])
                if dimensions_ok :
                    self._image_informations_read = True
                    loader = l
                    break
        return loader
    
    def _find_saver(self) :
        """ Return an instance of a subclass of itk.ImageIOBase that can read
            the current filename.
        """
        
        saver = None
        for save_class in self._io_classes :
            s = save_class.New()
            if s.CanWriteFile(self._filename) :
                saver = s
                break
        return saver

def main():
    import sys
    
    filename = sys.argv[1]
    
    loader = ITK(filename)
    print loader.can_load()
    metadata = loader.load_metadata()
    data = loader.load_data()
    print metadata, data.shape, data.dtype

if __name__ == "__main__" :
    main()