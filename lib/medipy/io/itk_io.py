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

import medipy.io # add PyArrayIO to itk
from medipy.io.io_base import IOBase
from medipy.itk import itk_image_to_array, medipy_image_to_itk_image

class ITK(IOBase) :
    
    _io_classes = [itk.BMPImageIO, itk.JPEGImageIO, itk.MetaImageIO, 
                   itk.NiftiImageIO, itk.NrrdImageIO, itk.PNGImageIO, 
                   itk.TIFFImageIO, itk.VTKImageIO]
    
    _larger_type = {
        itk.UC : itk.US,
        itk.US : itk.UL,
        itk.UI : itk.UL,
        itk.SC : itk.SS,
        itk.SS : itk.SL,
        itk.SI : itk.SL,
        itk.F : itk.D
    }
    
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
        
        self._filter = None
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
        InstantiatedTypes = set([itk.template(x[0])[1][0] 
                                 for x in itk.NumpyBridge.__template__.keys()])
        PixelType = ITK._component_type_to_type[self._loader.GetComponentType()]
        while PixelType not in InstantiatedTypes :
            PixelType = ITK._larger_type[PixelType]
        Dimension = self._loader.GetNumberOfDimensions()
        
        if self._filter == itk.ImageFileReader :
            ImageType = itk.Image[PixelType, Dimension]
            reader = itk.ImageFileReader[ImageType].New()
        elif self._filter == itk.PyArrayFileReader :
            reader = itk.PyArrayFileReader[PixelType, Dimension].New()
        
        reader.SetImageIO(self._loader)
        reader.SetFileName(self._filename)
        reader.Update()
        
        if self._filter == itk.ImageFileReader :
            array = itk_image_to_array(reader.GetOutput(), True)
        else :
            array = reader.GetArray()
        
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
    
    def can_save(self, image):
        if not self._saver :
            # No ImageIO is available for the given filename
            return False
        
        # Check the instantiations
        return (image.ndim in [x[1] for x in itk.Image.__template__])
    
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
            self._filter, self._loader = self._find_loader()
        self._saver = self._find_saver()
    
    filename = property(_get_filename, _set_filename)
    
    #####################
    # Private interface #
    #####################
    
    def _find_loader(self) :
        """ Return an instance of a subclass of itk.ImageIOBase that can read
            the current filename.
        """
        
        filter = None
        loader = None
        for load_class in self._io_classes :
            l = load_class.New()
            if l.CanReadFile(self._filename) :
                l.SetFileName(self._filename)
                l.ReadImageInformation()
                
                # Check if we can use itk::ImageFileReader, i.e. check if 
                # itk.Image is wrapped for the dimension as stored in the file
                use_image_file_reader = (l.GetNumberOfDimensions() in 
                                        [x[1] for x in itk.Image.__template__])
                # Check if we can use PyArrayImageReader
                use_py_array_file_reader = (l.GetNumberOfDimensions() in 
                                            [x[1] for x in itk.PyArrayFileReader.__template__])
                
                if use_image_file_reader :
                    filter = itk.ImageFileReader
                elif use_py_array_file_reader :
                    filter = itk.PyArrayFileReader
                    
                if use_image_file_reader or use_py_array_file_reader :
                    self._image_informations_read = True
                    loader = l
                    break
        return filter, loader
    
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
