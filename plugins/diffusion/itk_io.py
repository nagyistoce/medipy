##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import os

import itk

import numpy

import medipy.itk
import medipy.io

class Tensor2IO(medipy.io.IOBase):
    _io_classes = [itk.BMPImageIO, itk.JPEGImageIO, itk.MetaImageIO, 
                   itk.NiftiImageIO, itk.NrrdImageIO, itk.PNGImageIO, 
                   itk.TIFFImageIO, itk.VTKImageIO]
    
    def __init__(self, filename, *args, **kwargs) :
        medipy.io.IOBase.__init__(self, filename, *args, **kwargs)
        self._loader = None
        self._saver = None
        self._image_informations_read = False
        
        if filename is not None :
            self._set_filename(filename)
    
    def can_load(self):
        elements = self._filename.split(".")
        extension = elements[-1]
        if extension == "gz" :
            extension = ".".join(elements[-2:])
        
        # TODO : vtk
        if extension not in ["nii", "nii.gz"] :
            return False
        
        if extension in ["nii", "nii.gz"] :
            class_ = itk.NiftiImageIO
        
        itk_io = class_.New()
        if not itk_io.CanReadFile(self._filename) :
            return False
        
        itk_io.SetFileName(self._filename)
        itk_io.ReadImageInformation()
        
        pixel_type = itk_io.GetPixelTypeAsString(itk_io.GetPixelType())
        
        if pixel_type != "symmetric_second_rank_tensor" :
            return False
        
        return self._loader is not None
    
    def number_of_images(self):
        return 1
    
    def load_data(self, index=0) :
        if not self._image_informations_read :
            self._loader.SetFileName(self._filename)
            self._loader.ReadImageInformation()
            self._image_informations_read = True
        
        BridgedTypes = set([itk.template(x[0])[1][0] for x in itk.NumpyBridge])
        PixelType = medipy.itk.types.io_component_type_to_type[self._loader.GetComponentType()]
        while PixelType not in BridgedTypes :
            PixelType = medipy.itk.types.larger_type[PixelType]
        Dimension = self._loader.GetNumberOfDimensions()
        
        VectorImageType = itk.VectorImage[PixelType, Dimension]
        
        reader = itk.Tensor2ImageFileReader[VectorImageType].New()
        
        reader.SetImageIO(self._loader)
        reader.SetFileName(self._filename)
        reader.Update()
        
        array = medipy.itk.itk_image_to_array(reader.GetOutput(), True)
        return array
    
    def load_metadata(self, index=0) :
        if not self._image_informations_read :
            self._loader.SetFileName(self._filename)
            self._loader.ReadImageInformation()
            self._image_informations_read = True
            
        ndim = self._loader.GetNumberOfDimensions()
        
        metadata = {
            "data_type" : "vector",
            "image_type" : "tensor_2"
        }
        
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
        
        return metadata
    
    def can_save(self, image):
        if image.image_type != "tensor_2" :
            return False
        
        BridgedTypes = set([itk.template(x[0])[1][0] for x in itk.NumpyBridge])
        PixelType = medipy.itk.dtype_to_itk[image.dtype.type]
        while PixelType not in BridgedTypes :
            PixelType = medipy.itk.types.larger_type[PixelType]
        Dimension = image.ndim
        
        VectorImageType = itk.VectorImage[PixelType, Dimension]
        
        return VectorImageType in itk.Tensor2ImageFileWriter
    
    def save(self, image) :
        BridgedTypes = set([itk.template(x[0])[1][0] for x in itk.NumpyBridge])
        PixelType = medipy.itk.dtype_to_itk[image.dtype.type]
        while PixelType not in BridgedTypes :
            PixelType = medipy.itk.types.larger_type[PixelType]
        Dimension = image.ndim
        
        VectorImageType = itk.VectorImage[PixelType, Dimension]
        
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        
        writer = itk.Tensor2ImageFileWriter[VectorImageType].New(
            ImageIO = self._saver, Input = itk_image)
        
        writer.SetFileName(self._filename)
        writer.Update()
        
    ##############
    # Properties #
    ##############
    
    def _set_filename(self, filename):
        self._filter = None
        self._image_informations_read = False
        
        self._filename = str(filename)
        if os.path.isfile(self._filename) :
            self._loader = self._find_loader()
        self._saver = self._find_saver()
    
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
                
                instantiated = (l.GetNumberOfDimensions() in 
                                [itk.template(x[0])[1][1] 
                                 for x in itk.Tensor2ImageFileReader.__template__])
                
                if instantiated :
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
