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

import medipy.io # add PyArrayIO to itk
from medipy.io.io_base import IOBase
from medipy.itk import itk_image_to_array, medipy_image_to_itk_image, dtype_to_itk

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
        self._saver = None
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
        
        # Load gradient direction file if loading NIfTI
        if isinstance(self._loader, itk.NiftiImageIO) :
            base_name = os.path.splitext(self._filename)[0]
            if base_name.endswith(".nii") :
                base_name = os.path.splitext(base_name)[0]
            gradient_candidates = [
                base_name+".bvecs", # /foo/bar/image.bvecs
                base_name+".bvec", # /foo/bar/image.bvec
                os.path.join(os.path.dirname(self._filename), "bvecs"), # /foo/bar/bvecs
                os.path.join(os.path.dirname(self._filename), "bvec") # /foo/bar/bvec
            ]
            gradient_file = None
            for candidate in gradient_candidates :
                if os.path.isfile(candidate) :
                    gradient_file = candidate
                    break
            
            # Load b-values file
            bvalue_candidates = [
                base_name+".bvals", # /foo/bar/image.bvals
                base_name+".bval", # /foo/bar/image.bval
                os.path.join(os.path.dirname(self._filename), "bval"), # /foo/bar/bvals
                os.path.join(os.path.dirname(self._filename), "bvals") # /foo/bar/bval
            ]
            bvalue_file = None
            for candidate in bvalue_candidates :
                if os.path.isfile(candidate) :
                    bvalue_file = candidate
                    break
            
            if None not in [gradient_file, bvalue_file] :
                gradients = numpy.loadtxt(gradient_file, dtype=numpy.single)
                bvalues = numpy.loadtxt(bvalue_file, dtype=numpy.single)
                
                gradients = gradients.T
                
                mr_diffusion_sequence = []
                for index, gradient in enumerate(gradients) :
                    dataset = medipy.io.dicom.DataSet()
                    dataset.diffusion_directionality = "DIRECTIONAL"
                    
                    dataset.diffusion_bvalue = bvalues[index]
                    
                    gradient_dataset = medipy.io.dicom.DataSet()
                    gradient_dataset.diffusion_gradient_orientation = gradient
                    dataset.diffusion_gradient_direction_sequence = [gradient_dataset]
                    
                    mr_diffusion_sequence.append(dataset)
                
                metadata["mr_diffusion_sequence"] = mr_diffusion_sequence
        
        return metadata
    
    def can_save(self, image):
        if not self._saver :
            # No ImageIO is available for the given filename
            return False
        
        # Check the instantiations
        InstantiatedTypes = set([itk.template(x[0])[1][0] 
                                 for x in itk.NumpyBridge.__template__.keys()])
        PixelType = dtype_to_itk[image.dtype.type]
        while PixelType not in InstantiatedTypes :
            PixelType = ITK._larger_type[PixelType]
        Dimension = image.ndim
        
        return ((PixelType, Dimension) in itk.Image.__template__ or
                (PixelType, Dimension) in itk.PyArrayFileWriter.__template__)
    
    def save(self, image) :
        InstantiatedTypes = set([itk.template(x[0])[1][0] 
                                 for x in itk.NumpyBridge.__template__.keys()])
        PixelType = dtype_to_itk[image.dtype.type]
        while PixelType not in InstantiatedTypes :
            PixelType = ITK._larger_type[PixelType]
        Dimension = image.ndim
        
        if (PixelType, Dimension) in itk.Image.__template__ :
            itk_image = medipy_image_to_itk_image(image, False)
            writer = itk.ImageFileWriter[itk_image].New(
                ImageIO = self._saver, Input = itk_image)
        else :
            writer = itk.PyArrayFileWriter[PixelType, Dimension].New()
            writer.SetArray(image.data)
            writer.SetOrigin(image.origin)
            writer.SetSpacing(image.spacing)
            writer.SetDirection(image.direction)
        
        writer.SetFileName(self._filename)
        writer.Update()
        
        # Save gradient direction file if saving NIfTI
        if isinstance(self._saver, itk.NiftiImageIO) and "mr_diffusion_sequence" in image.metadata :
            gradients = [[], [], []]
            b_values = []
            for diffusion in image.metadata["mr_diffusion_sequence"] :
                gradient = diffusion.diffusion_gradient_direction_sequence[0].diffusion_gradient_orientation
                b_value = diffusion.diffusion_bvalue
                
                for index, value in enumerate(gradient) :
                    gradients[index].append(str(value))
                b_values.append(str(b_value))
            
            gradients = "\n".join([" ".join(direction) for direction in gradients])
            b_values = " ".join(b_values)
            
            base_name = os.path.splitext(self._filename)[0]
            if base_name.endswith(".nii") :
                base_name = os.path.splitext(base_name)[0]
            
            gradients_file = open("{0}.bvec".format(base_name), "w")
            gradients_file.write(gradients)
            gradients_file.close()
            
            bvalues_file = open("{0}.bval".format(base_name), "w")
            bvalues_file.write(b_values)
            bvalues_file.close()
    
    ##############
    # Properties #
    ##############
    
    def _set_filename(self, filename):
        self._filename = str(filename)
        if os.path.isfile(self._filename) :
            self._filter, self._loader = self._find_loader()
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
