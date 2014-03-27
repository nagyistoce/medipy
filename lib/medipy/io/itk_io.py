##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import logging
import os

import itk
import numpy

import medipy.io # add PyArrayIO to itk
import medipy.io.dicom
from medipy.io.io_base import IOBase
import medipy.itk
from medipy.itk import itk_image_to_array, medipy_image_to_itk_image, dtype_to_itk

class ITK(IOBase) :
    
    _io_classes = [itk.AnalyzeImageIO, itk.BMPImageIO, itk.JPEGImageIO, 
                   itk.MetaImageIO, itk.NiftiImageIO, itk.NrrdImageIO, 
                   itk.PNGImageIO, itk.TIFFImageIO, itk.VTKImageIO]
    
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

        # Loader informations        
        self._data_type = None
        self._image_type = None
        self._pixel_type_up_to_date = False
        self._image_informations_read = False
        
        IOBase.__init__(self, filename, report_progress)
        
        if filename is not None :
            self._set_filename(filename)
        
    def __deepcopy__(self, memo):
        # deepcopy does not like WrapITK objects
        new_object = ITK(self._filename, self._report_progress)
        return new_object
    
    def can_load(self):
        return self._loader is not None
    
    def number_of_images(self) :
        result = 1
        if self._loader.GetNameOfClass() == "NiftiImageIO":
            self._update_pixel_type()
            if self._loader.GetNumberOfDimensions() == 4 and self._data_type == "scalar":
                result = self._loader.GetDimensions(3)
        return result
    
    def _update_pixel_type(self):
        if self._pixel_type_up_to_date :
            # Nothing to do
            return
        
        # Update the loader
        if not self._image_informations_read :
            self._loader.SetFileName(self._filename)
            self._loader.ReadImageInformation()
            self._image_informations_read = True
        
        # Get the component type
        ComponentType = medipy.itk.types.io_component_type_to_type[
            self._loader.GetComponentType()]
        # Get the pixel type
        pixel_type = self._loader.GetPixelTypeAsString(self._loader.GetPixelType())
        if pixel_type == "scalar" :
            self._data_type = "scalar"
            self._image_type = "normal"
        elif pixel_type == "rgb" :
            self._data_type = "vector"
            self._image_type = "rgb"
        else :
            raise NotImplementedError(pixel_type)
        
        self._pixel_type_up_to_date = True
    
    def load_data(self, index=0) :
        self._update_pixel_type()
        
        InstantiatedTypes = set([itk.template(x[0])[1][0] 
                                 for x in itk.NumpyBridge.__template__.keys()])
        OriginalPixelType = medipy.itk.types.io_component_type_to_type[self._loader.GetComponentType()]
        PixelType = OriginalPixelType
        try_smaller_types = False
        while PixelType not in InstantiatedTypes :
            if PixelType not in medipy.itk.types.larger_type :
                # No instantiated pixel type larger that the pixel type in the
                # file ; try smaller types
                try_smaller_types = True
                break
            PixelType = medipy.itk.types.larger_type[PixelType]
        if try_smaller_types :
            PixelType = OriginalPixelType
            while PixelType not in InstantiatedTypes :
                if PixelType not in medipy.itk.types.smaller_type :
                    raise medipy.base.Exception(
                        "Cannot find an instantiated pixel type for {0}".format(
                            OriginalPixelType))
                PixelType = medipy.itk.types.smaller_type[PixelType]
            logging.warn("No instantiated type larger than {0}, using {1}".format(
                OriginalPixelType, PixelType))

        Dimension = self._loader.GetNumberOfDimensions()

        if self._filter == itk.ImageFileReader :
            ImageType = itk.Image[PixelType, Dimension]
            reader = itk.ImageFileReader[ImageType].New()
        elif self._filter == itk.PyArrayFileReader :
            reader = itk.PyArrayFileReader[PixelType, Dimension].New()
        
        # Re-using the same ImageIO causes a memory leak !
        loader_clone = itk.down_cast(self._loader.CreateAnother())
        reader.SetImageIO(loader_clone)
        reader.SetFileName(self._filename)
        reader.Update()
        
        if self._filter == itk.ImageFileReader :
            array = itk_image_to_array(reader.GetOutput(), True)
        else :
            array = reader.GetArray()
        
        if self._loader.GetNameOfClass() == "NiftiImageIO":
            if array.ndim == 4 and self._data_type == "scalar":
                # Make sure to /NOT/ use a view, which uses too much memory
                array = array[index].copy()
        
        return array
    
    def load_metadata(self, index=0) :
        
        self._update_pixel_type()
        
        metadata = {}

        metadata["data_type"] = self._data_type
        metadata["image_type"] = self._image_type
    
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
        
        if self._loader.GetNameOfClass() == "NiftiImageIO":
            if self._loader.GetNumberOfDimensions() == 4 and self._data_type == "scalar":
                metadata["direction"] = metadata["direction"][1:,1:]
                metadata["origin"] = metadata["origin"][1:]
                metadata["spacing"] = metadata["spacing"][1:]
        
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
                # Specify ndmin, since files with only 1 value will return
                # [a,b,c] instead of [[a,b,c]]. Same remark applies for bvalues.
                gradients = numpy.loadtxt(gradient_file, dtype=numpy.single)
                gradients = gradients.T
                
                if gradients.ndim == 1 :
                    gradients = numpy.asarray([gradients])
                bvalues = numpy.loadtxt(bvalue_file, dtype=numpy.single)
                if bvalues.ndim == 0 :
                    bvalues = numpy.asarray([bvalues])
                
                mr_diffusion_sequence = []
                for index, gradient in enumerate(gradients) :
                    dataset = medipy.io.dicom.DataSet()
                    dataset.diffusion_directionality = medipy.io.dicom.CS("DIRECTIONAL")
                    
                    dataset.diffusion_bvalue = medipy.io.dicom.FD(bvalues[index])
                    
                    gradient_dataset = medipy.io.dicom.DataSet()
                    gradient_dataset.diffusion_gradient_orientation = medipy.io.dicom.FD(gradient)
                    dataset.diffusion_gradient_direction_sequence = medipy.io.dicom.SQ([gradient_dataset])
                    
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
        try_smaller_types = False
        while PixelType not in InstantiatedTypes :
            if PixelType not in medipy.itk.types.larger_type :
                # No instantiated pixel type larger that the pixel type in the
                # file ; try smaller types
                try_smaller_types = True
                break
            PixelType = medipy.itk.types.larger_type[PixelType]
        if try_smaller_types :
            PixelType = dtype_to_itk[image.dtype.type]
            while PixelType not in InstantiatedTypes :
                if PixelType not in medipy.itk.types.smaller_type :
                    raise medipy.base.Exception(
                        "Cannot find an instantiated pixel type for {0}".format(
                            dtype_to_itk[image.dtype.type]))
                PixelType = medipy.itk.types.smaller_type[PixelType]
            logging.warn("No instantiated type larger than {0}, using {1}".format(
                dtype_to_itk[image.dtype.type], PixelType))
        
        Dimension = image.ndim
        
        return ((PixelType, Dimension) in itk.Image.__template__ or
                (PixelType, Dimension) in itk.PyArrayFileWriter.__template__)
    
    def save(self, image) :
        # Check the instantiations
        InstantiatedTypes = set([itk.template(x[0])[1][0] 
                                 for x in itk.NumpyBridge.__template__.keys()])
        PixelType = dtype_to_itk[image.dtype.type]
        try_smaller_types = False
        while PixelType not in InstantiatedTypes :
            if PixelType not in medipy.itk.types.larger_type :
                # No instantiated pixel type larger that the pixel type in the
                # file ; try smaller types
                try_smaller_types = True
                break
            PixelType = medipy.itk.types.larger_type[PixelType]
        if try_smaller_types :
            PixelType = dtype_to_itk[image.dtype.type]
            while PixelType not in InstantiatedTypes :
                if PixelType not in medipy.itk.types.smaller_type :
                    raise medipy.base.Exception(
                        "Cannot find an instantiated pixel type for {0}".format(
                            dtype_to_itk[image.dtype.type]))
                PixelType = medipy.itk.types.smaller_type[PixelType]
            logging.warn("No instantiated type larger than {0}, using {1}".format(
                dtype_to_itk[image.dtype.type], PixelType))
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
                # Make sure we have all the required elements
                if "diffusion_gradient_direction_sequence" not in diffusion :
                    continue
                elif not diffusion.diffusion_gradient_direction_sequence.value :
                    continue
                elif "diffusion_gradient_orientation" not in diffusion.diffusion_gradient_direction_sequence.value[0] :
                    continue
                elif "diffusion_bvalue" not in diffusion :
                    continue
                
                gradient = diffusion.diffusion_gradient_direction_sequence.value[0].diffusion_gradient_orientation
                b_value = diffusion.diffusion_bvalue.value
                
                for index, value in enumerate(gradient.value) :
                    gradients[index].append(str(value))
                b_values.append(str(b_value))
            
            if gradients[0] and b_values :
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
            self._pixel_type_up_to_date = False
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
                
                filter = None
                
                if l.GetPixelTypeAsString(l.GetPixelType()) in ["vector", "rgb"] :
                    if 1+l.GetNumberOfDimensions() in [x[1] for x in itk.PyArrayFileReader] :
                        filter = itk.PyArrayFileReader
                elif l.GetPixelTypeAsString(l.GetPixelType()) == "scalar" :
                    if l.GetNumberOfDimensions() in [x[1] for x in itk.Image] :
                        filter = itk.ImageFileReader
                    elif l.GetNumberOfDimensions() in [x[1] for x in itk.PyArrayFileReader] :
                        filter = itk.PyArrayFileReader
                    # Otherwise no filter can load this data
                    
                if filter is not None :
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
