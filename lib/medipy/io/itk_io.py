##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itertools
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
    
    _io_classes = [itk.BMPImageIO, itk.JPEGImageIO, 
                   itk.MetaImageIO, itk.NiftiImageIO, itk.NrrdImageIO, 
                   itk.PNGImageIO, itk.TIFFImageIO, itk.VTKImageIO]
    _maxdim = 4
    
    # Merge all supported read extensions, add a "*" before each of them       
    filenames = ["*"+str(x) for x in 
                 reduce(lambda l1, l2:list(l1)+list(l2), 
                        [x.New().GetSupportedReadExtensions() 
                         for x in _io_classes], 
                        [])]
    
    def __init__(self, filename=None, report_progress=None):
        IOBase.__init__(self, filename, report_progress)
        
        if filename is not None :
            self._set_filename(filename)
        
    def __deepcopy__(self, memo):
        # deepcopy does not like WrapITK objects
        new_object = ITK(self._filename, self._report_progress)
        return new_object
    
    def can_load(self):
        return all(self._find_loader())
    
    def number_of_images(self) :
        _, image_io = self._find_loader()
        
        result = 1
        if image_io.GetNameOfClass() == "NiftiImageIO":
            pixel_type = image_io.GetPixelTypeAsString(image_io.GetPixelType())
            if image_io.GetNumberOfDimensions() == 4 and pixel_type == "scalar":
                result = image_io.GetDimensions(3)
        return result
    
    def load_data(self, index=0) :
        ReaderClass, image_io = self._find_loader()
        
        reader = ReaderClass.New()
        
        # Re-using the same ImageIO causes a memory leak !
        image_io_clone = itk.down_cast(image_io.CreateAnother())
        reader.SetImageIO(image_io_clone)
        reader.SetFileName(self._filename)
        reader.Update()
        
        if reader.GetNameOfClass() == "ImageFileReader":
            array = itk_image_to_array(reader.GetOutput(), True)
        else :
            array = reader.GetArray()
        
        if image_io.GetNameOfClass() == "NiftiImageIO":
            data_type = image_io.GetPixelTypeAsString(image_io.GetPixelType())
            if array.ndim == 4 and data_type == "scalar":
                # Make sure to /NOT/ use a view, which uses too much memory
                array = array[index].copy()
        
        return array
    
    def load_metadata(self, index=0) :
        
        _, image_io = self._find_loader()
        
        pixel_type = image_io.GetPixelTypeAsString(image_io.GetPixelType())
        if pixel_type == "scalar" :
            data_type = "scalar"
            image_type = "normal"
        elif pixel_type == "rgb" :
            data_type = "vector"
            image_type = "rgb"
        else :
            raise NotImplementedError(pixel_type)
        
        metadata = {}

        metadata["data_type"] = data_type
        metadata["image_type"] = image_type
    
        ndim = image_io.GetNumberOfDimensions()
        # GetDirection returns columns of the direction matrix 
        # (cf. itkImageFileReader.txx), so we need to transpose the numpy array.
        # We then need to change the axes order from ITK to numpy
        metadata["direction"] = numpy.asarray([image_io.GetDirection(i) 
                                               for i in range(ndim)]).T
        metadata["direction"] = numpy.fliplr(numpy.flipud(metadata["direction"]))
        
        metadata["origin"] = numpy.asarray([image_io.GetOrigin(ndim-i-1) 
                                            for i in range(ndim)])
        metadata["spacing"] = numpy.asarray([image_io.GetSpacing(ndim-i-1) 
                                             for i in range(ndim)])
        
        if image_io.GetNameOfClass() == "NiftiImageIO":
            if image_io.GetNumberOfDimensions() == 4 and data_type == "scalar":
                metadata["direction"] = metadata["direction"][1:,1:]
                metadata["origin"] = metadata["origin"][1:]
                metadata["spacing"] = metadata["spacing"][1:]
        
        # TODO : other metadata from dictionary
        
        # Load gradient direction file if loading NIfTI
        if image_io.GetNameOfClass() == "NiftiImageIO":
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
                
                mr_diffusion = medipy.io.dicom.DataSet(
                    diffusion_directionality="DIRECTIONAL",
                    diffusion_bvalue=bvalues[index],
                    diffusion_gradient_direction_sequence=[
                        medipy.io.dicom.DataSet(
                            diffusion_gradient_orientation=gradients[index])]
                )
                
                metadata["mr_diffusion_sequence"] = [mr_diffusion]
        
        return metadata
    
    def can_save(self, image):
        return all(self._find_saver(image))
    
    def save(self, image) :
        WriterClass, image_io = self._find_saver(image)
        
        writer = WriterClass.New(FileName=self._filename)
        if writer.GetNameOfClass() == "ImageFileWriter":
            itk_image = medipy_image_to_itk_image(image, False)
            writer.SetImageIO(image_io)
            writer.SetInput(itk_image)
        elif writer.GetNameOfClass() == "PyArrayFileWriter":
            writer.SetArray(image.data)
            writer.SetOrigin(image.origin)
            writer.SetSpacing(image.spacing)
            writer.SetDirection(image.direction)
        else:
            raise medipy.base.Exception(
                "Unknown writer class: {}".format(writer.GetNameOfClass()))
        
        writer.Update()
        
        # Save gradient direction file if saving NIfTI
        if (image_io.GetNameOfClass() == "NiftiImageIO" and 
                "mr_diffusion_sequence" in image.metadata):
            gradients = [[], [], []]
            b_values = []
            for diffusion in image.metadata["mr_diffusion_sequence"] :
                # Make sure we have all the required elements
                try:
                    gradient = diffusion.diffusion_gradient_direction_sequence.value[0].diffusion_gradient_orientation
                    b_value = diffusion.diffusion_bvalue.value
                except:
                    continue
                
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
    
    #####################
    # Private interface #
    #####################

    def _find_image_io(self, mode):
        """ Return an instance of a subclass of itk.ImageIOBase that can read
            or write the current filename. The mode parameters decides whether
            reading ("Read") or writing ("Write") is required. If no ImageIO can
            be found, an exception is raised.
        """
        
        image_io = None
        
        for image_io_class in self._io_classes:
            candidate = image_io_class.New()
            
            can_read = False
            try:
                function = getattr(candidate, "Can{}File".format(mode))
                is_able = function(self._filename)
            except:
                # Do nothing, just avoid propagating the error
                pass
            
            if is_able:
                image_io = candidate
                break
        
        if image_io is None:
            raise medipy.base.Exception(
                "No ImageIO can {} file {}".format(mode, self._filename))
        
        return image_io

    def _find_filter(self, mode, data_type, component_type, dimensions_count):
        if data_type == "scalar":
            # Nothing to adjust
            pass
        elif data_type in ["vector", "rgb"]:
            dimensions_count += 1
        else:
            raise medipy.base.Exception(
                "Unknown pixel type {} in file {}".format(
                    pixel_type, self._filename))
        
        larger_types = medipy.itk.types.get_larger_types(component_type)
        
        filter_ = None
        
        dimensions = range(dimensions_count, 1+self._maxdim)
        
        if mode == "Read":
            ImageFilter = itk.ImageFileReader
            PyArrayFilter = itk.PyArrayFileReader
        elif mode == "Write":
            ImageFilter = itk.ImageFileWriter
            PyArrayFilter = itk.PyArrayFileWriter
        else:
            raise medipy.base.Exception("Unknown mode: {}".format(mode))
        
        filter_ = None
        for p in itertools.product([component_type]+larger_types, dimensions):
            filter_ = (
                ImageFilter.get(itk.Image.get(p)) or PyArrayFilter.get(p))
            if filter_:
                break
        
        if not filter_:
            # Try smaller types, with a warning
            smaller_types = medipy.itk.types.get_smaller_types(component_type)
            for p in itertools.product(smaller_types, dimensions):
                filter_ = (
                    ImageFilter.get(itk.Image.get(p)) or PyArrayFilter.get(p))
                if filter_:
                    logging.warn(
                        "No instantiated type larger than {}, using {}".format(
                            component_type.name, p[0].name))
                    break
        
        if not filter_:
            raise medipy.base.Exception(
                "No {}er for ({},{})".format(
                    mode, component_type.name, dimensions_count))
        
        return filter_
    
    def _find_loader(self) :
        """ Return an ITK filter class and an instance of a subclass of 
            itk.ImageIOBase that can read the current filename.
        """
        
        image_io = self._find_image_io("Read")
        image_io.SetFileName(self._filename)
        image_io.ReadImageInformation()
        
        reader = self._find_filter(
            "Read", image_io.GetPixelTypeAsString(image_io.GetPixelType()), 
            medipy.itk.types.io_component_type_to_type[
                image_io.GetComponentType()], image_io.GetNumberOfDimensions())
        
        return reader, image_io
    
    def _find_saver(self, image) :
        """ Return an instance of a subclass of itk.ImageIOBase that can write
            the current filename.
        """
        
        image_io = self._find_image_io("Write")
        writer = self._find_filter("Write", 
            image.data_type, medipy.itk.dtype_to_itk[image.dtype.type], 
            image.data.ndim)
        return writer, image_io
