##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import gzip
import logging
import math
import os
import re
import string
from struct import pack

import itk
import numpy

from medipy.base import ImageAnnotation, coordinate_system
from medipy.io.io_base import IOBase
import medipy.itk

encoding_char = "+"
encoding_max_length = 2
encoding_basis = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_\\-:\\/\\?\\.\\$"

class IPBColorMaps(object) :
    # Regular colormaps
    red = 1
    green = 2
    blue = 4
    
    yellow = 3
    purple = 5
    cyan = 6
    
    gray = 7
    gray_invert = 11
    gray_invert2 = 12
    
    rainbow = 8
    rainbow2 = 9
    hot_iron = 16
    cool = 20
    copper = 21
    
    # Stage colormaps
    fmri = 38
    epi = 39
    epi2 = 40
    sep = 41
    stats=42
    freehand = 44
    label128 = 45
    fmri_retard = 46
    p_value = 47
    z_score = 48
    stats2 = 49
    
    def get_name(cls, index):
        attributes = [ attr for attr in dir(cls) 
            if not (attr.startswith("_") or attr == "get_name") ]
        result = None
        for attribute in attributes :
            if getattr(cls, attribute) == index :
                result = attribute
        return result 
    get_name = classmethod(get_name)

def decode(s):
    """ Decode a string stored in IPB-encoded format. 
        
        This was recoded in Python from Medimax, cf. the function 
        IPBHeaderUnquote in noyau/annotations/IPBFileEntryList.cpp
    """
     
    pattern = "(\\%s[%s]{%d})"%(encoding_char, encoding_basis, encoding_max_length)
    pattern = re.compile(pattern)
    components = re.split(pattern, s)
    
    result = ""
    for component in components : 
        if component[0] == "+" :
            value = 0
            basis = 1
            
            for c in component[1:] :
                value += (basis*encoding_basis.index(c))
                basis *= len(encoding_basis)
            
            result += chr(value)
        else : 
            result += component
    
    return result 

class IPB(IOBase) :
    """I/O class for IPB format"""
    
    _magic = "IPB3"
    
    filenames = ["*.ipb"]
    
    _ipb_to_metadata = {
        "patient name" : "patients_name",
        "date of birth" : "patients_birth_date",
        "date of exam" : "series_date",
        "doctor name" : "referring_physicians_name",
    }
    
    def __init__(self, *args, **kwargs) :
        IOBase.__init__(self, *args, **kwargs)
        
        self.header = None
        self._header_filename = None
        self._data_filename = None
        self._is_ipb = False
        self._is_gz = False
        
        if self._filename is not None :
            self._set_filename(self._filename)
    
    def _set_filename(self, filename) :
        self._filename = filename
        self._header_filename, self._data_filename = \
            self._get_filenames(filename)
        if self._header_filename is not None :
            try :
                self.header = self._header_as_dict()
            except Exception, e:
                logging.debug("Exception caught while converting IPB header %s"%e)
                pass
    
    filename = property(lambda self:self._filename, _set_filename)
        
    def can_load(self) :
        return self._is_ipb
        
    def number_of_images(self) :
        return self.header["number of images"]
    
    def load_data(self, index=0) :
        buffer = ""
        if self._is_gz :
            file = gzip.open(self._data_filename, "rb")
        else :
            file = open(self._data_filename, "rb")
        
        if self._report_progress is not None :
            self._report_progress(0.)
            
        file_size = os.path.getsize(self._data_filename)
        nb_chunks = 10
        chunk_size = file_size/nb_chunks
        done = False
        chunk_index = 0
        while not done :
            chunk = file.read(chunk_size)
            
            if self._report_progress is not None :
                self._report_progress(float(chunk_index)/float(nb_chunks+1))
                
            if len(chunk) == 0 :
                done = True
            
            chunk_index += 1 
            buffer += chunk

        offset = self.header["offset in raw"][index]
        
        # Fetch dimensions and compute number of items
        dim = [ self.header[key][index] for key in ["depth", "height", "width"] ]
        number_of_items = reduce(lambda x,y:x*y, dim, 1)
        
        # Create the numpy array
        data = numpy.frombuffer(buffer, dtype="h", count=number_of_items, offset=offset)
        data.shape = dim
        
        # Fetch endianness of file
        if not self.header.has_key("endian") :
            file_endianness = ">"
        elif self.header["endian"].lower() == "bigendian" :
            file_endianness = ">"
        else:
            file_endianness = "<"
        
        if pack("<h", 1) == pack("=h",1):
            machine_endianness = "<"
        elif pack(">h", 1) == pack("=h",1):
            machine_endianness = ">"
        
        # Swap bytes if file and machine endianness differ
        # Note that numpy.ndarray.byteswap just calls array.byteswap and that the
        # endianness of the dtype does not get updated
        # The above comment is valid as of 2008-07-11
        if machine_endianness != file_endianness :
            data = data.byteswap()
        
        # Swap y and z
        data = data.swapaxes(1,0)
        # Mirror on y and z
        data = numpy.fliplr(data)
        data = numpy.flipud(data)
        
        file.close()
        
        return data
    
    def load_metadata(self, index=0) :
        metadata = {}
        
        # Normalize rescaling
        metadata["slope"] = 2.**self.header["icomp"][index]
        
        # Normalize spacing
        metadata["spacing"] = numpy.asarray(
            (self.header["dy"][index], 
                self.header["dz"][index], 
                self.header["dx"][index]))
        if (metadata["spacing"] == (0,0,0)).all() :
            metadata["spacing"] = numpy.ones(3)
        
        # Normalize Metadata
        for ipb_name, dicom_name in self._ipb_to_metadata.items() :
            if ipb_name in self.header :
                metadata[dicom_name] = self.header[ipb_name]
        
        # Normalize annotations
        metadata["annotations"] = self._get_annotations(index)
        
        # Normalize colormap
        if "colormap" in self.header :
            colormap = self.header["colormap"][index]
            colormap = IPBColorMaps.get_name(colormap)
        else :
            colormap = "gray"
        metadata["colormap"] = colormap
        
        metadata["filename"] = self.filename
        
        metadata["header"] = self.header
        
        return metadata
    
    def can_save(self, image) :
        return (self._header_filename is not None)
    
    def save(self, image) :
        self._write_header(image)
        self._write_data(image)
    
    def _get_annotation(self, image_index, element_index):
        element = "annotation_%04d_%04d"%(image_index+1, element_index)
        annotation = ImageAnnotation()
        
        postion_string = decode(self.header[element + "_position"])
        
        ipb_position = [float(x) for x in postion_string.split(" ")]
        annotation.position = [self.header["height"][image_index]-ipb_position[1], 
            self.header["width"][image_index]-ipb_position[2], 
            ipb_position[0]]
        
        annotation.color[0] = float(self.header[element + "_defaultColor.red"])/65535.
        annotation.color[1] = float(self.header[element + "_defaultColor.green"])/65535.
        annotation.color[2] = float(self.header[element + "_defaultColor.blue"])/65535.
        
        annotation.label = decode(self.header[element + "_textLabel"])
        
        shapes = {
            0 : ImageAnnotation.Shape.sphere,
            1 : ImageAnnotation.Shape.cube,
            2 : ImageAnnotation.Shape.cross,
            3 : ImageAnnotation.Shape.point,
        }
        
        if element + "_drawType" in self.header :
            annotation.shape = shapes[int(self.header[element + "_drawType"])]
        if element + "_radius" in self.header :
            annotation.size = float(self.header[element + "_radius"])
        if element + "_filled" in self.header :
            annotation.filled = (int(self.header[element + "_filled"]) == 1)
        if element + "_comment" in self.header :
            annotation.comment = decode(self.header[element + "_comment"])        
        
        return annotation
        
    def _get_annotations(self, index):
        list_name = "annotation_%04d"%(index+1)
        
        if list_name not in self.header.keys() : 
            return []
        
        list_size = int(self.header[list_name])
        
        l = []
        
        for i in range(list_size) :
            l.append(self._get_annotation(index, i))
        
        return l
    
    def _get_filenames(self, filename) :
        components = filename.split(".")
        radix = None
        extension = None
        if components[-1] == "gz" :
            radix = string.join(components[:-2], ".")
            extension = string.join(components[-2:], ".")
            self._is_gz = True 
        else :
            radix = string.join(components[:-1], ".")
            extension = string.join(components[-1:], ".")
            self._is_gz = False
        
        header_filename = None
        data_filename = None
        if extension.lower() == "ipb" :
            header_filename = filename
            if os.path.isfile(radix + ".img") :
                data_filename = radix + ".img"
                self._is_gz = False
            elif os.path.isfile(radix + ".img.gz") :
                data_filename = radix + ".img.gz"
                self._is_gz = True
        elif extension.lower() == "img" or extension.lower() == "img.gz":
            header_filename = radix + ".ipb"
            data_filename = filename
        
        return (header_filename, data_filename)
    
    _type_of_key = {
        # Data size and layout
        "number of images" : int,
        "width" : int,
        "height" : int,
        "depth" : int,
        "bits per pixel" : int,
        "offset in raw" : int,
        "endian" : str,
        # Space transformation
        "dx" : float,
        "dy" : float,
        "dz" : float,
        # Data transformation
        "icomp" : float,
        "rcoeff" : float,
        "min pixel" : int,
        "max pixel" : int,
        # Patient / physician identity
        "patient name" : str,
        "date of birth" : str,
        "doctor name" : str,
        "date of exam" : str,
        # Display
        "colormap" : int,
        "display style" : int,
        "x3dc" : int,
        "y3dc" : int,
        "z3dc" : int,
        "zoom" : int,
        "zoom_x" : int,
        "zoom_y" : int,
        "zoom_z" : int,
        # Talairach
        "AC_x" : float,
        "AC_y" : float,
        "AC_z" : float,
        "PC_x" : float,
        "PC_y" : float,
        "PC_z" : float,
        "AP_x" : float,
        "AP_y" : float,
        "AP_z" : float,
        "PP_x" : float,
        "PP_y" : float,
        "PP_z" : float,
        "RM_x" : float,
        "RM_y" : float,
        "RM_z" : float,
        "LM_x" : float,
        "LM_y" : float,
        "LM_z" : float,
        "SP_x" : float,
        "SP_y" : float,
        "SP_z" : float,
        "IP_x" : float,
        "IP_y" : float,
        "IP_z" : float,
        "IC_x" : float,
        "IC_y" : float,
        "IC_z" : float,
        "Tu_x" : float,
        "Tu_y" : float,
        "Tu_z" : float,
        "Tv_x" : float,
        "Tv_y" : float,
        "Tv_z" : float,
        "Tw_x" : float,
        "Tw_y" : float,
        "Tw_z" : float,
        "x_AP" : float,
        "x_PC" : float,
        "x_PP" : float,
        "z_IP" : float,
        "z_SP": float,
        "y_RM" : float,
        "y_LM" : float,
        # Coupes curvilignes
        "CURV_ISINIT" : int,
        "CURV_OFFSET" : int,
        "CURV_THETA_MIN" : float,
        "CURV_THETA_MAX" : float,
        "CURV_PHI_MIN" : float,
        "CURV_PHI_MAX" : float,
        "CURV_NB_THETA" : int,
        "CURV_NB_PHI" : int,
        "CURV_RAYON" : float,
        "CURV_OLDRAYON" : float,
        "CURV_THETA" : float,
        "CURV_PHI" : float,
        "CURV_ZOOM" : float,
        "CURV_XG" : int,
        "CURV_YG" : int,
        "CURV_ZG" : int,
        "CURV_NB_PARAMS" : int,
        # DTI
        "is_dti" : int,
        "nb_directions" : int,
        "gx" : float,
        "gy" : float,
        "gz" : float,
        "b" : float,
        # Other
        "overblancking" : int,
        "cutoff max" : int,
        "cutoff min" : int,
        "ORIGINAL_OFFSET" : int,
        "MASK_OFFSET" : int,
    }
    
    def _header_as_dict(self) :
        file = open(self._header_filename)
        lines = file.readlines()
        file.close()
        
        if lines[0].rstrip() == self._magic :
            self._is_ipb = True
        
        dictionary = {}
        for line in lines[1:] :
            line = line.rstrip()
            if len(line)==0 or line[0] == "-":
                continue
            else :
                key, value = line.split("=", 1)
                if value != "" :
                    dictionary[key] = value.split(" ")
                    if key in self._type_of_key : 
                        dictionary[key] = [ 
                            self._type_of_key[key](value) 
                            for value in dictionary[key] 
                        ]
                    else : 
                        dictionary[key] = dictionary[key][0]
        
        # Process global entries
        if "patient name" in dictionary :
            dictionary["patient name"] = " ".join(dictionary["patient name"])
        if "date of birth" in dictionary :
            dictionary["date of birth"] = " ".join(dictionary["date of birth"])
        if "doctor name" in dictionary :
            dictionary["doctor name"] = " ".join(dictionary["doctor name"])
        if "date of exam" in dictionary :
            dictionary["date of exam"] = " ".join(dictionary["date of exam"])
        dictionary["number of images"] = dictionary["number of images"][0]
        if dictionary.has_key("endian") :
            dictionary["endian"] = " ".join(dictionary["endian"])
        # Ignore ORIGINAL_OFFSET and MASK_OFFSET
        if dictionary.has_key("ORIGINAL_OFFSET") :
            del dictionary["ORIGINAL_OFFSET"]
        if dictionary.has_key("MASK_OFFSET") :
            del dictionary["MASK_OFFSET"]
        
        return dictionary

    def _get_icomp_and_rcoeff(self, data_min, data_max):
        epsilon = 1e-6
        
        max_abs = max(abs(data_min), data_max)
        
        if max_abs == 0 : # image has 0 everywhere
            rcoeff = 1
            icomp = 0
        else : 
            # Save as 16 bits signed
            rcoeff = max_abs/32767.
            icomp = math.log(rcoeff, 2.)
        
        # Get the closest integer icomp
        if icomp <= 0.0 :
            if abs(icomp-math.floor(icomp+0.5)) > epsilon : 
                icomp = math.ceil(icomp)
            else : 
                icomp = math.floor(icomp+0.5)
        else : 
            if abs(icomp-math.floor(icomp+0.5)) > epsilon : 
                icomp = math.ceil(icomp+1)
            else : 
                icomp = math.floor(icomp+1.0+0.5)
        
        rcoeff = 2**icomp
        
        return (icomp, rcoeff)

    def _get_format(self, value):
        format = ""
        if isinstance(value, (str, unicode)) :
            format = "%s" 
        elif type(value) == int : 
            format = "%i"
        elif type(value) == float : 
            format = "%.6f"
        elif type(value) == numpy.float32 : 
            format = "%.6f"
        return format

    def _write_header_entry(self, file, name, values):
        file.write(name + "=")
        for value in values[:-1] :
            if value is None : 
                value = ""
            file.write(self._get_format(value)%(value) + " ")
        
        if values[-1] is None : 
            value = ""
        else : 
            value = values[-1]
        file.write(self._get_format(value)%(value))
        file.write("\n")

    def _write_header(self, image):
        file = open(self._header_filename, "w")
        
        file.write(self._magic+"\n")
        
        if "patient" in image.metadata :
            patient = image.metadata["patient"]
            self._write_header_entry(file, "patient name", 
                                     [patient.get("patients_name", "")])
            self._write_header_entry(file, "date of birth", 
                                     [patient.get("patients_birth_date", "")])
        else : 
            self._write_header_entry(file, "patient name",  [""]) 
            self._write_header_entry(file, "date of birth", [""])
        
        if "study" in image.metadata :
            study = image.metadata["study"]
            self._write_header_entry(file, "doctor name", 
                                     [study.get("referring_physicians_name", "")])
        else :
            self._write_header_entry(file, "doctor name", [""])
        
        if "image" in image.metadata :
            image_metadata = image.metadata["image"]
            self._write_header_entry(file, "date of exam", 
                                     [image_metadata.get("series_date", "")])
        else :
            self._write_header_entry(file, "date of exam", [""])
        
        self._write_header_entry(file, "number of images", [1])
        
        self._write_header_entry(file, "colormap", [7])
        self._write_header_entry(file, "display style", [1])
        self._write_header_entry(file, "overblancking", [0])
        
        # Special axis order in IPB file format
        self._write_header_entry(file, "width",  [image.shape[-1]])
        self._write_header_entry(file, "height", [image.shape[-3]])
        self._write_header_entry(file, "depth",  [image.shape[-2]])
        
        # Find the icomp
        data_min = image.data.min()
        data_max = image.data.max()
        icomp, rcoeff = self._get_icomp_and_rcoeff(data_min, data_max)
        self._write_header_entry(file, "rcoeff",  [rcoeff])
        self._write_header_entry(file, "icomp",  [icomp])
        
        # Min and max value, in array units
        array_min = int(data_min/rcoeff)
        array_max = int(data_max/rcoeff)
        self._write_header_entry(file, "max pixel",  [array_max])
        self._write_header_entry(file, "min pixel",  [array_min])
        
        self._write_header_entry(file, "cutoff max", [array_max])
        self._write_header_entry(file, "cutoff min", [0]) 

        self._write_header_entry(file, "bits per pixel", [str("16")])
        
        self._write_header_entry(file, "dx", [image.spacing[-1]])
        self._write_header_entry(file, "dy", [image.spacing[-3]])
        self._write_header_entry(file, "dz", [image.spacing[-2]])
        
        self._write_header_entry(file, "x3dc", [str(128)])
        self._write_header_entry(file, "y3dc", [str(128)])
        self._write_header_entry(file, "z3dc", [str(128)])
        
        self._write_header_entry(file, "zoom", [str(1)])
        self._write_header_entry(file, "zoom_x", [str(0)])
        self._write_header_entry(file, "zoom_y", [str(0)])
        self._write_header_entry(file, "zoom_z", [str(0)])
 
        self._write_header_entry(file, "offset in raw" , [str(0)])
        
        if "dti" in image.metadata : 
            file.write("--- DTI ---\n")
            dti = image.metadata["dti"]
            self._write_header_entry(file, "is_dti" , [dti.type])
            self._write_header_entry(file, "nb_directions" , [dti.number_of_directions])
            self._write_header_entry(file, "gx" , [dti.gradient[2]])
            self._write_header_entry(file, "gy" , [dti.gradient[1]])
            self._write_header_entry(file, "gz" , [dti.gradient[0]])
            self._write_header_entry(file, "b" , [dti.b])
        
        if "talairach" in image.metadata :
            file.write("""--- Atlas de Talairach ---
AC_x=0.000000 
AC_y=0.000000 
AC_z=0.000000 
PC_x=0.000000 
PC_y=0.000000 
PC_z=0.000000 
AP_x=0.000000 
AP_y=0.000000 
AP_z=0.000000 
PP_x=0.000000 
PP_y=0.000000 
PP_z=0.000000 
RM_x=0.000000 
RM_y=0.000000 
RM_z=0.000000 
LM_x=0.000000 
LM_y=0.000000 
LM_z=0.000000 
SP_x=0.000000 
SP_y=0.000000 
SP_z=0.000000 
IP_x=0.000000 
IP_y=0.000000 
IP_z=0.000000 
IC_x=0.000000 
IC_y=0.000000 
IC_z=0.000000 
Tu_x=0.000000 
Tu_y=0.000000 
Tu_z=0.000000 
Tv_x=0.000000 
Tv_y=0.000000 
Tv_z=0.000000 
Tw_x=0.000000 
Tw_y=0.000000 
Tw_z=0.000000 
x_AP=0.000000 
x_PC=0.000000 
x_PP=0.000000 
z_IP=0.000000 
z_SP=0.000000 
y_RM=0.000000 
y_LM=0.000000""")
        #--- Coupes curvilignes --- 
        #CURV_ISINIT=0 
        #CURV_OFFSET=0 
        #CURV_THETA_MIN=0.000000 
        #CURV_THETA_MAX=0.000000 
        #CURV_PHI_MIN=0.000000 
        #CURV_PHI_MAX=0.000000 
        #CURV_NB_THETA=0 
        #CURV_NB_PHI=0 
        #CURV_RAYON=0.000000 
        #CURV_OLDRAYON=0.000000 
        #CURV_THETA=0.000000 
        #CURV_PHI=0.000000 
        #CURV_ZOOM=0.000000 
        #CURV_XG=0 
        #CURV_YG=0 
        #CURV_ZG=0 
        #CURV_NB_PARAMS=0
        if pack("<h", 1) == pack("=h",1):
            endianness = "littleEndian"
        elif pack(">h", 1) == pack("=h",1):
            endianness = "bigEndian"
        self._write_header_entry(file, "endian", [endianness])
        
        file.close()
    
    def _write_data(self, image):
        
        # Transform the image to the closest axis-aligned approximation of RAS space
        MatrixType = itk.Matrix[itk.D, 3, 3]
        MatrixBridge = itk.MatrixBridge[MatrixType]
        direction = coordinate_system.best_fitting_axes_aligned_matrix(image.direction)
        direction = numpy.dot(coordinate_system.RAS, direction)
        itk_direction = numpy.fliplr(numpy.flipud(direction))
        itk_direction = MatrixBridge.GetMatrixFromArray(itk_direction)
        
        itk_image = medipy.itk.medipy_image_to_itk_image(image, False)
        itk_image.SetDirection(itk_direction)
        
        itk_RAS = numpy.fliplr(numpy.flipud(coordinate_system.RAS))
        itk_RAS = MatrixBridge.GetMatrixFromArray(itk_RAS)
        
        orienter = itk.OrientImageFilter[itk_image, itk_image].New()
        orienter.UseImageDirectionOn()
        orienter.SetInput(itk_image)
        orienter.SetDesiredCoordinateDirection(itk_RAS)
        orienter.Update()
        
        ras_image = medipy.itk.itk_image_to_medipy_image(orienter.GetOutput(), None, True)
        
        # Find the icomp
        data_min = ras_image.data.min()
        data_max = ras_image.data.max()
        icomp, rcoeff = self._get_icomp_and_rcoeff(data_min, data_max)
        
        data = ras_image.data/rcoeff
        data = data.astype(numpy.int16)
        
        # Apply the inverse transformations done in load
        # Mirror on z and y
        data = numpy.flipud(data)
        data = numpy.fliplr(data)
        # Swap y and z
        data = data.swapaxes(1,0)
        
        last_dot = self._header_filename.rfind(".")
        data_filename = self._header_filename[:last_dot]
        data_filename += ".img"
        file = open(data_filename, "wb")
        data.tofile(file)

if __name__ == "__main__" :
    import sys
    
    loader = IPB(sys.argv[1])
    metadata = loader.load_metadata()
    del metadata["header"]
    print metadata