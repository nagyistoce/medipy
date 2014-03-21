##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import operator
import os
import struct

import numpy

import medipy.base
import unittest

def can_load(filename, shape, dtype, endianness="@", offset=0, skip_end=False):
    """ Test if an image of given shape and dtype can be loaded from a file.
    """
    
    filesize = os.path.getsize(filename)
    
    format = get_struct_format(shape, dtype, endianness, offset)
    size = struct.calcsize(format)
    
    if size < filesize and skip_end :
        can_import = True
    elif size == filesize :
        can_import = True
    else :
        can_import = False
    
    return can_import

def load(filename, shape, dtype, endianness="@", offset=0, skip_end=False):
    """ Load an image from raw data using the given shape and dtype.
    """
    
    format = get_struct_format(shape, dtype, endianness, offset)
    fd = open(filename, "rb")
    data = struct.unpack(format, fd.read())
    array = numpy.fromstring(data[0], dtype, reduce(operator.mul, shape, 1))
    array = array.reshape(shape)
    
    metadata = {"loader" : {"url" : filename } }
    
    return medipy.base.Image(data=array, metadata=metadata)

def get_struct_format(shape, dtype, endianness="@", offset=0):
    # Compute the length in individual items
    length = reduce(operator.mul, shape, 1)
    
    # Compute the length in bytes
    length *= numpy.nbytes[dtype]
    
    format = "{endianness}{offset}x{length}s".format(
        endianness=endianness, offset=offset, length=length)
    return format
