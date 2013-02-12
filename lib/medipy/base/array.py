##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import functools
import numpy

def crop(array, shape, keep_axis_begin=True) :
    """ Crop the array to the minimum of the given shape and the input array 
        shape. If ``keep_axis_begin`` is ``True``, then elements at the begin 
        of each axis are preserved, otherwise elements at the end of each axis 
        are preserved.
    """
    
    if array.shape == shape :
        return array
    
    slices = ()
    for axis in range(len(shape)) :
        if keep_axis_begin :
            begin = 0
            end = min(shape[axis], numpy.shape(array)[axis])
        else :
            end = numpy.shape(array)[axis]
            begin = max(0, end-shape[axis])
        
        s = slice(begin, end)
        slices += (s,)
    
    cropped_array = numpy.asarray(array)[slices]
    
    return cropped_array

def pad(array, shape, mode, keep_axis_begin=True, **kwargs) :
    """ Pad the array to the maximum of the given shape and the input array 
        shape. If ``keep_axis_begin`` is ``True``, then padding will be 
        inserted at the end of each axis, otherwise padding will be inserted at
        the beginning of each axis. The ``mode`` argument and the extra keywords
        arguments determine how the padding is done :
        
        * ``"constant"`` : a constant value is used (keyword argument : ``value``)
    """
    
    if array.shape == shape :
        return array
    
    # Padding functions, passed to numpy.apply_along_axis. All padding functions
    # must have a signature similar to
    #   my_padding_function(a, size, at_axis_end, ...)
    # where :
    #   * a : 1D array to be padded
    #   * size : size of the new (padded) array. Must be greater than len(a)
    #   * keep_axis_begin : cf. documentation of pad 
    
    def constant(a, size, keep_axis_begin, value) :
        """ Pad the array with a constant value
        """
        
        padded = value*numpy.ones((size,))
        if keep_axis_begin :
            padded[:len(a)] = a
        else :
            padded[-len(a):] = a
        return padded
    
    # Create the padding functor
    function = functools.partial(locals()[mode], **kwargs)
    
    padded = array
    for axis in range(len(shape)) :
        size = max(shape[axis], numpy.shape(array)[axis])
        padded = numpy.apply_along_axis(function, axis, padded, size, keep_axis_begin)
    
    return padded

def reshape(array, shape, mode, keep_axis_begin=True, **kwargs) :
    """ Crop and then pad the array so that the returned array has the given 
        shape. Cf. the documentation of :func:`crop` and :func:`pad` for the 
        description of the parameters.
    """ 
    
    cropped = crop(array, shape, keep_axis_begin)
    padded = pad(cropped, shape, mode, keep_axis_begin, **kwargs)
    
    return padded