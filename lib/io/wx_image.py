##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import logging
from cStringIO import StringIO
import os.path
import numpy
import wx

from medipy.io.io_base import IOBase

class WXImage(IOBase):
    """ Loader for the following common bitmap formats, based on wx.Image :
          * BMP
          * GIF
          * JPEG
          * PNG
          * PCX
          * PNM
          * TIFF
          * XPM
          * ICO
          * CUR
          * ANI
    """
    
    filenames = ["*.jpg", "*.jpeg", "*.png", "*.tif", "*.tiff"]  
    
    def __init__(self, *args, **kwargs):
        IOBase.__init__(self, *args, **kwargs)
        
        self.__image = None
        
    def can_load(self):
        try :
            result = wx.Image.CanRead(self._filename)
        except Exception, e :
            logging.debug("CanRead failed : %s"%e)
            return False
        else : 
            logging.debug("CanRead succeeded : %s"%result)
            return result
    
    def number_of_images(self):
        return wx.Image.GetImageCount(self._filename)
    
    def load_data(self, index=0):
        grey_image = self._image.ConvertToGreyscale()
        data = numpy.fromstring(grey_image.GetData(), dtype=numpy.uint8)
        data = data.reshape(grey_image.GetWidth(), grey_image.GetHeight(), 3)
        return data[:,:,0]
    
    def load_metadata(self, index=0):
        return {}
    
    def can_save(self):
        return False
    
    def _get__image(self, index=0):
        if self.__image is None :
            file_size = os.path.getsize(self._filename)
            nb_chunks = 10
            chunk_size = file_size/nb_chunks
            done = False
            chunk_index = 0
            
            buffer = ""
            file = open(self._filename, "rb")
            
            while not done :
                chunk = file.read(chunk_size)
                
                if self._report_progress is not None :
                    self._report_progress(float(chunk_index)/float(nb_chunks+1))
                    
                if len(chunk) == 0 :
                    done = True
                else :
                    chunk_index += 1 
                    buffer += chunk
            
            buffer = StringIO(buffer)
                
            if not wx.Image.CanReadStream(buffer) :
                raise Exception("Cannot find a handler for file %s"%self._filename)
            self.__image = wx.ImageFromStream(buffer, wx.BITMAP_TYPE_ANY, index)
        return self.__image
    
    _image = property(_get__image)