##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Load and save files from the local filesystem.
    
    The path is interpreted as an usual filesystem path.
    
    The fragment is of form [ tag "=" value { "&" tag "=" value } ], where tag
    is one of :
      * index : specifies the index of the image to load. Defaults to 0.
"""

import urlparse

import medipy.base
from medipy.base import Image

from medipy.io.dicom_io import Dicom
from medipy.io.ipb import IPB
from medipy.io.itk_io import ITK
from medipy.io.nifti_io import Nifti
from medipy.io.nmr2D import Nmr2D
from medipy.io.wx_image import WXImage
io_classes = [Dicom, ITK, IPB, WXImage, Nmr2D, Nifti]

from os.path import splitext
import numpy as np
import os
from dataset import DataSet
from collections import defaultdict


def load_serie(path, fragment=None, loader=None) :
    """ Load a volume as a serie of N images.
    
        If specified, the "loader" must be an image loader. If not specified,
        a suitable loader will be found automatically.
    """

    image = load(path,fragment,loader)

    # Create serie
    limages = []
    norigin = image.origin[1:]
    nspacing = image.spacing[1:]
    ndirection = image.direction[1:,1:]
    for ndata in list(image.data) :
        args = { "loader" : image.metadata["loader"] }
        limages.append(Image(data=ndata,origin=norigin,spacing=nspacing,direction=ndirection,metadata=args))

    # DWI nifti
    name = path.split('/')[-1].split('.')[0]
    base = "/".join(path.split('/')[:-1])+"/"
    ext = "."+".".join(path.split('/')[-1].split('.')[1:])
    print base,ext,name
    N = len(limages)
    if ext==".nii" or ext==".nii.gz" :
        lfbvec = [base+name+".bvecs",base+name+".bvec",base+"bvecs",base+"bvec"]
        lfbval = [base+name+".bvals",base+name+".bval",base+"bvals",base+"bval"]
        kbvec = [os.path.isfile(f) for f in lfbvec]
        kbval = [os.path.isfile(f) for f in lfbval]
        k1 = None
        if True in kbvec :
            k1 = kbvec.index(True)
        k2 = None
        if True in kbval :
            k2 = kbval.index(True)
        if k1!=None and k2!=None :
            fbvec = lfbvec[k1]
            fbval = lfbval[k2]
            bval = np.array(np.loadtxt(fbval))
            grad = np.array(np.loadtxt(fbvec),dtype=np.single)
            grad = grad.T
            if N==bval.shape[0] and N==grad.shape[0] :
                for cnt in range(N) :
                    dwi_dataset = DataSet()
                    dwi_dataset.diffusion_directionality = "DIRECTIONAL"
                    dwi_dataset.diffusion_bvalue = bval[cnt]
                    gradient_dataset = DataSet()
                    gradient_dataset.diffusion_gradient_orientation = grad[cnt]
                    dwi_dataset.diffusion_gradient_direction_sequence = [gradient_dataset]
                    limages[cnt].metadata["mr_diffusion_sequence"] = [dwi_dataset]

    return limages

def save_serie(limages, path, saver=None) :
    """ Save a serie of N 3D images as a volume.
    
        If specified, the "saver" must be an image saver. If not specified,
        a suitable saver will be found automatically.
    """
    
    N = len(limages)
    name = path.split('/')[-1].split('.')[0]
    base = "/".join(path.split('/')[:-1])+"/"
    ext = "."+".".join(path.split('/')[-1].split('.')[1:])
    if N>0 :
        nshape = [tuple(image.shape) for image in limages]
        nspacing = [tuple(image.spacing) for image in limages]
        dshape = {}
        dspacing = {}
        for t1,t2 in zip(nshape,nspacing) :
            dshape.setdefault(t1, [])
            dspacing.setdefault(t2, [])
        if len(dshape.keys())==1 and len(dspacing.keys())==1 :
            ndata = np.zeros((N,)+dshape.keys()[0],dtype=np.single)
            for cnt in range(N) :
                ndata[cnt] = limages[cnt].data
            args = { "loader" : limages[0].metadata["loader"] }
            ndirection = np.diag([1]*4)
            ndirection[1:,1:] = limages[0].direction
            image = Image(data=ndata,spacing=(1,)+dspacing.keys()[0],origin=(0,)+limages[0].origin,direction=ndirection,metadata=args)

            # Save volume nifti
            print image
            if ext==".nii" or ext==".nii.gz" :
                save(image,path,saver)
            else :
                save(image,base+name+".nii.gz",saver)

            # Save nifti DWI info 

        else :
            raise medipy.base.Exception("All images must have the same shape and spacing.")
    else :
        raise medipy.base.Exception("A serie of images must contain at least one image.")




def load(path, fragment=None, loader=None) :
    """ Load an image.
    
        If specified, the "loader" must be an image loader. If not specified,
        a suitable loader will be found automatically.
    """
    
    loader = loader or get_loader(path)
    
    index = _parse_fragment(fragment) if fragment else 0
    
    data = loader.load_data(index)
    metadata = loader.load_metadata(index)
    metadata["loader"] = {
        "index" : index,
        "loader" : loader
    }
    
    args = {}
    names = ["direction", "origin", "spacing", "data_type", "annotations"]
    if not isinstance(loader, Dicom) :
        names.append("image_type")
    for name in names :
        if name in metadata :
            args[name] = metadata[name]
            del metadata[name]
    args["metadata"] = metadata
    
    return Image(data=data, **args)

def save(image, path, saver=None) :
    """ Save an image.
    
        If specified, the "saver" must be an image saver. If not specified,
        a suitable saver will be found automatically.
    """
    
    saver = saver or get_saver(image, path)
    saver.save(image)

def number_of_images(path, fragment=None, loader=None):
    """ Return the number of images at given path.
    
        If specified, the "loader" must be an image loader. If not specified,
        a suitable loader will be found automatically.
    """
    
    loader = loader or get_loader(path)
    return loader.number_of_images()

def get_loader(path, report_progress=None) :
    """ Search for a suitable loader in io_classes.
    """
    
    for loader_class in io_classes : 
        loader = loader_class(path, report_progress=report_progress)
        if loader.can_load() :
            return loader
    
    # If we get here, no loader was found
    raise medipy.base.Exception("No loader available for {0}".format(path))

def get_saver(image, path, report_progress=None) :
    """ Search for a suitable saver in io_classes.
    """
    
    for saver_class in io_classes : 
        saver = saver_class(path, report_progress=report_progress)
        if saver.can_save(image) :
            return saver
    
    # If we get here, no loader was found
    raise medipy.base.Exception("No saver available for {0}".format(path))

def _parse_fragment(fragment) :
    """ Return the index value in the fragment
    """
    
    for tag, value in urlparse.parse_qsl(fragment) :
        if tag == "index" :
            try :
                index = int(value)
            except ValueError :
                raise medipy.base.Exception("Invalid index : \"{0}\"".format(value))
            else :
                return index
        else :
            raise medipy.base.Exception("Invalid fragment element : \"{0}\"".format(tag))
