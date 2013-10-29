##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import numpy
import medipy.medimax.recalage
import medipy.base
import medipy.io

def load_trf(ftrf,fmodel_wrap) :
    """ Load a .trf transformation field

    <gui>
        <item name="ftrf" type="File" label=".trf deformation field" />
        <item name="fmodel_wrap" type="File" label="Wrapped Diffusion Image"/>
        <item name="filed" type="Image" initializer="output=True" role="return" label="Deformation field"/>
    </gui>
    """

    images_wrap = medipy.io.io.load_serie(fmodel_wrap)
    model_wrap = images_wrap[0]

    ux = medipy.base.Image(shape=(1), dtype=numpy.float32)
    uy = medipy.base.Image(shape=(1), dtype=numpy.float32)
    uz = medipy.base.Image(shape=(1), dtype=numpy.float32)

    imSource = medipy.base.Image(shape=(1), spacing=model_wrap.spacing, origin=model_wrap.origin, direction=model_wrap.direction, dtype=numpy.float32)
    medipy.medimax.recalage.LoadTrfFile(str(ftrf),ux,uy,uz,imSource)

    field = medipy.base.Image(data=numpy.zeros(ux.shape+(3,),dtype=numpy.single),spacing=ux.spacing,data_type="vector")
    field[...,0] = ux.data
    field[...,1] = uy.data
    field[...,2] = uz.data

    return field

if __name__ == '__main__':

    fdata1 = "/home/grigis/data/SOM/01/data_fsl_eddy.nii.gz"
    fdata2 = "/home/grigis/data/SOM/02/data_fsl_eddy.nii.gz"
    ftrf = "/home/grigis/data/SOM/affnl.trf"

    print load_trf(ftrf,fdata2).shape
