##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


""" Basic tools for DTI
"""

import numpy as np
import medipy.base
import medipy.diffusion
import medipy.itk

def spectral_decomposition(slice_tensor):
    shape = slice_tensor.shape

    eigVal = medipy.base.Image(data=np.zeros(shape[:3]+(3,),dtype=np.single),data_type="vector")
    eigVec = medipy.base.Image(data=np.zeros(shape[:3]+(9,),dtype=np.single),data_type="vector")
    medipy.diffusion.dtiEigItk(slice_tensor,eigVal,eigVec)

    return eigVal,eigVec



def decompose_tensor(tensor):
    #outputs multiplicity as well so need to unique
    eigenvals, eigenvecs = np.linalg.eig(tensor)

    #need to sort the eigenvalues and associated eigenvectors
    order = eigenvals.argsort()[::-1]
    eigenvecs = eigenvecs[:, order]
    eigenvals = eigenvals[order]

    #Forcing negative eigenvalues to 0
    eigenvals = np.maximum(eigenvals, 0)

    return eigenvals, eigenvecs

def dti6to33(dt6):
    """ Full second order symmetric tensor from the six independent components
    """
    dt33 = np.zeros(dt6.shape[:3]+(3,3),dtype=np.single)
    dt33[:,:,:,0,0] = dt6[:,:,:,0]
    dt33[:,:,:,0,1] = dt6[:,:,:,1]
    dt33[:,:,:,0,2] = dt6[:,:,:,2]
    dt33[:,:,:,1,0] = dt6[:,:,:,1]
    dt33[:,:,:,1,1] = dt6[:,:,:,3]
    dt33[:,:,:,1,2] = dt6[:,:,:,4]
    dt33[:,:,:,2,0] = dt6[:,:,:,2]
    dt33[:,:,:,2,1] = dt6[:,:,:,4]
    dt33[:,:,:,2,2] = dt6[:,:,:,5]

    return dt33

def rotation33todt6(dt6,R):
    """ Apply a rotation R """

    dt6_r = np.zeros(dt6.shape,dtype=np.single)

    dt6_r[:,:,:,0] = R[0,0]*(R[0,0]*dt6[:,:,:,0] + R[0,1]*dt6[:,:,:,1] + R[0,2]*dt6[:,:,:,2])\
                   + R[0,1]*(R[0,0]*dt6[:,:,:,1] + R[0,1]*dt6[:,:,:,3] + R[0,2]*dt6[:,:,:,4])\
                   + R[0,2]*(R[0,0]*dt6[:,:,:,2] + R[0,1]*dt6[:,:,:,4] + R[0,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,1] = R[1,0]*(R[0,0]*dt6[:,:,:,0] + R[0,1]*dt6[:,:,:,1] + R[0,2]*dt6[:,:,:,2])\
                   + R[1,1]*(R[0,0]*dt6[:,:,:,1] + R[0,1]*dt6[:,:,:,3] + R[0,2]*dt6[:,:,:,4])\
                   + R[1,2]*(R[0,0]*dt6[:,:,:,2] + R[0,1]*dt6[:,:,:,4] + R[0,2]*dt6[:,:,:,5])
     	
    dt6_r[:,:,:,2] = R[2,0]*(R[0,0]*dt6[:,:,:,0] + R[0,1]*dt6[:,:,:,1] + R[0,2]*dt6[:,:,:,2])\
                   + R[2,1]*(R[0,0]*dt6[:,:,:,1] + R[0,1]*dt6[:,:,:,3] + R[0,2]*dt6[:,:,:,4])\
                   + R[2,2]*(R[0,0]*dt6[:,:,:,2] + R[0,1]*dt6[:,:,:,4] + R[0,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,3] = R[1,0]*(R[1,0]*dt6[:,:,:,0] + R[1,1]*dt6[:,:,:,1] + R[1,2]*dt6[:,:,:,2])\
                   + R[1,1]*(R[1,0]*dt6[:,:,:,1] + R[1,1]*dt6[:,:,:,3] + R[1,2]*dt6[:,:,:,4])\
                   + R[1,2]*(R[1,0]*dt6[:,:,:,2] + R[1,1]*dt6[:,:,:,4] + R[1,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,4] = R[2,0]*(R[1,0]*dt6[:,:,:,0] + R[1,1]*dt6[:,:,:,1] + R[1,2]*dt6[:,:,:,2])\
                   + R[2,1]*(R[1,0]*dt6[:,:,:,1] + R[1,1]*dt6[:,:,:,3] + R[1,2]*dt6[:,:,:,4])\
                   + R[2,2]*(R[1,0]*dt6[:,:,:,2] + R[1,1]*dt6[:,:,:,4] + R[1,2]*dt6[:,:,:,5])

    dt6_r[:,:,:,5] = R[2,0]*(R[2,0]*dt6[:,:,:,0] + R[2,1]*dt6[:,:,:,1] + R[2,2]*dt6[:,:,:,2])\
                   + R[2,1]*(R[2,0]*dt6[:,:,:,1] + R[2,1]*dt6[:,:,:,3] + R[2,2]*dt6[:,:,:,4])\
                   + R[2,2]*(R[2,0]*dt6[:,:,:,2] + R[2,1]*dt6[:,:,:,4] + R[2,2]*dt6[:,:,:,5])

    return dt6_r
