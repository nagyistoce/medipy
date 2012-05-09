"""
Functions used in peak initialization step
"""
import numpy as np
import os
import wx
from decimal import *
import fun
import conv
import fonc_util
def apic(im):
    M=np.size(im[:,1])
    N=np.size(im[1,:])
    mask=np.ones([M,N])
    pic=np.zeros([M,N])
    for i in range(10,M-10):
        for j in range(10,N-10):
            z=im[i-2:i+3,j-2:j+3]
            ind=fonc_util.find(z,z==np.amax(z))
            if ind[0][0]==2 and ind[1][0]==2 and np.amax(z)>10000:
                pic[i][j]=1
                mask[i-3:i+4,j-3:j+4]=0
    return pic            
    