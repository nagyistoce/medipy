"""
Test
"""

import numpy as np
import os
import wx
from medipy.gui import colormaps
from medipy.base import Image
from numpy import linalg
from medipy.components.io import load
import medipy.io.rbnmr as rbnmr
from decimal import *
import fun
import conv
import inipeak
import matplotlib.pyplot as plt
x=0
im={}
iH={}
iC={}
#for root,dirs,files in os.walk(os.getcwd()):
for root,dirs,files in os.walk('/base_image/nmr/colon-HSQC'):
    if '2rr'in files:
        x=x+1
        print root
        z= load(os.path.join(root,"2rr"))
        im[x]=z.data
        c=conv.convc(root)
        h=conv.convh(root)
        iH[x]=h
        iC[x]=c
        r=inipeak.apic(z.data)
       

#im[1][1][1]=100
#print iH[1][0]
#a=inipeak.apic(im[1])
