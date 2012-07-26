import numpy as np
import os
import wx
from medipy.gui import colormaps
from medipy.base import Image
from numpy import linalg
from medipy.io import load
#import segmentation as seg
import fonc_util
import corps
import scipy
import fark
from scipy.signal import sepfir2d
import main_annotate
import main_change
def Changedet(input1,input2) :
    """ multi-spectra disply
    
        <gui>
            <item name="input1" type="Directory" label="GROUP1"/>
            <item name="input2" type="Directory" label="GROUP2"/>
        </gui>
    """ 
 
    #main_annotate.mainanno(input1)
    main_change.runmainc(input1,input2)
    #print a


def annotate(input1) :
    """ Annotation de spectres
    
        <gui>
            <item name="input1" type="Directory" label="Path"/>
        </gui>
    """ 
 
    base_directory = input1
    i=0
    k={}
    for root,dirs,files in os.walk(base_directory):
        if '2rr'in files:
	   #print str(root)
	    k[i]=str(root)
	    i=i+1
    if len(k)>0:
	main_annotate.mainanno(k)
    else:
	 print 'No found spectrum '  
      
       
