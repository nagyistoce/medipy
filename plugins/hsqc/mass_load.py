import numpy as np
from medipy.gui import colormaps
from medipy.base import Image
from numpy import linalg
from medipy.io import load

def massive_load(input1) :
    """
    

    

       input1 : Directory
           path

       
      
       
     """
    
    base_directory=input1
    for root,dirs,files in os.walk(base_directory=input1):

        if 'wcmac.py'in files:
            print root
