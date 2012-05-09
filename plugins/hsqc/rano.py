"""
Test
"""
import corps
import logging
import sys
import warnings
import numpy as np
from medipy.io import load
import os
import fonc_util as f
#import segmentation as seg
from medipy.base import Image
import gen_in
import lawij
import ser_xml
import MCMC as m
import xml_gen
import locale
if __name__=="__main__":
    #da=np.zeros([21,21])
    #da[10,10]=1
    #db=np.zeros([21,21])
    #db[14,14]=1
    #db[12,13]=1
    #db[14,15]=1
    #db[14,17]=1
    #db[14,18]=1
    #db[8,13]=1.5
    #db[12,11]=1
    root='/home/miv/belghith/Bureau/test2/colon200/6/pdata/1'
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    #root='/home/miv/belghith/Bureau/KAROM/Akram/nmr/RP/5/pdata/1'
    tabp,tabc=gen_in.ini(root)
    xml_gen.gen_s_chggen(tabc,root,"peak1.xml")
    #lolo= load(os.path.join(root,"2rr"))
    #H=lolo.data
    #M=np.size(H[:,1])
    #N=np.size(H[1,:])
    #im1=np.zeros([M,N],dtype=np.single)
    #im2=1*n.ones([2000,2000],dtype=n.single)
    #im2[9,9]=2;
    #Ima1=Image(data=im1)
    #Ima2=Image(data=H)
    #print im2
    #ro=0.7413*corps.interquartile_range(H)
    #seg.threshold(Ima2,ro,1,5,7,Ima1)
    #im=Ima1.data
    #root='/base_image/nmr/colon-HSQC/57lanmic-HSQC/4/pdata/1'
    #lolo= load(os.path.join(root,"2rr"))
    #H1=lolo.data
    #ind=f.find(im,im==1)
    #for m in range(np.size(ind[0])):
        #r,indi=f.alig(ind[0][m],ind[1][m],H1,H)
        #if r[0]==0 and r[1]==0:
            #print r,indi,ind[0][m],ind[1][m]
            #print H1[ind[0][m]-2:ind[0][m]+3,ind[1][m]-2:ind[1][m]+3]
            #print H1[ind[0][m]-2:ind[0][m]+3,ind[1][m]-2:ind[1][m]+3]
    
    locale.setlocale(locale.LC_ALL, old_locale)