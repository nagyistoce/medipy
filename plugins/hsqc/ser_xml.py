"""
Corpus reading
"""
#import corps
#import logging
#import sys
#import warnings
#import numpy as np
#from pylab import *
#from medipy.components.io import load
#import os
#import fonc_util as f
#import medipy.components.belghith.segmentation as seg
#from medipy.base import Image
import numpy as np
def read(nom):
    z={}
    f=open(nom, 'r')
    l=f.readlines()
    #print np.size(l[:])
    dem=-1
    for h in range(np.size(l[:])):
        C=l[h][:]
        H="F1="
        toul= C.find(H)
        if toul>0:
            dem= dem+1
            touli=toul+4
            for h1 in range(touli,touli+20):
                if (l[h][h1]=='"'):
                    break
            z[dem,0]= l[h][toul+4:h1]
            H="F2="
            toul= C.find(H)
            touli=toul+4
            for h1 in range(touli,touli+20):
                if (l[h][h1]=='"'):
                    break
            z[dem,1]= l[h][toul+4:h1]
            H="annotation="
            toul= C.find(H)
            touli=toul+12
            for h1 in range(touli,touli+30):
                if (l[h][h1]=='"'):
                    break
            z[dem,2]= l[h][toul+12:h1]
            H="intensity="
            toul= C.find(H)
            touli=toul+11
            for h1 in range(touli,touli+30):
                if (l[h][h1]=='"'):
                    break
            z[dem,3]= l[h][toul+11:h1]
    return z
    f.close()
def read_s(nom):
    z={}
    f=open(nom, 'r')
    l=f.readlines()
    #print np.size(l[:])
    dem=-1
    for h in range(np.size(l[:])):
        C=l[h][:]
        H="F1="
        toul= C.find(H)
        if toul>0:
            dem= dem+1
            touli=toul+4
            for h1 in range(touli,touli+20):
                if (l[h][h1]=='"'):
                    break
            z[dem,0]= l[h][toul+4:h1]
            H="F2="
            toul= C.find(H)
            touli=toul+4
            for h1 in range(touli,touli+20):
                if (l[h][h1]=='"'):
                    break
            z[dem,1]= l[h][toul+4:h1]
            H="annotation="
            toul= C.find(H)
            touli=toul+12
            for h1 in range(touli,touli+30):
                if (l[h][h1]=='_'):
                    break
                elif (l[h][h1]=='"'):
                    break
            z[dem,2]= l[h][toul+12:h1]
            H="intensity="
            toul= C.find(H)
            touli=toul+11
            for h1 in range(touli,touli+30):
                if (l[h][h1]=='"'):
                    break
            z[dem,3]= l[h][toul+11:h1]
    return z
    f.close()    
def exceptions(nom):
    z={}
    nom=nom+'/except.txt'
    f=open(nom, 'r')
    l=f.readlines()
    #print np.size(l[:])
    dem=0
    for h in range(np.size(l[:])):
        for h1 in range(1,200):
                if (l[h][h1]=='"'):
                    break
        z[dem]=l[h][1:h1]
        #print z[dem]
        dem=dem+1
    f.close()
    return z