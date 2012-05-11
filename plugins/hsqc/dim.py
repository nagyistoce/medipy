"""
alignment initialization
"""
import numpy as np
import os
import wx
from medipy.gui import colormaps
from medipy.base import Image
from numpy import linalg
from medipy.io import load
import medipy.io.rbnmr as rbnmr
from scipy.interpolate import interpolate
from decimal import *
import fun
import conv
import inipeak
import MCMC
import warnings
import scipy
from scipy.signal import sepfir2d
import fonc_util as f
def alig(ic,jc,data2,amplitref):
    amplitref=float(amplitref)*1.00000000000000000000000001
    peak=[]
    amp=[]
    pos=[]
    #r1.append((boc))
    #peak=np.array([[0,0]])
    #amp= np.array([0])
    #pos=np.array([[0,0]])
    indi=0
    z=data2[ic-3:ic+4,jc-5:jc+6]
    m=np.size(z[:,1])
    n=np.size(z[1,:])
    for i in range(1,m-2):
        for j in range(1,n-1):
            zi=z[i-1:i+2,j-1:j+2]
            im=f.find(zi,zi==np.amax(zi))
            if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>0:
                #print i,j
                #t1=np.array([[3-i,5-j]])*1.
                peak.append((3-i,5-j))
                #peak= np.concatenate((peak, t1), axis=0)
                #t2=np.array([np.amax(zi)])*1.
                #amp= np.concatenate((amp, t2), axis=0)
                amp.append((np.amax(zi)))
    print peak
    #raw_input()
    if len(peak)>0:
                print 'd5al'
                #peak=np.delete(peak,0,0)
                #amp=np.delete(amp,0)
                if len(peak)==1:
                    pos=peak
                    #print 'oooook'
                else:
                    av=f.find(np.abs(peak[:][0]),np.abs(peak[:][0])==np.amin(np.abs(peak[:][0])))
                    if np.size(av)==1:
                        pos=peak[av]
                        print 'oooook'
                    else:
                        
                        pook=peak[av,:][0]
                        bv=f.find(np.abs(pook[:,1]),np.abs(pook[:,1])==np.amin(np.abs(pook[:,1])))
                        if np.size(bv)==1:
                            pos=pook[bv,:][0]
                        else:
                            piik=pook[bv,:][0]
                            amp1=amp[av]
                            amp2=amp1[bv]
                            dif=np.abs(amp[bv]-amplitref)
                            difa=f.find(dif,dif==np.amin(dif))
                            pos= piik[difa,:][0]
                pos=pos[0]                        
    else:
        amplitref=float(amplitref)*1.00000000000000000000000001
        peak=np.array([[0,0]])
        amp= z=np.array([0])
        pos=np.array([[0,0]])
        indi=0
        z=data2[ic-4:ic+5,jc-8:jc+9]
        m=np.size(z[:,1])
        n=np.size(z[1,:])
        for i in range(1,m-2):
            for j in range(1,n-1):
                zi=z[i-1:i+2,j-1:j+2]
                im=f.find(zi,zi==np.amax(zi))
                if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>0:
                #print i,j
                    t1=np.array([[4-i,8-j]])*1.
                    peak= np.concatenate((peak, t1), axis=0)
                    t2=np.array([np.amax(zi)])*1.
                    amp= np.concatenate((amp, t2), axis=0)
        #print peak
        if np.size(peak[:,1])>1:
                print 'd5al'
                peak=np.delete(peak,0,0)
                amp=np.delete(amp,0)
                if np.size(peak[:,1])==1:
                    pos=peak
                else:
                    av=f.find(np.abs(peak[:,0]),np.abs(peak[:,0])==np.amin(np.abs(peak[:,0])))
                    if np.size(av)==1:
                        pos=peak[av,:][0]
                    else:
                        
                        pook=peak[av,:][0]
                        bv=f.find(np.abs(pook[:,1]),np.abs(pook[:,1])==np.amin(np.abs(pook[:,1])))
                        if np.size(bv)==1:
                            pos=pook[bv,:][0]
                        else:
                            piik=pook[bv,:][0]
                            amp1=amp[av]
                            amp2=amp1[bv]
                            dif=np.abs(amp[bv]-amplitref)
                            difa=f.find(dif,dif==np.amin(dif))
                            pos= piik[difa,:][0]
                pos=pos[0]
        else:
                    pos=np.array([100,100])
                    indi=1
                    nz=data2[ic-4:ic+5,jc-5:jc+6]
                    val=amplitref
                    nm=np.size(nz[:,1])
                    nn=np.size(nz[1,:])
                    zios=np.zeros([nm,nn])
                    for ni in range(1,nm-2):
                        for nj in range(1,nn-2):
                            nzi=nz[ni-1:ni+2,nj-1:nj+2]
                            #print nzi[1,1]
                            if (float(nzi[1,1])>float(val/10)):
                                    difo=f.find(nzi,nzi[1,1]>nzi)
                                    zios[ni,nj]=np.size(difo[1])
                    #print zios    
                    if np.amax(zios)>5:
                        hv=f.find(zios,zios>5)
                        #print hv[0][0],hv[1][0]
                        crit=np.abs(nz-val)
                        pos=[hv[0][0],hv[1][0]]
                        for ll in range(1,np.size(hv[1])-1):
                            if crit[hv[0][ll],hv[1][ll]]<crit[pos[0],pos[1]]:
                                pos=[hv[0][ll],hv[1][ll]]   
                        tt=nz[pos[0]-1:pos[0]+2,pos[1]-1:pos[1]+2]
                        im=f.find(tt,tt==np.amax(tt))
                        pos[0]=pos[0]+im[0][0]-1
                        pos[1]=pos[1]+im[1][0]-1
                        pos=np.array([4,5])-pos
                         
    return pos,indi