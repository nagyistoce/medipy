"""
Functions used in Metabolite annotation
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
import warnings
import scipy
from scipy.signal import sepfir2d
from scipy.optimize import leastsq
def find(a,cond):
    b=np.nonzero(cond)
    return b
def safar(data):
    data[data<0]=0
    return data
    
def cauchy(p2,p3):
	siz=(p2-1)/2*(np.ones(2))
	xxrange = np.arange(-siz[1],siz[1]+1)
	yyrange = np.arange(-siz[1],siz[1]+1)
	X,Y = np.meshgrid(xxrange,yyrange)
	arg=((1/(p3[0]*3.14159))/((1/p3[0])**2+(X*X)))*((1/(p3[1]*3.14159))/((1/p3[1])**2+(Y*Y)))
	eps=2.2204*10**(-16)
	h=arg
	h[h<(eps*np.amax(h))]=0
	sumh=np.sum(h)
	if sumh!=0:
		h=h/sumh
	h=h/np.amax(h)
	return h
    
def expon(p2,p3):
    siz=(p2-1)/2*(np.ones(2))
    xxrange = np.arange(-siz[1],siz[1]+1)
    yyrange = np.arange(-siz[1],siz[1]+1)
    X,Y = np.meshgrid(xxrange,yyrange)
    #arg=((1/p3[0])/((1/p3[0])**2+(X*X)))*((1/p3[1])/((1/p3[1])**2+(Y*Y)))
    arg=(1/6.28)*(1/p3[0])*np.exp(-X*X/(2*p3[0]**2))*(1/p3[1])*np.exp(-X*X/(2*p3[1]**2))
    eps=2.2204*10**(-16)
    h=arg
    h[h<(eps*np.amax(h))]=0
    sumh=np.sum(h)
    if sumh!=0:
        h=h/sumh
    h=h/np.amax(h)
    return h

def exponcauchy(p2,p3):
    siz=(p2-1)/2*(np.ones(2))
    xxrange = np.arange(-siz[1],siz[1]+1)
    yyrange = np.arange(-siz[1],siz[1]+1)
    X,Y = np.meshgrid(xxrange,yyrange)
    #arg=((1/p3[0])/((1/p3[0])**2+(X*X)))*((1/p3[1])/((1/p3[1])**2+(Y*Y)))
    arg=(1/3.14159)*(1/p3[0])*np.exp(-X*X/(2*p3[0]**2))*((1/(p3[1]*3.14159))/((1/p3[1])**2+(Y*Y)))
    eps=2.2204*10**(-16)
    h=arg
    h[h<(eps*np.amax(h))]=0
    sumh=np.sum(h)
    if sumh!=0:
        h=h/sumh
    h=h/np.amax(h)
    return h

def subpix(z,ii,jj):
    trange = np.arange(11)
    ttrange = np.arange(11)
    X,Y = np.meshgrid(trange,ttrange)
    outgrid = interpolate.interp2d(X,Y,z,kind='quintic')
    xx=yy=arange(101)/10.
    l=outgrid(xx,yy)
    l=l[50-9:50+10,50-9:50+10]
    ind=find(l,l==np.amax(l))
    #print l
    #print ind[0][0],ind[1][0]
    ni=ii+(ind[0][0]-9.)/10
    nj=jj+(ind[1][0]-9.)/10
    #print ii,jj
    #print ni,nj
    return[ni,nj] 
def alig(ic,jc,data2,data1):
    peak=np.array([[0,0]])
    amp= z=np.array([0])
    pos=np.array([[0,0]])
    indi=0
    z=data2[ic-5:ic+6,jc-7:jc+8]
    m=np.size(z[:,1])
    n=np.size(z[1,:])
    for i in range(1,m-2):
        for j in range(1,n-1):
            zi=z[i-1:i+2,j-1:j+2]
            im=find(zi,zi==np.amax(zi))
            if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>0:
                #print i,j
                t1=np.array([[5-i,7-j]])*1.
                peak= np.concatenate((peak, t1), axis=0)
                t2=np.array([np.amax(zi)])*1.
                amp= np.concatenate((amp, t2), axis=0)
    if np.size(peak[:,1])>1:
                peak=np.delete(peak,0,0)
                amp=np.delete(amp,0)
                if np.size(peak[:,1])==1:
                    pos=peak
                else:
                    av=find(np.abs(peak[:,0]),np.abs(peak[:,0])==np.amin(np.abs(peak[:,0])))
                    if np.size(av)==1:
                        pos=peak[av,:][0]
                    else:
                        pook=peak[av,:][0]
                        bv=find(np.abs(pook[:,1]),np.abs(pook[:,1])==np.amin(np.abs(pook[:,1])))
                        if np.size(bv)==1:
                            pos=pook[bv,:][0]
                        else:
                            piik=pook[bv,:][0]
                            amp1=amp[av]
                            amp2=amp1[bv]
                            dif=np.abs(amp[bv]-data1[ic,jc])
                            difa=find(dif,dif==np.amin(dif))
                            pos= piik[difa,:][0]
                pos=pos[0]                        
    else:
        pos=np.array([100,100])
        indi=1
        nz=data2[ic-4:ic+5,jc-5:jc+6]
        val=data1[ic,jc]
        nm=np.size(nz[:,1])
        nn=np.size(nz[1,:])
        zios=np.zeros([nm,nn])
        for ni in range(1,nm-2):
            for nj in range(1,nn-2):
                nzi=nz[ni-1:ni+2,nj-1:nj+2]
                if (nzi[1,1]>val/10):
                        difo=find(nzi,nzi[1,1]>nzi)
                        zios[ni,nj]=np.size(difo[1])
        #print zios    
        if np.amax(zios)>5:
            hv=find(zios,zios>5)
            #print hv[0][0],hv[1][0]
            crit=np.abs(nz-val)
            pos=[hv[0][0],hv[1][0]]
            for ll in range(1,np.size(hv[1])-1):
                if crit[hv[0][ll],hv[1][ll]]<crit[pos[0],pos[1]]:
                    pos=[hv[0][ll],hv[1][ll]]   
            tt=nz[pos[0]-1:pos[0]+2,pos[1]-1:pos[1]+2]
            im=find(tt,tt==np.amax(tt))
            pos[0]=pos[0]+im[0][0]-1
            pos[1]=pos[1]+im[1][0]-1
            pos=np.array([4,5])-pos
                         
    return pos,indi
def subpix2(z,ii,jj):
    trange = np.arange(7)
    ttrange = np.arange(7)
    X,Y = np.meshgrid(trange,ttrange)
    outgrid = interpolate.interp2d(X,Y,z,kind='quintic')
    xx=yy=np.arange(61)/10.
    l=outgrid(xx,yy)
    l=l[30-9:30+10,30-9:30+10]
    ind=find(l,l==np.amax(l))
    #print l
    #print ind[0][0],ind[1][0]
    ni=ii+(ind[0][0]-9.)/10
    nj=jj+(ind[1][0]-9.)/10
    #print ii,jj
    #print ni,nj
    return[ni,nj]
def dephc(z):
    e = lambda v,z,: np.sum(np.abs(z-z[5,5]*expon(11,v)),1)
    vi=[1,1]
    #z[z<0]=0
    v, success = leastsq(e, vi, args=(z), maxfev=1000)
    cond='g'
    if v[0]<0.1 or v[0]>4 or v[1]<0.1 or v[1]>4 :
        cond='l'
        e = lambda v,z,: np.sum(np.abs(z-z[9,9]*cauchy(19,v)),1)
        vi=[0.3,0.3]    
        v, success = leastsq(e, vi, args=(z), maxfev=1000)
        if v[0]<0.001 or v[0]>2 or v[1]<0.001 or v[1]>2 :
            v[0]=v[1]=0.3+np.random.normal(0, 0.05, 1)
    return v,cond

def dephcl(z):
    e = lambda v,z,: np.sum(np.abs(z-z[9,9]*cauchy(19,v)),1)
    vi=[0.3,0.3]    
    v, success = leastsq(e, vi, args=(z), maxfev=1000)
    if v[0]<0.001 or v[0]>2 or v[1]<0.001 or v[1]>2 :
        v[0]=v[1]=0.3+np.random.normal(0, 0.05, 1)
    return v

def dephcg(z):
    e = lambda v,z,: np.sum(np.abs(z-z[9,9]*expon(19,v)),1)
    vi=[1,1]
    #z[z<0]=0
    v, success = leastsq(e, vi, args=(z), maxfev=1000)
    if v[0]<0.1 or v[0]>4 or v[1]<0.1 or v[1]>4 :
        v[0]=v[1]=2+np.random.normal(0, 0.05, 1)

    return v

def dephcaprio(z,a,b,c):
    if c=='g':
        e = lambda v,z,: np.sum(np.abs(z-z[9,9]*expon(19,v)),1)
        vi=[a,b]
        #z[z<0]=0
        v, success = leastsq(e, vi, args=(z), maxfev=1000)
        if np.abs(float(v[0]-a))>1: 
            v[0]=a+np.random.normal(0, 0.05, 1)
        if np.abs(float(v[1]-b))>1:    
            v[1]=b+np.random.normal(0, 0.05, 1)
    else:
        e = lambda v,z,: np.sum(np.abs(z-z[9,9]*cauchy(19,v)),1)
        vi=[a,b]    
        v, success = leastsq(e, vi, args=(z), maxfev=1000)
        if np.abs(float(v[0]-a))>0.5 or v[0]<0.08 or v[0]>8:
            v[0]=a+np.random.normal(0, 0.05, 1)
        if  np.abs(float(v[1]-b))>0.5 or v[1]<0.08:
            v[1]=b+np.random.normal(0, 0.05, 1)
    #print c                                        
    return v,c

def dephcaprio1(z,a,b,c):
    if c=='g':
        e = lambda v,z,: np.sum(np.abs(z-z[9,9]*expon(19,v)),1)
        vi=[a,b]
        #z[z<0]=0
        v, success = leastsq(e, vi, args=(z), maxfev=1000)
        if np.abs(float(v[0]-a))>1: 
            v[0]=a+np.random.normal(0, 0.05, 1)
        if np.abs(float(v[1]-b))>1:    
            v[1]=b+np.random.normal(0, 0.05, 1)
    else:
        z[z<0]=0
        e = lambda v,z,: np.sum(np.abs(z-z[9,9]*exponcauchy(19,v)),1)
        vi=[2,0.3]    
        v, success = leastsq(e, vi, args=(z), maxfev=1000)
        #if np.abs(float(v[0]-a))>0.5 or v[0]<0.08:
            #v[0]=a+np.random.normal(0, 0.05, 1)
        #if  np.abs(float(v[1]-b))>0.5 or v[1]<0.08:
            #v[1]=b+np.random.normal(0, 0.05, 1)
        v[1]=1/v[1]
    #print c                                        
    return v,c


def exp_hand(z,newn):
    c=0
    for i in range(len(z)):
        if z[i]==newn:
            c=1
            break
    return c
def exp_yn(amp,ampref,test):
    artest=np.array(test)
    ol=np.nonzero(artest==1)
    #print np.size(ol)
    #print len(amp)/2
    #print amp
    if np.size(ol)>len(amp)/2:
        ver=0
    else:
        o=np.nonzero(artest==0)
        vamp=np.array(amp)
        vampref=np.array(ampref)
        ivamref=vampref[o]
        ivam=vamp[o]
        ii=ivamref[ivamref==np.amax(ivamref)]
        jj=ivam[ivamref==np.amax(ivamref)]
        ver=1
        if np.size(ol)>1:
            ol=ol[0]
            #print np.size(ol)
        #print vampref
        #print vamp
        #print ol
        if np.size(ol)>1:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj))) 
                if (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))>50 or (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))<0.0200:
                    #print 'lela'
                    ver=0
        else:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0]) <0.001
                if (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])>50 or (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])<0.0200:
                    #print 'lela'
                    ver=0
            
    return ver
def exp_ync(amp,ampref,test):
    artest=np.array(test)
    ol=np.nonzero(artest==1)
    #print np.size(ol)
    #print len(amp)/2
    if np.size(ol)>len(amp)/2:
        ver=0
    else:
        o=np.nonzero(artest==0)
        vamp=np.array(amp)
        vampref=np.array(ampref)
        ivamref=vampref[o]
        ivam=vamp[o]
        ii=ivamref[ivamref==np.amax(ivamref)]
        jj=ivam[ivamref==np.amax(ivamref)]
        ver=1
        if np.size(ol)>1:
            ol=ol[0]
            #print np.size(ol)
        #print vampref
        #print vamp
        #print ol
        if np.size(ol)>1:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj))) 
                if (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))>10 or (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))<0.065:
                    #print 'lela'
                    ver=0
        else:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0]) 
                if (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])>10 or (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])<0.065:
                    #print 'lela'
                    ver=0
            
    return ver
        
def ser_mat(tabp):    
    boc=-1
    inp=[]
    inm=[]
    while boc<len(tabp)-1:
        boc+=1
        a=str(tabp[boc][2])
        newn=''
        for j in range(len(a)):
            if (a[j]=='_')==1:
                break
            newn+=a[j]
        
        r1=[]
        r1.append((boc))
        for jj in range(boc+1,len(tabp)):
            nomn=str(tabp[jj][2])
            #print nomn[0:j]
            try:
                if nomn[0:j+1]==newn+'_':
                    r1.append((jj))
                    #print 'ok'
                else:
                    break
            except:
                    break
        boc=boc+len(r1)-1
        inm.append(newn)
        inp.append(r1)
    return inm,inp

def tab2tab(X):
    Z=[]
    for i in range(len(X)):
        Z.append([X[i]])
    Z=np.array(Z)
    return Z
    