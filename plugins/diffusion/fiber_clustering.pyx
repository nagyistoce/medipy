# from dipy

cimport cython

import time
import numpy as np
cimport numpy as cnp

cdef extern from "math.h" nogil:
    float sqrt(float x)  
      
cdef extern from "stdlib.h" nogil:
    ctypedef unsigned long size_t
    void free(void *ptr)
    void *malloc(size_t size)
    void *calloc(size_t nelem, size_t elsize)
    void *realloc (void *ptr, size_t size)
    void *memcpy(void *str1, void *str2, size_t n)

cdef inline float* asfp(cnp.ndarray pt):
    return <float *>pt.data

# float 32 dtype for casting
cdef cnp.dtype f32_dt = np.dtype(np.float32)


def local_skeleton_clustering(tracks, d_thr=10):
    ''' Used in the HBM2010 abstract
    "Fast Dimensionality Reduction for Brain Tractography Clustering" by E.Garyfallidis et.al

    Parameters
    -----------
    tracks: sequence
        of tracks as arrays, shape (N1,3) .. (Nm,3)

    d_thr: float, average euclidean distance threshold

    Returns
    --------
    C: dict '''

    cdef :
        cnp.ndarray[cnp.float32_t, ndim=2] track
        cnp.ndarray[cnp.float32_t, ndim=2] h
        int lent,k,it
        float d[2]

    tracks_3 = [downsample(t,3) for t in tracks]    
    lent = len(tracks_3)

    C = {0:{'indices':[0],'hidden':tracks_3[0].copy(),'N':1}}
    ts=np.zeros((3,3),dtype=np.float32)
    
    for it from 1 <= it < lent by 1:      
        track = np.ascontiguousarray(tracks_3[it],dtype=f32_dt)            
        lenC = len(C.keys())       
        #if it%1000==0:
        #    print it,lenC        
        alld = np.zeros(lenC)
        flip = np.zeros(lenC)
        
        for k from 0 <= k < lenC by 1:         
            h = np.ascontiguousarray(C[k]['hidden']/C[k]['N'],dtype=f32_dt)
            track_direct_flip_3dist(asfp(track[0]),asfp(track[1]),asfp(track[2]),asfp(h[0]),asfp(h[1]),asfp(h[2]),d)            
            if d[1]<d[0]:                
                d[0] = d[1]
                flip[k] = 1               
            alld[k] = d[0]

        m_k = np.min(alld)
        i_k = np.argmin(alld)
        
        if m_k<d_thr:                       
            if flip[i_k]==1:                
                ts[0] = track[-1];ts[1]=track[1];ts[-1]=track[0]
                C[i_k]['hidden']+=ts
            else:                
                C[i_k]['hidden']+=track               
            C[i_k]['N']+=1
            C[i_k]['indices'].append(it)           
        else:
            C[lenC] = {}
            C[lenC]['hidden']= track.copy()
            C[lenC]['N'] = 1
            C[lenC]['indices'] = [it]
    
    return C

cdef inline void track_direct_flip_3dist(float *a1, float *b1,float  *c1,float *a2, float *b2, float *c2, float *out):
    ''' Calculate the euclidean distance between two 3pt tracks
    both direct and flip are given as output
    
    
    Parameters
    ----------------
    a1,b1,c1: 3 float[3] arrays representing the first track
    a2,b2,c2: 3 float[3] arrays representing the second track
    
    Returns
    -----------
    out: a float[2] array having the euclidean distance and the fliped euclidean distance '''
    
    cdef:
        int i
        float tmp1=0,tmp2=0,tmp3=0,tmp1f=0,tmp3f=0
           
    for i from 1 <= i < 3 by 1:
        tmp1 = tmp1 + (a1[i]-a2[i])*(a1[i]-a2[i])
        tmp2 = tmp2 + (b1[i]-b2[i])*(b1[i]-b2[i])
        tmp3 = tmp3 + (c1[i]-c2[i])*(c1[i]-c2[i])
        tmp1f = tmp1f + (a1[i]-c2[i])*(a1[i]-c2[i])
        tmp3f = tmp3f + (c1[i]-a2[i])*(c1[i]-a2[i])
                
    out[0] = (sqrt(tmp1)+sqrt(tmp2)+sqrt(tmp3))/3.0
    out[1] = (sqrt(tmp1f)+sqrt(tmp2)+sqrt(tmp3f))/3.0

def downsample(xyz,n_pols=3):
    ''' downsample for a specific number of points along the curve

    Uses the length of the curve. It works in as similar fashion to
    midpoint and arbitrarypoint.
    
    Parameters
    ----------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track
    n_pol : int
       integer representing number of points (poles) we need along the curve.

    Returns
    -------
    xyz2 : array shape (M,3)
       array representing x,z,z of M points that where extrapolated. M
       should be equal to n_pols  '''

    xyz = np.asarray(xyz)
    n_pts = xyz.shape[0]
    if n_pts == 0:
        raise ValueError('xyz array cannot be empty')
    if n_pts == 1:
        return xyz.copy().squeeze()
    cumlen = np.zeros(n_pts)
    cumlen[1:] = length(xyz, along=True)    
    step=cumlen[-1]/(n_pols-1)
    if cumlen[-1]<step:
        raise ValueError('Given number of points n_pols is incorrect. ')
    if n_pols<=2:
        raise ValueError('Given number of points n_pols needs to be'
                         ' higher than 2. ')
    xyz2=[_extrap(xyz,cumlen,distance) for distance in np.arange(0,cumlen[-1],step)]
    return np.vstack((np.array(xyz2),xyz[-1]))

def _extrap(xyz,cumlen,distance):
    ''' Helper function for extrapolate    
    '''    
    ind=np.where((cumlen-distance)>0)[0][0]
    len0=cumlen[ind-1]        
    len1=cumlen[ind]
    Ds=distance-len0
    Lambda = Ds/(len1-len0)
    return Lambda*xyz[ind]+(1-Lambda)*xyz[ind-1]

def length(xyz, along=False):
    ''' Euclidean length of track line

    Parameters
    ----------
    xyz : array-like shape (N,3)
       array representing x,y,z of N points in a track
    along : bool, optional
       If True, return array giving cumulative length along track,
       otherwise (default) return scalar giving total length.

    Returns
    -------
    L : scalar or array shape (N-1,)
       scalar in case of `along` == False, giving total length, array if
       `along` == True, giving cumulative lengths. '''

    xyz = np.asarray(xyz)
    if xyz.shape[0] < 2:
        if along:
            return np.array([0])
        return 0
    dists = np.sqrt((np.diff(xyz, axis=0)**2).sum(axis=1))
    if along:
        return np.cumsum(dists)
    return np.sum(dists)



        
