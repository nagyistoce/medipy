# This file contains code from DiPy (http://nipy.sourceforge.net/dipy/index.html),
# released under the following license.
#
# Copyright (c) 2009-2010, dipy developers
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.
# 
#     * Neither the name of the dipy developers nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cimport cython

import time
import numpy as np
cimport numpy as cnp

cdef extern from "math.h" nogil:
    float sqrt(float x) 
    float fabs(float x) 
      
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
cdef cnp.float32_t inf = np.inf


def local_skeleton_clustering(tracks, d_thr=10):
    ''' Used in the HBM2010 abstract
    "Fast Dimensionality Reduction for Brain Tractography Clustering" by E.Garyfallidis et.al

    Parameters :
    
    * tracks: sequence of tracks as arrays, shape (N1,3) .. (Nm,3)
    * d_thr: float, average euclidean distance threshold

    Return a dictionary with one entry per cluster. Each value contains the 
    indices of the fibers in the cluster (``"indices"``) and the number of fibers
    in the cluster (``"indices"``).
    '''

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

def most_similar_track_mam(tracks,metric='avg'):    
    ''' Find the most similar track in a bundle 
    using distances calculated from Zhang et. al 2008. 
    
    Parameters
    ----------
    tracks : sequence 
       of tracks as arrays, shape (N1,3) .. (Nm,3)
    metric : str
       'avg', 'min', 'max'
        
    Returns
    -------
    si : int
       index of the most similar track in tracks. This can be used as a
       reference track for a bundle.
    s : array, shape (len(tracks),)
        similarities between tracks[si] and the rest of the tracks in
        the bundle
    
    Notes
    -----
    A vague description of this function is given below:
    
    for (i,j) in tracks_combinations_of_2:

        calculate the mean_closest_distance from i to j  (mcd_i)  
        calculate the mean_closest_distance from j to i  (mcd_j)
        
        if 'avg':
            s holds the average similarities
        if 'min':
            s holds the minimum similarities
        if 'max':
            s holds the maximum similarities
        
    si holds the index of the track with min {avg,min,max} average metric
    '''
    cdef:
        size_t i, j, lent
        int metric_type
    if metric=='avg':
        metric_type = 0
    elif metric == 'min':
        metric_type = 1
    elif metric == 'max':
        metric_type = 2
    else:
        raise ValueError('Metric should be one of avg, min, max')
    # preprocess tracks
    cdef:
        size_t longest_track_len = 0, track_len
        cnp.ndarray[object, ndim=1] tracks32
    lent = len(tracks)
    tracks32 = np.zeros((lent,), dtype=object)
    # process tracks to predictable memory layout, find longest track
    for i in range(lent):
        tracks32[i] = np.ascontiguousarray(tracks[i], dtype=f32_dt)
        track_len = tracks32[i].shape[0]
        if track_len > longest_track_len:
            longest_track_len = track_len
    # buffer for distances of found track to other tracks
    cdef:
        cnp.ndarray[cnp.double_t, ndim=1] track2others
    track2others = np.zeros((lent,), dtype=np.double)
    # use this buffer also for working space containing summed distances
    # of candidate track to all other tracks
    cdef cnp.double_t *sum_track2others = <cnp.double_t *>track2others.data
    # preallocate buffer array for track distance calculations
    cdef:
        cnp.ndarray [cnp.float32_t, ndim=1] distances_buffer
        cnp.float32_t *t1_ptr, *t2_ptr, *min_buffer, distance
    distances_buffer = np.zeros((longest_track_len*2,), dtype=np.float32)
    min_buffer = <cnp.float32_t *> distances_buffer.data
    # cycle over tracks
    cdef:
        cnp.ndarray [cnp.float32_t, ndim=2] t1, t2
        size_t t1_len, t2_len
    for i from 0 <= i < lent-1:
        t1 = tracks32[i]
        t1_len = t1.shape[0]
        t1_ptr = <cnp.float32_t *>t1.data
        for j from i+1 <= j < lent:
            t2 = tracks32[j]
            t2_len = t2.shape[0]
            t2_ptr = <cnp.float32_t *>t2.data
            distance = czhang(t1_len, t1_ptr, t2_len, t2_ptr, min_buffer, metric_type)
            # get metric
            sum_track2others[i]+=distance
            sum_track2others[j]+=distance
    # find track with smallest summed metric with other tracks
    cdef double mn = sum_track2others[0]
    cdef size_t si = 0
    for i in range(lent):
        if sum_track2others[i] < mn:
            si = i
            mn = sum_track2others[i]
    # recalculate distance of this track from the others
    t1 = tracks32[si]
    t1_len = t1.shape[0]
    t1_ptr = <cnp.float32_t *>t1.data
    for j from 0 <= j < lent:
        t2 = tracks32[j]
        t2_len = t2.shape[0]
        t2_ptr = <cnp.float32_t *>t2.data
        track2others[j] = czhang(t1_len, t1_ptr, t2_len, t2_ptr, min_buffer, metric_type)
    return si, track2others
	
	
	
@cython.cdivision(True)
cdef inline cnp.float32_t czhang(size_t t1_len,
                                 cnp.float32_t *track1_ptr,
                                 size_t t2_len,
                                 cnp.float32_t *track2_ptr,
                                 cnp.float32_t *min_buffer,
                                 int metric_type) nogil:
    ''' Note ``nogil`` - no python calls allowed in this function '''
    cdef:
        cnp.float32_t *min_t2t1, *min_t1t2
    min_t2t1 = min_buffer
    min_t1t2 = min_buffer + t2_len
    min_distances(t1_len, track1_ptr,
                  t2_len, track2_ptr,
                  min_t2t1,
                  min_t1t2)
    cdef:
        size_t t1_pi, t2_pi
        cnp.float32_t mean_t2t1 = 0, mean_t1t2 = 0, dist_val
    for t1_pi from 0<= t1_pi < t1_len:
        mean_t1t2+=min_t1t2[t1_pi]
    mean_t1t2=mean_t1t2/t1_len
    for t2_pi from 0<= t2_pi < t2_len:
        mean_t2t1+=min_t2t1[t2_pi]
    mean_t2t1=mean_t2t1/t2_len
    if metric_type == 0:                
        dist_val=(mean_t2t1+mean_t1t2)/2.0
    elif metric_type == 1:        
        if mean_t2t1 < mean_t1t2:
            dist_val=mean_t2t1
        else:
            dist_val=mean_t1t2
    elif metric_type == 2:                
        if mean_t2t1 > mean_t1t2:
            dist_val=mean_t2t1
        else:
            dist_val=mean_t1t2                    
    return dist_val

@cython.cdivision(True)
cdef inline void min_distances(size_t t1_len,
                               cnp.float32_t *track1_ptr,
                               size_t t2_len,
                               cnp.float32_t *track2_ptr,
                               cnp.float32_t *min_t2t1,
                               cnp.float32_t *min_t1t2) nogil:
    cdef:
        cnp.float32_t *t1_pt, *t2_pt, d0, d1, d2
        cnp.float32_t delta2
        int t1_pi, t2_pi
    for t2_pi from 0<= t2_pi < t2_len:
        min_t2t1[t2_pi] = inf
    for t1_pi from 0<= t1_pi < t1_len:
        min_t1t2[t1_pi] = inf
    # pointer to current point in track 1
    t1_pt = track1_ptr
    # calculate min squared distance between each point in the two
    # lines.  Squared distance to delay doing the sqrt until after this
    # speed-critical loop
    for t1_pi from 0<= t1_pi < t1_len:
        # pointer to current point in track 2
        t2_pt = track2_ptr
        for t2_pi from 0<= t2_pi < t2_len:
            d0 = t1_pt[0] - t2_pt[0]
            d1 = t1_pt[1] - t2_pt[1]
            d2 = t1_pt[2] - t2_pt[2]
            delta2 = d0*d0 + d1*d1 + d2*d2
            if delta2 < min_t2t1[t2_pi]:
                min_t2t1[t2_pi]=delta2
            if delta2 < min_t1t2[t1_pi]:
                min_t1t2[t1_pi]=delta2
            t2_pt += 3 # to next point in track 2
        t1_pt += 3 # to next point in track 1
    # sqrt to get Euclidean distance from squared distance
    for t1_pi from 0<= t1_pi < t1_len:
        min_t1t2[t1_pi]=sqrt(min_t1t2[t1_pi])
    for t2_pi from 0<= t2_pi < t2_len:
        min_t2t1[t2_pi]=sqrt(min_t2t1[t2_pi])

def cut_plane(tracks,ref):
    ''' Extract divergence vectors and points of intersection 
    between planes normal to the reference fiber and other tracks
    
    Parameters
    ----------
    tracks: sequence
        of tracks as arrays, shape (N1,3) .. (Nm,3)
    ref: array, shape (N,3)
        reference track
        
    Returns
    -------
    hits: sequence
       list of points and rcds (radial coefficient of divergence)``
    
    Examples
    --------
    >>> from tools import cut_plane
    >>> refx = np.array([[0,0,0],[1,0,0],[2,0,0],[3,0,0]],dtype='float32')
    >>> bundlex = [np.array([[0.5,1,0],[1.5,2,0],[2.5,3,0]],dtype='float32')]
    >>> cut_plane(bundlex,refx)
    [array([[ 1.        ,  1.5       ,  0.        ,  0.70710683, 0,]], dtype=float32),
     array([[ 2.        ,  2.5       ,  0.        ,  0.70710677, 0.]], dtype=float32)]
        
    The orthogonality relationship
    np.inner(hits[p][q][0:3]-ref[p+1],ref[p+2]-ref[r][p+1])
    will hold throughout for every point q in the hits plane
    at point (p+1) on the reference track.
    '''
    cdef:
        size_t n_hits, hit_no, max_hit_len
        float alpha,beta,lrq,rcd,lhp,ld
        cnp.ndarray[cnp.float32_t, ndim=2] ref32
        cnp.ndarray[cnp.float32_t, ndim=2] track
        object hits
        cnp.ndarray[cnp.float32_t, ndim=1] one_hit
        float *hit_ptr
        cnp.ndarray[cnp.float32_t, ndim=2] hit_arr
        object Hit=[]
    # make reference fiber usable type
    ref32 = np.ascontiguousarray(ref, f32_dt)
    # convert all the tracks to something we can work with.  Get track
    # lengths
    cdef:
        size_t N_tracks=len(tracks)
        cnp.ndarray[cnp.uint64_t, ndim=1] track_lengths
        size_t t_no, N_track
    cdef object tracks32 = []
    track_lengths = np.empty((N_tracks,), dtype=np.uint64)
    for t_no in range(N_tracks):
        track = np.ascontiguousarray(tracks[t_no], f32_dt)
        track_lengths[t_no] = track.shape[0]
        tracks32.append(track)
    # set up loop across reference fiber points
    cdef:
        size_t N_ref = ref32.shape[0]
        size_t p_no, q_no
        float *this_ref_p, *next_ref_p, *this_trk_p, *next_trk_p
        float along[3], normal[3]
        float qMp[3], rMp[3], rMq[3], pMq[3]
        float hit[3], hitMp[3], *delta
    # List used for storage of hits.  We will fill this with lots of
    # small numpy arrays, and reuse them over the reference track point
    # loops.
    max_hit_len = 0
    hits = []
    # for every point along the reference track
    next_ref_p = asfp(ref32[0])
    for p_no in range(N_ref-1):
        # extract point to point vector into `along`
        this_ref_p = next_ref_p
        next_ref_p = asfp(ref32[p_no+1])
        csub_3vecs(next_ref_p, this_ref_p, along)
        # normalize
        cnormalized_3vec(along, normal)
        # initialize index for hits
        hit_no = 0
        # for every track
        for t_no in range(N_tracks):
            track=tracks32[t_no]
            N_track = track_lengths[t_no]
            # for every point on the track
            next_trk_p = asfp(track[0])
            for q_no in range(N_track-1):
                # p = ref32[p_no]
                # q = track[q_no]
                # r = track[q_no+1]
                # float* versions of above: p == this_ref_p
                this_trk_p = next_trk_p # q
                next_trk_p = asfp(track[q_no+1]) # r
                #if np.inner(normal,q-p)*np.inner(normal,r-p) <= 0:
                csub_3vecs(this_trk_p, this_ref_p, qMp) # q-p
                csub_3vecs(next_trk_p, this_ref_p, rMp) # r-p
                if (cinner_3vecs(normal, qMp) * cinner_3vecs(normal, rMp)) <=0:
                    #if np.inner((r-q),normal) != 0:
                    csub_3vecs(next_trk_p, this_trk_p, rMq)
                    beta = cinner_3vecs(rMq, normal)
                    if beta !=0:
                        #alpha = np.inner((p-q),normal)/np.inner((r-q),normal)
                        csub_3vecs(this_ref_p, this_trk_p, pMq)
                        alpha = (cinner_3vecs(pMq, normal) /
                                  cinner_3vecs(rMq, normal))
                        if alpha < 1:
                            # hit = q+alpha*(r-q)
                            hit[0] = this_trk_p[0]+alpha*rMq[0]
                            hit[1] = this_trk_p[1]+alpha*rMq[1]
                            hit[2] = this_trk_p[2]+alpha*rMq[2]
                            # h-p
                            csub_3vecs(hit, this_ref_p, hitMp)
                            # |h-p|
                            lhp = cnorm_3vec(hitMp)
                            delta = rMq # just renaming
                            # |r-q| == |delta|
                            ld = cnorm_3vec(delta)
                            ''' # Summary of stuff in comments
                            # divergence =((r-q)-inner(r-q,normal)*normal)/|r-q|
                            div[0] = (rMq[0]-beta*normal[0]) / ld
                            div[1] = (rMq[1]-beta*normal[1]) / ld
                            div[2] = (rMq[2]-beta*normal[2]) / ld
                            # radial coefficient of divergence d.(h-p)/|h-p|
                            '''
                            # radial divergence
                            # np.inner(delta, (hit-p)) / (ld * lhp)
                            if lhp > 0:
                                rcd = fabs(cinner_3vecs(delta, hitMp)
                                           / (ld*lhp))
                            else:
                                rcd=0
                            # hit data into array
                            if hit_no >= max_hit_len:
                                one_hit = np.empty((5,), dtype=f32_dt)
                                hits.append(one_hit)
                            else:
                                one_hit = hits[hit_no]
                            hit_ptr = <float *>one_hit.data
                            hit_ptr[0] = hit[0]
                            hit_ptr[1] = hit[1]
                            hit_ptr[2] = hit[2]
                            hit_ptr[3] = rcd
                            hit_ptr[4] = t_no
                            hit_no += 1
        # convert hits list to hits array
        n_hits = hit_no
        if n_hits > max_hit_len:
            max_hit_len = n_hits
        hit_arr = np.empty((n_hits,5), dtype=f32_dt)
        for hit_no in range(n_hits):
            hit_arr[hit_no] = hits[hit_no]
        Hit.append(hit_arr)
        #Div.append(divs[1:])
    return Hit[1:] 

cdef inline float cnorm_3vec(float *vec):
    ''' Calculate Euclidean norm of input vector

    Parameters
    ----------
    vec : float *
       length 3 float vector

    Returns
    -------
    norm : float
       Euclidean norm
    '''
    cdef float v0, v1, v2
    v0 = vec[0]
    v1 = vec[1]
    v2 = vec[2]
    return sqrt(v0 * v0 + v1*v1 + v2*v2) 

cdef inline float cinner_3vecs(float *vec1, float *vec2) nogil:
    cdef int i
    cdef float ip = 0
    for i from 0<=i<3:
        ip += vec1[i]*vec2[i]
    return ip

cdef inline void csub_3vecs(float *vec1, float *vec2, float *vec_out) nogil:
    cdef int i
    for i from 0<=i<3:
        vec_out[i] = vec1[i]-vec2[i]

cdef inline void cnormalized_3vec(float *vec_in, float *vec_out):
    ''' Calculate and fill normalized 3D vector 

    Parameters
    ----------
    vec_in : float *
       Length 3 vector to normalize
    vec_out : float *
       Memory into which to write normalized length 3 vector

    Returns
    -------
    void
    '''
    cdef float norm = cnorm_3vec(vec_in)
    cdef int i
    norm=norm+(norm==0)
    for i in range(3):
        vec_out[i] = vec_in[i] / norm

def mam_distances(xyz1,xyz2,t,metric='all'):
    ''' Min/Max/Mean Average Minimume Distance between tracks xyz1 and xyz2
    
    Based on the metrics in Zhang, Correia, Laidlaw 2008
    http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=4479455
    which in turn are based on those of Corouge et al. 2004
    
    Parameters
    ----------
    xyz1 : array, shape (N1,3), dtype float32
    xyz2 : array, shape (N2,3), dtype float32
       arrays representing x,y,z of the N1 and N2 points of two tracks
	t : float, minimum threshold so that distances below it are not considered
    metrics : {'avg','min','max','all'}
       Metric to calculate.  {'avg','min','max'} return a scalar. 'all'
       returns a tuple
       
    Returns
    -------
    avg_mcd: float
       average_mean_closest_distance
    min_mcd: float
       minimum_mean_closest_distance
    max_mcd: float
       maximum_mean_closest_distance
                    
    Notes
    -----
    Algorithmic description
    
    Lets say we have curves A and B.
    
    For every point in A calculate the minimum distance from every point
    in B stored in minAB
    
    For every point in B calculate the minimum distance from every point
    in A stored in minBA
    
    find average of minAB stored as avg_minAB
    find average of minBA stored as avg_minBA
    
    if metric is 'avg' then return (avg_minAB + avg_minBA)/2.0
    if metric is 'min' then return min(avg_minAB,avg_minBA)
    if metric is 'max' then return max(avg_minAB,avg_minBA)
    '''
    cdef:
        cnp.ndarray[cnp.float32_t, ndim=2] track1 
        cnp.ndarray[cnp.float32_t, ndim=2] track2
        size_t t1_len, t2_len
    track1 = np.ascontiguousarray(xyz1, dtype=f32_dt)
    t1_len = track1.shape[0]
    track2 = np.ascontiguousarray(xyz2, dtype=f32_dt)
    t2_len = track2.shape[0]
    # preallocate buffer array for track distance calculations
    cdef:
        cnp.float32_t *min_t2t1, *min_t1t2
        cnp.ndarray [cnp.float32_t, ndim=1] distances_buffer
    distances_buffer = np.zeros((t1_len + t2_len,), dtype=np.float32)
    min_t2t1 = <cnp.float32_t *> distances_buffer.data
    min_t1t2 = min_t2t1 + t2_len
    min_distances(t1_len, <cnp.float32_t *>track1.data,
                  t2_len, <cnp.float32_t *>track2.data,
                  min_t2t1,
                  min_t1t2)
    cdef:
        size_t t1_pi, t2_pi
        cnp.float32_t mean_t2t1 = 0, mean_t1t2 = 0
    for t1_pi from 0<= t1_pi < t1_len:
        if min_t1t2[t1_pi]>t:
            mean_t1t2+=min_t1t2[t1_pi]
    mean_t1t2=mean_t1t2/t1_len
    for t2_pi from 0<= t2_pi < t2_len:
        if min_t2t1[t2_pi]>t:
            mean_t2t1+=min_t2t1[t2_pi]
    mean_t2t1=mean_t2t1/t2_len
    if metric=='all':
        return ((mean_t2t1+mean_t1t2)/2.0,
                np.min((mean_t2t1,mean_t1t2)),
                np.max((mean_t2t1,mean_t1t2)))
    elif metric=='avg':
        return (mean_t2t1+mean_t1t2)/2.0
    elif metric=='min':            
        return np.min((mean_t2t1,mean_t1t2))
    elif metric =='max':
        return np.max((mean_t2t1,mean_t1t2))
    else :
        ValueError('Wrong argument for metric')
