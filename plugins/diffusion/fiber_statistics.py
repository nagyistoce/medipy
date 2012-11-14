
##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import medipy.base
import numpy as np
import medipy.diffusion.fiber_clustering
import sys


def compute_statistic_section(dt1,dt2,C_hits_1,C_hits_2,nb_plane=1,nbClust=27):
    ''' Statistic in bundle cross sections
    '''

    spm = []
    dim1 = dt1.shape
    dim2 = dt2.shape
    spacing1 = dt1.spacing
    spacing2 = dt2.spacing

    for c1,c2 in zip(C_hits_1.keys(),C_hits_2.keys()):
        sys.stdout.write("\riteration= %d, %d / %d, %d" %(c1,c2,len(C_hits_1.keys()),len(C_hits_2.keys())))
        sys.stdout.flush()

        population1 = C_hits_1[c1]['indices']
        population2 = C_hits_2[c2]['indices']
        skel1 = C_hits_1[c1]['skel']
        skel2 = C_hits_2[c2]['skel']

        T = []
        shift_plane = nb_plane-1
        if len(population1)==len(population2):
            l_pop = len(population1)
        else:
            medipy.base.Exception("Can only compare two populations of the same size")    

        # impossible to get neighborhood
        for cnt1 in range(shift_plane):
            T.append( (0.,3.,1.) )

        # iteration
        for cnt1 in range(shift_plane,l_pop-shift_plane):

            # get neighborhood
            points1 = np.zeros((0,3))
            weights1 = np.zeros((0,))
            points2 = np.zeros((0,3))
            weights2 = np.zeros((0,))
            for cnt2 in range(cnt1-shift_plane,cnt1+nb_plane):
	            pop1_tmp = []
	            for point in population1[cnt2]:
		            k = (point<dim1).all()
		            if k:
			            pop1_tmp.append( point )
	            if pop1_tmp!=[]:
		            pop1_tmp = np.asarray( pop1_tmp )
		            points1 = np.concatenate( (points1,pop1_tmp) )
		            weights1 = np.concatenate( (weights1,np.repeat(np.abs(cnt2-cnt1),pop1_tmp.shape[0])) )

	            pop2_tmp = []
	            for point in population2[cnt2]:
		            k = (point<dim2).all()
		            if k:
			            pop2_tmp.append( point )
	            if pop2_tmp!=[]:
		            pop2_tmp = np.asarray( pop2_tmp )
		            points2 = np.concatenate( (points2,pop2_tmp) )
		            weights2 = np.concatenate( (weights2,np.repeat(np.abs(cnt2-cnt1),pop2_tmp.shape[0])) )

            # reduce data set
            if points1.shape[0]>0:
	            center1 = skel1[cnt1]
	            points1 = cluster_voxel_2(points1,center1,weights1,spacing1)
            if points2.shape[0]>0:
	            center2 = skel2[cnt1]
	            points2 = cluster_voxel_2(points2,center2,weights2,spacing2)

            # statistic
            N1 = float(points1.shape[0])
            N2 = float(points2.shape[0])
            if N1>=nbClust and N2>=nbClust:
	            N1 = float(nbClust)
	            points1 = points1[:nbClust]
	            N2 = float(nbClust)
	            points2 = points2[:nbClust]
            else:
	            N1 = 1.
	            N2 = 1.

            N = N1
            if N==nbClust:

	            # interpolation 
	            t1 = []
	            for point in points1:
		            tmp = point/spacing1
		            #print tmp, shape
		            t1.append( dt1[int(tmp[0]),int(tmp[1]),int(tmp[2])] )
	            t1 = np.asarray(t1).T
	            t2 = []
	            for point in points2:
		            tmp = point/spacing2
		            t2.append( dt2[int(tmp[0]),int(tmp[1]),int(tmp[2])] )
	            t2 = np.asarray(t2).T

	            # compute parameters
	            m1 = np.mean(t1,1)
	            m2 = np.mean(t2,1)
	            d1 = t1 - np.reshape( np.repeat(m1,N1),(6,N1) )
	            d2 = t2 - np.reshape( np.repeat(m2,N2),(6,N2) )
	            s1 = np.sum( np.sum(d1[:3,:]**2,0) + 2*np.sum(d1[3:,:]**2,0),0 )
	            s2 = np.sum( np.sum(d2[:3,:]**2,0) + 2*np.sum(d2[3:,:]**2,0),0 )
	            s1 = np.sqrt( s1/(6*(N1-1)) )
	            s2 = np.sqrt( s2/(6*(N2-1)) )

	            # statistic
	            m1_33 = np.array([[m1[0],m1[3],m1[4]],[m1[5],m1[1],m1[3]],[m1[4],m1[5],m1[2]]])
	            m2_33 = np.array([[m2[0],m2[3],m2[4]],[m2[5],m2[1],m2[3]],[m2[4],m2[5],m2[2]]])
	            L1,V1 = np.linalg.eigh(m1_33)
	            L2,V2 = np.linalg.eigh(m2_33)
	            stat = np.sum( (L1 - L2)**2 )
	            s = 0.5 * ( s1**2+s2**2 )
	            if s!=0 and (N-1)>0:
		            stat = ((N-1.)*stat)/(6.*s)
	            else:
		            stat = 0.
            else:
	            stat = 0.

            # normalize stat
            df1 = 3.
            df2 = 12.*(N-1.)

            T.append( (stat,df1,df2) )

        # impossible to get neighborhood
        for cnt1 in range(shift_plane):
            T.append( (0.,3.,1.) )

        spm.append(np.asarray(T))

    return spm



def compute_statistic_bundle(dt1,dt2,C_hits_1,C_hits_2,d_thr=4.):
    ''' Interpolate tensors (NN,L,...) and consider a paired statistic on bundles
    '''

    spm = []
    dim1 = dt1.shape
    dim2 = dt2.shape
    spacing1 = dt1.spacing
    spacing2 = dt2.spacing

    for c1,c2 in zip(C_hits_1.keys(),C_hits_2.keys()):
        sys.stdout.write("\riteration= %d, %d / %d, %d" %(c1,c2,len(C_hits_1.keys()),len(C_hits_2.keys())))
        sys.stdout.flush()

        # ref1 used to discard registration errors
        ref1 = C_hits_1[c1]['skel']
        population1 = C_hits_1[c1]['indices']
        population2 = C_hits_2[c2]['indices']
        T = []

        if len(population1)==len(population2):
            N = len(population1)
        else:
            raise medipy.base.Exception("Can only compare two populations of the same size")   

        # density threshold
        density = 0.
        for pnts in population2:
            density+=len(pnts)/N

        # cluster & interpolation
        skel1 = []
        skel2 = []
        if density>d_thr:
            kk = True
        else:
            kk = False

        if kk:
            for cnt in range(N):
	            pnts1 = population1[cnt]
	            pnts2 = population2[cnt]
	            pntref = np.cast[int](np.round(ref1[cnt]/spacing1))     
	            pop1_tmp = []
	            pop2_tmp = []
	            if not (pntref<dim).all():
		            kk = False
		            break
	            for point in pnts1:
		            if (point<dim).all():
			            pop1_tmp.append( dt1[point[0],point[1],point[2]] )
	            for point in pnts2:
		            if (point<dim).all():
			            pop2_tmp.append( dt2[point[0],point[1],point[2]] )

	            pop1_tmp = np.asarray(pop1_tmp).T
	            pop2_tmp = np.asarray(pop2_tmp).T
	            skel1.append( np.mean(pop1_tmp,1) )
	            skel2.append( np.mean(pop2_tmp,1) )

        if kk:
            t1 = np.asarray(skel1).T
            t2 = np.asarray(skel2).T
            t1 = keep_eigenmatrice(t1)
            t2 = keep_eigenmatrice(t2)

            # compute parameters
            diff = t2-t1
            m = np.mean(diff,1)
            d = diff - np.reshape( np.repeat(m,N),(6,N) )
            s = np.sum( np.sum(d[:3,:]**2,0) + 2*np.sum(d[3:,:]**2,0),0 )
            s = np.sqrt( s/(6.*(float(N)-1.)) )

            # statistic
            m_33 = np.array([[m[0],m[3],m[4]],[m[5],m[1],m[3]],[m[4],m[5],m[2]]])
            L1,V1 = np.linalg.eigh(m_33)
            L2 = np.array([0.,0.,0.])
            stat = np.sum( (L1 - L2)**2 )
            s = s**2
            if s!=0 and (N-1)>0:
	            stat = ((float(N)-1.)*stat)/(3*s)
            else:
	            stat = 0.
        else:
            stat = 0.

        df1 = 3.
        df2 = 6.*(N-1.)
        for cnt in range(N):
            T.append( (stat,df1,df2) )

        spm.append(np.asarray(T))

    return spm


def keep_eigenmatrice(t):
    n1,n2 = t.shape
    res = []
    anti = np.zeros((3,))
    for cnt in range(n2):
        m = t[:,cnt]
        m_33 = np.array([[m[0],m[3],m[4]],[m[5],m[1],m[3]],[m[4],m[5],m[2]]])
        L,V = np.linalg.eigh(m_33)
        L = np.concatenate((L,anti))
        res.append( L )
    return np.asarray(res).T



def project_fiber_stat_nn(spm,C_hits,model):

    P = np.ones(model.shape,dtype=np.single)*(-1)
    dim = model.shape

    for c,t in zip(C_hits.keys(),spm):
        sys.stdout.write("\riteration= %d / %d" %(c,len(C_hits.keys())))
        sys.stdout.flush()

        population = C_hits[c]['indices']
        if t.ndim==1:
            T = t[:]
        else:
            T = t[:,0]

        for nb,points in enumerate(population):        
            for p in points:
                z,y,x = p
                if x>=0 and y>=0 and z>=0 and x<dim[2] and y<dim[1] and z<dim[0]:
                    cv = P[z,y,x]
                    nv = T[nb]
                    if nv>cv:
                        P[z,y,x] = nv

    return medipy.base.Image(data=P,spacing=model.spacing,origin=model.origin,direction=model.direction)




def bundle_cross_sections_voxel(hits,spacing,cutoff,pv_thr) :
    """ Compute bundle cross sections at a voxel level and reduce partial volume effects by thresholding the neighborhood
    """
    #spacing = spacing[::-1]

    N_hits = {}
    for cnt,key in enumerate(hits.keys()) :
	    sys.stdout.write("\riteration= %d / %d" %(cnt,len(hits.keys())))
	    sys.stdout.flush()

	    c_hits = hits[key]['indices']
	    skel = hits[key]['skel']
	    new_hits = []

	    for hit,center in zip(c_hits,skel) :
            # generate voxel neighborhood
		    d = np.sqrt( np.sum( (hit-center)**2,1 ) )
		    hit = hit[d<cutoff]
		    c = cluster_voxel(hit,spacing)

		    # remove partial volume effect voxels based on histogram
		    hist = np.zeros((c.shape[0],),dtype=int)
		    for point in hit:
			    point = np.round(point/spacing)
			    l = list(point==c)
			    l = [l[i].all() for i in range(len(l))]
			    indx = l.index(True)
			    hist[indx] += 1
		    dynamic = hist.max()-hist.min()
		    thr = np.floor(float(hist.min()) + float(dynamic)*(1.0-pv_thr))
		    new_hits.append( c[hist>=thr] )

	    N_hits[key]={'indices':new_hits,'skel':skel}

    return N_hits

def cluster_voxel(data,spacing):
    data_nn = np.round(data/spacing)
    data_nn = np.cast[int](data_nn)
    list1 = list(data_nn)
    list2 = [ list(y) for y in set([tuple(x) for x in list1]) ]
    list2 = [ np.asarray(x) for x in list2 ]
    return np.asarray(list2)

def cluster_voxel_2(data,center,weights,spacing):

    # sort by plane
    argmax = np.argsort(weights)
    weights = np.cast[int](weights[argmax])
    data = data[argmax]
    keys = np.unique(weights)
    points = np.zeros((0,3))
    for key in keys:
        idx = np.where(weights==key)
        tmp = data[idx]*spacing
        d = np.sqrt( np.sum( (tmp-center)**2,1 ) )
        argmax = np.argsort(d)
        tmp = tmp[argmax]
        points = np.concatenate( (points,tmp) )

    points = np.cast[int]( np.round(points/spacing) )

    list1 = list(points)
    list2 = [ list(y) for y in set([tuple(x) for x in list1]) ]
    list2 = [ np.asarray(x) for x in list2 ]
    list2 = np.asarray(list2)*spacing

    d = np.sqrt( np.sum( (list2-center)**2,1 ) )
    argmax = np.argsort(d)
    list2 = list2[argmax]

    return list2






def bundle_cross_sections(T,C,nb_samples=None):
    """ Compute bundle cross sections
    """

    C_hits={}
    for it,key in enumerate(C.keys()):
        sys.stdout.write("\riteration= %d / %d nb fibers= %d" %(it,len(C),C[key]['N']))
        sys.stdout.flush()

        N = C[key]['N']
        bundle_indices = C[key]['indices']
        indice_skel = C[key]['skel']

        B = []
        for indice in bundle_indices :
            B.append(T[indice])
        if nb_samples==None :
            skel = T[indice_skel]
        else :
            skel = medipy.diffusion.fiber_clustering.downsample(T[indice_skel],nb_samples)

        hits = cut_bundle(B,skel)

        C_hits[key] = {'indices':hits,'skel':skel}

    return C_hits

def cut_bundle(bundle,skel):
    """ Compute a bundle cut neighborhood
    """

    hits = medipy.diffusion.fiber_clustering.cut_plane(bundle,skel)
    hits_single = cut_bundle_single(hits,skel)

    for i in range(len(hits_single)):
        hits_single[i] = hits_single[i][:,:3]
    hits_single.insert(0,np.asarray((skel[0],)))
    hits_single.append(np.asarray((skel[-1],)))

    return hits_single

def cut_bundle_single(hits,skel):
    """ Remove multiple fiber crossing points using L2-norm
    """
      
    sk = skel[1:-1].copy()
    Nskel = sk.shape[0]
    Nhits = len(hits)
    if Nskel!=Nhits :
        raise IOError('shape mismatched')
	
    hits_single = []
    for hit,ref in zip(hits,sk) :
        # L2 intersection point distances to skeleton
        dt = np.sum( (hit[:,:3]-ref)**2,1 )
        n = np.where(dt==0)
        d = np.sqrt(dt)
        d[n] = 0
	
        # find duplicate
        C = {}
        Npoints = hit.shape[0]
        if Npoints> 0:
            C[hit[0,4]]={'indices':[0],'N':1,'v_min':d[0],'a_min':0}
            for i in range(1,Npoints) :
                if hit[i,4] in C.keys() :
                    C[hit[i,4]]['N'] += 1
                    C[hit[i,4]]['indices'].append(i)
                    if C[hit[i,4]]['v_min']>d[i] :
                        C[hit[i,4]]['v_min'] = d[i]
                        C[hit[i,4]]['a_min'] = i
                else :
                    C[hit[i,4]] = {'N': 1, 'indices': [i], 'v_min': d[i], 'a_min': i}

        # remove duplicate
        hit_single = []
        tmp = np.zeros((5,))
        tmp[:3] = ref
        hit_single.append(tmp)
        for c in C.keys():
            hit_single.append(hit[C[c]['a_min']])
		
        hits_single.append(np.asarray(hit_single))
	
    return hits_single




def cut_polydata(hits) : 
    """ Construct the vtk polydata 
    """

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(0.1)
    sphere.SetThetaResolution(8)
    sphere.SetPhiResolution(8)

    points= vtk.vtkPoints()
    scalars = vtk.vtkFloatArray()
    scalars.SetName("cut")

    i = 0
    for key in hits.keys()[:2] :
        fiber_hits = hits[key]['indices']
        for hit in fiber_hits :
            color = np.random.rand()
            for point in hit :
                z,y,x = point
                points.InsertPoint(i,(x,y,z))
                scalars.InsertTuple1(i,color)
                i += 1

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(scalars)

    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetInput(polydata)
    glyph.SetColorModeToColorByScalar()
    glyph.SetScaleModeToDataScalingOff() 
  
    glyph.Update()
    return glyph.GetOutput()



if __name__ == '__main__':
    import itk
    import vtk
    import medipy.diffusion
    import tensor
    from medipy.base import Image

    cutoff = 14
    pv_thr = 0.90
    spacing = (1,1,1)
    data = np.zeros((10,10,10,6),dtype=np.single)
    data[:,:,:] = (10,0,0,1,0,1)
    model = Image(data=data,spacing=spacing,data_type="vector")
    T,C = tensor.streamline_tractography(model,step=1.0,thr_fa=0.2,thr_angle=np.pi,thr_length=0,propagation_type="Euler",clustering_type="Fast",thr_clustering=1)
    tensor.skeleton(T,C)
    print "nb fibre =", len(T)
    print "nb cluster =", len(C)

    C_hits = bundle_cross_sections(T,C,nb_samples=None)
    N_hits = bundle_cross_sections_voxel(C_hits,spacing,cutoff,pv_thr)
    print "\nnb hits =", len(C_hits), len(N_hits)
    print "nb cuts =", len(C_hits[0]['indices']), len(N_hits[0]['indices'])

    dataset = tensor.fiber_polydata(T,C)
    cut_dataset = cut_polydata(C_hits)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName("/home/grigis/Bureau/fibers_test.vtk")
    writer.SetInput(cut_dataset)
    writer.Update()

