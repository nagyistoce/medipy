
##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################


import numpy as np
import medipy.diffusion.fiber_clustering
import sys


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

