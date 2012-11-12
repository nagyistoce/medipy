##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import itk
import numpy
import vtk

from medipy.base import Object3D 
import medipy.itk

import fiber_clustering
from utils import length, generate_image_sampling

def streamline_tractography_gui(model,*args,**kwargs) :
    """ Streamline 2nd order tensor tractography 

    <gui>
        <item name="model" type="Image" label="Input"/>
        <item name="step" type="Float" initializer="1.0" label="Propagation step"/>
        <item name="thr_fa" type="Float" initializer="0.2" label="FA threshold"/>
        <item name="thr_angle" type="Float" initializer="1.04" label="Angle threshold"/>
        <item name="thr_length" type="Float" initializer="50.0" label="Length threshold"/>
        <item name="thr_clustering" type="Float" initializer="20.0" label="Threshold in clustering"/>
        <item name="propagation_type" type="Enum" initializer="('Euler', 'Runge Kutta 4')" label="Choose propagation order"/>
        <item name="clustering_type" type="Enum" initializer="('None', 'Fast')" label="Choose the clustering strategy"/>
        <item name="output" type="Object3D" role="return" label="Output"/>
    </gui>
    """

    step = kwargs['step']
    thr_fa = kwargs['thr_fa']
    thr_angle = kwargs['thr_angle']
    thr_length = kwargs['thr_length']
    propagation_type = kwargs['propagation_type']
    clustering_type = kwargs['clustering_type']
    thr_clustering = kwargs['thr_clustering']

    T,C = streamline_tractography(model,step,thr_fa,thr_angle,thr_length,propagation_type,clustering_type,thr_clustering)
    vtk_polydata = fiber_polydata(T,C)
    output = Object3D(vtk_polydata,"Streamline Tractography")

    return output

def skeleton_gui(model,*args,**kwargs) :
    """ Skeleton of the white matter fibersStreamline 2nd order tensor tractography 

    <gui>
        <item name="model" type="Image" label="Input"/>
        <item name="step" type="Float" initializer="1.0" label="Propagation step"/>
        <item name="thr_fa" type="Float" initializer="0.2" label="FA threshold"/>
        <item name="thr_angle" type="Float" initializer="1.04" label="Angle threshold"/>
        <item name="thr_length" type="Float" initializer="50.0" label="Length threshold"/>
        <item name="thr_clustering" type="Float" initializer="20.0" label="Threshold in clustering"/>
        <item name="propagation_type" type="Enum" initializer="('Euler', 'Runge Kutta 4')" label="Choose propagation order"/>
        <item name="clustering_type" type="Enum" initializer="('Fast',)" label="Choose the clustering strategy"/>
        <item name="output" type="Object3D" role="return" label="Output"/>
    </gui>
    """

    step = kwargs['step']
    thr_fa = kwargs['thr_fa']
    thr_angle = kwargs['thr_angle']
    thr_length = kwargs['thr_length']
    propagation_type = kwargs['propagation_type']
    clustering_type = kwargs['clustering_type']
    thr_clustering = kwargs['thr_clustering']

    T,C = streamline_tractography(model,step,thr_fa,thr_angle,thr_length,propagation_type,clustering_type,thr_clustering)
    skeleton(T,C)
    vtk_polydata = fiber_polydata(T,C)
    output = Object3D(vtk_polydata,"Fiber Skeleton")

    return output

def skeleton(T,C) :
    """ Skeleton of the white matter fibersStreamline 2nd order tensor tractography 
    """

    for key in C.keys():
        N = C[key]['N']
        bundle_indices = C[key]['indices']
        B = []
        for indice in bundle_indices:
            B.append(T[indice])
        s,si = medipy.diffusion.fiber_clustering.most_similar_track_mam(B,metric='avg') 
        C[key]['skel'] = bundle_indices[s]  

def fiber_polydata(T,C) : 
    """ Construct the vtk polydata 
    """

    clusters_skel = np.ones((len(T),))
    clusters = np.ones((len(T),))
    nb_clusters = len(C.keys())
    flag = 'skel' in C[C.keys()[0]].keys()
    for cnt,c in enumerate(C.keys()) :
        indices = C[c]['indices']
        for i in indices :
            clusters[i] = float(cnt)/ float(nb_clusters)  
    if flag :
        for cnt,c in enumerate(C.keys()) :
            indice = C[c]['skel']
            clusters_skel[indice] = 0.0

    fusion = vtk.vtkAppendPolyData()

    for cluster,cluster_skel,fiber in zip(clusters,clusters_skel,T) :
        nb_points = fiber.shape[0]
        scalars = vtk.vtkFloatArray()
        scalars.SetName("cluster")
        scalars_skel = vtk.vtkFloatArray()
        scalars_skel.SetName("skeleton")
        points= vtk.vtkPoints()
        line = vtk.vtkCellArray()
        line.InsertNextCell(nb_points)

        for i in range(nb_points):
            z,y,x = fiber[i]
            points.InsertPoint(i,(x,y,z))
            scalars.InsertTuple1(i,cluster)
            scalars_skel.InsertTuple1(i,cluster_skel)
            line.InsertCellPoint(i)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        if flag :
            polydata.GetPointData().AddArray(scalars)
            polydata.GetPointData().AddArray(scalars_skel)
        else :
            polydata.GetPointData().SetScalars(scalars)  
        polydata.SetLines(line) 

        fusion.AddInput(polydata)
  
    fusion.Update()
    return fusion.GetOutput()

def streamline_tractography(model,step,thr_fa,thr_angle,thr_length,propagation_type,clustering_type,thr_clustering) :
    """ Streamline 2nd order tensor tractography 
    """

    seeds = generate_image_sampling(model,step=2.0*model.spacing)
    ndim = len(model.shape)

    tractography_filter = itk.StreamlineTractographyAlgorithm[itk.VectorImage[itk.F,ndim], itk.Image[itk.F,ndim]].New()

    itk_model = medipy.itk.medipy_image_to_itk_image(model, False)
    tractography_filter.SetInputModel(itk_model)
    for seed in seeds:
        tractography_filter.AppendSeed(seed)
    tractography_filter.SetStepSize(step)
    if propagation_type=="Euler" :
        tractography_filter.SetUseRungeKuttaOrder4(False)
    else :
        tractography_filter.SetUseRungeKuttaOrder4(True)
    tractography_filter.SetThresholdAngle(thr_angle)
    tractography_filter.SetThresholdFA(thr_fa)
    tractography_filter.Update()

    nb_fiber = tractography_filter.GetNumberOfFibers()
    fibers = []
    for i in range(nb_fiber) :
        fibers.append(tractography_filter.GetOutputFiberAsPyArray(i))

    final_fibers = []
    for fiber in fibers :
        if length(fiber,step)>=thr_length :
            final_fibers.append(fiber)

    C = {}
    if clustering_type=='Fast' :
        C = medipy.diffusion.fiber_clustering.local_skeleton_clustering(final_fibers, d_thr=thr_clustering)

    return final_fibers,C


def streamline_tractography_(model,seeds=None,step=1.0,thr_fa=0.2,thr_angle=np.pi,rk4=False):
    """ Streamline 2nd order tensor tractography """

    if seeds==None :
            seeds = generate_image_sampling(model,step=(1,1,1))
    print len(seeds)

    tractography_filter = itk.StreamlineTractographyAlgorithm[itk.VectorImage[itk.F,3], itk.Image[itk.F,3]].New()

    itk_model = medipy.itk.medipy_image_to_itk_image(model, False)
    tractography_filter.SetInputModel(itk_model)
    for seed in seeds:
        tractography_filter.AppendSeed(seed)
    tractography_filter.SetStepSize(step)
    tractography_filter.SetUseRungeKuttaOrder4(rk4)
    tractography_filter.SetThresholdAngle(thr_angle)
    tractography_filter.SetThresholdFA(thr_fa)
    tractography_filter.Update()

    nb_fiber = tractography_filter.GetNumberOfFibers()
    fibers = []
    for i in range(nb_fiber) :
        fibers.append(tractography_filter.GetOutputFiberAsPyArray(i))
    print [ fiber.shape for fiber in fibers]

    save_fibers(fibers,fout="/home/grigis/Bureau/fibers_test.vtk",step=step,thr_length=50.0)


def save_fibers(fibers,fout,step=None,thr_length=0.0,ensemble=1) :

    fusion = vtk.vtkAppendPolyData()

    for fiber in fibers :

        if length(fiber,step)>=thr_length :
            nb_points = fiber.shape[0]
            scalars = vtk.vtkFloatArray()
            points= vtk.vtkPoints()
            line = vtk.vtkCellArray()
            line.InsertNextCell(nb_points)

            for i in range(nb_points):
                z,y,x = fiber[i]
                points.InsertPoint(i,(x,y,z))
                scalars.InsertTuple1(i,ensemble)
                line.InsertCellPoint(i)

            polydata = vtk.vtkPolyData()
            polydata.SetPoints(points)
            polydata.GetPointData().SetScalars(scalars)
            polydata.SetLines(line) 

            fusion.AddInput(polydata)

    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(fout)
    writer.SetInput(fusion.GetOutput())
    writer.Update()
