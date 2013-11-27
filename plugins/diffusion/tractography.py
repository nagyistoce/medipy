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

from fiber_clustering import local_skeleton_clustering, most_similar_track_mam

def streamline(model, step=0.5, minimum_fa=0.2, maximum_angle=numpy.pi/3, 
               minimum_length=50, propagation_type="Euler", seed_spacing=None, 
               mask=None) :
    """ Deterministic streamline propagation algorithm, return a list of tracks,
        where a track is a list of points in physical coordinates.
    
        * ``model`` : tensor field
        * ``step`` : propagation step
        * ``minimum_fa`` : minimum fractional anisotropy value allowed for
          propagation
        * ``maximum_angle`` : minimum angle value (in radians) allowed for
          propagation
        * ``minimum_length`` : minimum fiber length, in physical units
        * ``propagation_type`` : propagation criterion, may be either ``"Euler"``
          or ``"RungeKuttaOrder4"``
        * ``seed_spacing`` : spacing between seeds in physical units, defaults
          to ``2*model.spacing``
        * ``mask`` : optional mask to restrict seeding
    """

    if seed_spacing is None : 
        seed_spacing = model.spacing*(2.0,)
    seeds = _generate_image_sampling(model, seed_spacing/model.spacing, mask)

    itk_model = medipy.itk.medipy_image_to_itk_image(model, False)

    ScalarImage = itk.Image[itk.template(itk_model)[1]]
    tractography_filter = itk.StreamlineTractographyAlgorithm[
        itk_model, ScalarImage].New()
    
    tractography_filter.SetInputModel(itk_model)
    for seed in seeds:
        tractography_filter.AppendSeed(seed[::-1])
    tractography_filter.SetStepSize(step)
    tractography_filter.SetUseRungeKuttaOrder4(propagation_type=="RungeKuttaOrder4")
    tractography_filter.SetMaximumAngle(maximum_angle)
    tractography_filter.SetMinimumFA(minimum_fa)
    
    mask_itk = None
    if mask :
        mask_itk = medipy.itk.medipy_image_to_itk_image(mask, False)
        tractography_filter.SetMask(mask_itk)
    
    tractography_filter.Update()

    fibers = []
    for i in range(tractography_filter.GetNumberOfFibers()) :
        fiber = tractography_filter.GetOutputFiberAsPyArray(i)
        if _length(fiber,step)>=minimum_length :
            fibers.append(fiber)
    
    return fibers

def skeleton(fibers, clusters, metric="avg") :
    """ Update each cluster with its skeleton, i.e. the most similar fiber it
        contains. 
    """

    for cluster in clusters.values():
        bundle = [fibers[index] for index in cluster['indices']]
        skeleton_index = most_similar_track_mam(bundle, metric)[0]
        cluster['skel'] = cluster['indices'][skeleton_index]

def create_object_3d(fibers, clusters=None) :
    """ Build an :class:`~medipy.base.Object3D` containing the ``fibers``. If 
        ``clusters`` is not ``None``, it must be a dictionary mapping the index
        of a cluster to the cluster information (TODO : describe this structure).
        If the clusters information is available the point data of the Object3D
        will have scalar field called ``"cluster"`` ; this field will contain 
        the index of the cluster which the fiber belongs to. 
    """
    
    if clusters :
        fiber_to_cluster_id = {}
        skeleton_fibers = set()
        have_skeleton = False
        for cluster_id, cluster_info in clusters.items() :
            for fiber_index in cluster_info["indices"] :
                fiber_to_cluster_id[fiber_index] = cluster_id
            if "skel" in cluster_info :
                have_skeleton = True
                skeleton_fibers.add(cluster_info["skel"])

    fusion = vtk.vtkAppendPolyData()

    for fiber_index, fiber in enumerate(fibers) :
        nb_points = len(fiber)
        
        if clusters :
            cluster = vtk.vtkFloatArray()
            cluster.SetName("cluster")
            
            skeleton = vtk.vtkFloatArray()
            skeleton.SetName("skeleton")
        
        points= vtk.vtkPoints()
        line = vtk.vtkCellArray()
        line.InsertNextCell(nb_points)
        for i in range(nb_points):
            points.InsertPoint(i, fiber[i][::-1])
            line.InsertCellPoint(i)
            if clusters :
                cluster.InsertTuple1(i, fiber_to_cluster_id[fiber_index])
                skeleton.InsertTuple1(i, 0 if fiber_index in skeleton_fibers else 1)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        if clusters :
            if have_skeleton :
                polydata.GetPointData().AddArray(cluster)
                polydata.GetPointData().AddArray(skeleton)
            else :
                polydata.GetPointData().SetScalars(cluster)
        polydata.SetLines(line) 

        fusion.AddInput(polydata)
  
    fusion.Update()
    
    return Object3D(fusion.GetOutput(), "Fibers")

def streamline_gui(model, step, minimum_fa, maximum_angle, minimum_length, 
                   propagation_type, seed_spacing, mask, 
                   clustering, maximum_cluster_distance) :
    """ Deterministic streamline propagation algorithm.
    
        <gui>
            <item name="model" type="Image" label="Input"/>
            <item name="step" type="Float" initializer="0.5" 
                  label="Propagation step"/>
            <item name="minimum_fa" type="Float" initializer="0.2" 
                  label="Minimum FA"/>
            <item name="maximum_angle" type="Float"
                  initializer="numpy.pi/3" label="Maximum angle"/>
            <item name="minimum_length" type="Float" initializer="50" 
                  label="Minimum length"/>
            <item name="propagation_type" type="Enum" 
                  initializer="('Euler', 'Runge Kutta 4')" 
                  label="Propagation type"/>
            <item name="seed_spacing" type="Array" 
                  initializer="float, 3, 3, (2,2,2)" label="Seed spacing"/>
            <item name="mask" type="Image" 
                  initializer="value=None, may_be_empty=True, may_be_empty_checked=True"
                  label="Mask"/>
            <item name="clustering" type="Enum" initializer="(None, 'Fast')"
                  label="Clustering"/>
            <item name="maximum_cluster_distance" type="Float" initializer="20"
                  label="Maximum cluster distance"/>
            <item name="output" type="Object3D" role="return" label="Output"/>
        </gui>
    """
    
    fibers = streamline(model, step, minimum_fa, maximum_angle, minimum_length, 
                        propagation_type, seed_spacing, mask)
    
    if clustering is None :
        clusters = {}
    elif clustering == "Fast" :
        clusters = local_skeleton_clustering(fibers, maximum_cluster_distance)
    else :
        raise medipy.base.Exception("Unknown clustering mode: {0!r}".format(clustering))
    
    return create_object_3d(fibers, clusters)

def skeleton_gui(model, step, minimum_fa, maximum_angle, minimum_length, 
                 propagation_type, seed_spacing, mask, 
                 clustering, maximum_cluster_distance) :
    """ Skeleton of the white matter fibersStreamline 2nd order tensor tractography 

        <gui>
            <item name="model" type="Image" label="Input"/>
            <item name="step" type="Float" initializer="0.5" 
                  label="Propagation step"/>
            <item name="minimum_fa" type="Float" initializer="0.2" 
                  label="Minimum FA"/>
            <item name="maximum_angle" type="Float"
                  initializer="numpy.pi/3" label="Maximum angle"/>
            <item name="minimum_length" type="Float" initializer="50" 
                  label="Minimum length"/>
            <item name="propagation_type" type="Enum" 
                  initializer="('Euler', 'Runge Kutta 4')" 
                  label="Propagation type"/>
            <item name="seed_spacing" type="Array" 
                  initializer="float, 3, 3, (2,2,2)" label="Seed spacing"/>
            <item name="mask" type="Image" 
                  initializer="value=None, may_be_empty=True, may_be_empty_checked=True"
                  label="Mask"/>
            <item name="clustering" type="Enum" initializer="('Fast',)"
                  label="Clustering"/>
            <item name="maximum_cluster_distance" type="Float" initializer="20"
                  label="Maximum cluster distance"/>
            <item name="output" type="Object3D" role="return" label="Output"/>
        </gui>
    """

    fibers = streamline(model, step, minimum_fa, maximum_angle, minimum_length, 
                        propagation_type, seed_spacing, mask)
    
    if clustering is None :
        clusters = {}
    elif clustering == "Fast" :
        clusters = local_skeleton_clustering(fibers, maximum_cluster_distance)
    else :
        raise medipy.base.Exception("Unknown clustering mode: {0!r}".format(clustering))
    
    skeleton(fibers, clusters)
    
    return create_object_3d(fibers, clusters)

def _generate_image_sampling(image, step, mask=None) :
    """ Generate seeds in physical space to initialize tractography. ``step`` is
        expressed in continuous index coordinates and in numpy order
    """

    begin = image.ndim*(0,)
    end = numpy.subtract(image.shape, 1)
    
    grid = numpy.mgrid[[slice(b,e,s) for b,e,s in zip(begin, end, step)]]
    
    if mask is None :
        seeds = [image.index_to_physical(i) for i in grid.reshape(3, -1).T]
    else :
        seeds = []
        for index_image in grid.reshape(3, -1).T :
            point_image = image.index_to_physical(index_image)
            index_mask = tuple(mask.physical_to_index(point_image, True))
            if mask.is_inside(index_mask) and mask[index_mask] != 0 :
                seeds.append(point_image)
    
    return seeds

def _length(xyz, constant_step=None):
    """ Euclidean length of track line in mm 
    """
    if constant_step==None :
        if xyz.shape[0] < 2 :
            return 0
        else :
            dists = numpy.sqrt((numpy.diff(xyz, axis=0)**2).sum(axis=1))
            return numpy.sum(dists)
    else :
        return (xyz.shape[0]-1)*constant_step
