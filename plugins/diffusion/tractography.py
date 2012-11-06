##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

def streamline_tractography(model,*args,**kwargs) :
    """ Streamline 2nd order tensor tractography 
 
    <gui>
        <item name="model" type="Image" label="Input"/>
        <item name="step" type="Float" initializer="1.0" label="Propagation step"/>
        <item name="thr_fa" type="Float" initializer="0.2" label="FA threshold"/>
        <item name="thr_angle" type="Float" initializer="1.04" label="Angle threshold"/>
        <item name="thr_length" type="Float" initializer="25.0" label="Length threshold"/>
        <item name="propagation_type" type="Enum" initializer="('Euler', 'Runge Kutta 4')" label="Choose propagation order"/>
        <item name="output" type="Object3D" role="return" label="Output"/>
    </gui>
    """
    step = kwargs['step']
    thr_fa = kwargs['thr_fa']
    thr_angle = kwargs['thr_angle']
    thr_length = kwargs['thr_length']
    propagation_type = kwargs['propagation_type']

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

    output = Object3D()
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
                scalars.InsertTuple1(i,1.0)
                line.InsertCellPoint(i)

            polydata = vtk.vtkPolyData()
            polydata.SetPoints(points)
            polydata.GetPointData().SetScalars(scalars)
            polydata.SetLines(line) 

            fusion.AddInput(polydata)
   
    fusion.Update()
    output.dataset = fusion.GetOutput()
    return output





def streamline_tractography_(model,seeds=None,step=1.0,thr_fa=0.2,thr_angle=np.pi,rk4=False):
    """ Streamline 2nd order tensor tractography """
 
    if seeds==None :
            seeds = generate_image_sampling(model,step=(1,1,1))
    print len(seeds)

    tractography_filter = itk.StreamlineTractographyAlgorithm[itk.VectorImage[itk.F,ndim], itk.Image[itk.F,ndim]].New()

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
    print "nb fibers=", nb_fiber
    fibers = []
    for i in range(nb_fiber) :
        fibers.append(tractography_filter.GetOutputFiberAsPyArray(i))
    print [ fiber.shape for fiber in fibers]
    print fibers[0]
    
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
