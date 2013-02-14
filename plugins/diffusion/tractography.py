##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

def streamline_tractography(model,seeds,step=1.0,thr_fa=0.2,thr_angle=np.pi,rk4=False):
    """ Streamline 2nd order tensor tractography """
 
    tractography_filter = itk.StreamlineTractographyAlgorithm[itk.VectorImage[itk.F,3], itk.Image[itk.F,3]].New()

    itk_model = medipy.itk.medipy_image_to_itk_image(model, False)
    print itk_model
    tractography_filter.SetInputModel(itk_model)
    for seed in seeds:
        tractography_filter.AppendSeed(seed)
    tractography_filter.SetStepSize(step)
    tractography_filter.SetUseRungeKuttaOrder4(rk4)
    tractography_filter.SetThresholdAngle(thr_angle)
    tractography_filter.SetThresholdFA(thr_fa)
    tractography_filter.Update()

    print tractography_filter.GetNumberOfFibers()