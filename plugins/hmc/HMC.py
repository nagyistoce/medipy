##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

import numpy

import itk

import medipy.itk
import medipy.hmc

def hmc(input_images, mask, atlas_bool, atlas_images, number_iter, number_classes, criterion_outliers, criterion_outliers_value, flair_bool, flair_position):
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input_images[0], False)
    nb_images = len(input_images)
    taille_cube = max_puissance(max(input_images[0].shape))[1]
    criterion_dico = {"percentage":0, "threshold":1}
    hmc_filter = itk.HiddenMarkovChainFilter[itk_input, itk_input].New(
        Number_images=nb_images, Atlas_bool=atlas_bool, Number_iter=number_iter, Number_classes=number_classes, Criterion_outliers=criterion_dico[criterion_outliers], Criterion_outliers_value=criterion_outliers_value, Taille_cube=taille_cube, Flair_bool=flair_bool, Position_Flair=flair_position)
 
    input_images_resized = []
    itk_inputs = []

    for idx, input_image in enumerate(input_images):
        input_images_resized.append(resize_input(input_image))
        itk_inputs.append(medipy.itk.medipy_image_to_itk_image(input_images_resized[len(input_images_resized)-1], False))
        hmc_filter.SetInputImage(idx, itk_inputs[len(itk_inputs)-1])

    if atlas_bool:
        for idx, atlas_image in enumerate(atlas_images):

            input_images_resized.append(resize_input(atlas_image))            
            itk_inputs.append(medipy.itk.medipy_image_to_itk_image(input_images_resized[len(input_images_resized)-1], False))
            hmc_filter.SetInputImage(idx+nb_images, itk_inputs[len(itk_inputs)-1])

        mask = resize_input(mask)
        itk_inputs.append( medipy.itk.medipy_image_to_itk_image(mask, False))
        hmc_filter.SetInputImage(idx+nb_images+1, itk_inputs[len(itk_inputs)-1])
    else:
        mask = resize_input(mask)
        itk_inputs.append(medipy.itk.medipy_image_to_itk_image(mask, False))
        hmc_filter.SetInputImage(nb_images, itk_inputs[len(itk_inputs)-1])
       

    hmc_filter.Update()

    itk_output_seg_image = hmc_filter.GetOuputSegImage()
    output_seg_image = medipy.itk.itk_image_to_medipy_image(itk_output_seg_image, None, True)
    itk_output_outliers_image = hmc_filter.GetOuputOutliersImage()
    output_outliers_image = medipy.itk.itk_image_to_medipy_image(itk_output_outliers_image, None, True)

    return [output_seg_image, output_outliers_image]


def resize_input(input_image):

    dim = input_image.shape
    dim_out = max_puissance(max(dim))[0]
    shape_out = [dim_out, dim_out, dim_out]
    data_out = numpy.zeros(shape_out, dtype=input_image.dtype)
    data_out[0:dim[0],0:dim[1],0:dim[2]] = input_image.data
    output = medipy.base.Image(data=data_out, shape=shape_out, dtype=input_image.dtype, data_type=input_image.data_type, image_type=input_image.image_type, annotations=input_image.annotations, metadata=input_image.metadata)
    output.copy_information(input_image)

    return output


def max_puissance(nb):

    res=1
    power=0
    while res<nb:
        res<<=1
        power+=1

    return [res,power]
 
