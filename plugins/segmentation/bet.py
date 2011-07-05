##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import itk

import medipy.itk

def bet(input, bt, output):
    """ Brain Extraction Tool.
        
        @TECHREPORT{SmithBET,
            author = Smith, Stevens M.,
            title = BET : Brain Extraction Tool,
            institution = FMRIB (Oxford Center for Functional Magnetic Resonance Imaging of the Brain),
            number = TR00SMS2b,
            url = http://www.fmrib.ox.ac.uk/analysis/research/bet/resolveuid/33c59904365e93773e533a5eb8e3b7d9
        
        <gui>
            <item name="input" type="Image" label="Input"/>
            <item name="bt" type="Float" initializer="0.5" label="b_t"
                tooltip="Main parameter (0 &lt;= b_t &lt;= 1)"/>
            <item name="output" type="Image" initializer="output=True" role="output"
                label="Output"/>
        </gui>
    """
    
    # We're in the same package as itkBetImageFilter, so it has already been
    # included in itk by __init__
    
    itk_input = medipy.itk.medipy_image_to_itk_image(input, False)
    bet_filter = itk.BETImageFilter[itk_input, itk_input].New(itk_input, BT=bt)
    itk_output = bet_filter()[0]
    medipy.itk.itk_image_to_medipy_image(itk_output, output, True)
