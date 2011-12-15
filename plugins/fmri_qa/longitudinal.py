##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

from medipy.fmri_qa import io

def longitudinal(summary_files, output_directory):
    
    snr = []
    sfnr = []
    fluctuation = []
    drift = []
    
    for summary_file in summary_files :
        summary = io.load_summary(summary_file)
        
        snr.append((summary["date"], summary["snr"]))
        sfnr.append((summary["date"], summary["sfnr"]))
        fluctuation.append((summary["date"], summary["fluctuation"]))
        drift.append((summary["date"], summary["drift"]))
    
    io.save_longitudinal(snr, sfnr, fluctuation, drift, output_directory)
    io.save_longitudinal_figures(snr, sfnr, fluctuation, drift, output_directory)