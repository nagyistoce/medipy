##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

from medipy.fmri_qa import io

def measure(summary_files, output_directory, baseline=None):
    
    snr = []
    sfnr = []
    fluctuation = []
    drift = []
    rdc = []
    
    for summary_file in summary_files :
        summary = io.load_summary(summary_file)
        
        snr.append((summary["date"], summary["snr"]))
        sfnr.append((summary["date"], summary["sfnr"]))
        fluctuation.append((summary["date"], summary["fluctuation"]))
        drift.append((summary["date"], summary["drift"]))
        rdc.append((summary["date"], summary["rdc"]))
    
    snr.sort(key=lambda x:x[0])
    sfnr.sort(key=lambda x:x[0])
    fluctuation.sort(key=lambda x:x[0])
    drift.sort(key=lambda x:x[0])
    
    io.save_longitudinal(snr, sfnr, fluctuation, drift, rdc, output_directory)
    io.save_longitudinal_figures(snr, sfnr, fluctuation, drift, output_directory, baseline)