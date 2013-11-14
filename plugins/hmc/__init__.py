from HMC import segmentation, post_process_outliers

import os.path
import medipy.itk

medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "MediPyHMC")
