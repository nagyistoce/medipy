import os.path
import medipy.itk

medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "BETImageFilter")
medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "StatisticalChangeDetectionImageFilters")

from api import *
