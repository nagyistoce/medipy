import os.path
import medipy.itk

medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "MediPySegmentation")

from api import *
