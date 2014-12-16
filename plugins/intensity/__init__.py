import os.path
import medipy.itk

medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "MediPyIntensity")

import normalization

from api import *
from replace import replace_labels
