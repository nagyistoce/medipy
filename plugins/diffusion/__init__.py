import os.path
import medipy.itk

from dtiEigItk import dtiEigItk

#medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "SymmetricSpectralAnalysisImageFilter")
medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "SecondOrderSymmetricTensorReconstructionFilter")
medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "MeanDiffusivityImageFilter")
medipy.itk.load_wrapitk_module(os.path.dirname(__file__), "FractionalAnisotropyImageFilter")
