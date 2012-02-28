import medipy.medimax.traitement
import numpy

#-------------------------------------------------------------
#  Correct Dark Bright Alternation of slices
#-------------------------------------------------------------

def CorrectDarkBrightCoronalAlternation_gui(im, imres, axis) :
	"""
	Correct Dark Bright Alternation of slices
	
			:gui:
				im : Image
					Image to correct
 	
				imres : Image : output=True
					Corrected image
		
				axis : Enum : ('x','y','z')
					Axis along which slices are corrected


	""" 
	#medipy.components.medimax3.traitement.imx_corr_sagittal_3d_p(im,imres)
	if str(axis)=='x':
		imtmp = numpy.swapaxes(im,2,1)
		medipy.medimax.traitement.imx_corr_sagittal_3d_p(imtmp,imtmp)
 		imres = numpy.swapaxes(imtmp,1,2)	

	elif str(axis)=='y':
		medipy.medimax.traitement.imx_corr_sagittal_3d_p(im,imres)
	else :
		imtmp = numpy.swapaxes(im,0,1)
		medipy.medimax.traitement.imx_corr_sagittal_3d_p(imtmp,imtmp)
 		imres = numpy.swapaxes(imtmp,1,0)			
                                              
