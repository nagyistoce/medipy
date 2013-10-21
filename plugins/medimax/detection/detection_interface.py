##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import medipy.medimax.detection
import numpy

#-------------------------------------------------------------
#  Compute a ROC Curve
#-------------------------------------------------------------

def ComputeRocCurve_GUI(imDetection, maskimDetection, maskBonnesDetec, nbPointsH, fichier) :
	"""
	Compute the ROC Curve to assess a detection result
		
		:gui:
			imDetection : Image
				Image of detection 
			
			maskimDetection : Image 
				ROI where to compute the ROC curve 
			
			maskBonnesDetec : Image 
				Ground truth Image 
			
			nbPointsH : Int : value=200
				Number of points
			
			fichier : File
				.dat file
	""" 

	print imDetection.flags, imDetection.shape
	print maskimDetection.flags, maskimDetection.shape
	print maskBonnesDetec.flags, maskBonnesDetec.shape
	print fichier
	print nbPointsH

	print 't'
	medipy.medimax.detection.ComputeRocCurve(imDetection, maskimDetection, maskBonnesDetec, nbPointsH, str(fichier))
	print 'tt'
