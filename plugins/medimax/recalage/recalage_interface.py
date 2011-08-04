##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011             
# Distributed under the terms of the CeCILL-B license, as published by 
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to            
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html       
# for details.                                                      
##########################################################################

import medipy.medimax.recalage
import numpy

#-------------------------------------------------------------
#  Interface entre type d'Interpolation (enum) et le numero dans medimax
#-------------------------------------------------------------

def InterpolationNumberInMedimax(inter_type) :
	""" 
	Convert an interpolation type (from ENUM list) to the corresponding number in medimax
	"""
	
	InterpolationConverter = {"Nearest" : 0,"Linear" : 1,"SinCard" : 2,"QuickSinCard2" : 3,"QuickSinCard3" : 4,"Bspline2" : 5,"Bspline3" : 6,"Bspline4" : 7,"Bspline5" : 8,	"Label" : 9}

	if not InterpolationConverter.has_key(str(inter_type)) :
		print("Interpolation type [ %s ] does not exist in Medimax ! Linear Interpolation is used by default."%str(inter_type))
		return 1
		
	return InterpolationConverter[str(inter_type)]



#-------------------------------------------------------------
#  Interface entre type de visu trf (enum) et le numero dans medimax
#-------------------------------------------------------------

def VisuTrfNumberInMedimax(visu_type) :
	""" 
	Convert a visu type (from ENUM list) to the corresponding number in medimax
	"""
	
	VisuConverter = {'Module' : 0,'Jacobian' : 1,'X-displacement' : 2,'Y-displacement' : 3,'Z-displacement' : 4}
	if not VisuConverter.has_key(str(visu_type)) :
		print("Visu type [ %s ] does not exist in Medimax ! Module is considered by default."%str(visu_type))
		return 0
		
	return VisuConverter[str(visu_type)]


#-------------------------------------------------------------
#  Multiresolution Linear registration based on MI
#-------------------------------------------------------------

def recalage_lineaire_multi_IM(registration_type, imreca, imref, imres, filename) :
	"""
	Multiresolution Linear registration based on Mutual Information, Simplex optimization and linear interpolation

		:gui:
			registration_type : Enum : ('Rigid','Rigid+Zoom','Affine')
				Registration type 

			imreca : Image
				Image to register
 	
			imref : Image
				Reference image
			
			imres : Image : output=True
				Registered image

			filename : String
				Name of saved .trf file
	""" 
	RegistrationTypeConverter = {"Rigid" : 0,"Rigid+Zoom" : 1,"Affine" : 2}
	medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,RegistrationTypeConverter[registration_type],5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Multiresolution Rigid registration based on MI
#-------------------------------------------------------------

def recalage_rigide_multi_IM(imreca, imref, imres, filename) :
	"""
	Multiresolution rigid registration based on Mutual Information, Simplex optimization and linear interpolation

		:gui:
			imreca : Image
				Image to register
 	
			imref : Image
				Reference image
			
			imres : Image : output=True
				Registered image

			filename : String
				Name of saved .trf file
	""" 
	medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,0,5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Multiresolution Rigid + Zoom registration based on MI
#-------------------------------------------------------------

def recalage_rigidezoom_multi_IM(imreca, imref, imres, filename) :
	"""
	Multiresolution rigid + zoom registration based on Mutual Information, Simplex optimization and linear interpolation

		:gui:
			imreca : Image
				Image to register
 	
			imref : Image
				Reference image
			
			imres : Image : output=True
				Registered image
		
			filename : String
				Name of saved .trf file

	""" 
	medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,1,5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Multiresolution Affine registration based on MI
#-------------------------------------------------------------

def recalage_affine_multi_IM(imreca, imref, imres, filename) :
	"""
	Multiresolution affine registration based on Mutual Information, Simplex optimization and linear interpolation

		:gui:
			imreca : Image
				Image to register
 	
			imref : Image
				Reference image
			
			imres : Image : output=True
				Registered image
		
			filename : String
				Name of saved .trf file

	""" 
	medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,2,5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Topology preserving Bspline-based registration 
#-------------------------------------------------------------

def recalage_Bspline_topo(imreca, imref, imres, resolf, symetrique, filename) :
	"""
	Topology preserving Bspline-based registration (default settings)

		:gui:
			imreca : Image
				Image to register
 	
			imref : Image
				Reference image
			
			imres : Image : output=True
				Registered image

			resolf : Int : value=1, range=(1,7)
				final resolution
			
			symetrique : Enum : ('Non Symmetric','Symmetric')
				use symmetric pairwise registration 
								
			filename : String
				Name of saved .trf file
	""" 
	dicoConverter = {"Non Symmetric" : 0,"Symmetric" : 1}
	medipy.medimax.recalage.BsplineRegistration3d(imref, imreca, imres,0,2,1.0,2,1,str(filename),resolf, 0.0,100000.0,12, 10,0, dicoConverter[str(symetrique)])


#-------------------------------------------------------------
#  Apply a transf_3d 
#-------------------------------------------------------------

def ApplyTransfo3d_GUI(imdeb, nomfichier, imres, inter_type) :
	"""
	Warp an image according to a .trf file
		
		:gui:
			imdeb : Image
				Image to warp 
			
			imres : Image : output=True
				Warped image
			
			inter_type : Enum : ('Nearest','Linear','SinCard','QuickSinCard2','QuickSinCard3','Bspline2','Bspline3','Bspline4','Bspline5','Label')
				Interpolation method 
				
			nomfichier : File
				trf file
	""" 
	medipy.medimax.recalage.ApplyTransfo3d(imdeb,str(nomfichier),imres,InterpolationNumberInMedimax(inter_type))

#-------------------------------------------------------------
#  Combine two transf_3d 
#-------------------------------------------------------------

def CombineTransfo3d_GUI(nomfichier1, nomfichier2, nomfichierres) :
	"""
	Combine two .trf files
		
		:gui:
			nomfichier1 : File
				first trf file
			
			nomfichier2 : File
				second trf file

			nomfichierres : String
				resulting trf file
	""" 
	medipy.medimax.recalage.CombineTransfo3d(str(nomfichier1), str(nomfichier2), str(nomfichierres), 5)


#-------------------------------------------------------------
#  Invert a transf_3d 
#-------------------------------------------------------------

def InvertTransfo3d_GUI(nomfichier, nomfichres,wdthref, hghtref, dpthref, dxref, dyref, dzref) :
	"""
	Invert a .trf file
		
		:gui:
			nomfichier : File
				.trf file to invert
			
			nomfichres : String
				resulting trf file

			wdthref : Int : value=-1
				width of the transformation  (optional)

			hghtref : Int : value=-1
				height of the transformation (optional)

			dpthref : Int : value=-1
				depth of the transformation (optional)

			dxref : Float : value=-1
				dx of the transformation (optional)

			dyref : Float : value=-1
				dy of the transformation (optional)

			dzref : Float : value=-1
				dz of the transformation (optional)
		""" 
	medipy.medimax.recalage.InvertTransfo3d(str(nomfichier), str(nomfichres), wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)


#-------------------------------------------------------------
#  MRI Info 3D
#-------------------------------------------------------------

def MriInfo3D_GUI(im) :
	"""
	Info MRI 3D
		
		:gui:
			im : Image
				Image to get info 
	""" 
	medipy.medimax.recalage.MriInfo3D(im)


#-------------------------------------------------------------
#  Visualisation of a transf_3d 
#-------------------------------------------------------------

def VisuTrfFile_GUI( nomfichier, output, visu_type) :
	"""
	Visualisation of a .trf file
		
		:gui:
			nomfichier : File
				trf file

			output : Image : output=True
				resulting image
			
			visu_type : Enum : ('Module','Jacobian','X-displacement','Y-displacement','Z-displacement')
				type of visualisation 
				
	""" 
	medipy.medimax.recalage.VisuTrfFile(str(nomfichier),output,VisuTrfNumberInMedimax(visu_type))

