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
    
    InterpolationConverter = {"Nearest" : 0,"Linear" : 1,"SinCard" : 2,"QuickSinCard2" : 3,"QuickSinCard3" : 4,"Bspline2" : 5,"Bspline3" : 6,"Bspline4" : 7,"Bspline5" : 8, "Label" : 9}

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

        <gui>
            <item name="registration_type" type="Enum" 
                  initializer="('Rigid','Rigid+Zoom','Affine')"
                  label="Registration type"/> 
            <item name="imreca" type="Image" label="Image to register" />
            <item name="imref" type="Image" label="Reference image" />
            <item name="imres" type="Image" initializer="output=True"
                role="output" label="Registered image" />
            <item name="filename" type="String" label="Name of saved .trf file"/>
        </gui>
    """ 
    RegistrationTypeConverter = {"Rigid" : 0,"Rigid+Zoom" : 1,"Affine" : 2}
    medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,RegistrationTypeConverter[registration_type],5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Multiresolution Rigid registration based on MI
#-------------------------------------------------------------

def recalage_rigide_multi_IM(imreca, imref, imres, filename) :
    """
    Multiresolution rigid registration based on Mutual Information, Simplex optimization and linear interpolation

        <gui>
            <item name="imreca" type="Image" label="Image to register"/>
            <item name="imref" type="Image" label="Reference image"/>
            <item name="imres" type="Image" role="output" 
                  initializer="output=True" label="Registered image"/>
            <item name="filename" type="String" label="Name of saved .trf file"/>
        </gui>
    """ 
    medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,0,5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Multiresolution Rigid + Zoom registration based on MI
#-------------------------------------------------------------

def recalage_rigidezoom_multi_IM(imreca, imref, imres, filename) :
    """
    Multiresolution rigid + zoom registration based on Mutual Information, Simplex optimization and linear interpolation

        <gui>
            <item name="imreca" type="Image" label="Image to register"/>
            <item name="imref" type="Image" label="Reference image"/>
            <item name="imres" type="Image" role="output" 
                  initializer="output=True" label="Registered image"/>
            <item name="filename" type="String" label="Name of saved .trf file"/>
        </gui>

    """ 
    medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,1,5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Multiresolution Affine registration based on MI
#-------------------------------------------------------------

def recalage_affine_multi_IM(imreca, imref, imres, filename) :
    """
    Multiresolution affine registration based on Mutual Information, Simplex optimization and linear interpolation

        <gui>
            <item name="imreca" type="Image" label="Image to register"/>
            <item name="imref" type="Image" label="Reference image"/>
            <item name="imres" type="Image" role="output" 
                  initializer="output=True" label="Registered image"/>
            <item name="filename" type="String" label="Name of saved .trf file"/>
        </gui>

    """ 
    medipy.medimax.recalage.LinearRegistration(imref,imreca, imres,2,5,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Topology preserving Bspline-based registration 
#-------------------------------------------------------------

def recalage_Bspline_topo(imreca, imref, imres, resolf, symetrique, filename) :
    """
    Topology preserving Bspline-based registration (default settings)

        <gui>
            <item name="imreca" type="Image" label="Image to register"/>
            <item name="imref" type="Image" label="Reference image"/>
            <item name="imres" type="Image" initializer="output=True" role="output" 
                label="Registered image"/>
            <item name="resolf" type="Int" initializer="value=1, range=(1,7)"
                label="final resolution"/>
            <item name="symetrique" type="Enum" initializer="('Non Symmetric','Symmetric')"
                label="use symmetric pairwise registration"/> 
            <item name="filename" type="String" label="Name of saved .trf file"/>
        </gui>
    """ 
    dicoConverter = {"Non Symmetric" : 0,"Symmetric" : 1}
    medipy.medimax.recalage.BsplineRegistration3d(imref, imreca, imres,0,2,1.0,4,1,str(filename),resolf, 0.0,100000.0,12, 10,0, dicoConverter[str(symetrique)], 0.5)


#-------------------------------------------------------------
#  Apply a transf_3d 
#-------------------------------------------------------------

def ApplyTransfo3d_GUI(imdeb, nomfichier, imres, inter_type) :
    """
    Warp an image according to a .trf file
        
        <gui>
            <item name="imdeb" type="Image" label="Image to warp"/> 
            <item name="imres" type="Image" role="output" initializer="output=True"
                label="Warped image"/>
            <item name="inter_type" type="Enum" initializer="('Nearest','Linear',
                'SinCard','QuickSinCard2','QuickSinCard3','Bspline2','Bspline3',
                'Bspline4','Bspline5','Label')" label="Interpolation method"/> 
            <item name="nomfichier" type="File" label="trf file"/>
        </gui>
    """ 
    medipy.medimax.recalage.ApplyTransfo3d(imdeb,str(nomfichier),imres,InterpolationNumberInMedimax(inter_type))

#-------------------------------------------------------------
#  Combine two transf_3d 
#-------------------------------------------------------------

def CombineTransfo3d_GUI(nomfichier1, nomfichier2, nomfichierres) :
    """
    Combine two .trf files
        
        <gui>
            <item name="nomfichier1" type="File" label="first trf file"/>
            <item name="nomfichier2" type="File" label="second trf file"/>
            <item name="nomfichierres" type="String" label="resulting trf file"/>
        </gui>
    """ 
    medipy.medimax.recalage.CombineTransfo3d(str(nomfichier1), str(nomfichier2), str(nomfichierres), 5)


#-------------------------------------------------------------
#  Invert a transf_3d 
#-------------------------------------------------------------

def InvertTransfo3d_GUI(nomfichier, nomfichres,wdthref, hghtref, dpthref, dxref, dyref, dzref) :
    """
    Invert a .trf file
        
        <gui>
            <item name="nomfichier" type="File" label=".trf file to invert"/>
            <item name="nomfichres" type="String" label="resulting trf file"/>
            <item name="wdthref" type="Int" initializer="value=-1" 
                  label="width of the transformation  (optional)"/>
            <item name="hghtref" type="Int" initializer="value=-1" 
                  label="height of the transformation (optional)"/>
            <item name="dpthref" type="Int" initializer="value=-1" 
                  label="depth of the transformation (optional)"/>
            <item name="dxref" type="Float" initializer="value=-1" 
                  label="dx of the transformation (optional)"/>
            <item name="dyref" type="Float" initializer="value=-1" 
                  label="dy of the transformation (optional)"/>
            <item name="dzref" type="Float" initializer="value=-1" 
                  label="dz of the transformation (optional)"/>
        </gui>
        """ 
    medipy.medimax.recalage.InvertTransfo3d(str(nomfichier), str(nomfichres), wdthref, hghtref, dpthref, dxref, dyref, dzref, 0.01)


#-------------------------------------------------------------
#  MRI Info 3D
#-------------------------------------------------------------

def MriInfo3D_GUI(im) :
    """
    Info MRI 3D
        
        <gui>
            <item name="im" type="Image" label="Image to get info"/>
        </gui> 
    """ 
    medipy.medimax.recalage.MriInfo3D(im)


#-------------------------------------------------------------
#  Visualisation of a transf_3d 
#-------------------------------------------------------------

def VisuTrfFile_GUI( nomfichier, output, visu_type) :
    """
    Visualisation of a .trf file
        
        <gui>
            <item name="nomfichier" type="File" label="trf file"/>
            <item name="output" type="Image" initializer="output=True" role="output"
                  label="resulting image"/>
            <item name="visu_type" type="Enum" initializer="('Module','Jacobian',
                  'X-displacement','Y-displacement','Z-displacement')"
                label="type of visualisation"/>
        </gui> 
                
    """ 
    medipy.medimax.recalage.VisuTrfFile(str(nomfichier),output,VisuTrfNumberInMedimax(visu_type))

#-------------------------------------------------------------
#  findSymmetryPlane
#-------------------------------------------------------------

def findSymmetryPlane(im, imres, filename) :
    """
    Find the symmetry plane of an image according to the sum of squared intensity differences
    The symmetry plane is x=width/2 in the resulting image
        <gui>
            <item name="im" type="Image" label="Input image" />
            <item name="imres" type="Image" initializer="output=True"
                role="output" label="resulting image" />
            <item name="filename" type="String" label="Name of saved .trf file"/>
        </gui>
    """ 
    medipy.medimax.recalage.LinearRegistration(im,im, imres,0,24,1,1,0,1,str(filename),0,2,0)


#-------------------------------------------------------------
#  Atrophy Simulation
#-------------------------------------------------------------

def AtrophySimulation(imref, mask, nomfichres, resolf, lambda_reg) :
    """
    Estimate a deformation field that simulates a give atrophy
        <gui>
            <item name="imref" type="Image" label="Desired level of atrophy" />
            <item name="mask" type="Image" initializer="may_be_empty=True" label="Mask of invariant point" />
            <item name="nomfichres" type="String" label="Name of saved .trf file"/>
            <item name="resolf" type="Int" initializer="value=1, range=(1,7)" label="final resolution"/>
            <item name="lambda_reg" type="Float" initializer="value=1" label="regularization weighting factor"/>
            
        </gui>
    """ 
    medipy.medimax.recalage.simulationAtrophie(imref, mask, str(nomfichres), resolf, lambda_reg)

#-------------------------------------------------------------
#  Apply a transf_3d  o a deformation field
#-------------------------------------------------------------

def warpDeformationField_GUI(nomfichierTrfSrc, nomfichierTrfApply, nomfichierTrfRes, inter_type) :
    """
    Warp a deformation field according to a .trf file
        
        <gui>
            <item name="nomfichierTrfSrc" type="File" label="trf file to warp"/>
            <item name="nomfichierTrfApply" type="File" label="trf file to apply"/>
            <item name="nomfichierTrfRes" type="String" label="resulting trf file"/>
            <item name="inter_type" type="Enum" initializer="('Nearest','Linear',
                'SinCard','QuickSinCard2','QuickSinCard3','Bspline2','Bspline3',
                'Bspline4','Bspline5','Label')" label="Interpolation method"/> 
            
        </gui>
    """ 
    medipy.medimax.recalage.warpDeformationField(str(nomfichierTrfSrc),str(nomfichierTrfApply),str(nomfichierTrfRes),InterpolationNumberInMedimax(inter_type))
