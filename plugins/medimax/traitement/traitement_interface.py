import medipy.medimax.traitement
import medipy.base
import numpy
import medipy.base.exception


#-------------------------------------------------------------
#  Correct Dark Bright Alternation of slices
#-------------------------------------------------------------

def CorrectDarkBrightCoronalAlternation_gui(im, axis) :
    """
    Correct Dark Bright Alternation of slices
    
            <gui>
                <item name ="im" type="Image" label="Image to correct"/>
                <item name ="imres" type="Image" label="Corrected image" initializer="output=True" role="return"/>
                <item name="axis" type="Enum" initializer="('x','y','z')" label="Axis along which slices are corrected"/>
            </gui>

    """ 
    #medipy.components.medimax3.traitement.imx_corr_sagittal_3d_p(im,imres)
    imres=medipy.base.Image(shape=im.shape, dtype=im.dtype)
    imres.copy_information(im)
     
    if str(axis)=='x':
        imres.data = numpy.ascontiguousarray(numpy.swapaxes(im,2,1), dtype=im.dtype)
        medipy.medimax.traitement.imx_corr_sagittal_3d_p(imres,imres)
        imres.data = numpy.ascontiguousarray(numpy.swapaxes(imres,1,2)  , dtype=im.dtype) 

    elif str(axis)=='y':
        medipy.medimax.traitement.imx_corr_sagittal_3d_p(im,imres)
    else :
        imres.data = numpy.ascontiguousarray(numpy.swapaxes(im,0,1), dtype=im.dtype)
        medipy.medimax.traitement.imx_corr_sagittal_3d_p(imres,imres)
        imres.data = numpy.ascontiguousarray(numpy.swapaxes(imres,1,0), dtype=im.dtype)            
      
    
    return imres
    

#-------------------------------------------------------------
#  IntensityNormalisationForRegisteredImages_gui
#-------------------------------------------------------------

def IntensityNormalisationForRegisteredImages_gui(imsrc, imref, mask_imsrc=None, mask_imref=None, method="oneParameterLinearRegression") :
    """
    Normalize intensity between two registered images
    
        <gui>
                <item name ="imsrc" type="Image" label="Image to correct"/>
                <item name ="imref" type="Image" label="Reference image"/>
                <item name ="mask_imsrc" type="Image" label="Mask of Image to correct" initializer="may_be_empty=True, may_be_empty_checked=True"/>
                <item name ="mask_imref" type="Image" label="Mask of Reference image" initializer="may_be_empty=True, may_be_empty_checked=True"/>
                <item name ="imres" type="Image" label="Corrected image" initializer="output=True" role="return"/>
                <item name="method" type="Enum" initializer="('oneParameterLinearRegression','twoParameterLinearRegression','totalLeastSquare')" label="method of correction"/>
        </gui>
    """
    
    if imsrc.shape != imref.shape :
        raise medipy.base.exception.Exception("Images must have the same size") 
    
    if mask_imsrc.shape != imsrc.shape or mask_imref.shape != imref.shape :
        raise medipy.base.exception.Exception("Mask must have the same size") 
    
    dicoMethodNormalisation = {"oneParameterLinearRegression":0, "twoParameterLinearRegression":1, "totalLeastSquare":2}
    
    if mask_imsrc == None :
        mask_imsrc= medipy.base.Image(shape=imsrc.shape, dtype=imsrc.dtype)
        mask_imsrc.data=numpy.copy(imsrc.data)
        mask_imsrc.copy_information(imsrc)
        
    if mask_imref == None :
        mask_imref= medipy.base.Image(shape=imref.shape, dtype=imref.dtype)
        mask_imref.data=numpy.copy(imref.data)
        mask_imref.copy_information(imref)
         
    imres=medipy.base.Image(shape=imsrc.shape, dtype=imsrc.dtype)
    imres.copy_information(imsrc)
    
    medipy.medimax.traitement.IntensityNormalisationForRegisteredImages(imsrc, imref, imres, mask_imsrc,mask_imref, dicoMethodNormalisation[method])
    
    return imres
