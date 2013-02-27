import medipy.medimax.recalage
import  medipy.io
import string
import os.path
import multiprocessing
import numpy as np
import tempfile
import medipy.intensity.normalization

#################################################
##   Some basics functions 
#################################################

def fileToListe(filename,prefixToAdd=None) :
    """
        renvoi une liste d'elements a partir du contenu d'un fichier
    """
    fichier = open(filename,'r')
    liste=fichier.read().splitlines()
    fichier.close()
    
    if prefixToAdd :
        newListe=[]
        for nomImage in liste :
            newListe.append(addPrefixToFilename(nomImage,prefixToAdd))    
        
        liste=newListe
        
    return(liste)
    

def addPrefixToFilename(filename,prefix) :
    base=os.path.basename(filename)
    dir=os.path.dirname(filename)
    if dir :
        return dir+"/"+prefix+base
    else :
        return prefix+base

def lanceMultithreadFunctionWithMultipleParam(param) :
    function=param[0]
    param.pop(0)
    function(*param)
    
#################################################
###  functions for groupwise registration
#################################################

###  Intialisation with affine registration ###

def RecalageAffineElementaire(filenameRef, filenameReca, filenameTrf) :
    imref = medipy.io.load(filenameRef) 
    imreca = medipy.io.load(filenameReca)
    imres = medipy.base.Image(shape=imref.shape, dtype=imreca.dtype)
    medipy.medimax.recalage.LinearRegistration(imref, imreca , imres,2,5,1,1,0,1,filenameTrf,0,2,0)
    

def groupwiseAffineRegistration(ImageTxtFile, ImrefFilename, TransfoIniTxtFile = None) :
    
    ImageList=fileToListe(ImageTxtFile)
    
    if TransfoIniTxtFile :
        TransfoIniList=fileToListe(TransfoIniTxtFile)
        if len(TransfoIniList) != len(ImageList) :
            print "The number of images and transfoIni should be the same !"
            exit(-1)
        else :
            TransfoIniList=fileToListe(TransfoIniTxtFile)
    else :
        TransfoIniList=[]
        for nom in fileToListe(ImageTxtFile) :
            TransfoIniList.append(string.rsplit(addPrefixToFilename(nom, 'aff_'),'.')[0]+".trf")
       
    param=[]
    for i in range(len(ImageList)) :
        param.append([RecalageAffineElementaire,ImrefFilename, ImageList[i], TransfoIniList[i] ])

    pool = multiprocessing.Pool() # Create a group of CPUs to run on
    pool.map(lanceMultithreadFunctionWithMultipleParam, [p for p in param])


### Computation of a series of warped images  ###
def createSeriesOfWarpedImages(ImageTxtFile,TransfoTxtFile, SerieFilename) :
    ImageList=fileToListe(ImageTxtFile)
    TransfoList=fileToListe(TransfoTxtFile)
        
    if len(TransfoList) != len(ImageList) :
        print "The number of images and transfo should be the same !"
        exit(-1)
    
    
    im4D=[] 
    for i in range(len(ImageList)) :
        im = medipy.io.load(ImageList[i])
        imres = medipy.base.Image(dtype=im.dtype)
        medipy.medimax.recalage.ApplyTransfo3d(im,TransfoList[i],imres,8)
        im4D.append(imres)
    
    medipy.io.save_serie(im4D, SerieFilename)   

### Verification that the sum of deformation fields vanishes ###
def verifSumOfDeformationFields(TransfoTxtFile,FilenameImref=None) :
    
    if FilenameImref :
        imref=medipy.io.load(FilenameImref)
    else :
        imref=None
    
    TransfoList=fileToListe(TransfoTxtFile)
    ux = medipy.base.Image()
    uy = medipy.base.Image()
    uz = medipy.base.Image()
    uxmean = medipy.base.Image()
    uymean = medipy.base.Image()
    uzmean = medipy.base.Image()
    
    medipy.medimax.recalage.LoadTrfFile(TransfoList[0], ux, uy,uz, imref)
    uxmean.data=ux.data
    uymean.data=uy.data
    uzmean.data=uz.data
   
        
    for i in range(1,len(TransfoList)) :
        ux = medipy.base.Image()
        uy = medipy.base.Image()
        uz = medipy.base.Image()
        medipy.medimax.recalage.LoadTrfFile(TransfoList[i], ux, uy,uz, imref)
        uxmean.data=uxmean.data+ux.data
        uymean.data=uymean.data+uy.data
        uzmean.data=uzmean.data+uz.data

    nb=float(len(TransfoList))
    uxmean.data=uxmean.data/nb
    uymean.data=uymean.data/nb
    uzmean.data=uzmean.data/nb
    
    print "ux "+str(np.amax(uxmean.data))+" / "+str(np.amin(uxmean.data))
    print "uy "+str(np.amax(uymean.data))+" / "+str(np.amin(uymean.data))
    print "uz "+str(np.amax(uzmean.data))+" / "+str(np.amin(uzmean.data))
    
    mean = medipy.base.Image()
    mean = np.sqrt(uxmean.data**2+uymean.data**2+uzmean.data**2)
    print "mean of all displacement field "+str(np.mean(mean)) 
        
###  Groupwise registration ###   
def groupwiseRegistrationFromTxtfile(ImageTxtFile, TransfoIniTxtFile=None, TransfoResTxtFile=None, resolution=6, regularisation=1,
wdth_template=256,hght_template=256, dpth_template=256,
dx_template=1, dy_template=1, dz_template=1, serieOfTemplate=None) :
    
    ImageList=fileToListe(ImageTxtFile)
    nb_image=len(ImageList)

    ## creation des noms de fichiers par defaut si les listes de transfo ne sont pas renseignees
    if TransfoIniTxtFile :
        TransfoIniList=fileToListe(TransfoIniTxtFile)
        if len(TransfoIniList) != nb_image :
            print "The number of images and transfoIni should be the same !"
            exit(-1)
        else :
            TransfoIniList=fileToListe(TransfoIniTxtFile)
    else :
        TransfoIniList=[]
        for nom in fileToListe(ImageTxtFile) :
            TransfoIniList.append(string.rsplit(addPrefixToFilename(nom, 'aff_'),'.')[0]+".trf")
    
    
    if TransfoResTxtFile :
        TransfoResList=fileToListe(TransfoResTxtFile)
        if len(TransfoResList) != nb_image :
            print "The number of images and TransfoRes should be the same !"
            exit(-1)
        else :
            TransfoResList=fileToListe(TransfoResTxtFile)
    else :
        TransfoResList=[]
        for nom in fileToListe(ImageTxtFile) :
            TransfoResList.append(string.rsplit(addPrefixToFilename(nom, 'groupwise_'),'.')[0]+".trf")
    
    
    #------------------------------------------------
    #---- Generation du fichier trf permettant de 
    #---- Regler les caracteristiques de l'espace du template
    #------------------------------------------------
    resample=tempfile.mkstemp(".trf")[1]
    f = open(resample, "w")
    f.write("TRANSF3D          \n")
    f.write("type=AFFINE3D        \n")
    f.write("nb_param=15       \n")
    f.write("width="+str(wdth_template)+"  \n")    
    f.write("height="+str(hght_template)+"    \n")
    f.write("depth="+str(dpth_template)+" \n")
    f.write("dx="+str(dx_template)+"       \n")
    f.write("dy="+str(dy_template)+"       \n")
    f.write("dz="+str(dz_template)+"       \n")
    f.write("a00=1             \n")
    f.write("a01=0             \n")
    f.write("a02=0             \n")
    f.write("a10=0             \n")
    f.write("a11=1             \n")
    f.write("a12=0             \n")
    f.write("a20=0             \n")
    f.write("a21=0            \n")
    f.write("a22=1            \n")
    f.write("Tx=0          \n")
    f.write("Ty=0          \n")
    f.write("Tz=0          \n")
    f.write("Xg=0.000000       \n")
    f.write("Yg=0.000000       \n")
    f.write("Zg=0.000000       \n")

    f.close()

    if serieOfTemplate  :
        im4D=[] 
   

    #------------------------------------------------
    #---- Initialisation du cerveau moyen -----------
    #------------------------------------------------
    print('#------------------------------------------------\n')
    print("Nombre total d'images : "+str(nb_image)+"\n")
    print('#------------------------------------------------\n')
    
    print('#------------------------------------------------\n')
    print('Initialisation du cerveau moyen\n')
    print('#------------------------------------------------\n')
    
    imfirst=medipy.io.load(ImageList[0])
    meanImage = medipy.base.Image(dtype=imfirst.dtype)
    medipy.medimax.recalage.CombineTransfo3d(TransfoIniList[0],resample, TransfoResList[0], 5)
    medipy.medimax.recalage.ApplyTransfo3d(imfirst,TransfoResList[0],meanImage,8)

    if serieOfTemplate :
        im4D.append(meanImage)

    bsplineTmp=tempfile.mkstemp(".trf")[1]
    
    for i in range(1,nb_image) :
        
        #------------------------------------------------
        #----       Recalage          -------
        #------------------------------------------------
        poids = float(i)/float(i+1) 
        
        print('#------------------------------------------------\n')
        print("Recalage de la "+str(i+1)+"eme image \n")
        print('#------------------------------------------------\n')
    
        medipy.medimax.recalage.CombineTransfo3d(TransfoIniList[i],resample, TransfoResList[i], 5)
        im=medipy.io.load(ImageList[i])
        imreca=medipy.base.Image(dtype=im.dtype)
        imres=medipy.base.Image(dtype=im.dtype)
        medipy.medimax.recalage.ApplyTransfo3d(im,TransfoResList[i],imreca,8)
        #print "MeanImage" + str(meanImage.shape)+ " "+str(meanImage.spacing)
        #print "imreca" + str(imreca.shape)+ " "+str(imreca.spacing)
        
        
        medipy.medimax.recalage.BsplineRegistration3d(meanImage, imreca, imres,0,2,regularisation,4,1,bsplineTmp,resolution, 0.0,100000.0,12, 10,0, 1, poids)
        medipy.medimax.recalage.CombineTransfo3d(TransfoResList[i],bsplineTmp.replace('.trf','_reca.trf'), TransfoResList[i], 5)
        
        #------------------------------------------------
        #----  Mise a jour des champs de deformation  ---
        #------------------------------------------------
        print('#------------------------------------------------\n')
        print("Mise a jour des champs de deformation \n")
        print('#------------------------------------------------\n')
    
        for j in range(i):
            medipy.medimax.recalage.CombineTransfo3d(TransfoResList[j],bsplineTmp.replace('.trf','_ref.trf'), TransfoResList[j], 5)
            
        #------------------------------------------------
        #----   Mise a jour de l'image moyenne  -------
        #------------------------------------------------
        print('#------------------------------------------------\n')
        print("Mise a jour de l'image moyenne \n")
        print('#------------------------------------------------\n')
    
        imref = medipy.base.Image(dtype=imfirst.dtype)
        meanImage = medipy.base.Image(dtype=np.float32)
        medipy.medimax.recalage.ApplyTransfo3d(imfirst,TransfoResList[0],imref,8)
        meanImage.data=np.copy(imref.data)
        meanImage.copy_information(imref)
        
        for j in range(1,i+1) :
            im=medipy.io.load(ImageList[j])
            imres=medipy.base.Image(dtype=im.dtype)       
            medipy.medimax.recalage.ApplyTransfo3d(im,TransfoResList[j],imres,8)
            meanImage.data=meanImage.data+medipy.intensity.normalization.mean_stdev_normalization(imref, imres, imref, imres).data
   
        meanImage.data=meanImage.data/float(i+1)
        
        if serieOfTemplate :
            im4D.append(meanImage)
    
    
    if serieOfTemplate :
        medipy.io.save_serie(im4D, serieOfTemplate) 
    
    # suppression fichiers temporaires
    os.remove(resample)
    os.remove(bsplineTmp)
    os.remove(bsplineTmp.replace('.trf', '.chp'))
    os.remove(bsplineTmp.replace('.trf', '_reca.trf'))
    os.remove(bsplineTmp.replace('.trf', '_reca.prm'))
    os.remove(bsplineTmp.replace('.trf', '_ref.trf'))
    os.remove(bsplineTmp.replace('.trf', '_ref.prm'))
    
    
    

if __name__ == "__main__" :
    #groupwiseRegistrationFromTxtfile("/home/miv/noblet/tmp/groupwise/listeImage.txt", "/home/miv/noblet/tmp/groupwise/listeAffID.txt", "/home/miv/noblet/tmp/groupwise/listeGroupwise.txt",resolution=4,regularisation=0,wdth_template=64,hght_template=64, dpth_template=64, dx_template=4, dy_template=4, dz_template=4)
    #groupwiseAffineRegistration("/home/miv/noblet/tmp/groupwise/listeImage.txt", "/home/miv/noblet/data/MNI_20patients_IPB_128/subject54.ipb", "/home/miv/noblet/tmp/groupwise/listeAff.txt")
    #createSeriesOfWarpedImages("/home/miv/noblet/tmp/groupwise/listeImage.txt", "/home/miv/noblet/tmp/groupwise/listeGroupwise.txt", "/home/miv/noblet/tmp/groupwise/serie.nii.gz")
    #verifSumOfDeformationFields("/home/miv/noblet/tmp/groupwise/listeGroupwise.txt")
    pass
