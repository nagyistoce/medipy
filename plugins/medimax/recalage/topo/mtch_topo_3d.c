/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/**-----------------------------------------------------------------------------
*** 
*** file:       mtch_topo_3d.c
***
*** project:    Imagix 1.01
***         
***
*** description:    Fichier pourle matching deformable des images 3D
***                     avec conservation de la topologie
*** 
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.
*** All rights are reserved.
***
***
***---------------------------------------------------------------------------*/
#include <config.h>
#include <stdio.h> 
#include <stdlib.h>     /* Needed when you use GET_FLOAT    */
#include <time.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/imx_lang.h"
#include "noyau/mani_3d.h"
#include "noyau/io/imx_log.h"
#include "noyau/io/imx_file.h"
#include "noyau/io/file_3d.h"
#include "noyau/io/imx_head.h"
#include "noyau/gui/imx_picture3d.h"

#include "math/imx_matrix.h"
#include "math/imx_bspline.h"   /* pour les fonctions d'interpolation */
#include "math/oper_3d.h"

#include "traitement/norm_3d.h"
#include "traitement/trai_3d.h"
#include "traitement/differential_bias_correction.h"

#include "morpho/morpho_3d.h"

#ifdef  HAS_ITK
#include "interface_itk/itk_algorithms.hpp"
#endif

#include "outils/imx_misc.h"

#include "recalage/mtch_3d.h"
#include "recalage/chps_3d.h"
#include "recalage/distance_3d.h"
#include "recalage/mtch_robust_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/validation_topo_3d.h"
#include "recalage/topo/distance_topo.h"
#include "recalage/topo/gradient_distance_topo.h"
#include "recalage/topo/hessien_distance_topo.h"
#include "recalage/topo/optimisation_topo.h"
#include "recalage/topo/analyse_intervalle.h"
#include "recalage/topo/normalisation_topo.h"

#ifndef COMPILE_FOR_MEDIPY
#include "dti/dti_recalage.h"
#endif //fin COMPILE_FOR_MEDIPY

int RENVERSEMENT_SEUIL=0; /* 0 par defaut (renversement du seuil) 1 on prend le seuil entre fm et Jm */ 

HistoJoint HISTOJOINT;

int LOG_ENERGY;  // variable globale pour identifier un fichier de log dans la fonction compare_critere
double  TOPO_ALPHA_ROBUST; /* parametre de la fonction de cout robuste de Geman Mclure*/
double TOPO_REGULARISATION_SIGMA_GAUSSIEN;
double TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
double TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE;
double TOPO_SEUIL_SIGNIFICATIVITE;
melange_gaussienne_2d MELANGE_GAUSSIENNE_2D;
double TOPO_PONDERATION_RECALAGE_SYMETRIQUE=1.0;
int PRECISION_SIMULATION;
ParamRecalageBspline _ParamRecalageBspline;

double PATCH_REF;
double PATCH_RECA;

double TOPO_PATCH_COEFF;
int TOPO_NLMhwnx;
int TOPO_NLMhwny;
int TOPO_NLMhwnz;
int TOPO_NLMhwvsx;
int TOPO_NLMhwvsy;
int TOPO_NLMhwvsz;


int TOPO_PATCH_SIZE=1;
int TOPO_PATCH_NEIGHBOR=2;

/*******************************************************************************
********************************************************************************
************************* RECALAGE BSPLINE *************************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**        matching_Bspline_topo_3d                                                
**                                                                        
**       recalage Bspline de deux images 3D avec conservation de la topologie                               
*******************************************************************************/

void matching_Bspline_topo_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[15],*nomfichres;
 int dist_type,inter_type,min_type,save_type,func_type,i,normalisation_type, reg_type=0;
 int nb_classe;
 double Jmin,Jmax;
 int resolf, e, err;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref= GET_PLACE3D(TEXT0031);
 im_res= GET_PLACE3D(TEXT0006);

 for(i=0;i<15;i++)
    quest[i]=CALLOC(80,char);
 
 /*type de fonction*/
  func_type=1;
  
 /*distance*/
  strcpy(quest[0],"L2");
  strcpy(quest[1],"L2 Sym");
  strcpy(quest[2],"Lp");
  strcpy(quest[3],"Lp Sym");
    strcpy(quest[4],"L1L2");
    strcpy(quest[5],"L1L2 Sym");
    strcpy(quest[6],"L1norm");
    strcpy(quest[7],"L2 Topo");
    strcpy(quest[8],"L2 Sym Topo");
    strcpy(quest[9],"Geman McLure");
    strcpy(quest[10],"Geman McLure Sym");
    strcpy(quest[11],"L2 pondere");
    strcpy(quest[12],"IM");
    strcpy(quest[13],"\0");

  dist_type=GETV_QCM("Distance",(const char **)quest);
  
    /*interpolation*/
    inter_type=4;

    /*minimisation*/
  strcpy(quest[0],"ICM");
  strcpy(quest[1],"Descente gradient");
if ((dist_type!=2)&&(dist_type!=3)&&(dist_type!=6))
    {
  strcpy(quest[2],"Levenberg-Marquardt");
  strcpy(quest[3],"Levenberg-Marquardt sans Topo");
  strcpy(quest[4],"Levenberg-Marquardt Parallel");

    strcpy(quest[5],"\0");
    }
else
    strcpy(quest[2],"\0");
    
  min_type=GETV_QCM("Minimisation",(const char **)quest);

   
    /*question sur l'enregistrement du champ resultat*/
    nomfichres=quest_save_result_3d(&save_type);
    

    /*resolution finale*/
    resolf= GET_INT("resolution", 4, &e);

    /*bornes sur le jacobien*/
    Jmin = 2;
    while (Jmin>=1)  
        Jmin = GET_DOUBLE("Jmin<1",0,&err);

    Jmax = 0;
    while (Jmax<=1)  
        Jmax = GET_DOUBLE("Jmax>1",100000,&err);

  nb_classe=-1;
    normalisation_type=2;

    for(i=0;i<15;i++)
    free(quest[i]);
        
imx_matching_Bspline_topo_3d(im_ref,im_reca,im_res,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax,normalisation_type,nb_classe,1,0);

 
 if (nomfichres)
  free(nomfichres);

  
 show_picture_3d(im_res);
 
}

/*******************************************************************************
**        matching_Bspline_topo_3d                                                
**                                                                        
**       recalage Bspline de deux images 3D avec conservation de la topologie                               
*******************************************************************************/

void matching_Bspline_topo_norm_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[20],*nomfichres;
 int dist_type,inter_type,min_type,save_type,func_type,i,normalisation_type;
 int nb_classe,adaptatif,reg_type,biais;
 double Jmin,Jmax;
 int resolf, e, err;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref= GET_PLACE3D(TEXT0031);
 im_res= GET_PLACE3D(TEXT0006);

 for(i=0;i<20;i++)
    quest[i]=CALLOC(80,char);
 
 /*type de fonction*/
  func_type=1;
  
 /*distance*/
  strcpy(quest[0],"L2");
  strcpy(quest[1],"L2 Sym");
  strcpy(quest[2],"Lp");
  strcpy(quest[3],"Lp Sym");
    strcpy(quest[4],"L1L2");
    strcpy(quest[5],"L1L2 Sym");
    strcpy(quest[6],"L1norm");
    strcpy(quest[7],"L2 Topo");
    strcpy(quest[8],"L2 Sym Topo");
    strcpy(quest[9],"Geman McLure");
    strcpy(quest[10],"Geman McLure Sym");
    strcpy(quest[11],"L2 pondere");
    strcpy(quest[12],"IM");
    strcpy(quest[13],"Patch-based");
    strcpy(quest[14],"\0");

  dist_type=GETV_QCM("Distance",(const char **)quest);
  
    /*interpolation*/
    inter_type=4;

    /*minimisation*/
  strcpy(quest[0],"ICM");
  strcpy(quest[1],"Descente gradient");
if ((dist_type!=2)&&(dist_type!=3)&&(dist_type!=6))
    {
  strcpy(quest[2],"Levenberg-Marquardt");
  strcpy(quest[3],"Levenberg-Marquardt sans Topo");
  strcpy(quest[4],"Levenberg-Marquardt Parallel");
    strcpy(quest[5],"\0");
    }
else
    strcpy(quest[2],"\0");
    
  min_type=GETV_QCM("Minimisation",(const char **)quest);

   
    /*question sur l'enregistrement du champ resultat*/
    nomfichres=quest_save_result_3d(&save_type);
    

    /*resolution finale*/
    resolf= GET_INT("resolution", 4, &e);

    /*bornes sur le jacobien*/
    Jmin = 2;
    while (Jmin>=1)  
        Jmin = GET_DOUBLE("Jmin<1",0,&err);

    Jmax = 0;
    while (Jmax<=1)  
        Jmax = GET_DOUBLE("Jmax>1",100000,&err);

  /* normalisation d'intensite */
  strcpy(quest[0],"Pas de normalisation");
  strcpy(quest[1],"ecart des moyennes");
  strcpy(quest[2],"Moyenne et ecart-type");
  strcpy(quest[3],"Regression lineaire a 2 parametres");
    strcpy(quest[4],"Rapport des moyennes");
    strcpy(quest[5],"Rapport Min/Max");
    strcpy(quest[6],"Mixture gaussiennes histogramme");
    strcpy(quest[7],"segmentation kmeans");
    strcpy(quest[8],"segmentation markov");
    strcpy(quest[9],"histogramme joint (lineaire par morceau)");
    strcpy(quest[10],"histogramme joint (spline)");
    strcpy(quest[11],"histogramme joint (mixture de gaussienne) standart");
    strcpy(quest[12],"histogramme joint (mixture de gaussienne) norme");
    strcpy(quest[13],"histogramme joint (mixture de gaussienne) reclasse");
    strcpy(quest[14],"Multimodal histogramme joint (mixture de gaussienne)");
    strcpy(quest[15],"EQUALISATION HISTO");
    strcpy(quest[16],"Moyenne et ecart-type robuste");
    strcpy(quest[17],"normalisation segmentation sous-jacente");
    #ifdef  HAS_ITK
    strcpy(quest[18],"normalisation histogramme par quantile");
    strcpy(quest[19],"\0");
    #else
    strcpy(quest[18],"\0");
    #endif
    
 normalisation_type=GETV_QCM("Normalisation",(const char **)quest);
 nb_classe=-1;
if ((normalisation_type>5)&&(normalisation_type!=10)&&(normalisation_type!=16))
    {
    if (normalisation_type==18)
    nb_classe = GET_INT("Nombre de quantile",7,&err);
    else
        if (normalisation_type==15)
            nb_classe = GET_INT("Nombre de bins",1024,&err);
        else
        nb_classe = GET_INT("Nombre de classe",10,&err);
    }
  

  /* correction de biais */
    strcpy(quest[0],"Non");
    strcpy(quest[1],"Oui");
    strcpy(quest[2],"\0");
    
    biais= GETV_QCM("Correction de biais ?",(const char **)quest);

    
  /* partition adaptative ? */
    strcpy(quest[0],"Oui");
    strcpy(quest[1],"Non");
    strcpy(quest[2],"\0");
    
    adaptatif= GETV_QCM("Adaptatif ?",(const char **)quest);


    /* regularisation */
    strcpy(quest[0],"Pas de regularisation");
  strcpy(quest[1],"Filtrage gaussien");
  strcpy(quest[2],"Regularisation membrane elastique");
  strcpy(quest[3],"Regularisation membrane elastique Lp");
  strcpy(quest[4],"Regularisation Log(J)^2");
  strcpy(quest[5],"Regularisation membrane elastique Log(J)");
  strcpy(quest[6],"Regularisation Log(J-Jmean)^2");
  strcpy(quest[7],"Regularisation identite");
  strcpy(quest[8],"Regularisation imagebased patch");
  strcpy(quest[9],"\0");
        
    reg_type= GETV_QCM("Regularisation ?",(const char **)quest);

    
    TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
    if (reg_type==1)
        TOPO_REGULARISATION_SIGMA_GAUSSIEN = GET_DOUBLE("sigma",1,&err);
    
    
    if (reg_type>1)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&err);
    
    
    if (reg_type==8)
        {
        TOPO_PATCH_SIZE=GET_INT("Taille patch (2n+1)",1,&err);
        TOPO_PATCH_NEIGHBOR=GET_INT("Taille voisinage (2n+1)",2,&err);
        }
    
    for(i=0;i<20;i++)
    free(quest[i]);
        
    


imx_matching_Bspline_topo_3d(im_ref,im_reca,im_res,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax,normalisation_type,nb_classe,adaptatif,biais);    
    
 
 if (nomfichres)
  free(nomfichres);

  
 show_picture_3d(im_res);
 
}


/*******************************************************************************
**        imx_matching_Bspline_topo_3d                                             
**                                                                       
**       recalage Bspline de 2 images 3D avec conservation de la topologie                                  
*******************************************************************************/
/*
  im_ref : numero de l'image sur laquelle on recale imreca
  im_reca : numero de l'image a recaler
  im_res : numero de l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  func_type : numero de la fonction d'interpolation du champ (1 pour Bspline
    degre 1)
  dist_type : numero de la fonction de distance (0 pour distance quadratique)
  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
*/
int imx_matching_Bspline_topo_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif,int biais)
{
int l,continu;
grphic3d *imref,*imreca,*imres;
imref=ptr_img_3d(im_ref);
imreca=ptr_img_3d(im_reca);
imres=ptr_img_3d(im_res);

 if ((imref->dx!=imreca->dx)||(imref->dy!=imreca->dy)||(imref->dz!=imreca->dz))
    {
      PUT_WARN("Les images n'ont pas les memes dx,dy,dz !");
      printf("Les images n'ont pas les memes dx,dy,dz !\n");
      
            return(-1) ;
    }

 if ((imref->width!=imreca->width)||(imref->height!=imreca->height)||(imref->depth!=imreca->depth))
    {
      PUT_WARN("Les images n'ont pas la meme taille !");
      printf("Les images n'ont pas la meme taille !\n");
      return(-1) ;
    }


continu=0;
for (l=0;l<10;l++)
    {
    if (imref->width==pow(2.0,l)) continu++;
    if (imref->height==pow(2.0,l)) continu++;
    if (imref->depth==pow(2.0,l)) continu++;
    if (imreca->width==pow(2.0,l)) continu++;
    if (imreca->height==pow(2.0,l)) continu++;
    if (imreca->depth==pow(2.0,l)) continu++;
    }

if (continu!=6)
    {
      PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
          printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
            //return(-1) ;
            save_type=3;
   }

 imx_matching_Bspline_topo_3d_p(imref,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax,normalisation_type, nb_classe,adaptatif,biais);
 
  return(1);
}


/*******************************************************************************
**        imx_matching_Bspline_topo_3d_p                                            
**                                                                        
**       recalage Bspline de 2 images 3D avec conservation de la topologie                                   
*******************************************************************************/
/*
  imref : pointeur sur l'image sur laquelle on recale imreca
  imreca : pointeur sur l'image a recaler
  imres : pointeur sur l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  func_type : numero de la fonction d'interpolation du champ (1 pour Bspline
    degre 1)
  dist_type : numero de la fonction de distance (0 pour distance quadratique)
  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
*/
int imx_matching_Bspline_topo_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif, int biais)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb;
 double *param,*param_norm;
 grphic3d *imtref,*imtreca,*imtres;
 field3d *champ;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 norm_func_locale_t normalisation;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;

#ifndef SAVE_INTERMEDIAIRE_2 
char nomfichier[255];
char temp[2];
char *ext2;
#endif

#ifndef TOPO_COMPARE_CRITERE 
    char nomfichier_nrj[255];
    char *ext;
#endif
#ifndef WIN32
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif


    wdth=imref->width;hght=imref->height;dpth=imref->depth;
    
    if (save_type==3) /* il faut transformer les images en puissance de 2 */
        {
        wdth_old=wdth; hght_old=hght; dpth_old=dpth;
        
        imx_realloc_mri_pow2_centered_p(imref);
        imx_realloc_mri_pow2_centered_p(imreca);
        //imx_realloc_mri_pow2_centered_p(imres);
        
        if (imreca->mask!=NULL)
            imx_realloc_mri_pow2_centered_p(imreca->mask);
    
        if (imref->mask!=NULL)
            imx_realloc_mri_pow2_centered_p(imref->mask);
        
        //if (imres->mask!=NULL)
         //   imx_realloc_mri_pow2_centered_p(imres->mask);
        
        wdth=imref->width;hght=imref->height;dpth=imref->depth;
        
        }
    
    tdebut=time(NULL);
    nbmax_param=3*pow((pow(2,resolf)-1),3);
 
   /* Allocation */
  if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
    if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
 init_log(NULL);
 aff_log("\nRECALAGE ");
 
 /*************mise en place des parametres du recalage******/
scal_func=topo_choose_bspline(func_type);
distance=topo_choose_distance(dist_type);
interpol=imx_choose_interpolation_fct(inter_type);  
minimisation=topo_choose_optimisation(min_type);
normalisation=topo_choose_normalisation(normalisation_type,nb_classe);
regularisation=topo_choose_regularisation(reg_type);
    
    aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
  /*creation des images temporaires*/
//  if (normalisation_type>0)
    {
    imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
    imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
        imx_norm_rcoeff_3d_p(imtref, imtreca);

        if ((imreca->mask!=NULL)&&(imref->mask!=NULL))
            {   imtreca->mask=cr_grphic3d(imreca->mask);imx_copie_3d_p(imreca->mask,imtreca->mask);
                imtref->mask=cr_grphic3d(imref->mask);imx_copie_3d_p(imref->mask,imtref->mask); }
        }
  /*else
    {
    imtref = imref;
    imtreca = imreca;
        imx_norm_rcoeff_3d_p(imtref, imtreca);
        }*/
        
  imtres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtres);

    if (imreca->mask!=NULL)
        {imtres->mask=cr_grphic3d(imreca->mask);imx_copie_3d_p(imreca->mask,imtres->mask);}
    

    /*allocation memoire de champ*/
   champ=cr_field3d(wdth,hght,dpth);
     

if (dist_type==12)
 init_HistoJoint(TOPO_SIZE_HISTOJOINT,TOPO_SIZE_HISTOJOINT);
    
  //calcul des variables globales pour la distance bas�e patch  
  if(dist_type==13)
  {
  
  int nx = imref->width;
  int ny = imref->height;
  int nz = imref->depth;
  int nbp = (nx-2)*(ny-2)*(nz-2);
  double ei_ref = 0;
  double ei_reca = 0;
  double sig_ref = 0;
  double sig_reca = 0;
  double beta = 0.5;
  int hwn = 3*3*3;
  int x,y,z;
 
  PATCH_REF = 0;
  PATCH_RECA = 0;
  
  for(x=1;x<nx-1;x++)
    for(y=1;y<ny-1;y++)
      for(z=1;z<nz-1;z++)
      {
      
       ei_ref = sqrt(6/7.0)*(imref->mri[x][y][z]  
-(imref->mri[x+1][y][z]+imref->mri[x-1][y][z]+imref->mri[x][y+1][z]+imref->mri[x][y-1][z]+imref->mri[x][y][z+1]+imref->mri[x][y][z-1])/6.0);

       ei_reca = sqrt(6/7.0)*(imreca->mri[x][y][z]  
-(imreca->mri[x+1][y][z]+imreca->mri[x-1][y][z]+imreca->mri[x][y+1][z]+imreca->mri[x][y-1][z]+imreca->mri[x][y][z+1]+imreca->mri[x][y][z-1])/6.0);

       sig_ref += ei_ref * ei_ref;
       sig_reca += ei_reca * ei_reca;
       
       

      }
      
  PATCH_REF = 2*beta*sig_ref*hwn / nbp;
  PATCH_RECA= 2*beta*sig_reca*hwn / nbp;
  
  } 
  
  if (reg_type == 8) // regularisation image-based patch
    {
    
  float minVoxSz = imref->dx;
  if(imref->dy < minVoxSz) minVoxSz = imref->dy;
  if(imref->dz < minVoxSz) minVoxSz = imref->dz;

 
    TOPO_NLMhwnx = (int)(0.5 + TOPO_PATCH_SIZE * minVoxSz / imref->dx);
    TOPO_NLMhwny = (int)(0.5 + TOPO_PATCH_SIZE * minVoxSz / imref->dy);;
    TOPO_NLMhwnz = (int)(0.5 + TOPO_PATCH_SIZE * minVoxSz / imref->dz);;
    
    TOPO_NLMhwvsx=(int)(0.5 + TOPO_PATCH_NEIGHBOR * minVoxSz / imref->dx);
    TOPO_NLMhwvsy=(int)(0.5 + TOPO_PATCH_NEIGHBOR * minVoxSz / imref->dy);
    TOPO_NLMhwvsz=(int)(0.5 + TOPO_PATCH_NEIGHBOR * minVoxSz / imref->dz);  
    
    TOPO_PATCH_COEFF=NLMSmoothComputation3D(imref, TOPO_NLMhwnx, TOPO_NLMhwny, TOPO_NLMhwnz, TOPO_PATCH_BETA, 0);
    
    _ParamRecalageBspline.ref_groupwise=imref;
    }
  


 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r,D;
  int ***masque_param;
  int *x0,*x1,*y0,*y1,*z0,*z1;
  int flag,compt;
  double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
  double sigma;
  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;
        
  /* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}

    if (nomfichres!=NULL)
        {
        #ifndef TOPO_COMPARE_CRITERE 
        strcpy(nomfichier_nrj,nomfichres);
    ext=strstr(nomfichier_nrj,".trf");
            if (ext!=NULL) *ext='\0';
        strcat(nomfichier_nrj,"_energy.txt");
        LOG_ENERGY=init_flog(nomfichier_nrj);      
        #endif
     }
     
      
    /* Debut traitement */
    for (r=resol;r<=resolf;r++)
  {


//for (z=0;z<2;z++)
//{

TOPO_SEUIL_SIGNIFICATIVITE=0.05;


    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    base_to_field_3d(nb_param,param,champ,NULL,NULL);
    interpol(imreca,champ,imtres); 
    
    if (imreca->mask!=NULL)
        inter_nearest_3d(imreca->mask,champ,imtres->mask);

    if (imtref->mask==NULL)
            {
            imtref->mask=cr_grphic3d(imtref);imx_copie_3d_p(imtref,imtref->mask);
            
            for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
            for (k=0;k<dpth;k++)
                    imtref->mask->mri[i][j][k]=1;
            }


    if (imtres->mask==NULL)
            {
            imtres->mask=cr_grphic3d(imtres);imx_copie_3d_p(imtres,imtres->mask);
            
            for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
            for (k=0;k<dpth;k++)
                    imtres->mask->mri[i][j][k]=1;
            }
            
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    topo_masque_param (imtref,imtreca,imtres,masque_param, nb_param);
    
    //annule_param_avec_masque_param (param, nb_param,masque_param);
    
    /* Dans le cas de l'IM, je choisi de conduire la derni�e r�olution avec ma normalisation d'intensit�*/
/*  if ((r>=(resolf-1))&&(distance==Energie_IM_locale_3d))
        {
        distance=Energie_quad_locale_3d;
        normalisation_type=14;
        normalisation=imx_normalisation_gaussienne2d_histo_joint_multimodal_p;
        nb_classe=20;
        }*/

    /*normalisation de imreca par rapport a imref: resultat dans imtreca*/ 
    if (normalisation_type>0)
    {
        printf("Debut normalisation\n");
        normalisation(imref,imtres,imtref,nb_classe);
        printf("Fin normalisation\n");
        imx_norm_rcoeff_3d_p(imtref, imtreca);
    }
    
    
    /*correction de biais de imreca par rapport a imref: resultat dans imtreca*/ 
    if (biais>0)
    {   
        grphic3d *im_tmp=cr_grphic3d(imtref);

        printf("Debut corrections de biais\n");
        imx_differential_bias_correction(imtref, imtres, im_tmp, 15, 0.5, 0.5);
        printf("Fin corrections de biais\n");
        imx_copie_3d_p(im_tmp,imtref);
        imx_norm_rcoeff_3d_p(imtref, imtreca);
        free_grphic3d(im_tmp);
    }
    
            
    if (imtref->rcoeff!=imtreca->rcoeff)
        printf("Attention les rcoeff ne sont pas identiques !!! Le recalage ne sera pas pertinent ....\n");
    
    if (adaptatif==0)
        test_correlation(imtref,imtres,nb_param,masque_param);

    if ((dist_type==9)||(dist_type==10))
        update_TOPO_ALPHA_ROBUST_global(imtref,imtreca, nb_param, param);

    if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
    else
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;

        
    if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)//&&(r==1))
        //update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE(imtref,imtreca, nb_param, param, param_norm, distance, regularisation,masque_param);
        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(imtref,imtreca,imtres,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);


    if (dist_type==12)
            update_maxnorm_HistoJoint(imtref,imtreca);

if (TOPO_REGULARISATION_SIGMA_GAUSSIEN != 0)
{
sigma=TOPO_REGULARISATION_SIGMA_GAUSSIEN/(x1[0]-x0[0]);
compt=0;
do
{
mini=minimisation(imtref,imtreca,imtres,champ,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);

imx_filtre_param_gaussien_3d_p(param, param_norm, nb_param,sigma);

flag=0;
    for(i=0;i<D;i++)
        for(j=0;j<D;j++)
            for(k=0;k<D;k++)
            if (masque_param[i][j][k]>0)
            {
            l=TOP_conv_ind(i,j,k,nb_param);
            if (1.0*fabs(param[l]-param_norm[l])/param[l]>0.05)
                flag++;
            else
                masque_param[i][j][k]=0;
            }
printf ("Il y a %f pourcent de param qui varie de plus de 10 pourcent\n",1.0*flag/nb_param);

for (i=0;i<nb_param;i++)
            param[i]=param_norm[i];
            
compt++;
            
}while ((flag>0.01*nb_param) && (compt<10));
}
else
 mini=minimisation(imtref,imtreca,imtres,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);

    #ifndef SAVE_INTERMEDIAIRE_2
    
if (nomfichres!=NULL)
   {
     
        /*calcul de la transformation*/
        base_to_field_3d(nb_param,param,champ,NULL,NULL);

     //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichres);
        ext2=strstr(nomfichier,".trf");
        if (ext2!=NULL) *ext2='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(r+48);
        temp[1]=0;
        strcat(nomfichier,temp);
      
    transf3d *transfo; 
    
    if ((save_type==0)||(save_type==3))
        {
            if (save_type==3) /* on tronque le champ */
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                imref->width=   wdth_old;
                imref->height=  hght_old;
                imref->depth=   dpth_old;
        
                imreca->width=  wdth_old;
                imreca->height= hght_old;
                imreca->depth=  dpth_old;
        
                
                } 
        
                    
         transfo=field_to_transf_3d(champ,imref,imreca); 
            
                if (save_type==3)
                    {
                champ->width=   wdth;
                champ->height=  hght;
                champ->depth=   dpth;
                
                imref->width=   wdth;
                imref->height=  hght;
                imref->depth=   dpth;
        
                imreca->width=  wdth;
                imreca->height= hght;
                imreca->depth=  dpth;
                    
                    }
        
        }
        else
    {  
     
     transfo=cr_transf3d(wdth,hght,dpth,NULL);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
         transfo->dx=imres->dx;
         transfo->dy=imres->dy;
         transfo->dz=imres->dz;
         
         for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*imreca->dx;
            transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
            transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
            } 
    }
    save_transf_3d(transfo,nomfichier);
    free_transf3d(transfo);
     
    

   } 
    #endif
 
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
    free_imatrix_3d(masque_param);
    }       

    if (dist_type==12)
        free_HistoJoint();
    
    /*calcul de la transformation*/
    base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
 if (imres != NULL)
    { 
            /*imtres->width=wdth_old;
            imtres->height=hght_old;
            imtres->depth=dpth_old;*/
            
            interpol=inter_qsinc3_3d;
        interpol(imreca,champ,imtres);
        imx_copie_3d_p(imtres,imres);
    }

  

    if (save_type==3) /* il faut transformer les images  dans leur taille d'origine*/
        {
    
        imx_undo_mri_pow2_centered_p(imref,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(imreca,wdth_old,hght_old,dpth_old);
        
        if ((imres != imreca) && (imres != imref))
            {imx_undo_mri_pow2_centered_p(imres,wdth_old,hght_old,dpth_old);
            if (imres->mask!=NULL)
                imx_undo_mri_pow2_centered_p(imres->mask,wdth_old,hght_old,dpth_old);
            }
    
        if (imreca->mask!=NULL)
        imx_undo_mri_pow2_centered_p(imreca->mask,wdth_old,hght_old,dpth_old);
        
        if (imref->mask!=NULL)
        imx_undo_mri_pow2_centered_p(imref->mask,wdth_old,hght_old,dpth_old);
        
        
        
        }
  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 
    transf3d *transfo; 
    
    if ((save_type==0)||(save_type==3))
        {
            if (save_type==3) /* on tronque le champ */
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                } 
        
        
         transfo=field_to_transf_3d(champ,imref,imreca); 
    }
        else
    {  
     
     transfo=cr_transf3d(wdth,hght,dpth,NULL);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
         transfo->dx=imres->dx;
         transfo->dy=imres->dy;
         transfo->dz=imres->dz;
         
         /* pour empecher des depalcement superieur a la taille de l'image */
         for (l=0;l<nb_param/3;l++) 
            {
            param[3*l]=MINI(param[3*l],wdth);
            param[3*l+1]=MINI(param[3*l+1],hght);
            param[3*l+2]=MINI(param[3*l+2],dpth);
            param[3*l]=MAXI(param[3*l],-1.0*wdth);
            param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
            param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
            } 
         
     for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*imreca->dx;
            transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
            transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
            } 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 
 }

    
  /*liberation de la base*/
   end_base_3d();

  /*caracteristique du champ*/
    aff_log("CHAMP: ");carac_field_3d(champ);

 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
  /*liberation memoire des images temporaires*/
  //if (normalisation_type>0)
  {
   free_grphic3d(imtref);free_grphic3d(imtreca); 
  }
   free_grphic3d(imtres);

 


 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 if (nomfichres!=NULL)
    {   
    #ifndef TOPO_COMPARE_CRITERE 
     end_flog(LOG_ENERGY);
    #endif
    }
 free(param);free(param_norm);
 
 return(1);  
}


/*******************************************************************************
**        topo_choose_bspline                                            
**                                                                        
**       Choix de la fonction Bspline                                
*******************************************************************************/
scal_func_t topo_choose_bspline(int func_type)
{
scal_func_t scal_func;  

    switch (func_type)
  {
    case 0: scal_func=Bsplined0;aff_log("BSPLINE DEGRE 0 ");break;
    case 1: scal_func=Bsplined1;aff_log("BSPLINE DEGRE 1 ");break;
    case 2: scal_func=Bsplined2;aff_log("BSPLINE DEGRE 2 ");break;
    case 3: scal_func=Bsplineod1;aff_log("BSPLINE O DEGRE 1 ");break;
    default: PUT_WARN("Use by default of BSPLINE DEGRE 1 "); scal_func=Bsplined1;aff_log("BSPLINE DEGRE 1 ");break;
  }
return scal_func;
}


/*******************************************************************************
**        topo_choose_distance                                            
**                                                                        
**       Choix du crit�e de similarite                                 
*******************************************************************************/
dist_func_locale_t topo_choose_distance(int dist_type)
{
dist_func_locale_t distance;
switch (dist_type)
  {
    case 0: distance=Energie_quad_locale_3d;aff_log("L2 ");break;
    case 1: distance=Energie_quad_sym_locale_3d;aff_log("L2 SYM ");break;
    case 2: distance=Energie_Lp_locale_3d;aff_log("Lp ");break;
    case 3: distance=Energie_Lp_sym_locale_3d;aff_log("LP SYM ");break;
        case 4: distance=Energie_L1L2_locale_3d;aff_log("L1L2 ");break;
        case 5: distance=Energie_L1L2_sym_locale_3d;aff_log("L1L2 SYM ");break;
        case 6: distance=Energie_L1norm_locale_3d;aff_log("L1NORM ");break;
        case 7: distance=Energie_quad_topo_locale_3d;aff_log("L2 TOPO ");break;
        case 8: distance=Energie_quad_sym_topo_locale_3d;aff_log("L2 SYM TOPO ");break;
        case 9: distance=Energie_geman_locale_3d;aff_log("GEMAN ");break;
    case 10: distance=Energie_geman_sym_locale_3d;aff_log("GEMAN SYM");break;
    case 11: distance=Energie_quad_locale_pondere_3d;aff_log("L2 PONDERE");break;
    case 12: distance=Energie_IM_locale_3d;aff_log("IM");break;
    case 13: distance=Energie_patch_locale_3d;aff_log("Patch");break;
    case 14: distance=Energie_DTI_quad_locale_3d;aff_log("DTI");break;
    default:PUT_WARN("Use by default of L2 ");distance=Energie_quad_locale_3d;aff_log("QUAD ");break;
  }
return(distance);
}

/*******************************************************************************
**        topo_choose_optimisation                                            
**                                                                        
**       Choix de la m�hode d'optimisation                                 
*******************************************************************************/
min_func_locale_t topo_choose_optimisation(int min_type)
{
min_func_locale_t minimisation;

  switch (min_type)
  {
    case 0: minimisation=TOP_min_icm_locale_3d;aff_log("ICM ");break;
    case 1: minimisation=TOP_min_desc_grad_locale_3d;aff_log("DESCGRAD ");break;
    case 2: minimisation=TOP_min_marquardt_locale_3d;aff_log("MARQUARDT ");break;
    case 3: minimisation=TOP_min_marquardt_penal_locale_3d;aff_log("MARQUARDT PENAL ");break;
    case 4: minimisation=TOP_min_marquardt_locale_3d_parallel;aff_log("MARQUARDT PARALLEL ");break;
    //case 4: minimisation=TOP_min_marquardt_penal_locale_gradient_3d;aff_log("MARQUARDT GRADIENT ");break;
    //case 5: minimisation=TOP_min_marquardt_adapt_penal_locale_3d;aff_log("MARQUARDT ADAPTATIF ");break;
    default:PUT_WARN("Use by default of ICM ");minimisation=TOP_min_icm_locale_3d;aff_log("ICM ");break;  
  }

return(minimisation);
}

/*******************************************************************************
**        topo_choose_normalisation                                            
**                                                                        
**       Choix de la m�hode de normalisation                                 
*******************************************************************************/
norm_func_locale_t topo_choose_normalisation(int normalisation_type, int nb_classe)
{
norm_func_locale_t normalisation;

 switch (normalisation_type)
  {
    case 0: normalisation=NULL;aff_log("SANS NORMALISATION ");break;
    case 1: normalisation=TOP_norm_mean;aff_log("NORMALISATION ECART MOYENNE ");break;
        case 2: normalisation=TOP_norm_meanecty;aff_log("NORMALISATION MOYENNE ECART-TYPE ");break;
    case 3: normalisation=TOP_norm_rl;aff_log("NORMALISATION REGRESSION LINEAIRE ");break;
    case 4: normalisation=TOP_norm_mean_ratio;aff_log("NORMALISATION RAPPORT MOYENNES ");break;
    case 5: normalisation=TOP_norm_minmax;aff_log("NORMALISATION RAPPORT MIN/MAX ");break;
    case 6: normalisation=TOP_norm_GM;aff_log("NORMALISATION HISTOGRAMME MIXTURE %d GAUSSIENNES ",nb_classe);break;
        case 7: normalisation=normalise_kmeans_p;aff_log("NORMALISATION SEGMENTATION KMEANS %d CLASSES ",nb_classe);break;
    case 8: normalisation=normalise_markov_p;aff_log("NORMALISATION SEGMENTATION MARKOV %d CLASSE ",nb_classe);break;
    case 9: normalisation=normalise_histo_joint;aff_log("NORMALISATION HISTOGRAMME JOINT (%d LINEAIRE PAR MORCEAUX)  ",nb_classe);break;
    case 10: normalisation=TOP_norm_spline;aff_log("NORMALISATION HISTOGRAMME JOINT (SPLINE) ",nb_classe);break;
    case 11: normalisation=imx_normalisation_gaussienne2d_histo_joint_standart_p;aff_log("NORMALISATION HISTOGRAMME JOINT (MIXTURE %d GAUSSIENNE) STANDART",nb_classe);break;
    case 12: normalisation=imx_normalisation_gaussienne2d_histo_joint_norm_p;aff_log("NORMALISATION HISTOGRAMME JOINT (MIXTURE %d GAUSSIENNE) NORME",nb_classe);break;
    case 13: normalisation=imx_normalisation_gaussienne2d_histo_joint_reclasse_p;aff_log("NORMALISATION HISTOGRAMME JOINT (MIXTURE %d GAUSSIENNE) RECLASSE",nb_classe);break;
    case 14: normalisation=imx_normalisation_gaussienne2d_histo_joint_multimodal_p;aff_log("NORMALISATION MULTIMODALE HISTOGRAMME JOINT (MIXTURE %d GAUSSIENNE) ",nb_classe);break;
    case 15: normalisation= normalisation_equalisation_histo_p;aff_log("EQUALISATION HISTO");break;
    case 16: normalisation= imx_norm_seuil_meanecty_robust_3d_p;aff_log("NORMALISATION MOYENNE ECART-TYPE ROBUSTE ");break;
        case 17: normalisation= imx_normalisation_segmentation_p;aff_log("NORMALISATION SEGMENTATION SOUS-JACENTE ( %d CLASSES) ",nb_classe);break;
        case 18: normalisation= TOP_HistogramMatchingImageFilter_3d_p;aff_log("NORMALISATION HISTOGRAMME ( %d QUANTILES) ",nb_classe);break;
     default: PUT_WARN("Use by default of NORMALISATION MOYENNE ECART-TYPE "); normalisation=TOP_norm_meanecty;aff_log("NORMALISATION MOYENNE ECART-TYPE ");break;  
  }
    
return(normalisation);
}

/*******************************************************************************
**        topo_choose_distance                                            
**                                                                        
**       Choix du crit�e de similarite                                 
*******************************************************************************/
reg_func_locale_t topo_choose_regularisation(int reg_type)
{
reg_func_locale_t regularisation;
switch (reg_type)
  {
    case 2: regularisation=regularisation_energie_membrane_local;aff_log("Membrane Elastique ");break;
    case 3: regularisation=regularisation_energie_membrane_Lp_local;aff_log("Membrane Elastique Lp ");break;
    case 4: regularisation=regularisation_log_jacobien_local;aff_log("Log(J)^2 ");break;
    case 5: regularisation=regularisation_energie_membrane_jacobien_local;aff_log("Membrane Elastique LogJ ");break;
    case 6: regularisation=regularisation_log_jacobien_centre_local;aff_log("Log(J-Jmean)^2 ");break;
    case 7: regularisation=regularisation_dist_identite_local;aff_log("Distance identite ");break;
    case 8: regularisation=regularisation_patch_imagebased_local;aff_log("image-based patch ");break;
    default:regularisation=NULL;break;
  }
return(regularisation);
}

/*******************************************************************************
**        topo_masque_param                                            
**                                                                        
**       Calcul le masque binaire sur lequel il faut optimiser les param�res                                
*******************************************************************************/
void topo_masque_param (grphic3d *imtref,grphic3d *imtreca,grphic3d *imtres,int ***masque_param, int nb_param)
{
int i,j,k,compt=0,ind,D;
int *x0,*x1,*y0,*y1,*z0,*z1;
int flag,topi,topj,topk;

x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
D=TOP_D(nb_param);

    for(i=0;i<D;i++)
        for(j=0;j<D;j++)
            for(k=0;k<D;k++)
                {
                ind=1.0*TOP_conv_ind(i,j,k,nb_param)/3;
                flag=0;
                topi=x0[ind];
                while ((topi<x1[ind])&&(flag==0))
                    {
                    topj=y0[ind];
                    while ((topj<y1[ind])&&(flag==0))
                        {
                        topk=z0[ind];
                        while ((topk<z1[ind])&&(flag==0))                   
                            {
                            if((imtref->mri[topi][topj][topk]!=0)||(imtres->mri[topi][topj][topk]!=0)/*||(imtreca->mri[topi][topj][topk]!=0)*/) flag=1;                 
                            topk++;
                            }
                        topj++;
                        }
                    topi++;
                    }
                if (flag==0)
                    {
                    masque_param[i][j][k]=0;
                    compt++;
                    }
                else
                    {
                    masque_param[i][j][k]=1;
                    }
                }
        printf("Reduction du nombre de parametres : %.2f %% \n",300.0*compt/nb_param);

}

/*******************************************************************************
**        annule_param_avec_masque_param                                            
**                                                                        
**       Met � zero les parametres en dehors du mask                                
*******************************************************************************/
void annule_param_avec_masque_param (double *param, int nb_param, int ***masque_param)
{
int i,j,k,ind,D;

D=TOP_D(nb_param);

for(i=0;i<D;i++)
for(j=0;j<D;j++)
for(k=0;k<D;k++)
    if(masque_param[i][j][k]==0)
    {
    ind=1.0*TOP_conv_ind(i,j,k,nb_param);
    param[ind]=0.0;
    param[ind+1]=0.0;
    param[ind+2]=0.0;
                
    }       
}
/*******************************************************************************
**        topo_masque_param_primitive                                            
**                                                                        
**       Calcul le masque binaire sur lequel il faut optimiser les param�res                                
*******************************************************************************/
void topo_masque_param_primitive (grphic3d *imtref,int ***masque_param, int nb_param)
{
int i,j,k,compt=0,ind,D;
int *x0,*x1,*y0,*y1,*z0,*z1;
int flag,topi,topj,topk;

x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
D=TOP_D(nb_param);

    for(i=0;i<D;i++)
        for(j=0;j<D;j++)
            for(k=0;k<D;k++)
                {
                ind=1.0*TOP_conv_ind(i,j,k,nb_param)/3;
                flag=0;
                topi=x0[ind];
                while ((topi<x1[ind])&&(flag==0))
                    {
                    topj=y0[ind];
                    while ((topj<y1[ind])&&(flag==0))
                        {
                        topk=z0[ind];
                        while ((topk<z1[ind])&&(flag==0))                   
                            {
                            if ((imtref->mri[topi][topj][topk]!=0)||(imtref->mask->mri[topi][topj][topk]!=0)) flag=1;                   
                            topk++;
                            }
                        topj++;
                        }
                    topi++;
                    }
                if (flag==0)
                    {
                    masque_param[i][j][k]=0;
                    compt++;
                    }
                else
                    {
                    masque_param[i][j][k]=1;
                    }
                }
        printf("Reduction du nombre de parametres : %.2f %% \n",300.0*compt/nb_param);

}


/*******************************************************************************
**        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE                                            
**                                                                        
*******************************************************************************/
void update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE(grphic3d *imtref,grphic3d *imtreca, int nb_param, double *param, double *param_norm, dist_func_locale_t distance, reg_func_locale_t regularisation, int ***masque_param)
{
double tmp,Ereg,Esim;
grphic3d *mask=NULL;
int wdth,hght,dpth;


if ((distance=Energie_ICP_sym_locale_3d)||(distance=Energie_ICP_locale_3d))
    mask=imtref;

wdth=imtref->width;
hght=imtref->height;
dpth=imtref->depth;

tmp=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;

/* calcul du crit�e de similarit�global sans la regularisation */
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
Esim=Energie_globale_3d(imtref,imtreca, nb_param, param,param_norm, distance,regularisation,masque_param);
Esim=Esim*Esim;

Ereg=regularisation_globale_3d(nb_param,param,param_norm, masque_param,mask,regularisation);
Ereg=Ereg*Ereg;

if (Ereg>0)
    //TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=pow(2.0,1.0*BASE3D.resol)*Esim/Ereg;
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=Esim/Ereg;
else
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=1;


//TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=1.0*Esim/(wdth*hght*dpth);


    
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=tmp;
printf("TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE : %f   TOPO_REGULARISATION_MEMBRANE_ELASTIQUE : %f \n",TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE,TOPO_REGULARISATION_MEMBRANE_ELASTIQUE);
}











/*******************************************************************************
**        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2                                            
**                                                                        
*******************************************************************************/
void update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin, min_func_locale_t minimisation,
                                            transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,
                                                                        int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
double tmp,Ereg,Esim;
double Ereg_min,Ereg_max,Esim_min,Esim_max;
double *p,*param_norm;
int wdth,hght,dpth,l;
grphic3d *mask=NULL;

if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d))
    mask=imref;


wdth=imref->width;
hght=imref->height;
dpth=imref->depth;

tmp=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;


p=(double *)malloc(nb_param*sizeof(double)); 
param_norm=(double *)malloc(nb_param*sizeof(double)); 

for (l=0;l<nb_param;l++)
    p[l]=param[l];


for (l=0;l<nb_param/3;l++)
{   param_norm[3*l]  =1.0*param[3*l]  /wdth;
    param_norm[3*l+1]=1.0*param[3*l+1]/hght;
    param_norm[3*l+2]=1.0*param[3*l+2]/dpth;}
    
/* calcul du crit�e de similarit�global sans la regularisation */
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
Esim=Energie_globale_3d(imref,imreca, nb_param, param,param_norm, dist,regularisation,masque_param);
Esim_max=Esim*Esim;

Ereg=regularisation_globale_3d(nb_param,param,param_norm, masque_param,mask,regularisation);
Ereg_min=Ereg*Ereg;



minimisation(imref,imreca,imres,champ_fin,the_transf,inter,dist, regularisation, nb_param,p,masque_param,Jmin,Jmax,nomfichiertrf);


for (l=0;l<nb_param/3;l++)
{   param_norm[3*l]  =1.0*p[3*l]  /wdth;
    param_norm[3*l+1]=1.0*p[3*l+1]/hght;
    param_norm[3*l+2]=1.0*p[3*l+2]/dpth;}

Esim=Energie_globale_3d(imref,imreca, nb_param, p,param_norm, dist,regularisation,masque_param);
Esim_min=Esim*Esim;

Ereg=regularisation_globale_3d(nb_param,p,param_norm, masque_param,mask,regularisation);
Ereg_max=Ereg*Ereg;


printf("Ereg_min : %f \t Ereg_max : %f \n",Ereg_min, Ereg_max);
if (fabs(Ereg_max-Ereg_min)>0)
    {//TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=pow(2.0,1.0*BASE3D.resol)*Esim/Ereg;
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=fabs((Esim_max-Esim_min)/(Ereg_max-Ereg_min));}
else
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=1;


//TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=1.0*Esim/(wdth*hght*dpth);

    
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=tmp;
printf("TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE : %f   TOPO_REGULARISATION_MEMBRANE_ELASTIQUE : %f \n",TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE,TOPO_REGULARISATION_MEMBRANE_ELASTIQUE);


//--- On prend comme initialisation, les valeurs estim�es
for (l=0;l<nb_param;l++)
    param[l]=p[l];




free(p);free(param_norm);
}
/*****************************************/

 /* toto = ptr_img_3d(4);

  imx_copie_3d_p(imtres,toto);
  toto->x3dc = toto->width/2;
  toto->y3dc = toto->height/2;
  toto->z3dc = toto->depth/2;
  imx_inimaxminpixel_3d_p(toto);
  Refresh_3d();
    XGET_OKCANCEL(&i);*/

/*****************************************/


/*******************************************************************************
**                           init_melange_gaussienne_2d                                            
**                                                                        
*******************************************************************************/
void init_melange_gaussienne_2d(int nb_classe)
{
Param_Gaussienne_2d *param;
grphic *histo;


MELANGE_GAUSSIENNE_2D.nb_classe =nb_classe;
param=(Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
MELANGE_GAUSSIENNE_2D.param=param;
histo=cr_grphic_modif(256,256,0,1,0);
MELANGE_GAUSSIENNE_2D.histo=histo;
}


/*******************************************************************************
**                           free_melange_gaussienne_2d                                            
**                                                                        
*******************************************************************************/
void free_melange_gaussienne_2d(void)
{
free(MELANGE_GAUSSIENNE_2D.param);
free_grphic(MELANGE_GAUSSIENNE_2D.histo);
}

/*******************************************************************************
**                           init_HistoJoint                                            
**                                                                        
*******************************************************************************/
void init_HistoJoint(int nbI,int nbJ)
{

HISTOJOINT.nbI=nbI;
HISTOJOINT.nbJ=nbJ;
HISTOJOINT.Ntot=0;
HISTOJOINT.maxI=0.0;
HISTOJOINT.maxJ=0.0;
HISTOJOINT.normI=0.0;
HISTOJOINT.normJ=0.0;
HISTOJOINT.IM=0.0;
HISTOJOINT.aux_IM=0.0;
HISTOJOINT.hij=alloc_dmatrix(nbI,nbJ);
HISTOJOINT.aux_hij=alloc_dmatrix(nbI,nbJ);
HISTOJOINT.hi=malloc(nbI*sizeof(double));
HISTOJOINT.hj=malloc(nbJ*sizeof(double));
HISTOJOINT.aux_hi=malloc(nbI*sizeof(double));
HISTOJOINT.aux_hj=malloc(nbJ*sizeof(double));
}


/*******************************************************************************
**                           update_maxnorm_HistoJoint                                            
**                                                                        
*******************************************************************************/
void update_maxnorm_HistoJoint(grphic3d *imref,grphic3d *imreca)
{

imx_inimaxminpixel_3d_p(imreca);
imx_inimaxminpixel_3d_p(imref);

HISTOJOINT.maxI=imreca->max_pixel;
HISTOJOINT.maxJ=imref->max_pixel;
HISTOJOINT.normI=1.0*(HISTOJOINT.nbI-1)/HISTOJOINT.maxI;
HISTOJOINT.normJ=1.0*(HISTOJOINT.nbJ-1)/HISTOJOINT.maxJ;

}


/*******************************************************************************
**                           free_HistoJoint                                            
**                                                                        
*******************************************************************************/
void free_HistoJoint(void)
{
free_dmatrix(HISTOJOINT.hij,HISTOJOINT.nbI,HISTOJOINT.nbJ);
free_dmatrix(HISTOJOINT.aux_hij,HISTOJOINT.nbI,HISTOJOINT.nbJ);
free(HISTOJOINT.hi);
free(HISTOJOINT.hj);
free(HISTOJOINT.aux_hi);
free(HISTOJOINT.aux_hj);
}


/*******************************************************************************
**                           update_HistoJoint_global                                            
**                                                                        
*******************************************************************************/

void update_HistoJoint_global(grphic3d *imtref,grphic3d *imtres)
{
int i,j,k,wdth,hght,dpth,Ntot,nbI,nbJ;
double x,y;
int x0,y0;
double IM;

wdth=imtref->width;
hght=imtref->height;
dpth=imtref->depth;

nbI=HISTOJOINT.nbI-1;
nbJ=HISTOJOINT.nbJ-1;



for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
    HISTOJOINT.hij[i][j]=0.0;


Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
    if((imtref->mri[i][j][k]>0)&&(imtres->mri[i][j][k]>0))
        {
        Ntot++;
        x=1.0*imtres->mri[i][j][k]/HISTOJOINT.maxI*nbI;
        y=1.0*imtref->mri[i][j][k]/HISTOJOINT.maxJ*nbJ;
        x0=(int)floor(x);
        y0=(int)floor(y);
        if ((x0<HISTOJOINT.nbI)&&(y0<HISTOJOINT.nbJ)) HISTOJOINT.hij[x0][y0]+=(x0-x+1.0)*(y0-y+1.0);
        
        if  (x0<nbI) HISTOJOINT.hij[x0+1][y0]+=(x-x0)*(y0-y+1.0);
        
        if  (y0<nbJ) HISTOJOINT.hij[x0][y0+1]+=(x0-x+1.0)*(y-y0);
        
        if  ((x0<nbI)&&(y0<nbJ)) HISTOJOINT.hij[x0+1][y0+1]+=(x-x0)*(y-y0);
        }


HISTOJOINT.Ntot=Ntot;

for (i=0;i<HISTOJOINT.nbI;i++)
    HISTOJOINT.hi[i]=0.0;


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
    HISTOJOINT.hi[i]+=HISTOJOINT.hij[i][j];



for (j=0;j<HISTOJOINT.nbJ;j++)
    HISTOJOINT.hj[j]=0.0;


for (j=0;j<HISTOJOINT.nbJ;j++)
for (i=0;i<HISTOJOINT.nbI;i++)
    HISTOJOINT.hj[j]+=HISTOJOINT.hij[i][j];
    

// calcul de l'information mutuelle
IM=0;


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
    if (HISTOJOINT.hij[i][j]>0)
    IM-=HISTOJOINT.hij[i][j]*log(Ntot*HISTOJOINT.hij[i][j]/HISTOJOINT.hi[i]/HISTOJOINT.hj[j]);

IM=IM/Ntot;
/*IM=-Ntot*log(Ntot);


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
    if (HISTOJOINT.hij[i][j]>0)
    IM-=HISTOJOINT.hij[i][j]*log(HISTOJOINT.hij[i][j]);

for (i=0;i<HISTOJOINT.nbI;i++)
    if (HISTOJOINT.hi[i]>0)
        IM+=HISTOJOINT.hi[i]*log(HISTOJOINT.hi[i]);

for (j=0;j<HISTOJOINT.nbJ;j++)
    if (HISTOJOINT.hj[j]>0)
    IM+=HISTOJOINT.hj[j]*log(HISTOJOINT.hj[j]);
    */

if (IM!=0)  
    IM=-1/IM;

HISTOJOINT.IM=IM;

}


/*******************************************************************************
**                           init_IM_locale_before_optimization_3d                                            
**                                                                        
*******************************************************************************/

void init_IM_locale_before_optimization_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, int topi, int topj, int topk)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3,Ntot,i0,j0;
double xrel,yrel,zrel,x,y;
double va,vb,vc,vd,ve,vf,vg,vh,Ptot;
field3d *champ;
vector3d ***data;
int nbI,nbJ;

nbI=HISTOJOINT.nbI-1;
nbJ=HISTOJOINT.nbJ-1;


Ptot=0;Ntot=0;
resol=BASE3D.resol;
width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
topD=(int)pow(2.0,1.0*resol)-1;


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;


champ=cr_field3d(topDx,topDy,topDz);
data=champ->raw;
for (i=0;i<topDx;i++)
    for(j=0;j<topDy;j++)
        for (k=0;k<topDz;k++)
            data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;
            
// Mise a jour du champ de deformation induit par les boites voisines

for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
    for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
        for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
        {
        l=TOP_conv_ind(i,j,k,nb_param)/3;

        x00=x0[l];x11=x1[l];y00=y0[l];y11=y1[l];z00=z0[l];z11=z1[l];
        px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
        for (auxi=MAXI(bx0,x00);auxi<MINI(bx1,x11);auxi++)
            for (auxj=MAXI(by0,y00);auxj<MINI(by1,y11);auxj++)
                for (auxk=MAXI(bz0,z00);auxk<MINI(bz1,z11);auxk++) 
                {
                f=fx[auxi-x00]*fy[auxj-y00]*fz[auxk-z00];
                ux=auxi-bx0;uy=auxj-by0;uz=auxk-bz0;
                data[ux][uy][uz].x=data[ux][uy][uz].x+px*f;
                data[ux][uy][uz].y=data[ux][uy][uz].y+py*f;
                data[ux][uy][uz].z=data[ux][uy][uz].z+pz*f;
                }
  
        }

// On balaye chacun des pixels de la boite
for (i=0;i<topDx;i++)
    for (j=0;j<topDy;j++)
        for (k=0;k<topDz;k++)
       if (imref->mri[i+bx0][j+by0][k+bz0]>0)
         {

    dx=(double)i+data[i][j][k].x+bx0;
    dy=(double)j+data[i][j][k].y+by0;
    dz=(double)k+data[i][j][k].z+bz0;
    i2=(int)dx;j2=(int)dy;k2=(int)dz;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
    
        va=(double)imreca->mri[i2][j2][k2];
        vb=(double)imreca->mri[i2][j2][k3];
        vc=(double)imreca->mri[i2][j3][k2];
        vd=(double)imreca->mri[i2][j3][k3];
        ve=(double)imreca->mri[i3][j2][k2];
        vf=(double)imreca->mri[i3][j2][k3];
        vg=(double)imreca->mri[i3][j3][k2];
        vh=(double)imreca->mri[i3][j3][k3];

        auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        
        if (auxE>0)
        {    
        x=1.0*auxE*HISTOJOINT.normI;
        y=1.0*imref->mri[i+bx0][j+by0][k+bz0]*HISTOJOINT.normJ;
        i0=(int)floor(x);
        j0=(int)floor(y);
        if ((i0<HISTOJOINT.nbI)&&(j0<HISTOJOINT.nbJ)) HISTOJOINT.hij[i0][j0]-=(i0-x+1.0)*(j0-y+1.0);
        if  (i0<nbI) HISTOJOINT.hij[i0+1][j0]-=(x-i0)*(j0-y+1.0);
        if  (j0<nbJ) HISTOJOINT.hij[i0][j0+1]-=(i0-x+1.0)*(y-j0);
        if  ((i0<nbI)&&(j0<nbJ)) HISTOJOINT.hij[i0+1][j0+1]-=(x-i0)*(y-j0);
        }
        
            }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
    
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
    if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imreca->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imreca->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imreca->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imreca->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imreca->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imreca->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imreca->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imreca->mri[i3][j3][k3];

        auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;

        if (auxE>0)
        {    
        x=1.0*auxE*HISTOJOINT.normI;
        y=1.0*imref->mri[i+bx0][j+by0][k+bz0]*HISTOJOINT.normJ;
        i0=(int)floor(x);
        j0=(int)floor(y);
        
        if ((i0<HISTOJOINT.nbI)&&(j0<HISTOJOINT.nbJ)) HISTOJOINT.hij[i0][j0]-=(i0-x+1.0)*(j0-y+1.0);
        if  (i0<nbI) HISTOJOINT.hij[i0+1][j0]-=(x-i0)*(j0-y+1.0);
        if  (j0<nbJ) HISTOJOINT.hij[i0][j0+1]-=(i0-x+1.0)*(y-j0);
        if  ((i0<nbI)&&(j0<nbJ)) HISTOJOINT.hij[i0+1][j0+1]-=(x-i0)*(y-j0);
        }
        
            }
       
      }
   }
   

free_field3d(champ);

    
}



/*******************************************************************************
**                           init_IM_locale_after_optimization_3d                                            
**                                                                        
*******************************************************************************/

void init_IM_locale_after_optimization_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, int topi, int topj, int topk)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3,Ntot,i0,j0;
double xrel,yrel,zrel,x,y;
double va,vb,vc,vd,ve,vf,vg,vh,Ptot;
field3d *champ;
vector3d ***data;
int nbI,nbJ;

nbI=HISTOJOINT.nbI-1;
nbJ=HISTOJOINT.nbJ-1;



Ptot=0;Ntot=0;
resol=BASE3D.resol;
width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
topD=(int)pow(2.0,1.0*resol)-1;


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;




champ=cr_field3d(topDx,topDy,topDz);
data=champ->raw;
for (i=0;i<topDx;i++)
    for(j=0;j<topDy;j++)
        for (k=0;k<topDz;k++)
            data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;
            
// Mise a jour du champ de deformation induit par les boites voisines

for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
    for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
        for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
        {
        l=TOP_conv_ind(i,j,k,nb_param)/3;

        x00=x0[l];x11=x1[l];y00=y0[l];y11=y1[l];z00=z0[l];z11=z1[l];
        px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
        for (auxi=MAXI(bx0,x00);auxi<MINI(bx1,x11);auxi++)
            for (auxj=MAXI(by0,y00);auxj<MINI(by1,y11);auxj++)
                for (auxk=MAXI(bz0,z00);auxk<MINI(bz1,z11);auxk++) 
                {
                f=fx[auxi-x00]*fy[auxj-y00]*fz[auxk-z00];
                ux=auxi-bx0;uy=auxj-by0;uz=auxk-bz0;
                data[ux][uy][uz].x=data[ux][uy][uz].x+px*f;
                data[ux][uy][uz].y=data[ux][uy][uz].y+py*f;
                data[ux][uy][uz].z=data[ux][uy][uz].z+pz*f;
                }
  
        }

// On balaye chacun des pixels de la boite
for (i=0;i<topDx;i++)
    for (j=0;j<topDy;j++)
        for (k=0;k<topDz;k++)
       if (imref->mri[i+bx0][j+by0][k+bz0]>0)
         {

    dx=(double)i+data[i][j][k].x+bx0;
    dy=(double)j+data[i][j][k].y+by0;
    dz=(double)k+data[i][j][k].z+bz0;
    i2=(int)dx;j2=(int)dy;k2=(int)dz;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
    
        va=(double)imreca->mri[i2][j2][k2];
        vb=(double)imreca->mri[i2][j2][k3];
        vc=(double)imreca->mri[i2][j3][k2];
        vd=(double)imreca->mri[i2][j3][k3];
        ve=(double)imreca->mri[i3][j2][k2];
        vf=(double)imreca->mri[i3][j2][k3];
        vg=(double)imreca->mri[i3][j3][k2];
        vh=(double)imreca->mri[i3][j3][k3];

        auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
           
                if (auxE>0)
        {    
            x=1.0*auxE*HISTOJOINT.normI;
        y=1.0*imref->mri[i+bx0][j+by0][k+bz0]*HISTOJOINT.normJ;
        i0=(int)floor(x);
        j0=(int)floor(y);
        if ((i0<HISTOJOINT.nbI)&&(j0<HISTOJOINT.nbJ)) HISTOJOINT.hij[i0][j0]+=(i0-x+1.0)*(j0-y+1.0);
        if (i0<nbI) HISTOJOINT.hij[i0+1][j0]+=(x-i0)*(j0-y+1.0);
        if (j0<nbJ) HISTOJOINT.hij[i0][j0+1]+=(i0-x+1.0)*(y-j0);
        if ((i0<nbI)&&(j0<nbJ)) HISTOJOINT.hij[i0+1][j0+1]+=(x-i0)*(y-j0);
        }
        
            }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
    
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
    if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imreca->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imreca->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imreca->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imreca->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imreca->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imreca->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imreca->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imreca->mri[i3][j3][k3];

        auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;

            if (auxE>0)
        {    
    
        x=1.0*auxE*HISTOJOINT.normI;
        y=1.0*imref->mri[i+bx0][j+by0][k+bz0]*HISTOJOINT.normJ;
        i0=(int)floor(x);
        j0=(int)floor(y);
        if ((i0<HISTOJOINT.nbI)&&(j0<HISTOJOINT.nbJ)) HISTOJOINT.hij[i0][j0]+=(i0-x+1.0)*(j0-y+1.0);
        if (i0<nbI) HISTOJOINT.hij[i0+1][j0]+=(x-i0)*(j0-y+1.0);
        if (j0<nbJ) HISTOJOINT.hij[i0][j0+1]+=(i0-x+1.0)*(y-j0);
        if ((i0<nbI)&&(j0<nbJ)) HISTOJOINT.hij[i0+1][j0+1]+=(x-i0)*(y-j0);
            }   
            }
            }
   }

// mise �jour de  HISTOJOINT.hi
for (i=0;i<HISTOJOINT.nbI;i++)
    HISTOJOINT.hi[i]=0.0;


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
    HISTOJOINT.hi[i]+=HISTOJOINT.hij[i][j];


free_field3d(champ);

}


/*******************************************************************************
**                           imx_realloc_mri_pow2_centered_p                                            
**      Modifie et r�lloue l'image pour avoir une taille en puissance de 2                                       
*******************************************************************************/


void imx_realloc_mri_pow2_centered_p (grphic3d *im)
{
int wdth, hght,dpth, nw, nh, nd, i,j,k,ideb,jdeb,kdeb,ifin,jfin,kfin;
TYPEMRI3D  ***mri;

wdth=im->width;
hght=im->height;
dpth=im->depth;

    /* calcul de la taille de l'image modifi� en puissance de 2 */
        nw=(int)pow(2.0,floor(log(wdth)/log(2)+.99999999));
        nh=(int)pow(2.0,floor(log(hght)/log(2)+.99999999));
        nd=(int)pow(2.0,floor(log(dpth)/log(2)+.99999999)); 
    
    /* alloction de la nouvelle image */    
        mri=cr_mri_3d(nw, nh, nd);
    
    if (mri==NULL)
        printf("Probleme d'allocation dans imx_realloc_mri_pow2_centered_p \n");
    
    /* calcul des coordon�s de d�art pour le remplissage afin que l'image soit centr� */
        ideb=(int)((nw-wdth)*0.5);
        jdeb=(int)((nh-hght)*0.5);
        kdeb=(int)((nd-dpth)*0.5);
        ifin=ideb+wdth;
        jfin=jdeb+hght;
        kfin=kdeb+dpth;
        
    
    /* remplissage le l'image allou� */
        for (i=0; i<nw; i++)
        for (j=0; j<nh; j++)
        for (k=0; k<nd; k++)
            if ((i>=ideb)&&(j>=jdeb)&&(k>=kdeb)&&(i<ifin)&&(j<jfin)&&(k<kfin))
            {mri[i][j][k]=im->mri[i-ideb][j-jdeb][k-kdeb];}
            else
            mri[i][j][k]=0; 
    
    
    /* mise �jour de im */
        free_mri_3d(im->mri);
        im->mri=mri;
        im->width=nw;
        im->height=nh;
        im->depth=nd;
        
}


/*******************************************************************************
**                           imx_undo_mri_pow2_centered_p                                            
**      Modifie et r�lloue l'image pour avoir une taille en puissance de 2                                       
*******************************************************************************/


void imx_undo_mri_pow2_centered_p (grphic3d *im, int wdth, int hght, int dpth)
{
int nw, nh, nd, i,j,k,ideb,jdeb,kdeb,ifin,jfin,kfin;
TYPEMRI3D  ***mri;


    /* calcul de la taille de l'image modifi� en puissance de 2 */
        nw=(int)pow(2.0,floor(log(wdth)/log(2)+.99999999));
        nh=(int)pow(2.0,floor(log(hght)/log(2)+.99999999));
        nd=(int)pow(2.0,floor(log(dpth)/log(2)+.99999999)); 
    
    /* alloction de la nouvelle image */    
        mri=cr_mri_3d(wdth, hght, dpth);
    
    /* calcul des coordon�s de d�art pour le remplissage afin que l'image soit centr� */
        ideb=(int)((nw-wdth)*0.5);
        jdeb=(int)((nh-hght)*0.5);
        kdeb=(int)((nd-dpth)*0.5);
        ifin=ideb+wdth;
        jfin=jdeb+hght;
        kfin=kdeb+dpth;
        
    
    /* remplissage le l'image allou� */
        for (i=0; i<wdth; i++)
        for (j=0; j<hght; j++)
        for (k=0; k<dpth; k++)
            mri[i][j][k]=im->mri[i+ideb][j+jdeb][k+kdeb];
            
    
    /* mise �jour de im */
        free_mri_3d(im->mri);
        im->mri=mri;
        im->width=wdth;
        im->height=hght;
        im->depth=dpth;
        
}


/*******************************************************************************
********************************************************************************
*************** RECALAGE BSPLINE  Primitives geometriques **********************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**        matching_Bspline_topo_primitive_3d                                                
**                                                                        
**       recalage Bspline de deux images binaires 3D avec conservation de la topologie                               
*******************************************************************************/

void matching_Bspline_topo_primitive_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[15],*nomfichres;
 int dist_type,inter_type,min_type,save_type,func_type,i, reg_type=0;
 double Jmin,Jmax;
 int resolf, e, err;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref= GET_PLACE3D(TEXT0031);
 im_res= GET_PLACE3D(TEXT0006);

 for(i=0;i<15;i++)
    quest[i]=CALLOC(80,char);
 
 /*type de fonction*/
  func_type=1;
  
     /*distance*/
    dist_type=0;
      
    /*interpolation*/
    inter_type=0;

    /*minimisation*/
    strcpy(quest[0],"ICM");
  strcpy(quest[1],"Descente gradient");
    strcpy(quest[2],"\0");
    
    /*
  strcpy(quest[2],"Levenberg-Marquardt");
  strcpy(quest[3],"Levenberg-Marquardt sans Topo");
    strcpy(quest[4],"\0");
    */
    
  min_type=GETV_QCM("Minimisation",(const char **)quest);

   
    /*question sur l'enregistrement du champ resultat*/
    nomfichres=quest_save_result_3d(&save_type);
    

    /*resolution finale*/
    resolf= GET_INT("resolution", 4, &e);

    /*bornes sur le jacobien*/
    Jmin = 2;
    while (Jmin>=1)  
        Jmin = GET_DOUBLE("Jmin<1",0,&err);

    Jmax = 0;
    while (Jmax<=1)  
        Jmax = GET_DOUBLE("Jmax>1",100000,&err);


/* regularisation */
    strcpy(quest[0],"Pas de regularisation");
  strcpy(quest[1],"Regularisation membrane elastique");
  strcpy(quest[3],"\0");
        
    reg_type= GETV_QCM("Regularisation ?",(const char **)quest);
//  reg_type=0;

    
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
    
    if (reg_type>0)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&err);
    

    for(i=0;i<15;i++)
    free(quest[i]);
        
imx_matching_Bspline_topo_primitive_3d(im_ref,im_reca,im_res,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax);

 
 if (nomfichres)
  free(nomfichres);

  
 show_picture_3d(im_res);
 
}

/*******************************************************************************
**        imx_matching_Bspline_topo_primitive_3d                                                
**                                                                        
**       recalage Bspline de deux images binaires 3D avec conservation de la topologie                               
*******************************************************************************/


int imx_matching_Bspline_topo_primitive_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax)
{
int l,continu;
grphic3d *imref,*imreca,*imres;
imref=ptr_img_3d(im_ref);
imreca=ptr_img_3d(im_reca);
imres=ptr_img_3d(im_res);

 if ((imref->dx!=imreca->dx)||(imref->dy!=imreca->dy)||(imref->dz!=imreca->dz))
    {
      PUT_WARN("Les images n'ont pas les memes dx,dy,dz !");
      printf("Les images n'ont pas les memes dx,dy,dz !\n");
      
            return(-1) ;
    }

 if ((imref->width!=imreca->width)||(imref->height!=imreca->height)||(imref->depth!=imreca->depth))
    {
      PUT_WARN("Les images n'ont pas la meme taille !");
      printf("Les images n'ont pas la meme taille !\n");
      return(-1) ;
    }


continu=0;
for (l=0;l<10;l++)
    {
    if (imref->width==pow(2.0,l)) continu++;
    if (imref->height==pow(2.0,l)) continu++;
    if (imref->depth==pow(2.0,l)) continu++;
    if (imreca->width==pow(2.0,l)) continu++;
    if (imreca->height==pow(2.0,l)) continu++;
    if (imreca->depth==pow(2.0,l)) continu++;
    }

if (continu!=6)
    {
      PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
          printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
            //return(-1) ;
            save_type=3;
   }

 imx_matching_Bspline_topo_primitive_3d_p(imref,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax);
 
  return(1);
}


/*******************************************************************************
**        imx_matching_Bspline_topo_primitive_3d_p                                                
**                                                                        
**       recalage Bspline de deux images binaires 3D avec conservation de la topologie                               
*******************************************************************************/

int imx_matching_Bspline_topo_primitive_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb;
 double *param,*param_norm;
 grphic3d *imtref,*imtreca,*imtres,*imtmp;
 field3d *champ;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation=NULL;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;
#ifndef TOPO_COMPARE_CRITERE 
    char nomfichier_nrj[255];
    char *ext;
#endif
#ifndef WIN32
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif


    wdth=imref->width;hght=imref->height;dpth=imref->depth;
    
    if (save_type==3) /* il faut transformer les images en puissance de 2 */
        {
        wdth_old=wdth; hght_old=hght; dpth_old=dpth;
        
        imx_realloc_mri_pow2_centered_p(imref);
        imx_realloc_mri_pow2_centered_p(imreca);
        imx_realloc_mri_pow2_centered_p(imres);
    
        wdth=imref->width;hght=imref->height;dpth=imref->depth;
        
        }
    
    tdebut=time(NULL);
    nbmax_param=3*pow((pow(2,resolf)-1),3);
 
   /* Allocation */
  if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
    if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
 init_log(NULL);
 aff_log("\nRECALAGE ");
 
 /*************mise en place des parametres du recalage******/
  scal_func=topo_choose_bspline(func_type);
    //distance=topo_choose_distance(dist_type);
    distance=Energie_ICP_sym_locale_3d;
    //distance=Energie_ICP_locale_3d;
    
    interpol=imx_choose_interpolation_fct(inter_type);  
  minimisation=topo_choose_optimisation(min_type);
    
    if (reg_type>0)
        regularisation=regularisation_energie_membrane_local;
    
    aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
  /*creation des images temporaires*/
    imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
    imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
        imtres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtres);
        imtreca->mask=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca->mask);
        imtref->mask=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref->mask); 
        imtres->mask=cr_grphic3d(imres);imx_copie_3d_p(imres,imtres->mask); 
        imtmp=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtmp);
        

  /* calcul du mask de imref avec l'image compl�entaire */
        imx_inv_3d_p(imtref,imtref->mask);
        imx_inv_3d_p(imreca,imtmp);
    
    
    /* calcul des cartes de distance */
    imx_chamfer_distance_3d_p(imreca,imtreca);
    imx_chamfer_distance_3d_p(imtmp,imtreca->mask);

    imx_mul_coe_3d_p(imtreca,1000.0,imtreca);
    imx_mul_coe_3d_p(imtreca->mask,1000.0,imtreca->mask);


    /*allocation memoire de champ*/
   champ=cr_field3d(wdth,hght,dpth);
     

    
 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r,D,bande;
  int ***masque_param;
  int *x0,*x1,*y0,*y1,*z0,*z1;
  double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;
        
  /* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}

    if (nomfichres!=NULL)
        {
        #ifndef TOPO_COMPARE_CRITERE 
        strcpy(nomfichier_nrj,nomfichres);
    ext=strstr(nomfichier_nrj,".trf");
            if (ext!=NULL) *ext='\0';
        strcat(nomfichier_nrj,"_energy.txt");
        LOG_ENERGY=init_flog(nomfichier_nrj);      
        #endif
     }
     
      
    /* Debut traitement */
    for (r=resol;r<=resolf;r++)
  {

    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    base_to_field_3d(nb_param,param,champ,NULL,NULL);
    interpol(imreca,champ,imtres); 
    
    bande=(int)((x1[0]-x0[0])/4.0);
    
    if (bande<15)
    {
    imx_dilat_3d_p(imtref,imtref->mask,1,bande);
    imx_sub_3d_p(imtref->mask,imtref,imtref->mask);
    }
        
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    
    
    /*  Remplissage de masque_param */
     topo_masque_param_primitive (imtref, masque_param, nb_param);

    //annule_param_avec_masque_param (param, nb_param,masque_param);

    if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
    else
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
        
/* A modifier et a adpater !!!*/
    if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(imtref,imtreca,imtres,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);



 mini=minimisation(imtref,imtreca,imtres,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);

   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
    free_imatrix_3d(masque_param);
    }       

    
    /*calcul de la transformation*/
    base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
  if (imres != NULL)
  { 
        interpol=inter_qsinc3_3d;
    interpol(imreca,champ,imtres);
    imx_copie_3d_p(imtres,imres);
  }

    if (save_type==3) /* il faut transformer les images  dans leur taille d'origine*/
        {
    
        imx_undo_mri_pow2_centered_p(imref,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(imreca,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(imres,wdth_old,hght_old,dpth_old);
    
        }
  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 
    transf3d *transfo; 
    
    if ((save_type==0)||(save_type==3))
        {
            if (save_type==3) /* on tronque le champ */
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                } 
        
        
         transfo=field_to_transf_3d(champ,imref,imreca); 
    }
        else
    {  
     
     transfo=cr_transf3d(wdth,hght,dpth,NULL);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
         transfo->dx=imres->dx;
         transfo->dy=imres->dy;
         transfo->dz=imres->dz;
         
         /* pour empecher des depalcement superieur a la taille de l'image */
         for (l=0;l<nb_param/3;l++) 
            {
            param[3*l]=MINI(param[3*l],wdth);
            param[3*l+1]=MINI(param[3*l+1],hght);
            param[3*l+2]=MINI(param[3*l+2],dpth);
            param[3*l]=MAXI(param[3*l],-1.0*wdth);
            param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
            param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
            } 
         
     for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*imreca->dx;
            transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
            transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
            } 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 
 }

  /*liberation de la base*/
   end_base_3d();

  /*caracteristique du champ*/
    aff_log("CHAMP: ");carac_field_3d(champ);

 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
  /*liberation memoire des images temporaires*/

   free_grphic3d(imtmp);
     free_grphic3d(imtref);
     free_grphic3d(imtreca); 
     free_grphic3d(imtres);

 


 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 if (nomfichres!=NULL)
    {   
    #ifndef TOPO_COMPARE_CRITERE 
     end_flog(LOG_ENERGY);
    #endif
    }
 free(param);free(param_norm); 
 return(1);  
}



/*******************************************************************************
**        simul_atrophie_topo_3d                                                
**                                                                        
**      Simulation d'atrophie : generation d'un champ de deformation preservant
**      la topologie a partir d'une carte de jacobien                              
*******************************************************************************/

void simul_atrophie_topo_3d()
{
 int im_jac,im_invariant,i;
 char *nomfichres,s[250];
 int inter_type,min_type,save_type,func_type,dist_type,reg_type,adaptatif;
 int resolf, e,ans;
 double Jmin, Jmax;
 char *quest[15];
 /*question sur les images a recaler*/
 im_jac=GET_PLACE3D("Carte de jacobien");
 
 sprintf(s,"Voulez vous consider un mask de points invariants\n");
 ans=GET_EXIST(s,&e);
 
 /*question sur les images a recaler*/
 if (ans==1)
 im_invariant=GET_PLACE3D("Mask des points invariants");
 else im_invariant=-100;
 
 
 
 for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);

 /*type de fonction*/
  func_type=1;
    
    /*interpolation*/
    inter_type=4;

    /*minimisation*/
    min_type=1;
    
    /* distance */
    strcpy(quest[0],"L2");
    strcpy(quest[1],"Log");
    strcpy(quest[2],"\0");
    
    
    dist_type=GETV_QCM("Critere",(const char **)quest);

    reg_type = 2;
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    /*bornes sur le jacobien*/
    Jmin = 2;
    while (Jmin>=1)  
        Jmin = GET_DOUBLE("Jmin<1",0,&e);

    Jmax = 0;
    while (Jmax<=1)  
        Jmax = GET_DOUBLE("Jmax>1",100000,&e);

    
    if (reg_type>1)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&e);
    
    /*question sur l'enregistrement du champ resultat*/
    nomfichres=quest_save_result_3d(&save_type);
    

    /*resolution finale*/
    resolf= GET_INT("resolution", 6, &e);

    /* adaptatif ? */
     adaptatif =1;
 
 
    /* distance */
    strcpy(quest[0],"precis");
    strcpy(quest[1],"tres precis");
    strcpy(quest[2],"\0");
    
    
    PRECISION_SIMULATION=GETV_QCM("Precision",(const char **)quest);
    PRECISION_SIMULATION+=1;
    
    for(i=0;i<3;i++)
        free(quest[i]);
    
imx_simul_atrophie_topo_3d(im_jac,im_invariant,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,adaptatif,Jmin, Jmax);

 
 if (nomfichres)
  free(nomfichres);
 
}

/*******************************************************************************
**        imx_simul_atrophie_topo_3d                                                
**                                                                        
**      Simulation d'atrophie : generation d'un champ de deformation preservant
**      la topologie a partir d'une carte de jacobien                              
*******************************************************************************/

int imx_simul_atrophie_topo_3d(int im_ref, int im_mask, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf,int adaptatif,double Jmin, double Jmax)
{
int l,continu;
grphic3d *imref, *mask;

imref=ptr_img_3d(im_ref);

if (im_mask==-100)
mask=NULL;
else
mask=ptr_img_3d(im_mask);

if ((imref->dx==0.0)||(imref->dy==0)||(imref->dz==0))
     {PUT_WARN("Les dx,dy et dz sont nuls");
          printf("Les dx,dy et dz sont nuls\n");
            return(-1) ;
    }




if (mask!=NULL)
 if ((imref->dx!=mask->dx)||(imref->dy!=mask->dy)||(imref->dz!=mask->dz))
    {
      PUT_WARN("La carte de jacobien et le mask de points invariants n'ont pas les memes dx,dy,dz !");
      printf("La carte de jacobien et le mask de points invariants n'ont pas les memes dx,dy,dz !\n");
      
            return(-1) ;
    }

if (mask!=NULL)
 if ((imref->width!=mask->width)||(imref->height!=mask->height)||(imref->depth!=mask->depth))
    {
      PUT_WARN("La carte de jacobien et le mask de points invariants n'ont pas la meme taille !");
      printf("La carte de jacobien et le mask de points invariants n'ont pas la meme taille !\n");
      return(-1) ;
    }




continu=0;
for (l=0;l<10;l++)
    {
    if (imref->width==pow(2.0,l)) continu++;
    if (imref->height==pow(2.0,l)) continu++;
    if (imref->depth==pow(2.0,l)) continu++;
}

if (continu!=3)
    {
      PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
          printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
            //return(-1) ;
            save_type=3;
   }

 imx_simul_atrophie_topo_3d_p(imref,mask,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, adaptatif, Jmin, Jmax);
 
  return(1);
}

/*******************************************************************************
**        imx_simul_atrophie_topo_3d_p                                                
**                                                                        
**      Simulation d'atrophie : generation d'un champ de deformation preservant
**      la topologie a partir d'une carte de jacobien                              
*******************************************************************************/
int imx_simul_atrophie_topo_3d_p(grphic3d *imref,grphic3d *mask, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int resolf,int adaptatif, double Jmin, double Jmax)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb;
 double *param,*param_norm;
 field3d *champ;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation=NULL;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;

#ifndef WIN32
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif


    wdth=imref->width;hght=imref->height;dpth=imref->depth;
    
    if (save_type==3) /* il faut transformer les images en puissance de 2 */
        {
        wdth_old=wdth; hght_old=hght; dpth_old=dpth;
        
        imx_realloc_mri_pow2_centered_p(imref);
        
        if (mask!=NULL)
            imx_realloc_mri_pow2_centered_p(mask);
        
    
        wdth=imref->width;hght=imref->height;dpth=imref->depth;
        }
    
    tdebut=time(NULL);
    nbmax_param=3*pow((pow(2,resolf)-1),3);
 
   /* Allocation */
  if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
    if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
 init_log(NULL);
 
 /*************mise en place des parametres du recalage******/
  scal_func=topo_choose_bspline(func_type);
    
    distance=Energie_atrophie_log_jacobien_locale_3d;
    if (dist_type==0) distance=Energie_atrophie_jacobien_locale_3d;
    if (dist_type==1) distance=Energie_atrophie_log_jacobien_locale_3d;
        
    
    interpol=imx_choose_interpolation_fct(inter_type);  
  minimisation=topo_choose_optimisation(min_type);
    
    if (reg_type>0)
        regularisation=regularisation_energie_membrane_local;
    
    aff_log("\n");
 

    /*allocation memoire de champ*/
   champ=cr_field3d(wdth,hght,dpth);
     

    
 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r,D,ind,flag,topi,topj,topk;
  int ***masque_param;
  int *x0,*x1,*y0,*y1,*z0,*z1;
  double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;
        
  /* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}


      
    /* Debut traitement */
    for (r=resol;r<=resolf;r++)
  {

    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    
        
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    
    /*  Remplissage de masque_param */
    // topo_masque_param_primitive (imtref, masque_param, nb_param);

    for (i=0;i<D;i++)
    for (j=0;j<D;j++)
    for (k=0;k<D;k++)
            masque_param[i][j][k]=1;
            
        if (mask!=NULL)
        for(i=0;i<D;i++)
        for(j=0;j<D;j++)
            for(k=0;k<D;k++)
                {
                ind=1.0*TOP_conv_ind(i,j,k,nb_param)/3;
                flag=0;
                topi=x0[ind];
                while ((topi<x1[ind])&&(flag==0))
                    {
                    topj=y0[ind];
                    while ((topj<y1[ind])&&(flag==0))
                        {
                        topk=z0[ind];
                        while ((topk<z1[ind])&&(flag==0))                   
                            {
                            if((mask->mri[topi][topj][topk]!=0)) flag=1;                    
                            topk++;
                            }
                        topj++;
                        }
                    topi++;
                    }
                if (flag==0)
                    {
                    masque_param[i][j][k]=1;
                    }
                else
                    {
                    masque_param[i][j][k]=0;
                    }
                }



    //annule_param_avec_masque_param (param, nb_param,masque_param);

        
    if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
    else
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
        
    if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(imref,imref,imref,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);
    
        mini=minimisation(imref,imref,imref,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);
    
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
    free_imatrix_3d(masque_param);
    }       

    
    /*calcul de la transformation*/
    base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
 
    if (save_type==3) /* il faut transformer les images  dans leur taille d'origine*/
        {
    
        imx_undo_mri_pow2_centered_p(imref,wdth_old,hght_old,dpth_old);
        
        if (mask!=NULL)
        imx_undo_mri_pow2_centered_p(mask,wdth_old,hght_old,dpth_old);
        
        }
        


  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 
    transf3d *transfo; 
    
    if ((save_type==0)||(save_type==3))
        {
            if (save_type==3) /* on tronque le champ */
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                } 
        
        
         transfo=field_to_transf_3d(champ,imref,imref); 
    }
        else
    {  
     
     transfo=cr_transf3d(wdth,hght,dpth,NULL);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
         transfo->dx=imref->dx;
         transfo->dy=imref->dy;
         transfo->dz=imref->dz;
         
         /* pour empecher des depalcement superieur a la taille de l'image */
         for (l=0;l<nb_param/3;l++) 
            {
            param[3*l]=MINI(param[3*l],wdth);
            param[3*l+1]=MINI(param[3*l+1],hght);
            param[3*l+2]=MINI(param[3*l+2],dpth);
            param[3*l]=MAXI(param[3*l],-1.0*wdth);
            param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
            param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
            } 
         
     for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*imref->dx;
            transfo->param[3*l+1]=param[3*l+1]*imref->dy;
            transfo->param[3*l+2]=param[3*l+2]*imref->dz;
            } 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 
 }

  /*liberation de la base*/
   end_base_3d();

  /*caracteristique du champ*/
    aff_log("CHAMP: ");carac_field_3d(champ);

 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
  /*liberation memoire des images temporaires*/

 


 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 if (nomfichres!=NULL)
    {   
    #ifndef TOPO_COMPARE_CRITERE 
     end_flog(LOG_ENERGY);
    #endif
    }
 free(param);free(param_norm); 
 return(1);  
}



/*******************************************************************************
********************************************************************************
************************* RECALAGE BSPLINE Symetrique **************************
********************************************************************************
*******************************************************************************/



/*******************************************************************************
**        matching_Bspline_topo_norm_symetrique_3d                                                
**                                                                        
**       recalage Bspline de deux images 3D avec conservation de la topologie    
**                  Avec mise en correspondance des intensit� et sym�rique                           
*******************************************************************************/

void matching_Bspline_topo_norm_symetrique_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[20],*nomfichres;
 int dist_type,inter_type,min_type,save_type,func_type,i,normalisation_type;
 int nb_classe,adaptatif,reg_type;
 double Jmin,Jmax;
 int resolf, e, err;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref= GET_PLACE3D(TEXT0031);
 im_res= GET_PLACE3D(TEXT0006);

 for(i=0;i<20;i++)
    quest[i]=CALLOC(80,char);
 
 /*type de fonction*/
  func_type=1;
  
 /*distance*/
  strcpy(quest[0],"L2");
  strcpy(quest[1],"L2 (parties manquantes)");
  
 /* strcpy(quest[1],"L2 Sym");
  strcpy(quest[2],"Lp");
  strcpy(quest[3],"Lp Sym");
    strcpy(quest[4],"L1L2");
    strcpy(quest[5],"L1L2 Sym");
    strcpy(quest[6],"L1norm");
    strcpy(quest[7],"L2 Topo");
    strcpy(quest[8],"L2 Sym Topo");
    strcpy(quest[9],"Geman McLure");
    strcpy(quest[10],"Geman McLure Sym");
    strcpy(quest[11],"L2 pondere");
    strcpy(quest[12],"IM");*/
    strcpy(quest[2],"\0");

  dist_type=GETV_QCM("Distance",(const char **)quest);
  
    /*interpolation*/
    inter_type=4;

    /*minimisation*/
  strcpy(quest[0],"ICM");
  strcpy(quest[1],"Descente gradient");
    strcpy(quest[2],"Levenberg-Marquardt");
  strcpy(quest[3],"Levenberg-Marquardt sans Topo");
  strcpy(quest[4],"Levenberg-Marquardt Parallel");
    strcpy(quest[5],"\0");
    
  min_type=GETV_QCM("Minimisation",(const char **)quest);

   
    /*question sur l'enregistrement du champ resultat*/
    nomfichres=quest_save_result_3d(&save_type);
    

    /*resolution finale*/
    resolf= GET_INT("resolution", 4, &e);

    /*bornes sur le jacobien*/
    Jmin = 2;
    while (Jmin>=1)  
        Jmin = GET_DOUBLE("Jmin<1",0,&err);

    Jmax = 0;
    while (Jmax<=1)  
        Jmax = GET_DOUBLE("Jmax>1",100000,&err);

  /* normalisation d'intensite */
  strcpy(quest[0],"Pas de normalisation");
  strcpy(quest[1],"ecart des moyennes");
  strcpy(quest[2],"Moyenne et ecart-type");
  strcpy(quest[3],"Regression lineaire a 2 parametres");
    strcpy(quest[4],"Rapport des moyennes");
    strcpy(quest[5],"Rapport Min/Max");
    strcpy(quest[6],"Mixture gaussiennes histogramme");
    strcpy(quest[7],"segmentation kmeans");
    strcpy(quest[8],"segmentation markov");
    strcpy(quest[9],"histogramme joint (lineaire par morceau)");
    strcpy(quest[10],"histogramme joint (spline)");
    strcpy(quest[11],"histogramme joint (mixture de gaussienne) standart");
    strcpy(quest[12],"histogramme joint (mixture de gaussienne) norme");
    strcpy(quest[13],"histogramme joint (mixture de gaussienne) reclasse");
    strcpy(quest[14],"Multimodal histogramme joint (mixture de gaussienne)");
    strcpy(quest[15],"EQUALISATION HISTO");
    strcpy(quest[16],"Moyenne et ecart-type robuste");
    strcpy(quest[17],"normalisation segmentation sous-jacente");
    #ifdef  HAS_ITK
    strcpy(quest[18],"normalisation histogramme par quantile");
    strcpy(quest[19],"\0");
    #else
    strcpy(quest[18],"\0");
    #endif
    
 normalisation_type=GETV_QCM("Normalisation",(const char **)quest);
 nb_classe=-1;
if ((normalisation_type>5)&&(normalisation_type!=10)&&(normalisation_type!=16))
    {
    if (normalisation_type==18)
    nb_classe = GET_INT("Nombre de quantile",7,&err);
    else
        if (normalisation_type==15)
            nb_classe = GET_INT("Nombre de bins",1024,&err);
        else
        nb_classe = GET_INT("Nombre de classe",10,&err);
    }
  
    
  
    /* regularisation */
    strcpy(quest[0],"Pas de regularisation");
  strcpy(quest[1],"Regularisation membrane elastique");
  strcpy(quest[2],"Regularisation distance identite");
  strcpy(quest[3],"\0");
        
        
    reg_type= GETV_QCM("Regularisation ?",(const char **)quest);

    
    TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
    if (reg_type==2)
        {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&err);
        reg_type=7; // c'est pour �re compatible avec la fonction topo_choose_regularisation
        }
    
    
    if (reg_type==1)
        {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&err);
        reg_type=2; // c'est pour �re compatible avec la fonction topo_choose_regularisation
        }
    
    
    
    /* dans le cas du recalage symetrique, on accorde une ponderation identique aux deux images */
    TOPO_PONDERATION_RECALAGE_SYMETRIQUE = -1;
    while ((TOPO_PONDERATION_RECALAGE_SYMETRIQUE<=0.0)||(TOPO_PONDERATION_RECALAGE_SYMETRIQUE>1))  
        TOPO_PONDERATION_RECALAGE_SYMETRIQUE = GET_DOUBLE("Ponderation de image de reference = Nref/(Nref+Nreca) (entre 0 et 1)",0.5,&err);
    
    
    adaptatif=1;
    
    for(i=0;i<20;i++)
    free(quest[i]);
        
    

    imx_matching_Bspline_topo_symetrique_3d(im_ref,im_reca,im_res,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax,normalisation_type,nb_classe,adaptatif);   
    
 
 if (nomfichres)
  free(nomfichres);

  
 show_picture_3d(im_res);
 
}


/*******************************************************************************
**        imx_matching_Bspline_topo_symetrique_3d                                             
**                                                                       
**       recalage Bspline de 2 images 3D avec conservation de la topologie                                  
*******************************************************************************/
/*
  im_ref : numero de l'image sur laquelle on recale imreca
  im_reca : numero de l'image a recaler
  im_res : numero de l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  func_type : numero de la fonction d'interpolation du champ (1 pour Bspline
    degre 1)
  dist_type : numero de la fonction de distance (0 pour distance quadratique)
  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
*/
int imx_matching_Bspline_topo_symetrique_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif)
{
int l,continu;
grphic3d *imref,*imreca,*imres;
imref=ptr_img_3d(im_ref);
imreca=ptr_img_3d(im_reca);
imres=ptr_img_3d(im_res);

 if ((imref->dx!=imreca->dx)||(imref->dy!=imreca->dy)||(imref->dz!=imreca->dz))
    {
      PUT_WARN("Les images n'ont pas les memes dx,dy,dz !");
      printf("Les images n'ont pas les memes dx,dy,dz !\n");
      
            return(-1) ;
    }

 if ((imref->width!=imreca->width)||(imref->height!=imreca->height)||(imref->depth!=imreca->depth))
    {
      PUT_WARN("Les images n'ont pas la meme taille !");
      printf("Les images n'ont pas la meme taille !\n");
      return(-1) ;
    }


continu=0;
for (l=0;l<10;l++)
    {
    if (imref->width==pow(2.0,l)) continu++;
    if (imref->height==pow(2.0,l)) continu++;
    if (imref->depth==pow(2.0,l)) continu++;
    if (imreca->width==pow(2.0,l)) continu++;
    if (imreca->height==pow(2.0,l)) continu++;
    if (imreca->depth==pow(2.0,l)) continu++;
    }

if (continu!=6)
    {
      PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
          printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
            //return(-1) ;
            save_type=3;
   }

 imx_matching_Bspline_topo_symetrique_3d_p(imref,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax,normalisation_type, nb_classe,adaptatif);
 
  return(1);
}


/*******************************************************************************
**        imx_matching_Bspline_topo_symetrique_3d_p                                            
**                                                                        
**       recalage Bspline de 2 images 3D avec conservation de la topologie                                   
*******************************************************************************/
/*
  imref : pointeur sur l'image sur laquelle on recale imreca
  imreca : pointeur sur l'image a recaler
  imres : pointeur sur l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  func_type : numero de la fonction d'interpolation du champ (1 pour Bspline
    degre 1)
  dist_type : numero de la fonction de distance (0 pour distance quadratique)
  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
*/
int imx_matching_Bspline_topo_symetrique_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb;
 double *param,*param_norm;
 grphic3d *imtref,*imtreca,*imtres, *imrefd, *imrecad;
 field3d *champ;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 norm_func_locale_t normalisation;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;
 char nomfichier_ref[1000];
 char nomfichier_reca[1000];
 char *ext2;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
 transf3d *transfo,*transfo_opp,*transfo_inv, *transfo_comb, *transfo_tmp; 
   int l;


printf("TOPO_PONDERATION_RECALAGE_SYMETRIQUE = %f \n",TOPO_PONDERATION_RECALAGE_SYMETRIQUE);
printf("lambda_ref = %f \n",lambda_ref);

#ifndef TOPO_COMPARE_CRITERE 
    char nomfichier_nrj[255];
    char *ext;
#endif
#ifndef WIN32
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif


    wdth=imref->width;hght=imref->height;dpth=imref->depth;
    
    if (save_type==3) /* il faut transformer les images en puissance de 2 */
        {
        wdth_old=wdth; hght_old=hght; dpth_old=dpth;
        
        imx_realloc_mri_pow2_centered_p(imref);
        imx_realloc_mri_pow2_centered_p(imreca);
        imx_realloc_mri_pow2_centered_p(imres);
    
        wdth=imref->width;hght=imref->height;dpth=imref->depth;
        
        }
    
    tdebut=time(NULL);
    nbmax_param=3*pow((pow(2,resolf)-1),3);
 
   /* Allocation */
  if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
    if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
 init_log(NULL);
 aff_log("\nRECALAGE ");
 
 /*************mise en place des parametres du recalage******/
  scal_func=topo_choose_bspline(func_type);
    
    
    
    distance=Energie_quad_locale_symetrique_3d;
    if (dist_type==0) distance=Energie_quad_locale_symetrique_3d;
    if (dist_type==1) distance=Energie_quad_locale_symetrique_coupe_3d;
    
    interpol=imx_choose_interpolation_fct(inter_type);  
    minimisation=topo_choose_optimisation(min_type);
    normalisation=topo_choose_normalisation(normalisation_type,nb_classe);
    regularisation=topo_choose_regularisation(reg_type);
    
    aff_log("\n");
  aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
  /*creation des images temporaires*/
    {
    imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
    imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
    
    imrefd=cr_grphic3d(imref);imx_copie_3d_p(imref,imrefd);
    imrecad=cr_grphic3d(imref);imx_copie_3d_p(imref,imrecad);
        
        
        imx_norm_rcoeff_3d_p(imtref, imtreca);

        if ((imreca->mask!=NULL)&&(imref->mask!=NULL))
            {   imtreca->mask=cr_grphic3d(imreca->mask);imx_copie_3d_p(imreca->mask,imtreca->mask);
                imtref->mask=cr_grphic3d(imref->mask);imx_copie_3d_p(imref->mask,imtref->mask); }
        
        }

  imtres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtres);

    if (imreca->mask!=NULL)
        {imtres->mask=cr_grphic3d(imreca->mask);imx_copie_3d_p(imreca->mask,imtres->mask);}
    

    /*allocation memoire de champ*/
   champ=cr_field3d(wdth,hght,dpth);
     

if (dist_type==12)
 init_HistoJoint(TOPO_SIZE_HISTOJOINT,TOPO_SIZE_HISTOJOINT);
    
 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r,D;
  int ***masque_param;
  int *x0,*x1,*y0,*y1,*z0,*z1;
  double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
  
  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;
        
  /* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}

    
     
    /* nom du fichier dans lequel on enregistre la transfo oppose */
    if (nomfichres!=NULL)
        {
        
        strcpy(nomfichier_ref,nomfichres);
        ext2=strstr(nomfichier_ref,".trf");
        if (ext2!=NULL) *ext2='\0';

        strcat(nomfichier_ref,"_ref.trf");


        strcpy(nomfichier_reca,nomfichres);
        ext2=strstr(nomfichier_reca,".trf");
        if (ext2!=NULL) *ext2='\0';
        strcat(nomfichier_reca,"_reca.trf");

     }
     
     
     
     
     if (nomfichres!=NULL)
        {
        #ifndef TOPO_COMPARE_CRITERE 
        strcpy(nomfichier_nrj,nomfichres);
    ext=strstr(nomfichier_nrj,".trf");
            if (ext!=NULL) *ext='\0';
        strcat(nomfichier_nrj,"_energy.txt");
        LOG_ENERGY=init_flog(nomfichier_nrj);      
        #endif
     }
     
      
    /* Debut traitement */
    for (r=resol;r<=resolf;r++)
  {



    TOPO_SEUIL_SIGNIFICATIVITE=0.05;


    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    base_to_field_3d(nb_param,param,champ,NULL,NULL);
    
    /*for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
            {
            champ->raw[i][j][k].x=champ->raw[i][j][k].x;
            champ->raw[i][j][k].y=champ->raw[i][j][k].y;
            champ->raw[i][j][k].z=champ->raw[i][j][k].z;
            }*/
    
    interpol(imreca,champ,imrecad); 
    
    for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
            {
            champ->raw[i][j][k].x=-1.0*lambda_ref*champ->raw[i][j][k].x;
            champ->raw[i][j][k].y=-1.0*lambda_ref*champ->raw[i][j][k].y;
            champ->raw[i][j][k].z=-1.0*lambda_ref*champ->raw[i][j][k].z;
            }
    
    interpol(imref,champ,imrefd); 
    
    
    //if (imreca->mask!=NULL)
        //inter_nearest_3d(imreca->mask,champ,imtres->mask);

    if (imtref->mask==NULL)
            {
            imtref->mask=cr_grphic3d(imtref);imx_copie_3d_p(imtref,imtref->mask);
            
            for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
            for (k=0;k<dpth;k++)
                    imtref->mask->mri[i][j][k]=1;
            }


    if (imtres->mask==NULL)
            {
            imtres->mask=cr_grphic3d(imtres);imx_copie_3d_p(imtres,imtres->mask);
            
            for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
            for (k=0;k<dpth;k++)
                    imtres->mask->mri[i][j][k]=1;
            }
            
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    topo_masque_param (imrefd,NULL,imrecad,masque_param, nb_param);
    
    //annule_param_avec_masque_param (param, nb_param,masque_param);

    
    
    /*normalisation de imreca par rapport a imref: resultat dans imtreca*/ 
    if (normalisation_type>0)
    {
        if (normalisation_type==12)
            {
            
            
            {
               
               
                
           
               transf3d *transfo,*transfo_opp,*transfo_inv, *transfo_comb; 
                int l;
                double *param_tmp;
                
                /*liberation de la base*/
                end_base_3d();

                // transfo direct
                transfo=cr_transf3d(wdth,hght,dpth,NULL);
                transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
                transfo->typetrans=BSPLINE3D;
                transfo->resol=BASE3D.resol;transfo->degre=func_type;
                transfo->dx=imtres->dx;
                transfo->dy=imtres->dy;
                transfo->dz=imtres->dz;

     
                for (l=0;l<nb_param/3;l++) 
                    {
                    transfo->param[3*l]=param[3*l]*imreca->dx;
                    transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
                    transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
                    }
    
    
                // transfo opposee
                transfo_opp=cr_transf3d(wdth,hght,dpth,NULL);
                transfo_opp->nb_param=nb_param;transfo_opp->param=CALLOC(nb_param,double);
                transfo_opp->typetrans=BSPLINE3D;
                transfo_opp->resol=BASE3D.resol;transfo_opp->degre=func_type;
                transfo_opp->dx=imtres->dx;
                transfo_opp->dy=imtres->dy;
                transfo_opp->dz=imtres->dz;
        
                for (l=0;l<nb_param;l++)
                transfo_opp->param[l]=-1.0*lambda_ref*transfo->param[l];
             

                // transfo inverse
                transfo_inv=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);
                inverse_transf_3d_p(transfo_opp,transfo_inv,0.01);
                
                // transfo composee
                transfo_comb=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);
                comb_transf_3d(transfo,transfo_inv,transfo_comb,5);
         
                transf_to_field_3d_noalloc(transfo_comb, champ, imref, imreca);
 
                interpol(imtreca,champ,imrecad); 

        
                param_tmp=CALLOC(nb_param,double);
                
             nb_param=init_base_3d(wdth,hght,dpth,scal_func);
            
            for(l=1; l<r;l++)
                {/*on passe a la resolution superieure*/
                nb_param=base_resol_up_3d(param_tmp,nb_param);
                }
                
                x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

                free(param_tmp);
                
                
                free_transf3d(transfo);
                free_transf3d(transfo_inv);
                free_transf3d(transfo_comb);
                free_transf3d(transfo_opp);
            
            }
            
            
            imx_normalisation_gaussienne2d_histo_joint_norm_p(imref,imrecad,imtref, nb_classe);
            
            imx_norm_rcoeff_3d_p(imtref, imrecad);

            
            //imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_nonsymetrise_p(imref,imreca, imrefd, imrecad, imtref, imtreca, nb_classe);
            
            //imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_p(imref,imreca, imrefd, imrecad, imtref, imtreca, nb_classe);      
            
            
            //imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_iteratif_p(imref,imreca, imrefd, imrecad, imtref, imtreca, nb_classe);
            }
        else
            {
         
            printf("Par defaut, on fait une normalisation moyenne ecart-type symetrique\n");
            imx_norm_seuil_meanecty_symetrique_3d_p(imref, imreca,  imrefd, imrecad,  imtref , imtreca);

            }
        }
            
    if (imtref->rcoeff!=imtreca->rcoeff)
        printf("Attention les rcoeff ne sont pas identiques !!! Le recalage ne sera pas pertinent ....\n");
    
    if (adaptatif==0)
        test_correlation(imtref,imtres,nb_param,masque_param);

    if ((dist_type==9)||(dist_type==10))
        update_TOPO_ALPHA_ROBUST_global(imtref,imtreca, nb_param, param);

    if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
    else
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;


        
    if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)//&&(r==1))
        {
        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(imtref,imtreca,imtres,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);
        //TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*(lambda_reca*lambda_reca+lambda_ref*lambda_ref)/lambda_reca/lambda_reca;
         // pour corriger le terme de regularisation qui est modifie a cause de lambda_reca
        }

    if (dist_type==12)
            update_maxnorm_HistoJoint(imtref,imtreca);

 mini=minimisation(imtref,imtreca,imtres,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);
 
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
    free_imatrix_3d(masque_param);
    }       

    if (dist_type==12)
        free_HistoJoint();
    
    /*calcul de la transformation*/
    base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
 /*liberation de la base*/
   end_base_3d();

}
  /*enregistrement du champ resultat dans un fichier*/
    
  
  // transfo direct
  transfo=cr_transf3d(wdth,hght,dpth,NULL);
  transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
  transfo->typetrans=BSPLINE3D;
  transfo->resol=resolf;transfo->degre=func_type;
  transfo->dx=imtres->dx;
  transfo->dy=imtres->dy;
  transfo->dz=imtres->dz;

    // pour empecher des depalcement superieur a la taille de l'image 
     for (l=0;l<nb_param/3;l++) 
        {
        param[3*l]=MINI(param[3*l],wdth);
        param[3*l+1]=MINI(param[3*l+1],hght);
        param[3*l+2]=MINI(param[3*l+2],dpth);
        param[3*l]=MAXI(param[3*l],-1.0*wdth);
        param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
        param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
        } 
         
        for (l=0;l<nb_param/3;l++) 
        {
        transfo->param[3*l]=param[3*l]*imreca->dx;
        transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
        transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
        }
    
    
  // transfo opposee
  transfo_opp=cr_transf3d(wdth,hght,dpth,NULL);
  transfo_opp->nb_param=nb_param;transfo_opp->param=CALLOC(nb_param,double);
  transfo_opp->typetrans=BSPLINE3D;
  transfo_opp->resol=resolf;transfo_opp->degre=func_type;
  transfo_opp->dx=imtres->dx;
  transfo_opp->dy=imtres->dy;
  transfo_opp->dz=imtres->dz;
        
    for (l=0;l<nb_param;l++)
     transfo_opp->param[l]=-1.0*lambda_ref*transfo->param[l];
             

  // transfo inverse
 transfo_inv=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);
 inverse_transf_3d_p(transfo_opp,transfo_inv,0.01);
                
  // transfo composee
 transfo_comb=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);
 comb_transf_3d(transfo,transfo_inv,transfo_comb,5);
        
    
  
transf_to_field_3d_noalloc(transfo_comb, champ, imref, imreca);
 
 
  if (imres != NULL)
  { 
    
    interpol=inter_qsinc3_3d;
    interpol(imreca,champ,imtres);
    imx_copie_3d_p(imtres,imres);
  }




    if (save_type==3) // il faut transformer les images  dans leur taille d'origine
    {
        imx_undo_mri_pow2_centered_p(imref,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(imreca,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(imres,wdth_old,hght_old,dpth_old);
    }
    
  
  // Enregistrement de la transfo resultat
  if (nomfichres!=NULL)
    {
    if (save_type==3) // on tronque le champ 
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)

                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                
                } 
        
        
        
         transfo_tmp=field_to_transf_3d(champ,imref,imreca); 
    
        
        save_transf_3d(transfo_tmp,nomfichres);
    
        free_transf3d(transfo_tmp);


                champ->width=   wdth;
                champ->height=  hght;
                champ->depth=   dpth;
                
                for (i=0;i<wdth;i++)
                for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                    {champ->raw[i][j][k].x=champ->raw[i][j][k].y=champ->raw[i][j][k].z=0.0;}
                

                imref->width=   wdth;
                imref->height=  hght;
                imref->depth=   dpth;
            
                imreca->width=  wdth;
                imreca->height= hght;
                imreca->depth=  dpth;
            
    /* enreristrement de la transformation de imreca */ 
     
        if (save_type==3) // on tronque le champ 
                {
                transf_to_field_3d_noalloc(transfo, champ, imref, imreca);
                
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                imref->width=   wdth_old;
                imref->height=  hght_old;
                imref->depth=   dpth_old;
            
                imreca->width=  wdth_old;
                imreca->height= hght_old;
                imreca->depth=  dpth_old;
                
                transfo_tmp=field_to_transf_3d(champ,imref,imreca); 
                save_transf_3d(transfo_tmp,nomfichier_reca);
                free_transf3d(transfo_tmp);
                }
            else
                save_transf_3d(transfo,nomfichier_reca); 
        
        
                champ->width=   wdth;
                champ->height=  hght;
                champ->depth=   dpth;
                
                imref->width=   wdth;
                imref->height=  hght;
                imref->depth=   dpth;
            
                imreca->width=  wdth;
                imreca->height= hght;
                imreca->depth=  dpth;
        for (i=0;i<wdth;i++)
                for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                    {champ->raw[i][j][k].x=champ->raw[i][j][k].y=champ->raw[i][j][k].z=0.0;}
                

    /* enreristrement de la transformation de imreca */ 

        if (save_type==3) // on tronque le champ 
                {
                transf_to_field_3d_noalloc(transfo_opp, champ, imref, imreca);
            
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                    imref->width=   wdth_old;
                imref->height=  hght_old;
                imref->depth=   dpth_old;
            
                imreca->width=  wdth_old;
                imreca->height= hght_old;
                imreca->depth=  dpth_old;
                
                transfo_tmp=field_to_transf_3d(champ,imref,imreca); 
                save_transf_3d(transfo_tmp,nomfichier_ref);
                free_transf3d(transfo_tmp);
                }
            else 
                save_transf_3d(transfo_opp,nomfichier_ref); 
 
    }   
        
    

free_transf3d(transfo);
free_transf3d(transfo_inv);
free_transf3d(transfo_comb);
        
        
  /*liberation de la base*/
  // end_base_3d();

  /*caracteristique du champ*/
    aff_log("CHAMP: ");carac_field_3d(champ);

 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
  /*liberation memoire des images temporaires*/
  //if (normalisation_type>0)
  {
   free_grphic3d(imtref);free_grphic3d(imtreca); 
  }
   free_grphic3d(imtres);
     free_grphic3d(imrefd);
     free_grphic3d(imrecad);
     

 


 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 if (nomfichres!=NULL)
    {   
    #ifndef TOPO_COMPARE_CRITERE 
     end_flog(LOG_ENERGY);
    #endif
    }
 free(param);free(param_norm);
 
 return(1);  
}



/*******************************************************************************
**        Project_Field_Bspline_topo_3d                                                
**                                                                        
**       estime la transformation Bspline la plus proche possible d'un champ donn�                           
*******************************************************************************/

void Project_Field_Bspline_topo_3d()
{
 char *quest[6],nomfichres[500],nomfichier[500],str[500],s[250];
 int min_type,i,reg_type=0,im_invariant, im_ROI;
 double Jmin,Jmax,prec,lambda;
 int resolf, e,ans;
 
 
    /*lecture du premier champ*/
    strcpy(nomfichier,GET_FILE("Fichier trf a projeter",&e));
    if(e != 0)
        return ;
    put_file_extension(nomfichier,".trf",nomfichier);
    
     if ((test_fichier_transf_3d(nomfichier))!=1) 
      {PUT_ERROR("This file contains no 3D transformation");return;}

    
    /*nom du fichier resultat*/
    sprintf(str,"File ? in %s",_CurrentPath);
    strcpy(nomfichres,SAVE_FILE(str,NULL,&e));
    put_file_extension(nomfichres,".trf",nomfichres);

 
    sprintf(s,"Voulez vous consider un mask de points invariants\n");
    ans=GET_EXIST(s,&e);
 
    /*question sur les images a recaler*/
    if (ans==1)
     im_invariant=GET_PLACE3D("Mask des points invariants");
    else im_invariant=-100;


    sprintf(s,"Voulez vous consider une ROI\n");
    ans=GET_EXIST(s,&e);
 
    /*question sur les images a recaler*/
    if (ans==1)
     im_ROI=GET_PLACE3D("ROI");
    else im_ROI=-100;


    for(i=0;i<6;i++)
        quest[i]=CALLOC(80,char);
 
    /*minimisation*/
    strcpy(quest[0],"ICM");
    strcpy(quest[1],"Descente gradient");
    strcpy(quest[2],"Levenberg-Marquardt");
  strcpy(quest[3],"Levenberg-Marquardt sans topo");
  strcpy(quest[4],"Levenberg-Marquardt Parallel");

    strcpy(quest[5],"\0");
    
    min_type=GETV_QCM("Minimisation",(const char **)quest);

    
    /*resolution finale*/
    resolf= GET_INT("resolution", 4, &e);

    /*bornes sur le jacobien*/
    Jmin = 2;
    while (Jmin>=1)  
        Jmin = GET_DOUBLE("Jmin<1",0,&e);

    Jmax = 0;
    while (Jmax<=1)  
        Jmax = GET_DOUBLE("Jmax>1",100000,&e);

    
    
    /* regularisation */
    strcpy(quest[0],"Pas de regularisation");
    strcpy(quest[1],"Regularisation membrane elastique");
    strcpy(quest[2],"\0");
        
        
    reg_type= GETV_QCM("Regularisation ?",(const char **)quest);

    
    TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
    if (reg_type==1)
        {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&e);
        reg_type=2; // c'est pour la compatible avec la fonction topo_choose_regularisation
        }

    
    prec = GET_DOUBLE("precision requise ",0.1,&e);

    lambda = GET_DOUBLE("contrainte sur le champ -lambda*u ",0,&e);
    
    if (lambda>0)
        TOPO_PONDERATION_RECALAGE_SYMETRIQUE=1.0/(1.0+lambda);
    else
        TOPO_PONDERATION_RECALAGE_SYMETRIQUE=0.0;

    for(i=0;i<5;i++)
        free(quest[i]);
        
    imx_Project_Field_Bspline_topo_3d(nomfichier,nomfichres,min_type,resolf,Jmin,Jmax,reg_type,prec,im_invariant, im_ROI);

 
   

}


/*******************************************************************************
**        imx_Project_Field_Bspline_topo_3d                                                
**                                                                        
**       estime la transformation Bspline la plus proche possible d'un champ donn�                           
*******************************************************************************/

int imx_Project_Field_Bspline_topo_3d(char *nomfichier,char* nomfichres,int min_type,int resolf,double Jmin,double Jmax,int
reg_type,double prec,int im_mask, int im_ROI)
{

transf3d* transfo,*transfores, *transfo_opp;
field3d* ch;
int i,j,k,wdth,hght,dpth;
grphic3d *mask, *ROI;    


if (im_mask==-100)
mask=NULL;
else
mask=ptr_img_3d(im_mask);

 
if (im_ROI==-100)
ROI=NULL;
else
ROI=ptr_img_3d(im_ROI);



transfo=load_transf_3d(nomfichier);
ch=transf_to_field_3d(transfo,NULL,NULL);


transfores= cr_transf3d(ch->width,ch->height,ch->depth,NULL); 
transfores->typetrans=BSPLINE3D;
transfores->degre=1;
transfores->resol=resolf;

if (TOPO_PONDERATION_RECALAGE_SYMETRIQUE>0.0)
{
transfo_opp= cr_transf3d(ch->width,ch->height,ch->depth,NULL); 
transfo_opp->typetrans=BSPLINE3D;
transfo_opp->degre=1;
transfo_opp->resol=resolf;
}
else
transfo_opp=NULL;


wdth=ch->width;
hght=ch->height;
dpth=ch->depth;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
    {
    ch->raw[i][j][k].x = ch->raw[i][j][k].x/transfo->dx;
    ch->raw[i][j][k].y = ch->raw[i][j][k].y/transfo->dy;
    ch->raw[i][j][k].z = ch->raw[i][j][k].z/transfo->dz;
    }



if (mask!=NULL)
 if ((ch->width!=mask->width)||(ch->height!=mask->height)||(ch->depth!=mask->depth))
    {
      PUT_WARN("Le champ et le mask de points invariants n'ont pas la meme taille !");
      printf("Le champ et le mask de points invariants n'ont pas la meme taille !\n");
      free_transf3d(transfo); free_transf3d(transfores); 
      free_field3d(ch);
      return(-1) ;
    }



if (ROI!=NULL)
 if ((ch->width!=ROI->width)||(ch->height!=ROI->height)||(ch->depth!=ROI->depth))
    {
      PUT_WARN("Le champ et la ROI n'ont pas la meme taille !");
      printf("Le champ et la ROI n'ont pas la meme taille  !\n");
      free_transf3d(transfo); free_transf3d(transfores); 
      free_field3d(ch);
      return(-1) ;
    }


imx_Project_Field_Bspline_topo_3d_p(ch,transfores,transfo_opp, min_type, resolf, Jmin, Jmax, reg_type,prec,mask,ROI);

    
save_transf_3d(transfores,nomfichres);

free_transf3d(transfo); free_transf3d(transfores); 

if (TOPO_PONDERATION_RECALAGE_SYMETRIQUE>0.0)
    {
        char nomfichier_opp[1000];
        char *ext2;
        
        strcpy(nomfichier_opp,nomfichres);
        ext2=strstr(nomfichier_opp,".trf");
        if (ext2!=NULL) *ext2='\0';

        strcat(nomfichier_opp,"_opp.trf");

     
    save_transf_3d(transfo_opp,nomfichier_opp);
    free_transf3d(transfo_opp); 
    }
    
free_field3d(ch); 
return(1);
}



/*******************************************************************************
**        imx_Project_Field_Bspline_topo_3d_p                                                
**                                                                        
**       estime la transformation Bspline la plus proche possible d'un champ donn�                           
*******************************************************************************/
// Attention il faut donner des champs en voxels !!!!!!!!!!!!
void imx_Project_Field_Bspline_topo_3d_p(field3d* ch, transf3d* transfores,transf3d* transfo_opp, int min_type,int resolf,double Jmin,double Jmax,int reg_type,double prec, grphic3d *mask, grphic3d *ROI)
{
int inter_type=1;
int i,j,k;
int func_type=1;
int wdth,hght,dpth,nb_param;
int tdebut,tfin,continu,nbmax_param;
double *param,*param_norm;
double  mini;
grphic3d *toto; // Attention, ce n'est pas une image alloue, mais juste pour passer des infos sur la taille (pour des raisons de compatibilite)
field3d *champ;

InterpolationFct interpol;
dist_func_locale_t distance;
reg_func_locale_t regularisation;
min_func_locale_t minimisation;
scal_func_t scal_func;

/*************mise en place des parametres du recalage******/
scal_func=topo_choose_bspline(func_type);
interpol=imx_choose_interpolation_fct(inter_type);  
minimisation=topo_choose_optimisation(min_type);
    
distance=Energie_quad_sous_champ_locale_3d;

if (reg_type>0)
regularisation=regularisation_energie_membrane_local;


wdth=ch->width;hght=ch->height;dpth=ch->depth;


if (ROI != NULL)
    toto = ROI;
else
    {
    toto = cr_grphic3d_modif(wdth, hght, dpth, 0, 1, 1);
    for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
        toto->mri[i][j][k]=1;
    }


/*allocation memoire de champ*/
 champ=cr_field3d(wdth,hght,dpth);



continu=0;
for (i=0;i<10;i++)
    {
    if (wdth==pow(2.0,i)) continu++;
    if (hght==pow(2.0,i)) continu++;
    if (dpth==pow(2.0,i)) continu++;
    }

if (continu!=3)
    {
        PUT_WARN("La taille du champ n'est pas une puissance de 2 ");
    printf("La taille du champ n'est pas une puissance de 2 \n");
    exit(0);
        }




tdebut=time(NULL);
    nbmax_param=3*pow((pow(2,resolf)-1),3);
 
/* Allocation */
if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); exit(0); }
if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); exit(0); }

    
 init_log(NULL);

_ParamRecalageBspline.chref=ch;
    
 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r,D,ind,flag,topi,topj,topk,compt;
  int ***masque_param;
  int *x0,*x1,*y0,*y1,*z0,*z1;
  double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE,dx,dy,dz;
  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;
        
  /* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}


      
    /* Debut traitement */
    for (r=resol;r<=resolf;r++)
    {

    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
    D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    
        
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    
    /*  Remplissage de masque_param */   // prevoir peut etre de restreindre la zone d'optimisation
    
    base_to_field_3d(nb_param,param,champ,NULL,NULL);
    
    compt=0;
    for(i=0;i<D;i++)
    for(j=0;j<D;j++)
    for(k=0;k<D;k++)
        {
        ind=1.0*TOP_conv_ind(i,j,k,nb_param)/3;
        flag=0;topi=x0[ind];
                while ((topi<x1[ind])&&(flag==0))
                    {
                    topj=y0[ind];
                    while ((topj<y1[ind])&&(flag==0))
                        {
                        topk=z0[ind];
                        while ((topk<z1[ind])&&(flag==0))                   
                            {
                            dx = champ->raw[topi][topj][topk].x -_ParamRecalageBspline.chref->raw[topi][topj][topk].x; 
                            dy = champ->raw[topi][topj][topk].y -_ParamRecalageBspline.chref->raw[topi][topj][topk].y; 
                            dz = champ->raw[topi][topj][topk].z -_ParamRecalageBspline.chref->raw[topi][topj][topk].z; 
                            
                            
                            if(sqrt(dx*dx+dy*dy+dz*dz)>prec)
                                 flag=1;                    
                            
                                
                            topk++;
                            }
                        topj++;
                        }
                    topi++;
                    }
                if (flag==0)
                    {
                    masque_param[i][j][k]=0;
                    compt++;
                    }
                else
                    {
                    masque_param[i][j][k]=1;
                    }
                }
        
        
        if (mask!=NULL)
        for(i=0;i<D;i++)
        for(j=0;j<D;j++)
        for(k=0;k<D;k++)
        {
        ind=1.0*TOP_conv_ind(i,j,k,nb_param)/3;
        flag=0;topi=x0[ind];
                while ((topi<x1[ind])&&(flag==0))
                    {
                    topj=y0[ind];
                    while ((topj<y1[ind])&&(flag==0))
                        {
                        topk=z0[ind];
                        while ((topk<z1[ind])&&(flag==0))                   
                            {
                                
                            if(mask->mri[topi][topj][topk]>0)
                                 flag=1;                    
                            
                                
                            topk++;
                            }
                        topj++;
                        }
                    topi++;
                    }
                if ((flag==1)&&(masque_param[i][j][k]>0))
                    {
                    masque_param[i][j][k]=0;
                    compt++;
                    }
        }

        printf("Reduction du nombre de parametres : %.2f %% \n",300.0*compt/nb_param);



        
    if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
    else
        TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
        
    if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(toto,toto,toto,NULL,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,NULL);
    
        mini=minimisation(toto,toto,toto,NULL,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,NULL);
    
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
    free_imatrix_3d(masque_param);
    }       

    

transfores->width=ch->width;transfores->height=ch->height;transfores->depth=ch->depth;
transfores->dx=ch->dx;transfores->dy=ch->dy;transfores->dz=ch->dz;
transfores->typetrans=BSPLINE3D;




     
         /* pour empecher des depalcement superieur a la taille de l'image */
         for (l=0;l<nb_param/3;l++) 
            {
            param[3*l]=MINI(param[3*l],wdth);
            param[3*l+1]=MINI(param[3*l+1],hght);
            param[3*l+2]=MINI(param[3*l+2],dpth);
            param[3*l]=MAXI(param[3*l],-1.0*wdth);
            param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
            param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
            } 
         
  
     
        for (l=0;l<nb_param/3;l++) 
        {
        param[3*l]=param[3*l]*ch->dx;
        param[3*l+1]=param[3*l+1]*ch->dy;
        param[3*l+2]=param[3*l+2]*ch->dz;
        }

  
  /*liberation de la base*/
   end_base_3d();
}

 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 
 transfores->nb_param=nb_param;transfores->param=param;


if (TOPO_PONDERATION_RECALAGE_SYMETRIQUE>0.0)
    {
    double *opp_param;
    double lambda=(1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
    int l;
    opp_param=calloc(transfores->nb_param,sizeof(double));
    
        
    for (l=0;l<transfores->nb_param;l++)
         opp_param[l]=-1.0*lambda*transfores->param[l];
        
    
    transfo_opp->width= transfores->width;
    transfo_opp->height=transfores->height;
    transfo_opp->depth= transfores->depth;
    
    transfo_opp->dx=transfores->dx;
    transfo_opp->dy=transfores->dy;
    transfo_opp->dz=transfores->dz;
    transfo_opp->typetrans=BSPLINE3D;
    transfo_opp->nb_param=transfores->nb_param;transfo_opp->param=opp_param;
    }    

 //free(param);

 /*liberation memoire du champ*/
 free_field3d(champ);

if (ROI == NULL)
    free_grphic3d(toto);

 free(param_norm); 

}



/*******************************************************************************
**        dti_matching_Bspline_topo_3d                                                
**                                                                        
**       recalage Bspline de deux serie DTI 3D avec conservation de la topologie                               
*******************************************************************************/
#ifndef COMPILE_FOR_MEDIPY

void dti_matching_Bspline_topo_3d()
{
char *quest[15],*nomfichres;
int dist_type,inter_type,min_type,save_type,func_type,i,normalisation_type, reg_type=0;
int nb_classe;
double Jmin,Jmax;
int resolf, e, err;
 
_ParamRecalageBspline.item_reca=imx_select_dti_serie_3d_p();
_ParamRecalageBspline.item_ref=imx_select_dti_serie_3d_p();
_ParamRecalageBspline.item_res=dti_copy_3d(_ParamRecalageBspline.item_ref);

if((_ParamRecalageBspline.item_ref==NULL)||(_ParamRecalageBspline.item_reca==NULL))
    return;

if((_ParamRecalageBspline.item_ref->S0->width!=_ParamRecalageBspline.item_reca->S0->width)||
(_ParamRecalageBspline.item_ref->S0->height!=_ParamRecalageBspline.item_reca->S0->height)||
(_ParamRecalageBspline.item_ref->S0->depth!=_ParamRecalageBspline.item_reca->S0->depth))
    {
    printf("Attention!!! Les series DTI n'ont pas les memes dimensions ! \n");
    return;
    }

if((_ParamRecalageBspline.item_ref->S0->dx!=_ParamRecalageBspline.item_reca->S0->dx)||
(_ParamRecalageBspline.item_ref->S0->dy!=_ParamRecalageBspline.item_reca->S0->dy)||
(_ParamRecalageBspline.item_ref->S0->dz!=_ParamRecalageBspline.item_reca->S0->dz))
    {
    printf("Attention!!! Les series DTI n'ont pas les memes dx,dy ou dz ! \n");
    return;
    }

if(_ParamRecalageBspline.item_ref->nb_directions!=_ParamRecalageBspline.item_reca->nb_directions)
    {
    printf("Attention!!! Les series DTI n'ont pas le meme nombre de directions ! \n");
    return;
    }


if (_ParamRecalageBspline.item_ref==_ParamRecalageBspline.item_reca)
    {
    printf("Attention!!! Les series DTI source et cible sont identiques !! \n");
    return;
    }
for(i=0;i<5;i++)
    quest[i]=CALLOC(80,char);
 
//  type de fonction
func_type=1;
  

// distance
dist_type=13;
  
//interpolation
inter_type=4;

//minimisation
strcpy(quest[0],"ICM");
//strcpy(quest[1],"Descente gradient");
//strcpy(quest[2],"Levenberg-Marquardt");
//strcpy(quest[3],"Levenberg-Marquardt sans Topo");
strcpy(quest[2],"\0");

min_type=GETV_QCM("Minimisation",(const char **)quest);

   
// question sur l'enregistrement du champ resultat
nomfichres=quest_save_result_3d(&save_type);
    
// resolution finale
resolf= GET_INT("resolution", 4, &e);

// bornes sur le jacobien
Jmin = 2;
while (Jmin>=1)  
    Jmin = GET_DOUBLE("Jmin<1",0,&err);

Jmax = 0;
while (Jmax<=1)  
    Jmax = GET_DOUBLE("Jmax>1",100000,&err);

nb_classe=-1;
normalisation_type=2;


// regularisation 
strcpy(quest[0],"Pas de regularisation");
strcpy(quest[1],"Regularisation membrane elastique");
strcpy(quest[2],"\0");
        
        
reg_type= GETV_QCM("Regularisation ?",(const char **)quest);

    
TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
if (reg_type==1)
    {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&e);
    reg_type=2; // c'est pour la compatible avec la fonction topo_choose_regularisation
    }



for(i=0;i<5;i++)
    free(quest[i]);
        
imx_dti_matching_Bspline_topo_3d_p(func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax,normalisation_type,nb_classe);

 
 if (nomfichres)
  free(nomfichres);

}


/*******************************************************************************
**        imx_dti_matching_Bspline_topo_3d_p                                            
**                                                                        
**       recalage Bspline de 2 images 3D avec conservation de la topologie                                   
*******************************************************************************/
/*
  imref : pointeur sur l'image sur laquelle on recale imreca
  imreca : pointeur sur l'image a recaler
  imres : pointeur sur l'image ou l'on stocke le resultat du recalage de imreca
    (peut etre NULL)
  func_type : numero de la fonction d'interpolation du champ (1 pour Bspline
    degre 1)
  dist_type : numero de la fonction de distance (0 pour distance quadratique)
  inter_type : numero de la fonction d'interpolation des images (1 pour lineaire)
  min_type : numero de la fonction de minimisation utilisee aux plus faibles
    echelles de la transformation (5 pour quasi-newton)
  save_type : 0 si on ne veut rien sauver, 1 si on veut sauver les parametres
    de la transformation (fichier trf), 2 si on veut sauver le champ
  nomfichres : fichier ou le resultat est sauve, si (save_type != 0)
  resolf : resolution finale du recalage ; typiquement 4 pour des images
    128*128*128
  nomfichiertrf : nom du fichier trf a partir duquel on commence le recalage
    (peut etre NULL, mais ce n'est pas recommande ; il convient d'effectuer
    un recalage affine et de donner la transformation resultante comme etat
    initial)
  renormalisation : TRUE si l'on veut effectuer une renormalisation des
    images (recommande)
*/
int imx_dti_matching_Bspline_topo_3d_p(int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb;
 double *param,*param_norm;
 grphic3d *imtres;
 field3d *champ;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 norm_func_locale_t normalisation;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;

#ifndef SAVE_INTERMEDIAIRE_2 
char nomfichier[255];
char temp[2];
char *ext2;
#endif

#ifndef TOPO_COMPARE_CRITERE 
    char nomfichier_nrj[255];
    char *ext;
#endif
#ifndef WIN32
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif

save_type=3;

wdth=_ParamRecalageBspline.item_ref->S0->width;hght=_ParamRecalageBspline.item_ref->S0->height;dpth=_ParamRecalageBspline.item_ref->S0->depth;


if (save_type==3) /* il faut transformer les images en puissance de 2 */
    {
    wdth_old=wdth; hght_old=hght; dpth_old=dpth;
    
    imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.item_ref->S0);
    imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.item_reca->S0);
    imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.item_res->S0);
    
    for (i=0;i<_ParamRecalageBspline.item_ref->nb_directions;i++)
        {
        imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.item_ref->img_3d[i]);
        imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.item_reca->img_3d[i]);
        imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.item_res->img_3d[i]);
        }
    
    wdth=_ParamRecalageBspline.item_ref->S0->width;hght=_ParamRecalageBspline.item_ref->S0->height;dpth=_ParamRecalageBspline.item_ref->S0->depth;
    }
    
tdebut=time(NULL);
nbmax_param=3*pow((pow(2,resolf)-1),3);
 
/* Allocation */
if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
init_log(NULL);
aff_log("\nRECALAGE ");
 
/*************mise en place des parametres du recalage******/
scal_func=topo_choose_bspline(func_type);
distance=topo_choose_distance(dist_type);
interpol=imx_choose_interpolation_fct(inter_type);  
minimisation=topo_choose_optimisation(min_type);
normalisation=topo_choose_normalisation(normalisation_type,nb_classe);
regularisation=topo_choose_regularisation(reg_type);
    
aff_log("\n");
aff_log("recalage de %s sur %s \n",(_ParamRecalageBspline.item_reca->S0->patient).name,(_ParamRecalageBspline.item_ref->S0->patient).name);
 
    

/*allocation memoire de champ*/
 champ=cr_field3d(wdth,hght,dpth);
     
imtres=cr_grphic3d(_ParamRecalageBspline.item_reca->S0);imx_copie_param_3d_p(_ParamRecalageBspline.item_reca->S0,imtres);
 
 /********************RECALAGE ONDELETTE******************************************/
 {
int l,resol,r,D;
int ***masque_param;
int *x0,*x1,*y0,*y1,*z0,*z1;
double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
nb_param=init_base_3d(wdth,hght,dpth,scal_func);
resol=BASE3D.resol;
        
/* Initialisation tableau */
for (l=0;l<nb_param;l++)
   {param[l]=0.0;}

if (nomfichres!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    strcpy(nomfichier_nrj,nomfichres);
    ext=strstr(nomfichier_nrj,".trf");
    if (ext!=NULL) *ext='\0';
    strcat(nomfichier_nrj,"_energy.txt");
    LOG_ENERGY=init_flog(nomfichier_nrj);      
    #endif
    }
     
/* Debut traitement */
for (r=resol;r<=resolf;r++)
    {
    TOPO_SEUIL_SIGNIFICATIVITE=0.05;

    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
    D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    base_to_field_3d(nb_param,param,champ,NULL,NULL);
    interpol(_ParamRecalageBspline.item_reca->S0,champ,imtres); 

    
            
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    topo_masque_param (_ParamRecalageBspline.item_ref->S0,_ParamRecalageBspline.item_reca->S0,imtres,masque_param, nb_param);




if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
else
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;


if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(_ParamRecalageBspline.item_ref->S0,_ParamRecalageBspline.item_reca->S0,imtres,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);

mini=minimisation(_ParamRecalageBspline.item_ref->S0,_ParamRecalageBspline.item_reca->S0,imtres,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);

#ifndef SAVE_INTERMEDIAIRE_2
    
if (nomfichres!=NULL)
    {
    /*calcul de la transformation*/
    base_to_field_3d(nb_param,param,champ,NULL,NULL);

    //si l'extension .trf existe on la supprime
    strcpy(nomfichier,nomfichres);
    ext2=strstr(nomfichier,".trf");
    if (ext2!=NULL) *ext2='\0';

    strcat(nomfichier,"_");
    temp[0]=(char)(r+48);
    temp[1]=0;
    strcat(nomfichier,temp);
      
        transf3d *transfo; 
    
        if ((save_type==0)||(save_type==3))
        {
        if (save_type==3) /* on tronque le champ */
        {
        ideb=(int)((wdth-wdth_old)*0.5);
        jdeb=(int)((hght-hght_old)*0.5);
        kdeb=(int)((dpth-dpth_old)*0.5);

        for (i=0;i<wdth_old;i++)
        for (j=0;j<hght_old;j++)
        for (k=0;k<dpth_old;k++)
            {
            champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
            }
                
        champ->width=   wdth_old;
        champ->height=  hght_old;
        champ->depth=   dpth_old;
                
        _ParamRecalageBspline.item_ref->S0->width=  wdth_old;
        _ParamRecalageBspline.item_ref->S0->height= hght_old;
        _ParamRecalageBspline.item_ref->S0->depth=  dpth_old;
        
        _ParamRecalageBspline.item_reca->S0->width= wdth_old;
        _ParamRecalageBspline.item_reca->S0->height=    hght_old;
        _ParamRecalageBspline.item_reca->S0->depth= dpth_old;
        } 
        
        transfo=field_to_transf_3d(champ,imref,imreca); 
            
        if (save_type==3)
            {
            champ->width=   wdth;
            champ->height=  hght;
            champ->depth=   dpth;
                
            _ParamRecalageBspline.item_ref->S0->width=   wdth;
            _ParamRecalageBspline.item_ref->S0->height=  hght;
            _ParamRecalageBspline.item_ref->S0->depth=   dpth;
        
            _ParamRecalageBspline.item_reca->S0->width=  wdth;
            _ParamRecalageBspline.item_reca->S0->height= hght;
            _ParamRecalageBspline.item_reca->S0->depth=  dpth;
            }
        
        }
        else
            {  
            transfo=cr_transf3d(wdth,hght,dpth,NULL);
            transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
            transfo->typetrans=BSPLINE3D;
            transfo->resol=resolf;transfo->degre=func_type;
        transfo->dx=imtres->dx;
        transfo->dy=imtres->dy;
        transfo->dz=imtres->dz;
         
         for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*_ParamRecalageBspline.item_reca->S0->dx;
            transfo->param[3*l+1]=param[3*l+1]*_ParamRecalageBspline.item_reca->S0->dy;
            transfo->param[3*l+2]=param[3*l+2]*_ParamRecalageBspline.item_reca->S0->dz;
            } 
    }
    save_transf_3d(transfo,nomfichier);
    free_transf3d(transfo);
     
    

   } 
    #endif
 
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
free_imatrix_3d(masque_param);
}       

/*calcul de la transformation*/
base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
  

if (save_type==3) /* il faut transformer les images  dans leur taille d'origine*/
    {
    imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.item_ref->S0,wdth_old,hght_old,dpth_old);
    imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.item_reca->S0,wdth_old,hght_old,dpth_old);
    imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.item_res->S0,wdth_old,hght_old,dpth_old);
    
    
    for (i=0;i<_ParamRecalageBspline.item_ref->nb_directions;i++)
        {
        imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.item_ref->img_3d[i],wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.item_reca->img_3d[i],wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.item_res->img_3d[i],wdth_old,hght_old,dpth_old);
        }
    
    }
    
    
/*enregistrement du champ resultat dans un fichier*/
if (nomfichres!=NULL)
{ 
    transf3d *transfo; 
    
        if ((save_type==0)||(save_type==3))
        {
        if (save_type==3) /* on tronque le champ */
            {
            ideb=(int)((wdth-wdth_old)*0.5);
            jdeb=(int)((hght-hght_old)*0.5);
            kdeb=(int)((dpth-dpth_old)*0.5);

            for (i=0;i<wdth_old;i++)
            for (j=0;j<hght_old;j++)
            for (k=0;k<dpth_old;k++)
                {
                champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                }
                
            champ->width=   wdth_old;
            champ->height=  hght_old;
            champ->depth=   dpth_old;
            } 
                
         transfo=field_to_transf_3d(champ,_ParamRecalageBspline.item_ref->S0,_ParamRecalageBspline.item_reca->S0); 
    }
        else
    {  
     
        transfo=cr_transf3d(wdth,hght,dpth,NULL);
        transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
        transfo->typetrans=BSPLINE3D;
        transfo->resol=resolf;transfo->degre=func_type;
    transfo->dx=imtres->dx;
    transfo->dy=imtres->dy;
    transfo->dz=imtres->dz;
         
     /* pour empecher des depalcement superieur a la taille de l'image */
     for (l=0;l<nb_param/3;l++) 
        {
        param[3*l]=MINI(param[3*l],wdth);
        param[3*l+1]=MINI(param[3*l+1],hght);
        param[3*l+2]=MINI(param[3*l+2],dpth);
        param[3*l]=MAXI(param[3*l],-1.0*wdth);
        param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
        param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
        } 
         
     for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*_ParamRecalageBspline.item_reca->S0->dx;
            transfo->param[3*l+1]=param[3*l+1]*_ParamRecalageBspline.item_reca->S0->dy;
            transfo->param[3*l+2]=param[3*l+2]*_ParamRecalageBspline.item_reca->S0->dz;
            } 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 
 }

    
/*liberation de la base*/
end_base_3d();

/*caracteristique du champ*/
aff_log("CHAMP: ");carac_field_3d(champ);

/*liberation memoire des variables dynamiques*/ 
/*liberation memoire du champ*/
free_field3d(champ);

/*liberation memoire des images temporaires*/
free_grphic3d(imtres);

/*affichage du temps de calcul total*/
tfin=time(NULL);
{
int h,m,s,ttotal;
ttotal=tfin-tdebut;
h=ttotal/3600;
m=(ttotal-h*3600)/60;
s=(ttotal-h*3600-m*60);
aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
}
end_log();
if (nomfichres!=NULL)
    {   
    #ifndef TOPO_COMPARE_CRITERE 
     end_flog(LOG_ENERGY);
    #endif
    }
free(param);free(param_norm);
 
 return(1);  
}
#endif //fin COMPILE_FOR_MEDIPY

/*******************************************************************************
**        groupwise_matching_Bspline_topo_3d                                                
**                                                                        
**       recalage Bspline de N images 3D avec conservation de la topologie                               
*******************************************************************************/

void groupwise_matching_Bspline_topo_3d()
{

char * file;
char * nb_img_ini;
int e;
int nb_img=0,num_ref, ans;
char *quest[15],*nomfichres;
int dist_type=0,inter_type=1,min_type=0,save_type=0,func_type=0,i,normalisation_type=0, reg_type=0,resolf,nb_classe=0;
double Jmin,Jmax;
char    s[256];
char *nomfichiertrf;
 


file=strdup(GET_FILE("Serie d'images a recaler (*.ipb)", &e)); 
    if (e)
    {
        free(file);
        return;
    }    

nb_img_ini = getheader_interfile(file,"number of images=", ENTIER,0);   
nb_img=atoi(nb_img_ini);

_ParamRecalageBspline.nb_tot_groupwise=nb_img;


TOPO_PONDERATION_RECALAGE_SYMETRIQUE=((double)nb_img-1.0)/(double)nb_img; //ie ponderation de l'image de reference

num_ref= GET_INT("numero de l'image a recaler [1..N]", 1, &e);

_ParamRecalageBspline.nb_reca_groupwise=num_ref;

if ((num_ref<1)||(num_ref>nb_img))
    {
    printf("Le numero d'image n'est pas licite !!! \n");
    free(file);
    return;
    }
num_ref=num_ref-1; // pour pointer sur le bon element du tableau _ParamRecalageBspline.serie_groupwise

 for(i=0;i<15;i++)
    quest[i]=CALLOC(80,char);


/*distance*/
strcpy(quest[0],"variance nonsym");
strcpy(quest[1],"variance sym");
strcpy(quest[2],"\0");

dist_type=GETV_QCM("Critere",(const char **)quest);

 
if (dist_type==0)
    TOPO_PONDERATION_RECALAGE_SYMETRIQUE=1.0; 


 
/*minimisation*/
strcpy(quest[0],"ICM");
strcpy(quest[1],"Descente gradient");
strcpy(quest[2],"Levenberg-Marquardt");
strcpy(quest[3],"Levenberg-Marquardt sans Topo");
  strcpy(quest[4],"Levenberg-Marquardt Parallel");
strcpy(quest[5],"\0");

min_type=GETV_QCM("Minimisation",(const char **)quest);

  
/*question sur l'enregistrement du champ resultat*/
nomfichres=quest_save_result_3d(&save_type);
    
/*resolution finale*/
resolf= GET_INT("resolution", 4, &e);

/*bornes sur le jacobien*/
Jmin = 2;
while (Jmin>=1)  
    Jmin = GET_DOUBLE("Jmin<1",0,&e);

Jmax = 0;
while (Jmax<=1)  
    Jmax = GET_DOUBLE("Jmax>1",100000,&e);

/* regularisation */
strcpy(quest[0],"Pas de regularisation");
strcpy(quest[1],"Regularisation membrane elastique");
strcpy(quest[2],"Distance a l'identite");
strcpy(quest[3],"\0");
        
        
reg_type= GETV_QCM("Regularisation ?",(const char **)quest);

    
TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
    
if (reg_type==2)
    {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&e);
    reg_type=7; // c'est pour �re compatible avec la fonction topo_choose_regularisation
    }


if (reg_type==1)
    {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&e);
    reg_type=2; // c'est pour �re compatible avec la fonction topo_choose_regularisation
    }


 sprintf(s,"Use a transformation as initialisation?\n");
 ans=GET_EXIST(s,&e);
  if (ans==1)
  {
    char nomfichier[250];
    strcpy(nomfichier,GET_FILE("*.trf",&e));
    put_file_extension(nomfichier,".trf",nomfichier);
    nomfichiertrf=CALLOC(strlen(nomfichier),char);
    strcpy(nomfichiertrf,nomfichier);
  }
  else
    nomfichiertrf = NULL;
 


for(i=0;i<15;i++)
    free(quest[i]);

 /*type de fonction*/
  func_type=1;

        
imx_groupwise_matching_Bspline_topo_3d(file,nb_img,num_ref,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax,normalisation_type,nb_classe,1,nomfichiertrf);

 
 if (nomfichres)
  free(nomfichres);


if (nomfichiertrf)
  free(nomfichiertrf);

free(file);

}



/*******************************************************************************
**        imx_groupwise_matching_Bspline_topo_3d                                             
**                                                                       
**       recalage Bspline de N images 3D avec conservation de la topologie                                  
*******************************************************************************/
int imx_groupwise_matching_Bspline_topo_3d(char* file, int nb_img, int num_ref, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax, int normalisation_type, int nb_classe,int adaptatif, char* nomfichiertrf)
{
int l,continu,wdth,hght,dpth,max=0,step=0,i,j,k;
float dx,dy,dz;

// Mise a jour de la structure ParamRecalageBspline
_ParamRecalageBspline.serie_groupwise=CALLOC(nb_img,ptr_grphic3d);

wdth = getmri_width(file, 1, MODE_3D);
hght = getmri_height(file, 1, MODE_3D);
dpth = getmri_depth(file, 1, MODE_3D);

dx =  getmri_dx (file, 1, MODE_3D);
dy =  getmri_dy (file, 1, MODE_3D);
dz =  getmri_dz (file, 1, MODE_3D);


 for (l = 0; l<nb_img;l++) {
    _ParamRecalageBspline.serie_groupwise[l] = cr_grphic3d_modif(wdth, hght, dpth, 1., 1., 1);      

#ifndef COMPILE_FOR_MEDIPY
    load_mri_3d(file,  _ParamRecalageBspline.serie_groupwise[l], l+1, max, step);
#endif //fin COMPILE_FOR_MEDIPY

  }

// test si toutes les images de la serie ont la meme taille
 for (l = 0; l<nb_img;l++) 
    {
     if ((_ParamRecalageBspline.serie_groupwise[l]->width!=wdth)||(_ParamRecalageBspline.serie_groupwise[l]->height!=hght)||(_ParamRecalageBspline.serie_groupwise[l]->depth!=dpth))
        { PUT_WARN("Les images n'ont pas la meme taille !");
            printf("Les images n'ont pas la meme taille !\n");
            return(-1) ;
        }
    }


// test si toutes les images de la serie ont le meme dx,dy,dz
 for (l = 0; l<nb_img;l++) 
    {
     if ((_ParamRecalageBspline.serie_groupwise[l]->dx!=dx)||(_ParamRecalageBspline.serie_groupwise[l]->dy!=dy)||(_ParamRecalageBspline.serie_groupwise[l]->dz!=dz))
        {   
        PUT_WARN("Les images n'ont pas les memes dx,dy,dz !");
            printf("Les images n'ont pas les memes dx,dy,dz !\n");
        return(-1) ;
        }
    }

 
 
// Mise a jour de reca_groupwise,ref_groupwise et ref2_groupwise
_ParamRecalageBspline.reca_groupwise = cr_grphic3d_modif(wdth, hght, dpth, 1., 1., 1);      
imx_copie_3d_p( _ParamRecalageBspline.serie_groupwise[num_ref],_ParamRecalageBspline.reca_groupwise); 

_ParamRecalageBspline.dbl_ref_groupwise = alloc_dmatrix_3d(wdth, hght, dpth);

// cacul de la somme des N-1 images
_ParamRecalageBspline.ref_groupwise = cr_grphic3d_modif(wdth, hght, dpth, 1., 1., 1);   

for (i=0;i<wdth; i++)
for (j=0;j<hght; j++)
for (k=0;k<dpth; k++)
    _ParamRecalageBspline.dbl_ref_groupwise[i][j][k]=0.0;


for (l=0;l<nb_img;l++)
    if (l!=num_ref)
    {
    for (i=0;i<wdth; i++)
    for (j=0;j<hght; j++)
    for (k=0;k<dpth; k++)
        _ParamRecalageBspline.dbl_ref_groupwise[i][j][k]+=_ParamRecalageBspline.serie_groupwise[l]->mri[i][j][k]*_ParamRecalageBspline.serie_groupwise[l]->rcoeff;
    
    }   
    
imx_convert_dmatrix3d_grphic3d(_ParamRecalageBspline.dbl_ref_groupwise,_ParamRecalageBspline.ref_groupwise);



// cacul de la somme des N-1 images au carre
_ParamRecalageBspline.ref2_groupwise = cr_grphic3d_modif(wdth, hght, dpth, 1., 1., 1);      
_ParamRecalageBspline.dbl_ref2_groupwise = alloc_dmatrix_3d(wdth, hght, dpth);

for (i=0;i<wdth; i++)
for (j=0;j<hght; j++)
for (k=0;k<dpth; k++)
    _ParamRecalageBspline.dbl_ref2_groupwise[i][j][k]=0.0;

for (l=0;l<nb_img;l++)
    if (l!=num_ref)
    {
    for (i=0;i<wdth; i++)
    for (j=0;j<hght; j++)
    for (k=0;k<dpth; k++)
        _ParamRecalageBspline.dbl_ref2_groupwise[i][j][k]+=_ParamRecalageBspline.serie_groupwise[l]->mri[i][j][k]
                *_ParamRecalageBspline.serie_groupwise[l]->mri[i][j][k]
                *_ParamRecalageBspline.serie_groupwise[l]->rcoeff
                *_ParamRecalageBspline.serie_groupwise[l]->rcoeff;
        
    }   
    
imx_convert_dmatrix3d_grphic3d(_ParamRecalageBspline.dbl_ref2_groupwise,_ParamRecalageBspline.ref2_groupwise);


 
continu=0;
for (l=0;l<10;l++)
    {
    if (wdth==pow(2.0,l)) continu++;
    if (hght==pow(2.0,l)) continu++;
    if (dpth==pow(2.0,l)) continu++;
    }

if (continu!=3)
    {
        PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
    printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
    save_type=3;
    }

imx_groupwise_matching_Bspline_topo_3d_p(func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax,normalisation_type, nb_classe,adaptatif,nomfichiertrf);
 
 

// Liberation memoire des elements de la structure ParamRecalageBspline
for (l = 0; l<nb_img;l++) 
    {
    free_grphic3d(_ParamRecalageBspline.serie_groupwise[l]);
    }

free(_ParamRecalageBspline.serie_groupwise);
free_grphic3d(_ParamRecalageBspline.reca_groupwise);
free_grphic3d(_ParamRecalageBspline.ref_groupwise);
free_grphic3d(_ParamRecalageBspline.ref2_groupwise);
free_dmatrix_3d(_ParamRecalageBspline.dbl_ref_groupwise);
free_dmatrix_3d(_ParamRecalageBspline.dbl_ref2_groupwise);
  return(1);
}

/*******************************************************************************
**        imx_groupwise_matching_Bspline_topo_3d                                             
**                                                                       
**       recalage Bspline de N images 3D avec conservation de la topologie                                  
*******************************************************************************/
int imx_groupwise_matching_Bspline_topo_3d_p(int func_type,int dist_type, int reg_type, int inter_type,int min_type,int save_type, char* nomfichres,int resolf, double Jmin, double
Jmax, int normalisation_type, int nb_classe,int adaptatif, char* nomfichiertrf)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb;
 double *param,*param_norm;
 grphic3d *imtres,*imrefd,*imrecad;
 field3d *champ;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 norm_func_locale_t normalisation;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
 char nomfichier_ref[1000];
 char nomfichier_reca[1000];
  transf3d *transfo,*transfo_opp,*transfo_inv, *transfo_comb, *transfo_tmp; 
   int l;
 

char *ext2;

#ifndef SAVE_INTERMEDIAIRE_2 
char nomfichier[255];
char temp[2];
#endif

#ifndef TOPO_COMPARE_CRITERE 
    char nomfichier_nrj[255];
    char *ext;
#endif
#ifndef WIN32
    setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif

printf("TOPO_PONDERATION_RECALAGE_GROUPWISE = %f \n",TOPO_PONDERATION_RECALAGE_SYMETRIQUE);
printf("lambda_ref = %f \n",lambda_ref);

wdth=_ParamRecalageBspline.serie_groupwise[0]->width;hght=_ParamRecalageBspline.serie_groupwise[0]->height;dpth=_ParamRecalageBspline.serie_groupwise[0]->depth;


if (save_type==3) /* il faut transformer les images en puissance de 2 */
    {
    wdth_old=wdth; hght_old=hght; dpth_old=dpth;
    
    imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.ref_groupwise);
    imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.ref2_groupwise);
    imx_realloc_mri_pow2_centered_p(_ParamRecalageBspline.reca_groupwise);
    
    
    wdth=_ParamRecalageBspline.ref_groupwise->width;hght=_ParamRecalageBspline.ref_groupwise->height;dpth=_ParamRecalageBspline.ref_groupwise->depth;
    }
    
tdebut=time(NULL);
nbmax_param=3*pow((pow(2,resolf)-1),3);
 
/* Allocation */
if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
init_log(NULL);
aff_log("\nRECALAGE ");
 
/*************mise en place des parametres du recalage******/
scal_func=topo_choose_bspline(func_type);
interpol=imx_choose_interpolation_fct(inter_type);  
minimisation=topo_choose_optimisation(min_type);
normalisation=topo_choose_normalisation(normalisation_type,nb_classe);
regularisation=topo_choose_regularisation(reg_type);


if (dist_type==0)
distance=Energie_groupwise_variance_locale_nonsym_3d;
else
distance=Energie_groupwise_variance_locale_3d;

    
aff_log("\n");
aff_log("recalage Groupwise de l'image %d  sur %d \n",_ParamRecalageBspline.nb_reca_groupwise,_ParamRecalageBspline.nb_tot_groupwise);
 
    

/*allocation memoire de champ*/
 champ=cr_field3d(wdth,hght,dpth);
     
imtres=cr_grphic3d(_ParamRecalageBspline.reca_groupwise);imx_copie_param_3d_p(_ParamRecalageBspline.reca_groupwise,imtres);
imrefd=cr_grphic3d(_ParamRecalageBspline.ref_groupwise);imx_copie_3d_p(_ParamRecalageBspline.ref_groupwise,imrefd);
imrecad=cr_grphic3d(_ParamRecalageBspline.reca_groupwise);imx_copie_3d_p(_ParamRecalageBspline.reca_groupwise,imrecad);

 /********************RECALAGE ONDELETTE******************************************/
 {
int l,resol,r,D;
int ***masque_param;
int *x0,*x1,*y0,*y1,*z0,*z1;
double lambda=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
transf3d *transfo_ini;


nb_param=init_base_3d(wdth,hght,dpth,scal_func);
resol=BASE3D.resol;

/* Initialisation tableau */
for (l=0;l<nb_param;l++)
   {param[l]=0.0;}


if (nomfichiertrf!=NULL)
{
transfo_ini=load_transf_3d(nomfichiertrf);
    

if ((transfo_ini->width!=wdth)&&(transfo_ini->height!=hght)&&(transfo_ini->depth!=dpth))
    printf("La transfo ini n'a pas la bonne taille .... \n");
    
if ((transfo_ini->dx!=_ParamRecalageBspline.reca_groupwise->dx)&&(transfo_ini->dy!=_ParamRecalageBspline.reca_groupwise->dy)&&(transfo_ini->dz!=_ParamRecalageBspline.reca_groupwise->dz))
    printf("La transfo ini n'a pas les bons dx, dy, dz .... \n");


if (transfo_ini->typetrans!=BSPLINE3D)
    {
    printf("Il faudrait projeter sur une base de fonctions Bspline d'abord .... ! A FAIRE .... \n");
    return(0);
    
    }


if (transfo_ini->resol>resolf)
    {
    printf("La resolution souhaitee est plus faible que la resolution du champ passe en initialisation .... ! A FAIRE .... \n");
    return(0);
    
    }

resol= transfo_ini->resol;

for (l=1; l<resol; l++) 
  nb_param=base_resol_up_3d(param,nb_param);

         
        for (l=0;l<nb_param/3;l++) 
        {
        param[3*l]=transfo_ini->param[3*l]/_ParamRecalageBspline.reca_groupwise->dx;
        param[3*l+1]=transfo_ini->param[3*l+1]/_ParamRecalageBspline.reca_groupwise->dy;
        param[3*l+2]=transfo_ini->param[3*l+2]/_ParamRecalageBspline.reca_groupwise->dz;
        }
   
}

        

    /* nom du fichier dans lequel on enregistre la transfo oppose */
    if (nomfichres!=NULL)
        {
        
        strcpy(nomfichier_ref,nomfichres);
        ext2=strstr(nomfichier_ref,".trf");
        if (ext2!=NULL) *ext2='\0';

        strcat(nomfichier_ref,"_ref.trf");


        strcpy(nomfichier_reca,nomfichres);
        ext2=strstr(nomfichier_reca,".trf");
        if (ext2!=NULL) *ext2='\0';
        strcat(nomfichier_reca,"_reca.trf");

     }
     



if (nomfichres!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    strcpy(nomfichier_nrj,nomfichres);
    ext=strstr(nomfichier_nrj,".trf");
    if (ext!=NULL) *ext='\0';
    strcat(nomfichier_nrj,"_energy.txt");
    LOG_ENERGY=init_flog(nomfichier_nrj);      
    #endif
    }
     
/* Debut traitement */
for (r=resol;r<=resolf;r++)
    {
    TOPO_SEUIL_SIGNIFICATIVITE=0.05;

    x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
    D=TOP_D(nb_param);
    
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);

    base_to_field_3d(nb_param,param,champ,NULL,NULL);
    interpol(_ParamRecalageBspline.reca_groupwise,champ,imrecad); 

    
    for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
            {
            champ->raw[i][j][k].x=-1.0*lambda_ref*champ->raw[i][j][k].x;
            champ->raw[i][j][k].y=-1.0*lambda_ref*champ->raw[i][j][k].y;
            champ->raw[i][j][k].z=-1.0*lambda_ref*champ->raw[i][j][k].z;
            }
    
    interpol(_ParamRecalageBspline.ref_groupwise,champ,imrefd); 

    
    /* Creation du masque binaire r��encant les Slpqr sur lesquel l'image ne contient psa d'info */
    masque_param=alloc_imatrix_3d(D,D,D);
    //topo_masque_param(imrefd,NULL,imrecad,masque_param, nb_param); //a verifier si c'est pertinent
    
    for (i=0;i<D;i++)
    for (j=0;j<D;j++)
    for (k=0;k<D;k++)
        masque_param[i][j][k]=1;



if ((x1[0]-x0[0])<SIZE_MIN_REGULARISATION)
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
else
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;


if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.reca_groupwise,imtres,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);


mini=minimisation(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.reca_groupwise,imtres,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);


 
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }

        

    /*affichage du temps de calcul partiel*/
    tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
  
free_imatrix_3d(masque_param);
}       

    /*calcul de la transformation*/
    base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
 /*liberation de la base*/
   end_base_3d();

}
  /*enregistrement du champ resultat dans un fichier*/
    
  
  // transfo direct
  transfo=cr_transf3d(wdth,hght,dpth,NULL);
  transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
  transfo->typetrans=BSPLINE3D;
  transfo->resol=resolf;transfo->degre=func_type;
  transfo->dx=_ParamRecalageBspline.reca_groupwise->dx;
  transfo->dy=_ParamRecalageBspline.reca_groupwise->dy;
  transfo->dz=_ParamRecalageBspline.reca_groupwise->dz;


    // pour empecher des depalcement superieur a la taille de l'image 
     for (l=0;l<nb_param/3;l++) 
        {
        param[3*l]=MINI(param[3*l],wdth);
        param[3*l+1]=MINI(param[3*l+1],hght);
        param[3*l+2]=MINI(param[3*l+2],dpth);
        param[3*l]=MAXI(param[3*l],-1.0*wdth);
        param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
        param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
        } 
         
        for (l=0;l<nb_param/3;l++) 
        {
        transfo->param[3*l]=param[3*l]*_ParamRecalageBspline.reca_groupwise->dx;
        transfo->param[3*l+1]=param[3*l+1]*_ParamRecalageBspline.reca_groupwise->dy;
        transfo->param[3*l+2]=param[3*l+2]*_ParamRecalageBspline.reca_groupwise->dz;
        }
    
    
  // transfo opposee
  transfo_opp=cr_transf3d(wdth,hght,dpth,NULL);
  transfo_opp->nb_param=nb_param;transfo_opp->param=CALLOC(nb_param,double);
  transfo_opp->typetrans=BSPLINE3D;
  transfo_opp->resol=resolf;transfo_opp->degre=func_type;
  transfo_opp->dx=_ParamRecalageBspline.reca_groupwise->dx;
  transfo_opp->dy=_ParamRecalageBspline.reca_groupwise->dy;
  transfo_opp->dz=_ParamRecalageBspline.reca_groupwise->dz;
        
    for (l=0;l<nb_param;l++)
     transfo_opp->param[l]=-1.0*lambda_ref*transfo->param[l];
             

  // transfo inverse
 transfo_inv=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);
 inverse_transf_3d_p(transfo_opp,transfo_inv,0.01);
                
  // transfo composee
 transfo_comb=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);
 comb_transf_3d(transfo,transfo_inv,transfo_comb,5);
        
    
  
transf_to_field_3d_noalloc(transfo_comb, champ, _ParamRecalageBspline.ref_groupwise, _ParamRecalageBspline.reca_groupwise);
 
 




    if (save_type==3) // il faut transformer les images  dans leur taille d'origine
    {
        imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.reca_groupwise,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.ref_groupwise,wdth_old,hght_old,dpth_old);
        imx_undo_mri_pow2_centered_p(_ParamRecalageBspline.ref2_groupwise,wdth_old,hght_old,dpth_old);
    }
    
  
  // Enregistrement de la transfo resultat
  if (nomfichres!=NULL)
    {
    if (save_type==3) // on tronque le champ 
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)

                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                
                } 
        
        
        
         transfo_tmp=field_to_transf_3d(champ,_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.reca_groupwise); 
    
        transfo_tmp->dx=_ParamRecalageBspline.reca_groupwise->dx;
        transfo_tmp->dy=_ParamRecalageBspline.reca_groupwise->dy;
        transfo_tmp->dz=_ParamRecalageBspline.reca_groupwise->dz;

        if (dist_type!=0)
        save_transf_3d(transfo_tmp,nomfichres);
    
        free_transf3d(transfo_tmp);


                champ->width=   wdth;
                champ->height=  hght;
                champ->depth=   dpth;
                
                for (i=0;i<wdth;i++)
                for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                    {champ->raw[i][j][k].x=champ->raw[i][j][k].y=champ->raw[i][j][k].z=0.0;}
                

                _ParamRecalageBspline.ref_groupwise->width= wdth;
                _ParamRecalageBspline.ref_groupwise->height=    hght;
                _ParamRecalageBspline.ref_groupwise->depth= dpth;
            
                _ParamRecalageBspline.reca_groupwise->width=    wdth;
                _ParamRecalageBspline.reca_groupwise->height=   hght;
                _ParamRecalageBspline.reca_groupwise->depth=    dpth;
            
    /* enreristrement de la transformation de imreca */ 
     
        if (save_type==3) // on tronque le champ 
                {
                transf_to_field_3d_noalloc(transfo, champ, _ParamRecalageBspline.ref_groupwise, _ParamRecalageBspline.reca_groupwise);
                
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                _ParamRecalageBspline.ref_groupwise->width= wdth_old;
                _ParamRecalageBspline.ref_groupwise->height=    hght_old;
                _ParamRecalageBspline.ref_groupwise->depth= dpth_old;
            
                _ParamRecalageBspline.reca_groupwise->width=    wdth_old;
                _ParamRecalageBspline.reca_groupwise->height=   hght_old;
                _ParamRecalageBspline.reca_groupwise->depth=    dpth_old;
                
                transfo_tmp=field_to_transf_3d(champ,_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.reca_groupwise); 
                
                transfo_tmp->dx=_ParamRecalageBspline.reca_groupwise->dx;
                transfo_tmp->dy=_ParamRecalageBspline.reca_groupwise->dy;
                transfo_tmp->dz=_ParamRecalageBspline.reca_groupwise->dz;

                if (dist_type!=0)
                    save_transf_3d(transfo_tmp,nomfichier_reca);
                else
                    save_transf_3d(transfo_tmp,nomfichres);
                
                
                free_transf3d(transfo_tmp);
                }
            else
                {if (dist_type!=0)
                    save_transf_3d(transfo,nomfichier_reca); 
                else
                    save_transf_3d(transfo,nomfichres);
                }
        
                champ->width=   wdth;
                champ->height=  hght;
                champ->depth=   dpth;
                
                _ParamRecalageBspline.ref_groupwise->width= wdth;
                _ParamRecalageBspline.ref_groupwise->height=    hght;
                _ParamRecalageBspline.ref_groupwise->depth= dpth;
            
                _ParamRecalageBspline.reca_groupwise->width=    wdth;
                _ParamRecalageBspline.reca_groupwise->height=   hght;
                _ParamRecalageBspline.reca_groupwise->depth=    dpth;
        
                for (i=0;i<wdth;i++)
                for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                    {champ->raw[i][j][k].x=champ->raw[i][j][k].y=champ->raw[i][j][k].z=0.0;}
                

    /* enreristrement de la transformation de imreca */ 

        if (save_type==3) // on tronque le champ 
                {
                transf_to_field_3d_noalloc(transfo_opp, champ, _ParamRecalageBspline.ref_groupwise, _ParamRecalageBspline.reca_groupwise);
            
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                    _ParamRecalageBspline.ref_groupwise->width= wdth_old;
                _ParamRecalageBspline.ref_groupwise->height=    hght_old;
                _ParamRecalageBspline.ref_groupwise->depth= dpth_old;
            
                _ParamRecalageBspline.reca_groupwise->width=    wdth_old;
                _ParamRecalageBspline.reca_groupwise->height=   hght_old;
                _ParamRecalageBspline.reca_groupwise->depth=    dpth_old;
                
                transfo_tmp=field_to_transf_3d(champ,_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.reca_groupwise); 
                
                transfo_tmp->dx=_ParamRecalageBspline.reca_groupwise->dx;
                transfo_tmp->dy=_ParamRecalageBspline.reca_groupwise->dy;
                transfo_tmp->dz=_ParamRecalageBspline.reca_groupwise->dz;

                
                if (dist_type!=0)
                    save_transf_3d(transfo_tmp,nomfichier_ref);
                    
                free_transf3d(transfo_tmp);
                }
            else 
                if (dist_type!=0)
                    save_transf_3d(transfo_opp,nomfichier_ref); 
 
    }   
        
    

free_transf3d(transfo);
free_transf3d(transfo_inv);
free_transf3d(transfo_comb);
        
        
  /*liberation de la base*/
  // end_base_3d();

  /*caracteristique du champ*/
    aff_log("CHAMP: ");carac_field_3d(champ);

 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);
  /*liberation memoire des images temporaires*/

   free_grphic3d(imtres);
     free_grphic3d(imrefd);
     free_grphic3d(imrecad);
     

 


 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 if (nomfichres!=NULL)
    {   
    #ifndef TOPO_COMPARE_CRITERE 
     end_flog(LOG_ENERGY);
    #endif
    }
 free(param);free(param_norm);
 
 return(1);  
}




/*******************************************************************************
**        matching_Bspline_topo_momentgeometrique_3d                                                
**                                                                        
*******************************************************************************/

void matching_Bspline_topo_momentgeometrique_3d()
{
 int im_ref,im_reca,im_res;
 char *quest[3],*nomfichres;
 int dist_type,inter_type,min_type,save_type,func_type,i;
 int reg_type;
 double Jmin,Jmax;
 int resolf, e;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref= GET_PLACE3D(TEXT0031);
 im_res= GET_PLACE3D(TEXT0006);

 for(i=0;i<3;i++)
    quest[i]=CALLOC(80,char);
 
 /*type de fonction*/
  func_type=1;   // 1 = fonction BSPLINE de degr� 1
  
 /*distance*/
  dist_type=0;   // 0 = norme L2 
  
 /*interpolation*/
 inter_type=0;   // 0 = plus proche voisin

 /*minimisation*/
 min_type=2;    // 2 = Levenberg-Marquard
  
 /*question sur l'enregistrement du champ resultat*/
 nomfichres=quest_save_result_3d(&save_type);


 /*resolution finale*/
 resolf= GET_INT("resolution", 4, &e);

    
/*bornes sur le jacobien*/
Jmin=0;
Jmax=100000;


/* regularisation */
  strcpy(quest[0],"Pas de regularisation");
  strcpy(quest[1],"Regularisation membrane elastique");
  strcpy(quest[2],"\0");
        


reg_type= GETV_QCM("Regularisation ?",(const char **)quest);

    
TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
    
if (reg_type==1)
    {TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = GET_DOUBLE("facteur lambda",1,&e);
    reg_type=2; // c'est pour la compatible avec la fonction topo_choose_regularisation
    }

for(i=0;i<3;i++)
    free(quest[i]);
        
    imx_matching_Bspline_topo_momentgeometrique_3d(im_ref,im_reca,im_res,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax);   
    
 
 if (nomfichres)
  free(nomfichres);

  
 show_picture_3d(im_res);
 
}


/*******************************************************************************
**        imx_matching_Bspline_topo_momentgeometrique_3d                                             
**                                                                       
*******************************************************************************/
int imx_matching_Bspline_topo_momentgeometrique_3d(int im_ref, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax)
{
int l,continu;
grphic3d *imref,*imreca,*imres;
imref=ptr_img_3d(im_ref);
imreca=ptr_img_3d(im_reca);
imres=ptr_img_3d(im_res);

 if ((imref->dx!=imreca->dx)||(imref->dy!=imreca->dy)||(imref->dz!=imreca->dz))
    {
      PUT_WARN("Les images n'ont pas les memes dx,dy,dz !");
      printf("Les images n'ont pas les memes dx,dy,dz !\n");
      
            return(-1) ;
    }

 if ((imref->width!=imreca->width)||(imref->height!=imreca->height)||(imref->depth!=imreca->depth))
    {
      PUT_WARN("Les images n'ont pas la meme taille !");
      printf("Les images n'ont pas la meme taille !\n");
      return(-1) ;
    }


continu=0;
for (l=0;l<10;l++)
    {
    if (imref->width==pow(2.0,l)) continu++;
    if (imref->height==pow(2.0,l)) continu++;
    if (imref->depth==pow(2.0,l)) continu++;
    if (imreca->width==pow(2.0,l)) continu++;
    if (imreca->height==pow(2.0,l)) continu++;
    if (imreca->depth==pow(2.0,l)) continu++;
    }

if (continu!=6)
    {
      PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
          printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
            //return(-1) ;
            save_type=3;
   }

 imx_matching_Bspline_topo_momentgeometrique_3d_p(imref,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax);
 
  return(1);
}


/*******************************************************************************
**        imx_matching_Bspline_topo_momentgeometrique_3d_p                                            
**                                                                        
*******************************************************************************/

int imx_matching_Bspline_topo_momentgeometrique_3d_p(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax)
{
 int    wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,j,k,ideb,jdeb,kdeb,e;
 double *param,*param_norm;
 grphic3d *imtref,*imtreca,*imtres, *imrefmome, *imresmome, *imrecamome, *imtresmome;
 field3d *champ,*champinv;
 InterpolationFct interpol;
 dist_func_locale_t distance;
 reg_func_locale_t regularisation;
 min_func_locale_t minimisation;
 scal_func_t scal_func;
 double  mini;
 int tdebut,tfin;
 int nbmax_param;
 transf3d *transfinv; 
 bool regu = FALSE;
 int num=1;
 float rayon=3;
 int x,y,z;

 
wdth=imref->width;hght=imref->height;dpth=imref->depth;
    

//------- Pour gerer le cas ou les images ne sont pas en puissance de 2
if (save_type==3) /* il faut transformer les images en puissance de 2 */
    {
    wdth_old=wdth; hght_old=hght; dpth_old=dpth;

    imx_realloc_mri_pow2_centered_p(imref);
    imx_realloc_mri_pow2_centered_p(imreca);
    imx_realloc_mri_pow2_centered_p(imres);

    wdth=imref->width;hght=imref->height;dpth=imref->depth; 
    }
    

tdebut=time(NULL);
nbmax_param=3*pow((pow(2,resolf)-1),3);
 
/* Allocation */
if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }
if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); }

    
 init_log(NULL);
 aff_log("\nRECALAGE ");



 /*************mise en place des parametres du recalage******/
scal_func=topo_choose_bspline(func_type);
distance=topo_choose_distance(dist_type);
interpol=imx_choose_interpolation_fct(inter_type);  
minimisation=topo_choose_optimisation(min_type);
regularisation=topo_choose_regularisation(reg_type);
    
aff_log("\n");
aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);


//------- creation des images temporaires
imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
imx_norm_rcoeff_3d_p(imtref, imtreca);
imtres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtres);
//------- creation des images des moments
imrefmome=cr_grphic3d(imref);imx_copie_param_3d_p(imref,imrefmome);
imrecamome=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imrecamome);
imresmome=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imresmome);
imtresmome=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtresmome);
imx_norm_rcoeff_3d_p(imrefmome, imrecamome);
//------- calcul du moment de l'image de r�f�rence
imx_ComputeRotationInvariantMoment3D_p(imref, imrefmome, num, rayon);   

/*allocation memoire des champs*/
champ=cr_field3d(wdth,hght,dpth);
champinv=cr_field3d(wdth,hght,dpth);

 /********************RECALAGE ONDELETTE******************************************/
 {
  int l,resol,r,D,i;
  int ***masque_param;
  int *x0,*x1,*y0,*y1,*z0,*z1;
  nb_param=init_base_3d(wdth,hght,dpth,scal_func);
  resol=BASE3D.resol;
        
/* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}

/* Debut traitement : BOUCLE SUR LES DIFFERENTES ECHELLES  */
for (r=resol;r<=resolf;r++)    
{
    regu = FALSE;
    /* Debut traitement : BOUCLE au niveau des moments  */
    for (i=1;i<=5;i++)
    {
        //----------- a)Appliquer le champs sur imreca pour touver imtres dans l'espace de imref
        x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
        D=TOP_D(nb_param);
        aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);
        base_to_field_3d(nb_param,param,champ,NULL,NULL);
        interpol(imreca,champ,imtres); 
        //----------- b)calculer le moment g�om�trique de imtres
        imx_ComputeRotationInvariantMoment3D_p(imtres, imtresmome, num, rayon);
            
        //Pr�paration de transfo � partir du param 
            transfinv=cr_transf3d(wdth,hght,dpth,NULL);// a chaque iteration je fais cette creation et je la vide : voir apr�s ??????BON
        transfinv->nb_param=nb_param;transfinv->param=CALLOC(nb_param,double);
            transfinv->typetrans=BSPLINE3D;
            transfinv->resol=r;transfinv->degre=func_type;
        transfinv->dx=imtres->dx;
        transfinv->dy=imtres->dy;
        transfinv->dz=imtres->dz;
 
         /* pour empecher des depalcement superieur a la taille de l'image */ //cette boucle change le param?????????????????BON
         for (l=0;l<nb_param/3;l++) 
            {
            param[3*l]=MINI(param[3*l],wdth);
            param[3*l+1]=MINI(param[3*l+1],hght);
            param[3*l+2]=MINI(param[3*l+2],dpth);
            param[3*l]=MAXI(param[3*l],-1.0*wdth);
            param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
            param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
            } 

        for (l=0;l<nb_param/3;l++) 
            {
            transfinv->param[3*l]=param[3*l]*transfinv->dx;
            transfinv->param[3*l+1]=param[3*l+1]*transfinv->dy;
            transfinv->param[3*l+2]=param[3*l+2]*transfinv->dz;
            } 
        //----------- c) L inversion de la transform�e correspondant au champs
        e=inv_bspline_3d(transfinv, champinv,0.001);
         
                        
        // conversion du champ en voxel
        
        
        for (x=0;x<wdth; x++)
        for (y=0;y<hght; y++)
        for (z=0;z<dpth; z++)
            {
            champinv->raw[x][y][z].x=champinv->raw[x][y][z].x/transfinv->dx;
            champinv->raw[x][y][z].y=champinv->raw[x][y][z].y/transfinv->dy;
            champinv->raw[x][y][z].z=champinv->raw[x][y][z].z/transfinv->dz;
            }
        
        
        free_transf3d(transfinv);
        
        //----------- d) appliquer Le champs inverse au moment imresmome
        inter_linear_3d(imtresmome,champinv,imrecamome); 

        //save_mri_ipb_3d_p("/home/miv/noblet/medimax/imtresmome.ipb",imtresmome);
        
        /* Creation du masque binaire referencant les Slpqr sur lesquel l'image ne contient psa d'info */
        masque_param = alloc_imatrix_3d(D,D,D);// je le cr�e chaque it�ration???????????????????????????????????????????
        topo_masque_param (imrefmome,imrecamome,imresmome,masque_param, nb_param);
    
                    
        if (imrefmome->rcoeff!=imrecamome->rcoeff)
            printf("Attention les rcoeff ne sont pas identiques !!! Le recalage ne sera pas pertinent ....\n");

        //recherche du param�tre de r�gularisation convenable
        if (regu == FALSE)
        {
        regu = TRUE;
        
        
        if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE != 0)
        update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2           (imrefmome,imrecamome,imresmome,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);
        }

        //-----------e)  Optimisation des parametres du champ de deformation entre les moments imrecamome et imrefmome
        //save_mri_ipb_3d_p("/home/miv/noblet/medimax/imrefmome.ipb",imrefmome);
        //save_mri_ipb_3d_p("/home/miv/noblet/medimax/imrecamome.ipb",imrecamome);

        mini=minimisation(imrefmome,imrecamome,imresmome,champ,base_to_field_3d,interpol,distance, regularisation,      nb_param,param,masque_param,Jmin,Jmax,nomfichres);
    }
    /* FIN traitement : BOUCLE au niveau des moments  */




    if (r!=resolf)
    {/*on passe a la resolution superieure*/
        nb_param=base_resol_up_3d(param,nb_param);
        }   

        /*affichage du temps de calcul partiel*/
        tfin=time(NULL);
        {
         int h,m,s,ttotal;
         ttotal=tfin-tdebut;
         h=ttotal/3600;
         m=(ttotal-h*3600)/60;
         s=(ttotal-h*3600-m*60);
         aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
        }
  
    free_imatrix_3d(masque_param);

}
/* FIN traitement : FIN DE BOUCLE SUR LES DIFFERENTES ECHELLES  */    

    

/*calcul de la transformation*/
base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
 if (imres != NULL)
    { 
    interpol(imreca,champ,imtres);
        imx_copie_3d_p(imtres,imres);
    }

  

 
 //------- Pour gerer le cas ou les images ne sont pas en puissance de 2
 if (save_type==3) /* il faut transformer les images  dans leur taille d'origine*/
    {
    imx_undo_mri_pow2_centered_p(imref,wdth_old,hght_old,dpth_old);
    imx_undo_mri_pow2_centered_p(imreca,wdth_old,hght_old,dpth_old);
    imx_undo_mri_pow2_centered_p(imres,wdth_old,hght_old,dpth_old);
    }
  

  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 
    transf3d *transfo; 
    
    if ((save_type==0)||(save_type==3))
        {
            if (save_type==3) /* on tronque le champ */
                {
                ideb=(int)((wdth-wdth_old)*0.5);
                jdeb=(int)((hght-hght_old)*0.5);
                kdeb=(int)((dpth-dpth_old)*0.5);

                for (i=0;i<wdth_old;i++)
                for (j=0;j<hght_old;j++)
                for (k=0;k<dpth_old;k++)
                    {
                    champ->raw[i][j][k]=champ->raw[i+ideb][j+jdeb][k+kdeb];
                    }
                
                champ->width=   wdth_old;
                champ->height=  hght_old;
                champ->depth=   dpth_old;
                
                } 
        
        
         transfo=field_to_transf_3d(champ,imref,imreca); 
    }
        else
    {  
     /////////Ce code fait le passage d un param a un transfo
     transfo=cr_transf3d(wdth,hght,dpth,NULL);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=func_type;
         transfo->dx=imres->dx;
         transfo->dy=imres->dy;
         transfo->dz=imres->dz;
         
         /* pour empecher des depalcement superieur a la taille de l'image */
         for (l=0;l<nb_param/3;l++) 
            {
            param[3*l]=MINI(param[3*l],wdth);
            param[3*l+1]=MINI(param[3*l+1],hght);
            param[3*l+2]=MINI(param[3*l+2],dpth);
            param[3*l]=MAXI(param[3*l],-1.0*wdth);
            param[3*l+1]=MAXI(param[3*l+1],-1.0*hght);
            param[3*l+2]=MAXI(param[3*l+2],-1.0*dpth);
            } 
         
     for (l=0;l<nb_param/3;l++) 
            {
            transfo->param[3*l]=param[3*l]*imreca->dx;
            transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
            transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
            } 
    }
    save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 


 }
 /********************FIN RECALAGE ONDELETTE******************************************/


  /*liberation de la base*/
   end_base_3d();

  /*caracteristique du champ*/
    aff_log("CHAMP: ");carac_field_3d(champ);

 /*liberation memoire des variables dynamiques*/ 
  /*liberation memoire du champ*/
   free_field3d(champ);// et l autre inv??????????????????????
  /*liberation memoire des images temporaires*/
    free_grphic3d(imtref);free_grphic3d(imtreca); 
    free_grphic3d(imtres);
//free les autres???????????????????????????????????????
 


 
 /*affichage du temps de calcul total*/
 tfin=time(NULL);
 {
  int h,m,s,ttotal;
  ttotal=tfin-tdebut;
  h=ttotal/3600;
  m=(ttotal-h*3600)/60;
  s=(ttotal-h*3600-m*60);
  aff_log("TEMPS TOTAL: %dh %dm %ds \n",h,m,s);
 }
 end_log();
 
 free(param);free(param_norm);




 return(1);  
}



