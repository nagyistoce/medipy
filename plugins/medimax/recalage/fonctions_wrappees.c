/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "mtch_3d.h"
#include "chps_3d.h"
#include "noyau/io/imx_log.h"
#include "outils/imx_misc.h"
#include "fonctions_wrappees.h"
#include "noyau/gui/imx_picture3d.h"
#include "distance_3d.h"
#include "recalage/topo/mtch_topo_3d.h"


void SetGrphic3dDxDyDzPositive(grphic3d *im)
{
if ((im->dx < 0) || (im->dy < 0) || (im->dz < 0 ))
     {
         im->dx=fabs(im->dx);
         im->dy=fabs(im->dy);
         im->dz=fabs(im->dz);
         printf("Attention, les dx,dy,dz negatif sont transformes en positif pour respecter la compatibilite avec medimax\n");
    }
     
     
}

/*******************************************************************************
**  Linear registration
*******************************************************************************/
int LinearRegistration( grphic3d *imref, 
                        grphic3d *imreca, 
                        grphic3d *imres, 
                        int registration_type, 
                        int dist_type, 
                        int inter_type, 
                        int min_type, 
                        int multistart, 
                        int save_type, 
                        char *nomfichres, 
                        int precision,
                        int start_resol,
                        int end_resol
                        )
{
    

transf3d *inittrf=NULL;
transf3d *transfres=NULL;

char *nomfichtrf = NULL; // pas de transfo ini                                  


eMatchPrecision matchPrecision;

//intervalle de recherche
eResearchInterv FieldOfResearch = MEDIUM;

SetGrphic3dDxDyDzPositive(imref);
SetGrphic3dDxDyDzPositive(imreca);


if (save_type==2)
    nomfichres=NULL;

if (strlen(nomfichres)<1)
    {
    save_type=2;
    nomfichres=NULL;
    }
switch(precision)
    {
    case 0:
        matchPrecision=NORMAL;
        break;
    case 1:
        matchPrecision=PRECIS;
        break;
    case 2:
        matchPrecision=TRES_PRECIS;
        break;
    case 3:
        matchPrecision=PRECISION_MAX;
        break;
    default:        
        matchPrecision=NORMAL;
        break;
    }





   if (nomfichtrf !=NULL) // initialisation de la transfo par un fichier.
     inittrf=load_transf_3d(nomfichtrf);
   transfres=cr_transf3d_p(imref,NULL);

  imx_inimaxminpixel_3d_p(imref);
  imx_inimaxminpixel_3d_p(imreca);
  if ((imref->min_pixel<0)||(imreca->min_pixel<0))
  { PUT_WARN("L'image contient des valeurs negatives!!\nOn ajoute -val_min pour la recaler\n"); }

   init_log(NULL);
   aff_log("\n");
   aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);

   switch(registration_type)
    {
    case 0: // rigid
        
        if ((dist_type==11) || (dist_type==12)) 
            {
            aff_log("\n RECALAGE RIGIDE ICP");
            dist_type=dist_type-11;
            imx_matching_rigide_ICP_3d_p(imref, imreca, imres, dist_type, inter_type, min_type, save_type, transfres,inittrf, FieldOfResearch, matchPrecision);
            break;
            }
        
        if (multistart==1)
        {
         aff_log("\n RECALAGE RIGIDE EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
         imx_matching_rigide_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);
        }
        else
        {
        aff_log("\n RECALAGE RIGIDE SANS ZOOM EN MULTIRESOLUTION");
        imx_matching_rigide_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision,start_resol,end_resol);
        }       
        break;
    case 1: // rigid + zoom
         if ((dist_type==11) || (dist_type==12)) 
            {
            aff_log("\n RECALAGE RIGIDE ICP");
            dist_type=dist_type-11;
            imx_matching_rigide_zoom_ICP_3d_p(imref, imreca, imres, dist_type, inter_type, min_type, save_type, transfres,inittrf, FieldOfResearch, matchPrecision);
            break;
            }
 
        
        if (multistart==1)
        {
        aff_log("\n RECALAGE RIGIDE+ZOOM EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
        imx_matching_rigide_zoom_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);
        }       
        else
        {
        aff_log("\n RECALAGE RIGIDE AVEC ZOOM EN MULTIRESOLUTION");
        imx_matching_rigide_zoom_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision, start_resol, end_resol);
        }       
        break;
    case 2: // affine
        if ((dist_type==11) || (dist_type==12)) 
            {
            aff_log("\n RECALAGE RIGIDE ICP");
            dist_type=dist_type-11;
            imx_matching_affine_ICP_3d_p(imref, imreca, imres, dist_type, inter_type, min_type, save_type, transfres,inittrf, FieldOfResearch, matchPrecision);
            break;
            }
 
        
        if (multistart==1)
        {
        aff_log("\n RECALAGE AFFINE EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
        imx_matching_affine_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);
        }       
        else
        {
        aff_log("\n RECALAGE AFFINE EN MULTIRESOLUTION");
        imx_matching_affine_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision, start_resol, end_resol);
        }       
        break;
    default: // affine
        if (multistart==1)
        {
        aff_log("\n RECALAGE AFFINE EN MULTIRESOLUTION AVEC MULTISTART (methode Jenkison&Smith) ");
        imx_matching_affine_multistart_multires_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,FieldOfResearch,matchPrecision);
        }       
        else
        {
        aff_log("\n RECALAGE AFFINE EN MULTIRESOLUTION");
        imx_matching_affine_multi_3d_p(imref,imreca,imres,dist_type,inter_type,min_type,save_type,transfres,inittrf,matchPrecision, start_resol, end_resol);
        }       
        break;
    }
   
 
   /*enregistrement du champ resultat*/
   if (nomfichres!=NULL)
   { 
     save_transf_3d(transfres,nomfichres);
   } 

   if (inittrf)
     free_transf3d(inittrf);
   free_transf3d(transfres);

   end_log();
   return(1);
}


/*******************************************************************************
**  Apply a transf_3d
*******************************************************************************/

int ApplyTransfo3d(grphic3d *imdeb, char *nomfichier, grphic3d *imres, int inter_type)
{
 
    transf3d *transfo;
    char tmp[PATH_LEN];
    InterpolationFct inter_fct;
    int err=0;
      
    SetGrphic3dDxDyDzPositive(imdeb);
    put_file_extension(nomfichier,".trf",tmp);
  
    /*chargement de la transformation*/
    transfo = load_transf_3d(tmp);
    if (! transfo)
    {
      return(0);
    }
  
    inter_fct = imx_choose_interpolation_fct(inter_type);
  
    err = imx_apply_transf_3d_p(imres, imdeb, transfo, inter_fct);

    imx_reinitialise_visual_params_3d_p(imres);

    free_transf3d(transfo); 

    if (err==0)
    return(1);
    else 
    return (0);
    
}


/*******************************************************************************
**  Compute similarity measure between images
*******************************************************************************/
double SimilarityMeasure3d(grphic3d *im1, grphic3d *im2, int err_type)
{
    double res;
    dist_func_t  func_erreur;
    ptr_distance_param dist_param=CALLOC(1, distance_param);

    func_erreur = imx_choose_distance_fct(err_type);

    dist_param->imreca=im1;
    dist_param->imref=im2;
    res = func_erreur(dist_param);
    FREE(dist_param);
    return res;
}


/*******************************************************************************
**  Bspline Registration
*******************************************************************************/
int BsplineRegistration3d(  grphic3d *imref, 
                            grphic3d *imreca, 
                            grphic3d *imres, 
                            int dist_type,
                            int reg_type,
                            double reg_factor, 
                            int min_type,
                            int save_type,
                            char *nomfichres,
                            int resolf,
                            double Jmin,
                            double Jmax,
                            int normalisation_type,
                            int nb_classe,
                            int biais,
                            int symetrique)
{


int l,continu,func_type,inter_type,adaptatif;

SetGrphic3dDxDyDzPositive(imref);
SetGrphic3dDxDyDzPositive(imreca);

// test sur la sauvegarde du fichier
if (save_type==2) // pas d'enregistrement de transfo
    {
    nomfichres=NULL;
    }
else 
    {
    if (strlen(nomfichres)<1)
        {
        save_type=2;
        nomfichres=NULL;
        }
    }
 // type de fonction : pour le moment, ne fonctionne que pour des BSpline de degre 1
  func_type=1; 

// type d'interpolation pour la deformation finale de l'image (sinc3d)
    inter_type=4;

// test si la methode d'optimisation est compatible avec la fonction de cout choisie
if ((dist_type==2)||(dist_type==3)||(dist_type==6))
    {
    if ((min_type==2)||(min_type==3))
        {
        printf("The optimization method is not compatible with the cost function !\n");
        return(1);
        }
    }
    
// test sur la resolution
if ((resolf<1) || (resolf>7))
        {
        printf("The resolution should be set between 1 and 7 !\n");
        return(1);
        }

// test de la compatibilite sur les bornes du jacobien
if (Jmin>=1)
        {
        printf("Jmin should be lower than 1 !\n");
        return(1);
        }

if (Jmax<=1)
        {
        printf("Jmax should be higher than 1 !\n");
        return(1);
        }


// test sur la compatibilite du nombre de classe en fonction de la methode de normalisation d'intensite
if ((normalisation_type>5)&&(normalisation_type!=10)&&(normalisation_type!=16))
    {
    if (normalisation_type==18)
        printf("You chose to consider %d quantile to normalize histogramms \n",nb_classe);
    else
        if (normalisation_type==15)
            printf("You chose to consider %d bins for histogramm equalization\n",nb_classe);
        else
            printf("You chose to consider %d classes for intensity normalization \n",nb_classe);
    }
  

// On ne choisit pas le recalage adaptatif
adaptatif=1;

// parametre de la regularisation
TOPO_REGULARISATION_SIGMA_GAUSSIEN=0;
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;

if (reg_type==1)
    TOPO_REGULARISATION_SIGMA_GAUSSIEN = reg_factor;
    
if (reg_type>1)
    TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=reg_factor;


// test sur la regularisation
if (reg_type>7)
        {
        printf("You chose a wrong type of regularization !\n");
        return(1);
        }
        

// verification sur les caracteristiques de l'images
if ((imref->dx!=imreca->dx)||(imref->dy!=imreca->dy)||(imref->dz!=imreca->dz))
    {
      PUT_WARN("Images do not have the same dx,dy,dz !");
      printf("Images do not have the same dx,dy,dz ! \n");
      
            return(1) ;
    }

 if ((imref->width!=imreca->width)||(imref->height!=imreca->height)||(imref->depth!=imreca->depth))
    {
      PUT_WARN("Images do not have the same dimensions !");
      printf("Images do not have the same dimensions !\n");
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
      PUT_WARN("Images dimensions are not a power of 2 !!! Temporary modification of image dimensions. The result will be saved as a .chp file");
      printf("Images dimensions are not a power of 2 !!! Temporary modification of image dimensions. The result will be saved as a .chp file\n");
            save_type=3;
   }

// mise a jour des min-max de l'image car le wrapper ne met pas a jour cette info
imx_inimaxminpixel_3d_p(imref);
imx_inimaxminpixel_3d_p(imreca);
        
if (symetrique > 0)
    {
    TOPO_PONDERATION_RECALAGE_SYMETRIQUE=0.5;
    normalisation_type=2; //il y a encore un pb avec le melange de gaussienne
    imx_matching_Bspline_topo_symetrique_3d_p(imref,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax,normalisation_type, nb_classe,adaptatif);
    }
else
    {
    imx_matching_Bspline_topo_3d_p(imref,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax,normalisation_type, nb_classe,adaptatif,biais);    
    }


 

return(0);

}


/*******************************************************************************
**  Combine two trf file
*******************************************************************************/

int CombineTransfo3d(char *nomfichier1, char *nomfichier2, char *nomfichierres, int interpolation)
{
   
transf3d *transfo1,*transfo2,*transfores;
int wdth, hght, dpth;
  
/*lecture des champ*/
transfo1=load_transf_3d(nomfichier1);
transfo2=load_transf_3d(nomfichier2);

wdth=transfo1->width;hght=transfo1->height;dpth=transfo1->depth;

// si transfo pas bspline, alors la taille de transfo1 n'a pas d'importance
if (transfo1->typetrans==BSPLINE3D)
    {
    if (transfo2->width!=wdth || transfo2->height!=hght || transfo2->depth!=dpth)
       {
        PUT_ERROR("The two fields have not the same size");
        return(0);
       }
            
     }
transfores=cr_transf3d(transfo2->width,transfo2->height,transfo2->depth,NULL);
             
comb_transf_3d(transfo1,transfo2,transfores,interpolation);
     
/*enregistrement du champ resultat*/
save_transf_3d(transfores,nomfichierres);
   
/* Libere memoire allouee */
free_transf3d(transfo1); free_transf3d(transfo2);free_transf3d(transfores); 

return(1);
}

/*******************************************************************************
**  Invert a trf file                              
*******************************************************************************/
int InvertTransfo3d(char *nomfichier, char *nomfichres, int wdthref, int hghtref, int dpthref, float dxref, float dyref, float dzref, float prec)
{
char  nomfich1[500];
transf3d *transfo1,*transfores;
int wdth,hght,dpth,i, wdth_res, hght_res, dpth_res, e=1;
float dxres, dyres, dzres;
double *param;
field3d *chres;
 
 
put_file_extension(nomfichier,".trf",nomfich1);

if ((test_fichier_transf_3d(nomfich1))!=1) 
{PUT_ERROR("This file contains no 3D transformation");return(0);}


/*nom du fichier resultat*/
put_file_extension(nomfichres,".trf",nomfichres);

/*lecture des champ*/
transfo1=load_transf_3d(nomfich1);
        
wdth=transfo1->width;
hght=transfo1->height;
dpth=transfo1->depth;

if ((wdthref>0)&&(hghtref>0)&&(dpthref>0)&&(dxref>0)&&(dyref>0)&&(dzref>0))
    {
    wdth_res=wdthref;
    hght_res=hghtref;
    dpth_res=dpthref;

    dxres=dxref;
    dyres=dyref;
    dzres=dzref;
    }
else
    {
    wdth_res=transfo1->width;
    hght_res=transfo1->height;
    dpth_res=transfo1->depth;

    dxres=transfo1->dx;
    dyres=transfo1->dy;
    dzres=transfo1->dz;
    }


transfores=cr_transf3d(wdth_res,hght_res,dpth_res,NULL);
        
if ((transfo1->typetrans==RIGID3D)||    (transfo1->typetrans==RIGIDZOOM3D)||(transfo1->typetrans==AFFINE3D))
        {   
        /* copie des parametres */
        transfores->typetrans=AFFINE3D;
        transfores->dx=dxres;
        transfores->dy=dyres;
        transfores->dz=dzres;
        transfores->nb_param=15;
        transfores->param=CALLOC(transfores->nb_param,double);
        param=CALLOC(transfores->nb_param,double);
        
        for(i=0;i<transfo1->nb_param;i++)
            param[i]=transfo1->param[i];
        
        if (transfo1->typetrans==RIGID3D)
            {   rigid_to_rigidz_3d(param);rigidz_to_affine_3d(param);}
    
        if (transfo1->typetrans==RIGIDZOOM3D)
            {   rigidz_to_affine_3d(param);}
    
        transfo1->param=param;
                 
        e=inv_affine_3d(transfo1,transfores);
        }
else 
        {
        if (transfo1->typetrans==BSPLINE3D)
            {
            chres=cr_field3d(wdth,hght,dpth);
            e=inv_bspline_3d(transfo1,chres,prec);
            }
        else
            {
            // Je suis pas sur que cela fonctionne vraiment ....
            PUT_WARN("L'inversion de chp peut etre hasardeuse ...");
            chres=cr_field3d(wdth,hght,dpth);
            e=inv_field_3d(transfo1,chres,prec);
            }
    
        transfores=field_to_transf_3d(chres,NULL,NULL);
        transfores->typetrans=CHAMP3D;
        transfores->dx=dxres;
        transfores->dy=dyres;
        transfores->dz=dzres;
    
        free_field3d(chres);
        }   
        
/*enregistrement du champ resultat*/
save_transf_3d(transfores,nomfichres);

/* Libere memoire allouee */
free_transf3d(transfo1); free_transf3d(transfores);

return(e);    
}   

/*******************************************************************************
**  Load trf file                              
*******************************************************************************/
int LoadTrfFile(char *nomfichier, grphic3d *ux, grphic3d *uy,grphic3d *uz, grphic3d *imSource)
{
field3d *champ;
transf3d *transfo;
grphic3d imRes;
transfo=load_transf_3d(nomfichier);

imRes.width=transfo->width;
imRes.height=transfo->height;
imRes.depth=transfo->depth; 
imRes.dx=transfo->dx;
imRes.dy=transfo->dy;
imRes.dz=transfo->dz;


    
    
champ=transf_to_field_3d(transfo,&imRes,imSource);


// allocation des images
resize_grphic3d_buffer(ux, transfo->width,transfo->height,transfo->depth);
resize_grphic3d_buffer(uy, transfo->width,transfo->height,transfo->depth);
resize_grphic3d_buffer(uz, transfo->width,transfo->height,transfo->depth);


composante_field_to_image_3d(champ,ux,0);
composante_field_to_image_3d(champ,uy,1);
composante_field_to_image_3d(champ,uz,2);

// Mise a jour des dx,dy,dz de l'image resultat
ux->dx=transfo->dx;
ux->dy=transfo->dy;
ux->dz=transfo->dz;

uy->dx=transfo->dx;
uy->dy=transfo->dy;
uy->dz=transfo->dz;

uz->dx=transfo->dx;
uz->dy=transfo->dy;
uz->dz=transfo->dz; 

//printf("%d %d %d \n",ux->width, uy->height, uz->depth);
free_field3d(champ);
free_transf3d(transfo);

return (1);
}

/*******************************************************************************
**  Save trf file                              
*******************************************************************************/
int SaveTrfFile(char *nomfichier, grphic3d *ux, grphic3d *uy,grphic3d *uz, grphic3d *imSource)
{
field3d *champ;
transf3d *transfo;
int i,j,k;
grphic3d imRes;





if ((ux->width != uy->width)||(ux->width != uz->width)||(ux->height != uy->height)||(ux->height != uz->height)||(ux->depth != uy->depth)||(ux->depth != uz->depth))
     {
     PUT_ERROR("The three images do not have the same size");
     exit(-1);
     }


if ((ux->dx != uy->dx)||(ux->dx != uz->dx)||(ux->dy != uy->dy)||(ux->dy != uz->dy)||(ux->dz != uy->dz)||(ux->dz != uz->dz))
     {
     PUT_ERROR("The three images do not have the same dx, dy, dz");
     exit(-1);
     }
 
champ=cr_field3d(ux->width,ux->height,ux->depth);



// conversion en champs voxelique
for(i=0;i<ux->width;i++)
for(j=0;j<ux->height;j++)
for(k=0;k<ux->depth;k++)
    {
    champ->raw[i][j][k].x=ux->mri[i][j][k]*ux->rcoeff;//*ux->dx;
    champ->raw[i][j][k].y=uy->mri[i][j][k]*uy->rcoeff;//*ux->dy;
    champ->raw[i][j][k].z=uz->mri[i][j][k]*uz->rcoeff;//*ux->dz;     
    }

 
imRes.width=ux->width;
imRes.height=ux->height;
imRes.depth=ux->depth;  
imRes.dx=ux->dx;
imRes.dy=ux->dy;
imRes.dz=ux->dz;
   
transfo=field_to_transf_3d(champ,&imRes,imSource);

transfo->dx=ux->dx;
transfo->dy=ux->dy;
transfo->dz=ux->dz;

save_transf_3d(transfo,nomfichier);
free_transf3d(transfo);


free_field3d(champ);

return (1);
}
/*******************************************************************************
**  MRI info 3D
*******************************************************************************/
void MriInfo3D( grphic3d *p)
{
    imx_inimaxminpixel_3d_p(p);
    printf("Dimensions (pixels) : %dx%dx%d\n",  p->width, p->height, p->depth);
    printf("Bits par pixel : %d\n",  p->bitppixel);
    printf("Zoom : %d\n",  p->zoom);
    printf("Valeurs de MRI : min=%ld max=%ld\n",  p->min_pixel, p->max_pixel);
    printf("Coefficient de valeur reelle : rcoeff=%f icomp=%f\n",  p->rcoeff, p->icomp);
    printf("Valeurs de l'image : min=%f max=%f\n", p->min_pixel*p->rcoeff, p->max_pixel*p->rcoeff);
    printf("Resolution (mm) : dx=%f dy=%f dz=%f\n", p->dx, p->dy, p->dz);
} 


/*******************************************************************************
**  Visualisation of trf file                              
*******************************************************************************/
void VisuTrfFile(char *nomfichier, grphic3d *output, int type)
{
transf3d *transfo;
field3d *champ;

transfo=load_transf_3d(nomfichier);
champ=transf_to_field_3d(transfo,NULL,NULL);

resize_grphic3d_buffer(output, transfo->width,transfo->height,transfo->depth);


switch (type)
    {
    case 0 :/*module*/
        module_field_to_image_3d(champ,output);
        break;
    case 1 : /*jacobien*/
        jacobien_field_to_image_3d(champ,output);
        break;
    case 2 : /* composante suivant x*/
        composante_field_to_image_3d(champ,output,0);
        break;
    case 3 : /* composante suivant y*/
        composante_field_to_image_3d(champ,output,1);
        break;
    case 4 : /* composante suivant z*/
        composante_field_to_image_3d(champ,output,2);
        break;
    default: /* par defaut, on affiche le module */
        module_field_to_image_3d(champ,output);
        break;
    
    }           

free_field3d(champ);
free_transf3d(transfo);
 
}
