/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h> 	
#include <time.h>
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "recalage/chps_3d.h"
//#include "noyau/gui/imx_picture3d.h"
#include "recalage/reca_3d.h"
#include "morpho/morpho_3d.h"
#include "noyau/mani_3d.h"
#include "noyau/imx_3d.h"
#include "noyau/io/imx_log.h"
#include "outils/imx_list.h"
#include "outils/imx_misc.h"
#include "recalage/topo/recasegconjoint.h"
#include "math/imx_matrix.h"
#include "math/imx_bspline.h"
#include "noyau/imx_lang.h"
#include "math/oper_3d.h"
#include "traitement/trai_3d.h"
#include "noyau/mani_3d.h"


extern "C" 
{
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/divers_topo.h"
#include "recalage/topo/distance_topo.h"
#include "segmentation/otsu_3d.h"
}
#include<map>
#include<string>
#include<list>
#include<iostream>
#include<vector>
#include<set>
#include<algorithm>
#define SAVE_INTERMEDIAIRE_2
using namespace std;


int MINIMAL = 0;
int MAXIMAL = 0;
int RESOLUTION = 0;
int ***Proportion;
static double *histoOs = NULL;
static double *histoPartieMolle = NULL;
static grphic3d *imageSegmentee = NULL;
int tailleHisto;

void lancementrecaseg()
{
}		       

/*DEFINITION DES ALGORITHMES DE SEGMENTATION*/
//Algorithme de segmentation : seuillage !!!
//Cet algorithme determine la segmentation de l image ImageASegmenter a partir 
//de la connaissance de deux images de seuil (ImageSeuilHaut et ImageSeuilBas).
//Le resultat de la segmentation est mis dans l image ImageResultat

void Segmentation   (float ***ImageSeuilHaut     , //Image representant le seuil haut
                     float ***ImageSeuilBas      , //Image representant le seuil bas
                     grphic3d *ImageASegmenter    , //Image a segmenter
                     grphic3d *ImageResultat      //Image resultat
		     ) 
{
short ***ImSeuilHaut;
short ***ImSeuilBas;

ImSeuilHaut = alloc_matrix_3dSHORT(ImageASegmenter->width,ImageASegmenter->height,ImageASegmenter->depth);
ImSeuilBas = alloc_matrix_3dSHORT(ImageASegmenter->width,ImageASegmenter->height,ImageASegmenter->depth);

for (unsigned int i=0;i<ImageASegmenter->width;i++)
for (unsigned int j=0;j<ImageASegmenter->height;j++)
for (unsigned int k=0;k<ImageASegmenter->depth;k++)
	{
	ImSeuilHaut[i][j][k]=ImageSeuilHaut[i][j][k] ;
	}

for (unsigned int i=0;i<ImageASegmenter->width;i++)
for (unsigned int j=0;j<ImageASegmenter->height;j++)
for (unsigned int k=0;k<ImageASegmenter->depth;k++)
	{
	ImSeuilBas[i][j][k]=ImageSeuilBas[i][j][k] ;
	}

Points_3dInt pointBG;
Points_3dInt pointHD;
pointBG.x=pointBG.y=pointBG.z=0;
pointHD.x=pointHD.y=pointHD.z=255;


if (imageSegmentee == NULL)
	{
	
	imageSegmentee = cr_grphic3d(ImageResultat);
	for (unsigned int i=0;i<ImageASegmenter->width;i++)
	for (unsigned int j=0;j<ImageASegmenter->height;j++)
	for (unsigned int k=0;k<ImageASegmenter->depth;k++)
		imageSegmentee->mri[i][j][k]=0;
	
	segmentationNicolas(ImageASegmenter,ImSeuilHaut,ImSeuilBas,imageSegmentee);	
	
	for (unsigned int i=0;i<ImageASegmenter->width;i++)
	for (unsigned int j=0;j<ImageASegmenter->height;j++)
	for (unsigned int k=0;k<ImageASegmenter->depth;k++)
		{
	if (	(imageSegmentee->mri[i][j][k]== 1) || (imageSegmentee->mri[i][j][k]== 2) ) 
		ImageResultat->mri[i][j][k]= 1;
	else
		ImageResultat->mri[i][j][k]= 0;
		}

	}		
else
	{
	 segmentationNicolas2(ImSeuilHaut,ImSeuilBas,ImageASegmenter,imageSegmentee,ImageResultat,&pointBG,&pointHD,NULL,NULL);
	 for (unsigned int i=0;i<ImageASegmenter->width;i++)
	for (unsigned int j=0;j<ImageASegmenter->height;j++)
	for (unsigned int k=0;k<ImageASegmenter->depth;k++)
		{
		imageSegmentee->mri[i][j][k] = ImageResultat->mri[i][j][k];
		if (	(ImageResultat->mri[i][j][k]== 1) || (ImageResultat->mri[i][j][k]== 2) ) 
			ImageResultat->mri[i][j][k]= 1;
		else
			ImageResultat->mri[i][j][k]= 0;
		}
	}


separation_composantes(ImageASegmenter,ImSeuilHaut,ImSeuilBas,ImageResultat);
generation_tunnels(ImageASegmenter,ImSeuilHaut,ImSeuilBas,ImageResultat);

free_matrix_3dSHORT(ImSeuilHaut);
free_matrix_3dSHORT(ImSeuilBas);	
}	

/*LANCEMENT DE L ALGORITHME*/
//Cette fonction prend en parametre l image Template (ImageTemplate), l image a segmenter (ImageASegmenter) 
//Elle lance l algo d optimisation echelle par echelle. 
int imx_segmentation_p2(grphic3d *ImageTemplate , 
			grphic3d *ImageTemplateNonDecoupe , 
		       grphic3d *ImageASegmenter,
		       grphic3d *ImageResultat  ,
		       double   *param	        ,
		       int      nb_param        ,
		       double   seuilHaut       ,
		       double   seuilBas	,              
		       field3d *chres           ,
		       int optimisationEchelle0)
{
int bande;
double mini;
grphic3d *ImageTemplateOs,*ImageTemplateOsDilate,*temporaire;
scal_func_t scal_func;
InterpolationFct interpol2;
interpol2 = inter_labeled_3d;


for (unsigned int i=0;i<ImageTemplateNonDecoupe->width;i++)
for (unsigned int j=0;j<ImageTemplateNonDecoupe->height;j++)
for (unsigned int k=0;k<ImageTemplateNonDecoupe->depth;k++)
	{
	if (ImageTemplateNonDecoupe->mri[i][j][k])
		ImageTemplateNonDecoupe->mri[i][j][k]=1;	
	}
ImageTemplateNonDecoupe->max_pixel=1;

aff_log("\nSEGMENTATION RESOLUTION %d    nb param %d \n",BASE3D.resol,nb_param);

//creation des images temporaires segmentees
temporaire  	      = cr_grphic3d_modif(ImageResultat->width, ImageResultat->height, ImageResultat->depth, 0,1, 1); 
ImageTemplateOsDilate = cr_grphic3d_modif(ImageResultat->width, ImageResultat->height, ImageResultat->depth, 0,1, 1);
ImageTemplateOs	      = cr_grphic3d_modif(ImageResultat->width, ImageResultat->height, ImageResultat->depth, 0,1, 1); 

for (unsigned int i=0;i<ImageResultat->width;i++)
for (unsigned int j=0;j<ImageResultat->height;j++)
for (unsigned int k=0;k<ImageResultat->depth;k++)
	ImageResultat->mri[i][j][k]=ImageTemplateOs->mri[i][j][k] = 
	ImageTemplateOsDilate->mri[i][j][k] = temporaire->mri[i][j][k]=0;

//Calcul des images ou il est genant de mettre des 1, ou il est genant de mettre des 0
if (RESOLUTION == 0)
	bande = 20;
if (RESOLUTION == 1)
	bande = 15;
if (RESOLUTION == 2)
	bande = 10;	
if (RESOLUTION == 3)
	bande = 6;		
if (RESOLUTION == 4)
	bande = 4;	
if (RESOLUTION == 5)
	bande = 3;	
if (RESOLUTION == 6)
	bande = 2;		

//Calcul des images templates
char nomfichier_nrj[255];
if (chres!=NULL)
	{
	for (int i=0;i<ImageTemplate->width;i++)
	for (int j=0;j<ImageTemplate->height;j++)
	for (int k=0;k<ImageTemplate->depth;k++)
		if (ImageTemplate->mri[i][j][k]!=0)
			ImageTemplate->mri[i][j][k]=1;
	ImageTemplate->rcoeff = 1;
	ImageTemplate->icomp = 0;
	ImageTemplate->max_pixel=1;          
	for (int i=0;i<ImageTemplateNonDecoupe->width;i++)
	for (int j=0;j<ImageTemplateNonDecoupe->height;j++)
	for (int k=0;k<ImageTemplateNonDecoupe->depth;k++)
		if (ImageTemplateNonDecoupe->mri[i][j][k]!=0)
			ImageTemplateNonDecoupe->mri[i][j][k]=1;
	ImageTemplateNonDecoupe->rcoeff = 1;
	ImageTemplateNonDecoupe->icomp = 0;
	ImageTemplateNonDecoupe->max_pixel=1;
	
	interpol2(ImageTemplateNonDecoupe,chres,ImageTemplateOs); //J applique la transformation inverse a l image template : je mets le resultat dans imtres->mask 	 
	interpol2(ImageTemplate,chres,temporaire);                //J applique la transformation inverse a l image template : je mets le resultat dans imtres->mask 	 	
	
	
	imx_dilat_3d_p(temporaire,ImageTemplateOsDilate,1,bande);   
	for (int i=0;i<ImageTemplateOsDilate->width;i++)
	for (int j=0;j<ImageTemplateOsDilate->height;j++)
	for (int k=0;k<ImageTemplateOsDilate->depth;k++)
		if (ImageTemplateOsDilate->mri[i][j][k]!=0)
			ImageTemplateOsDilate->mri[i][j][k]=1;
	ImageTemplateOsDilate->rcoeff = 1;
	ImageTemplateOsDilate->icomp = 0;
	ImageTemplateOsDilate->max_pixel=1;          

	/*imageSegmentee = cr_grphic3d(temporaire);
	imx_dilat_3d_p(ImageTemplateOs,temporaire,1,bande);
	for (unsigned int i=0;i<temporaire->width;i++)
	for (unsigned int j=0;j<temporaire->height;j++)
	for (unsigned int k=0;k<temporaire->depth;k++)
		if (temporaire->mri[i][j][k])
			imageSegmentee->mri[i][j][k]=1;
		else
			imageSegmentee->mri[i][j][k]=0;
	
	for (int i=0;i<temporaire->width;i++)
	for (int j=0;j<temporaire->height;j++)
	for (int k=0;k<temporaire->depth;k++)
		{
		if (temporaire->mri[i][j][k] == 1)
			{
			for (int ii=MAXI(i-1,0);i<=MINI(i+1,temporaire->width-1);ii++)
			for (int jj=MAXI(j-1,0);j<=MINI(j+1,temporaire->height-1);jj++)
			for (int kk=MAXI(k-1,0);k<=MINI(k+1,temporaire->depth-1);kk++)
				{
				if ( (temporaire->mri[ii][jj][kk] == 3) || (temporaire->mri[ii][jj][kk] == 0))
					temporaire->mri[ii][jj][kk] = 2;
				}
			}
		if (temporaire->mri[i][j][k] == 0)
			{
			for (int ii=MAXI(i-1,0);i<=MINI(i+1,temporaire->width-1);ii++)
			for (int jj=MAXI(j-1,0);j<=MINI(j+1,temporaire->height-1);jj++)
			for (int kk=MAXI(k-1,0);k<=MINI(k+1,temporaire->depth-1);kk++)
				{
				if ( (temporaire->mri[ii][jj][kk] == 2) || (temporaire->mri[ii][jj][kk] == 1))
					temporaire->mri[ii][jj][kk] = 3;
				}
			}

		}*/


	/*imx_dilat_3d_p(ImageTemplate,temporaire,1,bande);
	for (int i=0;i<temporaire->width;i++)
	for (int j=0;j<temporaire->height;j++)
	for (int k=0;k<temporaire->depth;k++)
		if (temporaire->mri[i][j][k]!=0)
			temporaire->mri[i][j][k]=1;
	temporaire->rcoeff = 1;
	temporaire->icomp = 0;
	temporaire->max_pixel=1;
	interpol2(temporaire,chres,ImageTemplateOsDilate); //J applique la transformation inverse a l image template : je mets le resultat dans imtres->mask 	 */
	}
else
	{
	imx_copie_3d_p(ImageTemplateNonDecoupe,ImageTemplateOs);
	imx_dilat_3d_p(ImageTemplate,temporaire,1,bande);
	imx_copie_3d_p(temporaire,ImageTemplateOsDilate);
	}
scal_func=topo_choose_bspline(1);
end_base_3d();
nb_param=init_base_3d(ImageResultat->width,ImageResultat->height,ImageResultat->depth,scal_func);


int topD  = (int)pow(2.0,1.0*(BASE3D.resol))-1;
int *x0,*x1,*y0,*y1,*z0,*z1,l;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;

if (optimisationEchelle0==0)
	{
	Proportion=(int ***)alloc_imatrix_3d(topD,topD,topD);				
	for (int topi=0; topi<topD; topi++)
	for (int topj=0; topj<topD; topj++)
	for (int topk=0; topk<topD; topk++)
		{
		int compteur ;		   	
		compteur = 0;	
		l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
	
		for (int i=x0[l];i<x1[l];i++)
		for (int j=y0[l];j<y1[l];j++)
		for (int k=z0[l];k<z1[l];k++)
			{		
			if ((ImageTemplateOsDilate->mri[i][j][k]!=0)&&(ImageTemplateOs->mri[i][j][k]>0))			
				compteur ++;		
			}
		Proportion[topi][topj][topk]=compteur;
		}
	mini=minimisation22(ImageTemplateOsDilate  ,
			    ImageASegmenter    ,
			    ImageResultat      ,
			    nb_param	       ,
			    param              ,
			    seuilBas           ,
			    seuilHaut          ,
			    1,
			    ImageTemplateOs);
	seuilBas = param[0] ;
	for (l=0;l<nb_param;l++)
		param[l]=0.0;
	free_imatrix_3d(Proportion);
	}
else
	{
	for (int echelleSegmentation=1;echelleSegmentation<RESOLUTION;echelleSegmentation++)
		{	
		topD  = (int)pow(2.0,1.0*(BASE3D.resol))-1;	
		x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;

		Proportion=(int ***)alloc_imatrix_3d(topD,topD,topD);				
		for (int topi=0; topi<topD; topi++)
		for (int topj=0; topj<topD; topj++)
		for (int topk=0; topk<topD; topk++)
			{
			int compteur ;		   	
			compteur = 0;	
			l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
			
			for (int i=x0[l];i<x1[l];i++)
			for (int j=y0[l];j<y1[l];j++)
			for (int k=z0[l];k<z1[l];k++)
			{
		
			if ((ImageTemplateOsDilate->mri[i][j][k]!=0)&&(ImageTemplateOs->mri[i][j][k]>0))			
			compteur ++;		
			}
			Proportion[topi][topj][topk]=compteur;
			}
		
		//Lancement de l algorithme pour une echelle fixe. 	
		if (echelleSegmentation == RESOLUTION -1)
		mini=minimisation22(ImageTemplateOsDilate  ,
			  	  ImageASegmenter    ,
		  		  ImageResultat      ,
			  	  nb_param	     ,
			   	  param              ,
			  	  seuilBas           ,
		  		  seuilHaut          ,
			  	  0,
			  	  ImageTemplateOs);
	
		free_imatrix_3d(Proportion);
		nb_param=base_resol_up_3d(param,nb_param);
		}
	}
free_grphic3d(imageSegmentee);
imageSegmentee = NULL;

for (unsigned int i=0;i<ImageResultat->width;i++)
for (unsigned int j=0;j<ImageResultat->height;j++)
for (unsigned int k=0;k<ImageResultat->depth;k++)
	if (ImageTemplateOsDilate->mri[i][j][k]==0)
		ImageResultat->mri[i][j][k]=0;

//liberation memoire
free_grphic3d(ImageTemplateOs); 
free_grphic3d(ImageTemplateOsDilate);
free_grphic3d(temporaire);
return(1);  
}



//LANCEMENT DE LA SEGMENTATION ET DE L OPTIMISATION POUR UNE ECHELLE DONNEE : ON PASSE DONC EN REVUE TOUS LES BLOCS
double minimisation22(grphic3d *ImageTemplateOsDilate, 
		      grphic3d *ImageASegmenter  ,
		      grphic3d *ImageResultat    ,
                      int nb_param	       ,
		      double *param	       ,
		      double seuilBas            ,
		      double seuilHaut          ,
		      int optimisationEchelle0,
		      grphic3d *ImageTemplateOs)
{
int ***masque_param = NULL;

//Creation des images de seuil
float ***ImageSeuilBas  = (float ***)alloc_matrix_3d(ImageASegmenter->width,ImageASegmenter->height,ImageASegmenter->depth); 
float ***ImageSeuilHaut = (float ***)alloc_matrix_3d(ImageASegmenter->width,ImageASegmenter->height,ImageASegmenter->depth);

for (unsigned int i=0;i<ImageASegmenter->width;i++)
for (unsigned int j=0;j<ImageASegmenter->height;j++)
for (unsigned int k=0;k<ImageASegmenter->depth;k++)
	{
	ImageResultat->mri[i][j][k]=0;
	ImageSeuilBas[i][j][k] = 0;
	if (ImageTemplateOsDilate->mri[i][j][k]!=0)
		ImageSeuilHaut[i][j][k] = 30000;
	else
		ImageSeuilHaut[i][j][k] = -1;
	}
	

int topD  = (int)pow(2.0,1.0*BASE3D.resol)-1;
masque_param=(int ***)alloc_imatrix_3d(topD,topD,topD);


//On donne les parametres, l image a segmenter, les images template ImageTemplateOs et ImageTemplateVide ainsi que seuilBas et seuilHaut.
// l algo fournit  les images de seuil (ImageSeuilHaut et ImageSeuilBas), l image segmentee : ImageSegmenteeCourante, 
//ainsi que la fonction de cout.

//On regarde les blocs sur lesquels il faut optimiser.
for (int topi=0; topi<topD; topi++)
for (int topj=0; topj<topD; topj++)
for (int topk=0; topk<topD; topk++)
	{		
	int *x0,*x1,*y0,*y1,*z0,*z1,l;
	Points_3dInt pointBG       ; 
	Points_3dInt pointHD       ;
				   
	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
	l=(int)TOP_conv_ind(topi,topj,topk,nb_param)/3;
	pointBG.x=x0[l];pointHD.x=x1[l];
	pointBG.y=y0[l];pointHD.y=y1[l];
	pointBG.z=z0[l];pointHD.z=z1[l];	
	masque_param[topi][topj][topk] = -1; //On n optimise jamais sur le bloc 
	for (int i=pointBG.x;i<pointHD.x;i++)
	for (int j=pointBG.y;j<pointHD.y;j++)
	for (int k=pointBG.z;k<pointHD.z;k++)
		{
		if ((ImageTemplateOsDilate->mri[i][j][k]!=0)&&(ImageASegmenter->mri[i][j][k]!=0))
			masque_param[topi][topj][topk]=1; //On optimise y a peut etre de l os	
		}																									
	}


for (int topi=0; topi<topD; topi++)
for (int topj=0; topj<topD; topj++)
for (int topk=0; topk<topD; topk++)
	{
	if (masque_param[topi][topj][topk]!=1)
		{
		int li = (int)TOP_conv_ind(topi,topj,topk,nb_param);			
		param[li]=  MAXIMAL-seuilBas;
		}
	if (masque_param[topi][topj][topk]==1)
		{
		int *x0,*x1,*y0,*y1,*z0,*z1,l;
		Points_3dInt pointBG       ; 
		Points_3dInt pointHD       ;
				   
		x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
		l=(int)TOP_conv_ind(topi,topj,topk,nb_param)/3;
		pointBG.x=x0[l];pointHD.x=x1[l];
		pointBG.y=y0[l];pointHD.y=y1[l];
		pointBG.z=z0[l];pointHD.z=z1[l];	
		short int min,max;
		int nimportequoi = 0;
		min = 0;
		max = 0;
		int nbreVoxel = 0;
		for (int i=pointBG.x;i<pointHD.x;i++)
		for (int j=pointBG.y;j<pointHD.y;j++)
		for (int k=pointBG.z;k<pointHD.z;k++)
			{
			if ((ImageTemplateOsDilate->mri[i][j][k]!=0)&&(ImageASegmenter->mri[i][j][k]!=0))
					{
					if (nimportequoi ==0)
						{
						min = max = ImageASegmenter->mri[i][j][k];
						nimportequoi=1;						
						}
					
					if (ImageASegmenter->mri[i][j][k]<min)
						{
						min = ImageASegmenter->mri[i][j][k];
						}
					if (ImageASegmenter->mri[i][j][k]>max)
						max = ImageASegmenter->mri[i][j][k];	
					nbreVoxel++;		
					}
					
			}
		if (max == min)
			{
			masque_param[topi][topj][topk]=-1;
			int li = (int)TOP_conv_ind(topi,topj,topk,nb_param);			
			param[li] =  MAXIMAL - seuilBas ;		
			li = li/3;
			pointBG.x=x0[li];pointHD.x=x1[li];
			pointBG.y=y0[li];pointHD.y=y1[li];
			pointBG.z=z0[li];pointHD.z=z1[li];	

			for (int i=pointBG.x;i<pointHD.x;i++)
			for (int j=pointBG.y;j<pointHD.y;j++)
			for (int k=pointBG.z;k<pointHD.z;k++)
				{
				ImageSeuilHaut[i][j][k]=-1;
				}	
			}
		else					
			{
			vector<unsigned int> histoC ;	
			histoC.resize(max-min+1); 
			for (int i=0;i<(max-min+1);i++)
				histoC[i]=0;			
			for (int i=pointBG.x;i<pointHD.x;i++)
			for (int j=pointBG.y;j<pointHD.y;j++)
			for (int k=pointBG.z;k<pointHD.z;k++)
				{
				if ((ImageTemplateOsDilate->mri[i][j][k]!=0)&&(ImageASegmenter->mri[i][j][k]!=0))
					{									
					histoC[ImageASegmenter->mri[i][j][k]-min]++;				
					}
				}
			double seuil = 0;
			if (histoOs != NULL)
				{
				double p1 = (double)Proportion[topi][topj][topk] / (double)nbreVoxel; //Proportion d os
				if (p1<0.05)
					p1 = 0.05;
				if (p1>0.95)
					p1 = 0.95;
				int th = ChoixSeuil(histoC,min,p1);	
				if ( (th == histoC.size()+min) ) //Il n y a pas d os dans le bloc.
					{
					masque_param[topi][topj][topk]=-1;
					th =  MAXIMAL- seuilBas ;	
					int li = (int)TOP_conv_ind(topi,topj,topk,nb_param);			
					param[li] =  MAXIMAL - seuilBas ;		
					li = li/3;
					pointBG.x=x0[li];pointHD.x=x1[li];
					pointBG.y=y0[li];pointHD.y=y1[li];
					pointBG.z=z0[li];pointHD.z=z1[li];	

					for (int i=pointBG.x;i<pointHD.x;i++)
					for (int j=pointBG.y;j<pointHD.y;j++)
					for (int k=pointBG.z;k<pointHD.z;k++)
						{
						ImageSeuilHaut[i][j][k]=-1;
						}
					}
				seuil = th - 0.5 ; //cf fonction ChoixSeuil				
				}
			else
				{
				int itMilieu =0;
				int compteur = 0;
				for (int it=max-min; it >=0;it--)
					{
					compteur += histoC[it];				
					if (compteur>Proportion[topi][topj][topk])
						{					
						itMilieu = it;
						break;
						}
					}
				histoPartieMolle = (double*) malloc(sizeof(double)*(max-min+1));
				histoOs = (double*)malloc(sizeof(double)*(max-min+1));
				tailleHisto = max - min + 1;
				MINIMAL = min;
				MAXIMAL = max;
				seuil = itMilieu + min;
				}
			l = (int)TOP_conv_ind(topi,topj,topk,nb_param);
			seuil *= ImageASegmenter->rcoeff;			
			if (optimisationEchelle0 == 0)
				param[l] = seuil - seuilBas ;		
			else
				seuilBas = seuil;
			}			
		} 
	}

//Determination des parametres...

if (optimisationEchelle0 == 0)
	DetermineParameter(param	        ,
			   nb_param		,
			   ImageTemplateOsDilate,
			   ImageSeuilBas        ,
			   seuilBas             ,
			   masque_param	       );
else			   
	for (unsigned int i=0;i<ImageASegmenter->width;i++)
	for (unsigned int j=0;j<ImageASegmenter->height;j++)
	for (unsigned int k=0;k<ImageASegmenter->depth;k++)
		ImageSeuilBas[i][j][k] = seuilBas;

//Determination de la carte de seuil
/*TranslateBase3dToWholeImage(ImageSeuilBas ,
			    param         ,
			    1		  ,	 
			    nb_param    ,
			    seuilBas     );*/
			    			

//Calcul de la segmentation
Segmentation   (ImageSeuilHaut     ,
         	ImageSeuilBas      ,
               	ImageASegmenter    ,
               	ImageResultat     );


for (int i=0;i<tailleHisto;i++)
	{
	histoPartieMolle[i] = histoOs[i] = 0;
	}
for (int i=0;i<ImageResultat->width;i++)
for (int j=0;j<ImageResultat->height;j++)
for (int k=0;k<ImageResultat->depth;k++)
	{
	if ((ImageTemplateOsDilate->mri[i][j][k]!=0)&&(ImageASegmenter->mri[i][j][k]!=0))
		{
		if (ImageResultat->mri[i][j][k]!=0)
			histoOs[ImageASegmenter->mri[i][j][k]-MINIMAL]++;
		else
			histoPartieMolle[ImageASegmenter->mri[i][j][k]-MINIMAL]++;		
		}
	}

EstimeGauss(histoOs, tailleHisto);
EstimeGauss(histoPartieMolle, tailleHisto);

if (optimisationEchelle0 != 0)
	param[0]= seuilBas;
	
//liberation memoire
free_matrix_3d(ImageSeuilHaut)        ;
free_matrix_3d(ImageSeuilBas)         ;
free_imatrix_3d(masque_param);
return 0;
}		     
		    




//FONCTION PERMETTANT DE FAIRE LE PASSAGE ENTRE LES BASES DE SPLINE ET LES IMAGES DE SEUIL HAUT ET BAS
void TranslateBase3dToWholeImage(float ***ImageRes     ,
				 double *param          ,
			   	 int numImage           ,	 
			    	 int nb_param      	,
				 double seuil           )
{
int topD=(int)pow(2.0,1.0*BASE3D.resol)-1;

for (int topi=0; topi<topD; topi+=2)
for (int topj=0; topj<topD; topj+=2)
for (int topk=0; topk<topD; topk+=2)
	TranslateBase3dToImage(ImageRes,param,numImage,nb_param,topi,topj,topk,seuil);

}

		    
void TranslateBase3dToImage(float ***ImageRes     ,
			    double *param          ,
			    int numImage           ,	 
			    int nb_param           ,			    
			    int topi               ,
			    int topj               ,
			    int topk               ,
			    double seuil           )
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz,ux,uy,uz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py;
int i,j,k,auxi,auxj,auxk;
field3d *champ;
vector3d ***data;


resol=BASE3D.resol;
width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
topD=(int)pow(2.0,1.0*resol)-1;


l=(int)TOP_conv_ind(topi,topj,topk,nb_param)/3;
bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;


champ=cr_field3d(topDx,topDy,topDz);

data=champ->raw;
for (i=0;i<topDx;i++)
for (j=0;j<topDy;j++)
for (k=0;k<topDz;k++)
	data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;

			
// Determination des valeurs de l image dans la boite topi,topj,topk.
for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
	{
	l=(int)TOP_conv_ind(i,j,k,nb_param)/3;

	x00=x0[l];x11=x1[l];y00=y0[l];y11=y1[l];z00=z0[l];z11=z1[l];
   	px=param[3*l];py=param[3*l+1];//pz=param[3*l+2];
   	for (auxi=MAXI(bx0,x00);auxi<MINI(bx1,x11);auxi++)
    	for (auxj=MAXI(by0,y00);auxj<MINI(by1,y11);auxj++)
     	for (auxk=MAXI(bz0,z00);auxk<MINI(bz1,z11);auxk++) 
      		{
		f=fx[auxi-x00]*fy[auxj-y00]*fz[auxk-z00];
      		ux=auxi-bx0;uy=auxj-by0;uz=auxk-bz0;
      		data[ux][uy][uz].x=data[ux][uy][uz].x+px*f;
      		data[ux][uy][uz].y=data[ux][uy][uz].y+py*f;
      		//data[ux][uy][uz].z=data[ux][uy][uz].z+pz*f;
		}  
	}

// On met a jour l image
for (i=0;i<topDx;i++)
for (j=0;j<topDy;j++)
for (k=0;k<topDz;k++)
	{
	if (numImage == 1)
		ImageRes[i+bx0][j+by0][k+bz0] = (data[i][j][k].x+seuil) ;
	else
		ImageRes[i+bx0][j+by0][k+bz0] = (data[i][j][k].y+seuil) ;
	}
free_field3d(champ);
}	
	




/*******************************************************************************
********************************************************************************
*************** RECALAGE BSPLINE  Primitives geometriques non segmentees********
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**        matching_Bspline_topo_primitive_non_seg_3d                                                
**                                                                        
**       recalage Bspline de deux images binaires 3D (dont l une n est pas encore segmentee)avec conservation de la topologie                               
*******************************************************************************/

void matching_Bspline_topo_primitive_non_seg_3d()
{

char *quest[15],*nomfichres;
int dist_type,inter_type,min_type,save_type,func_type,i, reg_type=0;
double Jmin,Jmax;
int resolf, e, err;
 
//question sur les images a recaler
int im_reca=GET_PLACE3D(TEXT0030);
int im_ref = GET_PLACE3D("template decoupe");
int im_ref2= GET_PLACE3D("template non decoupe");
int im_res= GET_PLACE3D(TEXT0006);

for(i=0;i<15;i++)
   quest[i]=CALLOC(80,char);
 
//type de fonction
func_type=1;

//distance
dist_type=0;

//interpolation
inter_type=0;

//minimisation
strcpy(quest[0],"ICM");
strcpy(quest[1],"Descente gradient");
strcpy(quest[2],"\0");
//min_type=GETV_QCM("Minimisation",(char **)quest);
min_type = 1;
//question sur l'enregistrement du champ resultat*/
nomfichres=quest_save_result_3d(&save_type);

//resolution finale
resolf= GET_INT("resolution", 4, &e);

//bornes sur le jacobien
Jmin = 2;
while (Jmin>=1)  
	Jmin = GET_DOUBLE("Jmin<1",0,&err);
Jmax = 0;
while (Jmax<=1)  
	Jmax = GET_DOUBLE("Jmax>1",100000,&err);


//regularisation 
strcpy(quest[0],"Pas de regularisation");
strcpy(quest[1],"Regularisation membrane elastique");
strcpy(quest[3],"\0");
//reg_type= GETV_QCM("R�gularisation ?",(char **)quest);
reg_type= 1;
TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=0;
if (reg_type>0)
	TOPO_REGULARISATION_MEMBRANE_ELASTIQUE = 1;//GET_DOUBLE("facteur lambda",1,&err);
for(i=0;i<15;i++)
    free(quest[i]);
		
imx_matching_Bspline_topo_primitive_non_seg_3d(im_ref,im_ref2,im_reca,im_res,func_type,dist_type,reg_type,inter_type,min_type,save_type,nomfichres,resolf,Jmin,Jmax);

 
 if (nomfichres)
  free(nomfichres);

  free(histoPartieMolle);free(histoOs);
  histoOs = NULL;
  histoPartieMolle = NULL;
  show_picture_3d(im_res);
 
}

/*******************************************************************************
**        imx_matching_Bspline_topo_primitive_non_seg_3d                                                
**                                                                        
**       recalage Bspline de deux images binaires 3D (l une n est pas encore segmentee)avec conservation de la topologie                               
*******************************************************************************/


int imx_matching_Bspline_topo_primitive_non_seg_3d(int im_ref,int im_ref2, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax)
{
int l,continu;
grphic3d *imref,*imref2,*imreca,*imres;
imref=ptr_img_3d(im_ref);
imref2=ptr_img_3d(im_ref2);
imreca=ptr_img_3d(im_reca);
imres=ptr_img_3d(im_res);
imageSegmentee = NULL;
		
/*imx_Copie_imgmask_3d_p(imref,imref);
imx_Copie_imgmask_3d_p(imref2,imref2);*/

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
	if (imref->width==(unsigned int)pow(2.0,l)) continu++;
	if (imref->height==(unsigned int)pow(2.0,l)) continu++;
	if (imref->depth==(unsigned int)pow(2.0,l)) continu++;
	if (imreca->width==(unsigned int)pow(2.0,l)) continu++;
	if (imreca->height==(unsigned int)pow(2.0,l)) continu++;
	if (imreca->depth==(unsigned int)pow(2.0,l)) continu++;
	}

if (continu!=6)
	{
      PUT_WARN("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp");
		  printf("Les tailles des images ne sont pas des puissances de 2 !!! Modification temporaire de la taille des images. Le resultat sera stocke sous forme de fichier chp\n");
			//return(-1) ;
			save_type=3;
   }

imx_matching_Bspline_topo_primitive_non_seg_3d_p(imref,imref2,imreca,imres,func_type,dist_type, reg_type, inter_type,min_type,save_type,nomfichres,resolf, Jmin, Jmax);
 
return(1);
}




/*******************************************************************************
**        imx_matching_Bspline_topo_primitive_3d_p                                                
**                                                                        
*       recalage Bspline de deux images binaires 3D avec conservation de la topologie                               
*******************************************************************************/

int imx_matching_Bspline_topo_primitive_non_seg_3d_p(grphic3d *imref,grphic3d *imref2, grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax)
{
int	wdth,hght,dpth,nb_param, wdth_old=0, hght_old=0, dpth_old=0,i,j,k,ideb,jdeb,kdeb,trucMuche;
double *param,*param_norm,*param_segmentation,*paramBidon;
double seuilBas = 200;
double seuilHaut = 400000;
grphic3d *imtref,*imtref2,*imtreca,*imtres,*imtmp,*imtseg;
field3d *champ;
InterpolationFct interpol;
dist_func_locale_t distance;
reg_func_locale_t regularisation=NULL;
min_func_locale_t minimisation;
scal_func_t scal_func;
double  mini;
int tdebut,tfin;
int nbmax_param;
char nomfichier_nrj[255];

//#ifndef SAVE_INTERMEDIAIRE_2 
char nomfichier[255];
char nomfichier2[255];
char temp[5];
char *ext2;
//#endif
/*#ifndef TOPO_COMPARE_CRITERE 
	char nomfichier_nrj[255];
	char *ext;
#endif*/
#ifndef WIN32
 	setvbuf(stdout, (char *)NULL, _IONBF, 0); // setbuf(stdout,NULL);
#endif


wdth=imref->width;hght=imref->height;dpth=imref->depth;
	
if (save_type==3) // il faut transformer les images en puissance de 2 
	{
	wdth_old=wdth; hght_old=hght; dpth_old=dpth;
	
	imx_realloc_mri_pow2_centered_p(imref);
	imx_realloc_mri_pow2_centered_p(imreca);
	imx_realloc_mri_pow2_centered_p(imres);
	
	wdth=imref->width;hght=imref->height;dpth=imref->depth;
		
	}
	
tdebut=time(NULL);
nbmax_param=3*(int)pow((pow(2.0,resolf)-1.0),3.0);
 
// Allocation 
if((param = (double*)calloc(nbmax_param, sizeof(double))) == NULL) 
	{
	PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); 
	}
if((param_norm = (double*)calloc(nbmax_param, sizeof(double))) == NULL) 
	{
	PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); 
	}
if((param_segmentation = (double*)calloc(nbmax_param, sizeof(double))) == NULL) 
	{
	PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); 
	}

if((paramBidon = (double*)calloc(nbmax_param, sizeof(double))) == NULL) 
	{
	PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(0); 
	}
for (i=0;i<nbmax_param;i++)
	{
	param[i] = param_norm[i] = param_segmentation[i] = paramBidon[i]=0;
	}
	
init_log(NULL);
aff_log("\nRECALAGE ");
 
//mise en place des parametres du recalage
scal_func=topo_choose_bspline(func_type);
distance=Energie_ICP_sym_locale_3d; 
interpol=imx_choose_interpolation_fct(inter_type);  
minimisation=topo_choose_optimisation(min_type);
if (reg_type>0)
	regularisation=regularisation_energie_membrane_local;

aff_log("\n");
aff_log("recalage de %s sur %s \n",(imreca->patient).name,(imref->patient).name);
 
//creation des images temporaires
imtref=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref);
imtref2=cr_grphic3d(imref2);imx_copie_3d_p(imref2,imtref2);
imtreca=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca);
imtres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtres);
imtreca->mask=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtreca->mask);
imtref->mask=cr_grphic3d(imref);imx_copie_3d_p(imref,imtref->mask); 
imtref2->mask=cr_grphic3d(imref2);imx_copie_3d_p(imref2,imtref2->mask); 
imtres->mask=cr_grphic3d(imres);imx_copie_3d_p(imres,imtres->mask); 
imtmp=cr_grphic3d(imreca);imx_copie_3d_p(imreca,imtmp);
imtseg = cr_grphic3d(imreca);		

//allocation memoire de champ
champ=cr_field3d(wdth,hght,dpth);
imx_inv_3d_p(imtref,imtref->mask);
	 

transf3d *transfo; 
field3d *chres ;
trucMuche = 1;
/********************RECALAGE ONDELETTE******************************************/
{
int l,resol,r,D,bande;
int ***masque_param;
int *x0,*x1,*y0,*y1,*z0,*z1;
nb_param=init_base_3d(wdth,hght,dpth,scal_func);


resol=BASE3D.resol;
		
//Initialisation tableau 
for (l=0;l<nb_param;l++)
	param[l]=0.0;

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
	 

//SEGMENTATION A L ECHELLE R
r = 0;	
chres=NULL;
RESOLUTION = 1;
imx_segmentation_p2(imtref    , 	        //Image template
		    imtref2 	      , 	//IMAGE TEMPLATE NON DeCOUPE 
	       	    imreca	      , 	//Image a segmenter
	       	    imtseg  	      , 	//Image segmente
	       	    param_segmentation, 	//Parametre de segmentation
	       	    nb_param          ,
	       	    seuilHaut         ,
	       	    seuilBas	      ,                       	   
	       	    chres	      ,
		    1                );
		    
//char nomfichier_nrj[500];
//sprintf(nomfichier_nrj,"/home/faisan/yu.ipb");			  
//save_mri_ipb_3d_p(nomfichier_nrj,imtseg);

for (int ia=0;ia<nbmax_param;ia++)
	{
	param_segmentation[ia] = 0;
	}	

//int iterationFinale = 0;
int iterationFinale = 1;
// Debut traitement 
double lambda = TOPO_REGULARISATION_MEMBRANE_ELASTIQUE;
for (r=resol;r<=resolf;r++)
	{
	RESOLUTION = r;
	x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;	
	if (r>=3)
		TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
	else
		TOPO_REGULARISATION_MEMBRANE_ELASTIQUE=lambda;
	D=(int)TOP_D(nb_param);	
	 
	//calcul du mask de imref avec l'image compl�mentaire 
	imx_inv_3d_p(imtseg,imtmp);	
	
	//calcul des cartes de distance 
	imx_chamfer_distance_3d_p(imtseg,imtreca); //Sur l image de segmentation
	imx_chamfer_distance_3d_p(imtmp,imtreca->mask); 

	imx_mul_coe_3d_p(imtreca,1000.0,imtreca);
	imx_mul_coe_3d_p(imtreca->mask,1000.0,imtreca->mask);

	
	
	aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);
	base_to_field_3d(nb_param,param,champ,NULL,NULL);
	interpol(imtseg,champ,imtres); 	
	bande=(int)((x1[0]-x0[0])/4.0);	
	if (bande<10)
		{
		imx_dilat_3d_p(imtref,imtref->mask,1,bande);
		imx_sub_3d_p(imtref->mask,imtref,imtref->mask);
		}
	//Creation du masque binaire r�f�rencant les Slpqr sur lesquel l'image ne contient pas d'info 
	masque_param=alloc_imatrix_3d(D,D,D);
	// Remplissage de masque_param 
	topo_masque_param_primitive (imtref, masque_param, nb_param);

	// A modifier et a adpater
	if ((TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)&&(trucMuche==1))
		update_TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE2(imtref,imtreca,imtres,champ,minimisation,base_to_field_3d,interpol,distance,regularisation,nb_param,param,masque_param,Jmin,Jmax,nomfichres);

	mini=minimisation(imtref,imtreca,imtres,champ,base_to_field_3d,interpol,distance, regularisation, nb_param,param,masque_param,Jmin,Jmax,nomfichres);

	trucMuche=0;
	#ifdef SAVE_INTERMEDIAIRE_2
	
	if (nomfichres!=NULL)
		{
		int unite = iterationFinale %10;
		int decimal = iterationFinale / 10;
	 	/*calcul de la transformation*/
	 	base_to_field_3d(nb_param,param,champ,NULL,NULL);

	 	//si l'extension .trf existe on la supprime
   		strcpy(nomfichier,nomfichres);
  		ext2=strstr(nomfichier,".trf");
  		if (ext2!=NULL) *ext2='\0';

		strcat(nomfichier,"_");
 		temp[0]=(char)(r+48);
 		temp[1]='_';
		temp[2]=(char)(decimal+48);				
		temp[3]=(char)(unite+48);
 		temp[4]=0;
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
				
				champ->width=	wdth_old;
				champ->height=	hght_old;
				champ->depth=	dpth_old;
				
				imref->width=	wdth_old;
				imref->height=	hght_old;
				imref->depth=	dpth_old;
		
				imreca->width=	wdth_old;
				imreca->height=	hght_old;
				imreca->depth=	dpth_old;		
				} 
							
			 transfo=field_to_transf_3d(champ,imref,imreca); 
    		
			if (save_type==3)
				{
				champ->width=	wdth;
				champ->height=	hght;
				champ->depth=	dpth;
				
				imref->width=	wdth;
				imref->height=	hght;
				imref->depth=	dpth;
		
				imreca->width=	wdth;
				imreca->height=	hght;
				imreca->depth=	dpth;
					
				}		
			}
		else
			{  
     		     	transfo=cr_transf3d(wdth,hght,dpth,NULL);
		        transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
		        transfo->typetrans=BSPLINE3D;
		        transfo->resol=BASE3D.resol;transfo->degre=func_type;
			transfo->dx=imtref->dx;
			transfo->dy=imtref->dy;
			transfo->dz=imtref->dz;
			for (l=0;l<nb_param/3;l++) 
		 		{
				transfo->param[3*l]=param[3*l]*imreca->dx;
				transfo->param[3*l+1]=param[3*l+1]*imreca->dy;
				transfo->param[3*l+2]=param[3*l+2]*imreca->dz;
				}
			//On met a 0 les parametres qui nous interessent pas.
			int topDD  = (int)pow(2.0,1.0*BASE3D.resol)-1;
			for (int topi=0; topi<topDD; topi++)
			for (int topj=0; topj<topDD; topj++)
			for (int topk=0; topk<topDD; topk++)
				{		
				int *x0,*x1,*y0,*y1,*z0,*z1,l;
				Points_3dInt pointBG       ; 
				Points_3dInt pointHD       ;
				   
				x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
				l=(int)TOP_conv_ind(topi,topj,topk,nb_param)/3;
				pointBG.x=x0[l];pointHD.x=x1[l];
				pointBG.y=y0[l];pointHD.y=y1[l];
				pointBG.z=z0[l];pointHD.z=z1[l];	
				int tutu = -1;
				for (int i=pointBG.x;i<pointHD.x;i++)
				for (int j=pointBG.y;j<pointHD.y;j++)
				for (int k=pointBG.z;k<pointHD.z;k++)
					{
					if (imtref->mri[i][j][k]!=0)
						tutu=1;	
					}																									
				if (tutu==-1)
					{
					transfo->param[3*l] = transfo->param[3*l+1] = transfo->param[3*l+2] = 0;
					}
				}			
							
		        
			}
	save_transf_3d(transfo,nomfichier);
        free_transf3d(transfo);
	} 
	#endif	    	 
	
	chres=cr_field3d(wdth,hght,dpth)    ;	
	transfo=cr_transf3d(wdth,hght,dpth,NULL);
	transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
	transfo->typetrans=BSPLINE3D;
	transfo->resol=r;transfo->degre=func_type;
	transfo->dx=imref->dx;
	transfo->dy=imref->dy;
	transfo->dz=imref->dz;
		 
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
	//INVERSION DU CHAMPS
	end_base_3d(); //ON SUPPRIME LA BASE 3D CAR LA FONCTION INVERSION LA RECREE
	inv_bspline_3d(transfo,chres,0.3)  ; ///CCCCCCC
	nb_param=init_base_3d(wdth,hght,dpth,scal_func); //ON REMET COMME C ETAIT AVANT
	free_transf3d(transfo);
		
	for (int ia=0;ia<nbmax_param;ia++)
		{
		param_segmentation[ia] = 0;
		}
	imx_segmentation_p2(imtref    , 	//Image template
		    imtref2 	      , 	//IMAGE TEMPLATE NON DeCOUPE 
	       	    imreca	      , 	//Image a segmenter
	       	    imtseg  	      , 	//Image segmente
	       	    param_segmentation, 	//Parametre de segmentation
	       	    nb_param          ,
	       	    seuilHaut         ,
	       	    seuilBas	      ,                       	   
	       	    chres	      ,
		    1                );
		   
	if (chres!=NULL)	
		{
		free_field3d(chres);
		chres=NULL;
		}

	if (r!=resolf)
	    {
	    end_base_3d(); // On supprime BASE3D pour faire les changements d echelles avec les parametres de segmentation
	    nb_param=init_base_3d(wdth,hght,dpth,scal_func); //ON REMET COMME C ETAIT AVANT
	    while(BASE3D.resol!=r)
		nb_param=base_resol_up_3d(paramBidon,nb_param);
	    if (iterationFinale == 1)
	    	{
		imx_copie_3d_p(imtseg,imres);
		r --; //Pour faire une nouvelle iteration
		iterationFinale++;	
		/*imx_copie_3d_p(imtseg,imres);
		nb_param=base_resol_up_3d(param,nb_param);	   			
		iterationFinale = 1;*/
		}	
	    else
	    	{
		printf("%d KKKKKKKKKK\n",iterationFinale);
		int unite = iterationFinale %10;
		int decimal = iterationFinale / 10;
	 	/*calcul de la transformation*/
	 	
		strcpy(nomfichier,nomfichres);
  		ext2=strstr(nomfichier,".trf");
  		if (ext2!=NULL) *ext2='\0';
		strcat(nomfichier,"_");
 		temp[0]=(char)(r+48);
 		temp[1]='_';
		temp[2]=(char)(decimal+48);				
		temp[3]=(char)(unite+48);
 		temp[4]=0;
 		strcat(nomfichier,temp);	 
		unite = (iterationFinale-1) %10;
		decimal = (iterationFinale-1) / 10;
	 	/*calcul de la transformation*/
	 	
		strcpy(nomfichier2,nomfichres);
  		ext2=strstr(nomfichier2,".trf");
  		if (ext2!=NULL) *ext2='\0';

		strcat(nomfichier2,"_");
 		temp[0]=(char)(r+48);
 		temp[1]='_';
		temp[2]=(char)(decimal+48);				
		temp[3]=(char)(unite+48);
 		temp[4]=0;
 		strcat(nomfichier2,temp); 
		 
		char n2[200];
		char n1[200];
		 
		put_file_extension(nomfichier,".trf",n1);
		put_file_extension(nomfichier2,".trf",n2);
		 end_base_3d(); // On supprime BASE3D pour faire les changements d echelles avec les parametres de segmentation
		 transf3d *transfo1=load_transf_3d(n1);
		 transf3d *transfo2=load_transf_3d(n2);
		 field3d *champ1=transf_to_field_3d(transfo1,NULL,NULL);
		 field3d *champ2=transf_to_field_3d(transfo2,NULL,NULL);
		 double valueMax = 0;
		  for(i=0;i<transfo1->width;i++)
		  for(j=0;j<transfo1->height;j++)
		  for(k=0;k<transfo1->depth;k++)
			{
			double value = 0;
			value+= ( (champ2->raw[i][j][k].x - champ1->raw[i][j][k].x)*(champ2->raw[i][j][k].x - champ1->raw[i][j][k].x));
			value+= ( (champ2->raw[i][j][k].y - champ1->raw[i][j][k].y)*(champ2->raw[i][j][k].y - champ1->raw[i][j][k].y));		
			value+= ( (champ2->raw[i][j][k].z - champ1->raw[i][j][k].z)*(champ2->raw[i][j][k].z - champ1->raw[i][j][k].z));				 
			value = sqrt(value);
			if (value>valueMax)
				valueMax = value;
			}
		free_field3d(champ1);
		free_field3d(champ2);
		free_transf3d(transfo1);
		free_transf3d(transfo2);		
		
	    	nb_param=init_base_3d(wdth,hght,dpth,scal_func); //ON REMET COMME C ETAIT AVANT
	    	while(BASE3D.resol!=r)
			nb_param=base_resol_up_3d(paramBidon,nb_param);
		
		printf("valueMax %f iterationFinale %d \n",valueMax,iterationFinale);
		if ( (valueMax>0.25)&&(iterationFinale<3) )
			{
			imx_copie_3d_p(imtseg,imres);
			r --; //Pour faire une nouvelle iteration
			iterationFinale++;	
			}
		else
			{
			imx_copie_3d_p(imtseg,imres);
			nb_param=base_resol_up_3d(param,nb_param);	   			
			iterationFinale = 1;
			trucMuche = 1;
			}
		
		}

	   }
	else
	    {
	    end_base_3d(); // On supprime BASE3D pour faire les changements d echelles avec les parametres de segmentation
	    nb_param=init_base_3d(wdth,hght,dpth,scal_func); //ON REMET COMME C ETAIT AVANT
	    while(BASE3D.resol!=r)
		nb_param=base_resol_up_3d(paramBidon,nb_param);
	    if (iterationFinale ==1 )
	    	{		
		imx_copie_3d_p(imtseg,imres);
		r = resolf-1; //Pour faire une nouvelle iteration
		iterationFinale++;
		/*imx_copie_3d_p(imtseg,imres);*/
		}
	    else
	    	{
		int unite = iterationFinale %10;
		int decimal = iterationFinale / 10;
	 	/*calcul de la transformation*/
	 	
		strcpy(nomfichier,nomfichres);
  		ext2=strstr(nomfichier,".trf");
  		if (ext2!=NULL) *ext2='\0';

		strcat(nomfichier,"_");
 		temp[0]=(char)(r+48);
 		temp[1]='_';
		temp[2]=(char)(decimal+48);				
		temp[3]=(char)(unite+48);
 		temp[4]=0;
 		strcat(nomfichier,temp);
		 
		 
		unite = (iterationFinale-1) %10;
		decimal = (iterationFinale-1) / 10;
	 	/*calcul de la transformation*/
	 	
		strcpy(nomfichier2,nomfichres);
  		ext2=strstr(nomfichier2,".trf");
  		if (ext2!=NULL) *ext2='\0';

		strcat(nomfichier2,"_");
 		temp[0]=(char)(r+48);
 		temp[1]='_';
		temp[2]=(char)(decimal+48);				
		temp[3]=(char)(unite+48);
 		temp[4]=0;
 		strcat(nomfichier2,temp); 
		 
		char n2[200];
		char n1[200];
		 end_base_3d();
		put_file_extension(nomfichier,".trf",n1);
		put_file_extension(nomfichier2,".trf",n2);
		  
		 transf3d *transfo1=load_transf_3d(n1);
		 transf3d *transfo2=load_transf_3d(n2);
		 field3d *champ1=transf_to_field_3d(transfo1,NULL,NULL);
		 field3d *champ2=transf_to_field_3d(transfo2,NULL,NULL);
		 double valueMax = 0;
		  for(i=0;i<transfo1->width;i++)
		  for(j=0;j<transfo1->height;j++)
		  for(k=0;k<transfo1->depth;k++)
			{
			double value = 0;
			value+= ( (champ2->raw[i][j][k].x - champ2->raw[i][j][k].x)*(champ1->raw[i][j][k].x - champ1->raw[i][j][k].x));
			value+= ( (champ2->raw[i][j][k].y - champ2->raw[i][j][k].y)*(champ1->raw[i][j][k].y - champ1->raw[i][j][k].y));		
			value+= ( (champ2->raw[i][j][k].z - champ2->raw[i][j][k].z)*(champ1->raw[i][j][k].z - champ1->raw[i][j][k].z));				 
			value = sqrt(value);
			if (value>valueMax)
				valueMax = value;
			}
		free_field3d(champ1);//RAJOUT
		free_field3d(champ2);//RAJOUT
		free_transf3d(transfo1);//RAJOUT
		free_transf3d(transfo2);//RAJOUT				
	   
	    nb_param=init_base_3d(wdth,hght,dpth,scal_func); //ON REMET COMME C ETAIT AVANT
	    while(BASE3D.resol!=r)
		nb_param=base_resol_up_3d(paramBidon,nb_param);
		//printf("compteur %d iterationFinale %d \n",compteur,iterationFinale);
		if ( (valueMax>0.25)&&(iterationFinale<3) )
			{
			imx_copie_3d_p(imtseg,imres);
			r = resolf-1; //Pour faire une nouvelle iteration
			}
		iterationFinale++;			
		}
	    }	  
	free_imatrix_3d(masque_param);
	}       

	
//calcul de la transformation
base_to_field_3d(nb_param,param,champ,NULL,NULL);
 
if (imres != NULL)
	{ 
	imx_copie_3d_p(imtseg,imres);
	}
if (save_type==3) // il faut transformer les images  dans leur taille d'origine
	{	
	imx_undo_mri_pow2_centered_p(imref,wdth_old,hght_old,dpth_old);
	imx_undo_mri_pow2_centered_p(imreca,wdth_old,hght_old,dpth_old);
	imx_undo_mri_pow2_centered_p(imres,wdth_old,hght_old,dpth_old);
	}
//enregistrement du champ resultat dans un fichier
if (nomfichres!=NULL)
	{ 
	transf3d *transfo; 
    
	if ((save_type==0)||(save_type==3))
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
			champ->width=	wdth_old;
			champ->height=	hght_old;
			champ->depth=	dpth_old;
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
		//On met a 0 les parametres qui nous interessent pas.
		int topDD  = (int)pow(2.0,1.0*resolf)-1;
		for (int topi=0; topi<topDD; topi++)
		for (int topj=0; topj<topDD; topj++)
		for (int topk=0; topk<topDD; topk++)
			{		
			int *x0,*x1,*y0,*y1,*z0,*z1,l;
			Points_3dInt pointBG       ; 
			Points_3dInt pointHD       ;
				   
			x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
			l=(int)TOP_conv_ind(topi,topj,topk,nb_param)/3;
			pointBG.x=x0[l];pointHD.x=x1[l];
			pointBG.y=y0[l];pointHD.y=y1[l];
			pointBG.z=z0[l];pointHD.z=z1[l];	
			int tutu = -1;
			for (int i=pointBG.x;i<pointHD.x;i++)
			for (int j=pointBG.y;j<pointHD.y;j++)
			for (int k=pointBG.z;k<pointHD.z;k++)
				{
				if (imtref->mri[i][j][k]!=0)
					tutu=1;	
				}																									
			if (tutu==-1)
				{
				transfo->param[3*l] = transfo->param[3*l+1] = transfo->param[3*l+2] = 0;
				}
			}	
    		}
    	save_transf_3d(transfo,nomfichres);
    	free_transf3d(transfo);
   	} 
}//FIN RECALAGE ONDELLETE

//liberation de la base
end_base_3d();
if (chres!=NULL)	
	{
	free_field3d(chres);
	chres=NULL;
	}
//caracteristique du champ
aff_log("CHAMP: ");carac_field_3d(champ);


//liberation memoire du champ
free_field3d(champ);
//liberation memoire des images temporaires
free_grphic3d(imtmp);
free_grphic3d(imtref);
free_grphic3d(imtref2);
free_grphic3d(imtreca); 
free_grphic3d(imtres);
free_grphic3d(imtseg);



//affichage du temps de calcul total
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

free(param);
free(param_norm);
free(param_segmentation); 

return(1);  

}




/*******************************************
** 		  imx_convert_dmatrix3d_grphic3d() 
**
**  Conversion d'une matrice de double dans
**  une image grphic3d
**
** renvoie 0 si probl�me et 1 sinon
**  
********************************************/

int imx_convert_fmatrix3d_grphic3d(float ***mat,grphic3d* im)
{
int i,j,k;
int wdth=im->width,hght=im->height,dpth=im->depth;
double max,min,rcoeff;


/* recherche du min et du max */
min = HUGE_VAL; max	= -HUGE_VAL;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (mat[i][j][k]<min)
		min=mat[i][j][k];
	
	if (mat[i][j][k]>max)
		max=mat[i][j][k];
	}


imx_brukermax_3d(max, min, im);
rcoeff=im->rcoeff;

/* remplissage de l'image */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	im->mri[i][j][k]=(int)(mat[i][j][k]/rcoeff);
	}
	
return(1);
}

double gaussienneValue(double x,double mu,double sigma)
{
double	 prob;

if (sigma == 0)
	return 0;

if	((x < (mu-4.0*sigma)) || (x > (mu+4.0*sigma)))
	{
	return 0;
	}
else
	{		
	prob=( exp( -((x-mu)*(x-mu)) / (2.0*sigma*sigma) ) ) / ( sigma*sqrt( (double)(2.0*3.1415926539) ) );
	return prob;	
	}
}



double **alloc_matrice_2d_doubleT(int width, int height)
{
double	**t2;
int	i;
	
if	((t2=(double**)calloc((width),sizeof(double*))) == NULL)
	{PUT_WNDMSG( ERR_MEM_ALLOC);return NULL;}
	
if	((t2[0]=(double*)calloc((width*height),sizeof(double))) == NULL)
	{free(t2);PUT_WNDMSG( ERR_MEM_ALLOC);return NULL;}
	
for	(i=1 ; i<width ; i++) 
	t2[i]=t2[i-1]+height; 
	
return t2;
}


void free_matrice_2d_doubleT(double **D)
{
free(D[0]);
free(D);
}


int ChoixSeuil(vector<unsigned int> cumul, int min,double p1)
{
double mino;
int threshold;
double res1 = 0;
double res2 = 0;
double res = 0;
int entre = 0;
double pp1Opt;
//printf("on recherche sur l intervalle :%d %d\n",min,min+cumul.size()-1);

double total = 0;
for (int i=0;i<cumul.size();i++) //Le seuil est fixe a la valeur de i
	{
	total += cumul[i];
	}

for (int i=0;i<cumul.size()+1;i++) //Le seuil est fixe a la valeur de i // ce qui est egal a superieur ou egal a i est considere comme os...
	{
	//if ( (pp1<p1*1.3)&&(pp1>p1*0.5) )
		{
		res1=0;
		res2=0;
		res =0;
		for (int j=i;j<cumul.size();j++) //On compte combien de parties molles va apparaitre a tord dans la partie droite.
			{
			 double facteur = (1.-p1)*(double)histoPartieMolle[j+min-MINIMAL] /
		 	( (1.-p1)*(double)histoPartieMolle[j+min-MINIMAL]+p1*(double)histoOs[j+min-MINIMAL]);
			 res1 += (facteur * cumul[j]);
			}
		for (int j=0;j<i;j++) //On compte combien d os va apparaitre a tord dans la partie gauche.
			{
			 double facteur = p1*(double)histoOs[j+min-MINIMAL] /
		 		((1.-p1)*(double)histoPartieMolle[j+min-MINIMAL]+p1*(double)histoOs[j+min-MINIMAL]);
			 res2 += (facteur * cumul[j]);
			}
		res = res1+res2;
		if (entre == 0)
			{
			entre = 1;
			mino = res;			
			threshold = min + i;				
			}
		if (res<=mino)
			{
			mino = res;
			threshold = min + i;				
			}
		}
	}
//printf("FIN ALGO %f %f %d\n",p1,pp1Opt,threshold);
return threshold;
}
		 




void filtrage(double *y, long int nbpts)
{
  double sigma = 0.5;
  int wnd = 7;
  double filter[7];
  int swap = (wnd-1)/2;
  int i,l,p;
  
  for(i=0;i<wnd;i++)
  	filter[i]=exp(-1.0*((double)(i-swap))*((double)(i-swap))/((double)2.0*sigma*sigma));
	double norm = 0;
  for(i=0;i<wnd;i++)
    norm+=filter[i];
  for(i=0;i<wnd;i++)
    filter[i]/=norm;
  
  
  double *yres;


  /* Allocation du tableau de  yres */
  if( (yres = (double*)malloc(sizeof(double)*nbpts)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return;
  }


for (l=0;l<nbpts;l++)
	yres[l]=0.;
  
  for (l=0;l<nbpts;l++)
  for(p=-swap;p<=swap;p++)
      {
        if ( (l+p>=0)&&(l+p<nbpts))
        {
          yres[l]+=(double)y[l+p]*filter[p+swap];
  
        }
  
  	}
  
  for (l=0;l<nbpts;l++)
    y[l]=yres[l];


  free(yres);
}
void EstimeGauss(double *y, long int nbpts)
{
  double sigma = 0;
  double mu = 0; 
  double compteur = 0;
  int i;
  for(i=0;i<nbpts;i++)
  	{
	compteur+=y[i];
	mu+=(double)i*(double)y[i];
	sigma+=(double)i*(double)i*(double)y[i];;
	}

  mu/=compteur;
  sigma/=compteur;
  sigma-=(mu*mu);

  if (sigma<0)
  	sigma = 0;
  sigma=sqrt(sigma); 
  printf("SIGMA : %f MU : %f \n",sigma ,mu + MINIMAL );


  double compteur2=0;
  for(i=0;i<nbpts;i++)
  	{
	if (compteur2<(compteur/2))
		compteur2+=y[i];
	else 
		break;		
	}
  mu = i -0.5;

  vector<float> test;
  for(i=0;i<nbpts;i++)
  	{
	float res = i - mu;
	if (res<0)
		res= -res;
	for (int k=0;k<y[i];k++)
		{
		test.push_back(res);
		}
	}
   sort(test.begin(),test.end());
   sigma = test[test.size()/2]*1.4826;
printf("SIGMA : %f MU : %f \n",sigma ,mu + MINIMAL );

   for(i=0;i<nbpts;i++)
  	{
	y[i] = gaussienneValue(i,mu,sigma);
	}
   double sum = 0;
   for(i=0;i<nbpts;i++)
  	{
	sum += y[i];
	}
   if (sum!=0)
   	sum = 1./sum;

   for(i=0;i<nbpts;i++)
  	{
	y[i] = y[i]*sum;
	}
   
   double min = 1000;
   for(i=0;i<nbpts;i++)
  	{
	if ((y[i]!=0)&&(y[i]<min))
		min = y[i];
	}
min = min /2.;
   for(i=0;i<nbpts;i++)
  	{
	if (y[i]==0)
		{
		y[i]=min/2.;
		}
	}
}



/* CODE NICOLAS POUR LA SEGMENTATION */
  /****************************************************************************************************************/

  void MarqueCCObjet2(int x,int y,int z,char tampon[3][3][3]);

  /****************************************************************************************************************/

  void MarqueCCObjet2(int x,int y,int z,char tampon[3][3][3])
  {
    if ((x>=0)&&(y>=0)&&(z>=0)&&(x<3)&&(y<3)&&(z<3))
    {
      if (tampon[x][y][z]==2)
      {
	tampon[x][y][z]=-2;
	MarqueCCObjet2(x-1,y,z,tampon);
	MarqueCCObjet2(x+1,y,z,tampon);
	MarqueCCObjet2(x,y-1,z,tampon);
	MarqueCCObjet2(x,y+1,z,tampon);
	MarqueCCObjet2(x,y,z-1,tampon);
	MarqueCCObjet2(x,y,z+1,tampon);
      }  
    }  
  }  

  int CCObjet2(char I[3][3][3])
  {
    char tampon[3][3][3];
    for(int x=0;x<3;x++)
    {  
      for(int y=0;y<3;y++)
      {
	for(int z=0;z<3;z++)
	{
          if (I[x][y][z]==0) 
          {    
            tampon[x][y][z]=1;
          }      
          else 
          {    
            tampon[x][y][z]=2;
          }      
	}
      }    
    }    
    int nbCC=0;
    if ((I[0][1][1]!=0)&&(tampon[0][1][1]==2))
    {
      nbCC++;
      MarqueCCObjet2(0,1,1,tampon);
    }
    if ((I[1][0][1]!=0)&&(tampon[1][0][1]==2))
    {
      nbCC++;
      MarqueCCObjet2(1,0,1,tampon);
    }
    if ((I[1][1][0]!=0)&&(tampon[1][1][0]==2))
    {
      nbCC++;
      MarqueCCObjet2(1,1,0,tampon);
    }
    if ((I[2][1][1]!=0)&&(tampon[2][1][1]==2))
    {
      nbCC++;
      MarqueCCObjet2(2,1,1,tampon);
    }
    if ((I[1][2][1]!=0)&&(tampon[1][2][1]==2))
    {
      nbCC++;
      MarqueCCObjet2(1,2,1,tampon);
    }
    if ((I[1][1][2]!=0)&&(tampon[1][1][2]==2))
    {
      nbCC++;
      MarqueCCObjet2(1,1,2,tampon);
    }
    return nbCC;
  }

  /****************************************************************************************************************/

  void MarqueCCObjet(int x,int y,int z,char tampon[3][3][3])
  {
    if ((x>=0)&&(y>=0)&&(z>=0)&&(x<3)&&(y<3)&&(z<3))
    {
      if (tampon[x][y][z]==2)
      {
	tampon[x][y][z]=-2;
	for(int dx=-1;dx<=1;dx++)
	{    
          for(int dy=-1;dy<=1;dy++)
          {    
            for(int dz=-1;dz<=1;dz++)
            {      
              MarqueCCObjet(x+dx,y+dy,z+dz,tampon);
            }      
          }       
	}      
      }
    }  
  }

  int CCObjet(char I[3][3][3])
  {
    char tampon[3][3][3];
    for(int x=0;x<3;x++)
    {  
      for(int y=0;y<3;y++)
      {
	for(int z=0;z<3;z++)
	{
          if (I[x][y][z]==0) 
          {    
            tampon[x][y][z]=1;
          }      
          else 
          {    
            tampon[x][y][z]=2;
          }      
	}
      }    
    }    
    int nbCC=0;
    for(int x=0;x<3;x++)
    {  
      for(int y=0;y<3;y++)
      {  
	for(int z=0;z<3;z++)
	{      
          if ((I[x][y][z]!=0) && (tampon[x][y][z]==2))
          {
            nbCC++;
            MarqueCCObjet(x,y,z,tampon);
          }
	}    
      }    
    }    

    return nbCC;
  }

  /****************************************************************************************************************/

  int calcT1(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres)
  {
    char imX[3][3][3];
    char N26s[3][3][3];
    char imA1[3][3][3];
    
    for(int x=0;x<3;x++)
    {
      for(int y=0;y<3;y++)
      {
	for(int z=0;z<3;z++)
	{
	  imX[x][y][z]=0;
	  N26s[x][y][z]=1;
	  imA1[x][y][z]=0;
	}
      }
    }
    for(int dx=-1;dx<=1;dx++)
    {
      for(int dy=-1;dy<=1;dy++)
      {
	for(int dz=-1;dz<=1;dz++)
	{
          if ((x+dx>=0)&&(y+dy>=0)&&(z+dz>=0)&&(x+dx<256)&&(y+dy<256)&&(z+dz<256))
          {
            if (imres->mri[x+dx][y+dy][z+dz]==1)
            {
              imX[1+dx][1+dy][1+dz]=1;
            }
          }
          if ((dx==0)&&(dy==0)&&(dz==0))
          {
            N26s[1+dx][1+dy][1+dz]=0;
          }
          imA1[1+dx][1+dy][1+dz]=N26s[1+dx][1+dy][1+dz]*imX[1+dx][1+dy][1+dz];
	}
      }
    }  
    return CCObjet(imA1);
  }

  /****************************************************************************************************************/
  int calcT2(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres)
  {
    char imXb[3][3][3];
    char N18s[3][3][3];
    char imA2[3][3][3];
    
    for(int x=0;x<3;x++)
    {
      for(int y=0;y<3;y++)
      {
	for(int z=0;z<3;z++)
	{
	  imXb[x][y][z]=0;
	  N18s[x][y][z]=0;
	  imA2[x][y][z]=0;
	}
      }
    }
    
    for(int dx=-1;dx<=1;dx++)
    {
      for(int dy=-1;dy<=1;dy++)
      {
	for(int dz=-1;dz<=1;dz++)
	{
          if ((dx!=0)||(dy!=0)||(dz!=0))
          {
            if ((dx==0)||(dy==0)||(dz==0))
            {
              N18s[1+dx][1+dy][1+dz]=1;
            }
            if ((x+dx>=0)&&(y+dy>=0)&&(z+dz>=0)&&(x+dx<256)&&(y+dy<256)&&(z+dz<256))
            {
              if (imres->mri[x+dx][y+dy][z+dz]==0)
              {
        	imXb[1+dx][1+dy][1+dz]=1;
              }
            }
            else
            {
              imXb[1+dx][1+dy][1+dz]=1;
            }
          }
          imA2[1+dx][1+dy][1+dz]=N18s[1+dx][1+dy][1+dz]*imXb[1+dx][1+dy][1+dz];
	}
      }
    }  
    return CCObjet2(imA2);
  }

  /****************************************************************************************************************/
  bool simple(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres)
  {
    char imX[3][3][3];
    char imXb[3][3][3];
    char N26s[3][3][3];
    char N18s[3][3][3];
    char imA1[3][3][3];
    char imA2[3][3][3];
    
    for(int x=0;x<3;x++)
    {
      for(int y=0;y<3;y++)
      {
	for(int z=0;z<3;z++)
	{
	  imX[x][y][z]=0;
	  imXb[x][y][z]=0;
	  N26s[x][y][z]=0;
	  N18s[x][y][z]=0;
	  imA1[x][y][z]=0;
	  imA2[x][y][z]=0;
	}
      }
    }
    
    for(int dx=-1;dx<=1;dx++)
    {
      for(int dy=-1;dy<=1;dy++)
      {
	for(int dz=-1;dz<=1;dz++)
	{
          if ((dx==0)&&(dy==0)&&(dz==0))
          {
            imX[1+dx][1+dy][1+dz]=1;
          }
          else
          {
            N26s[1+dx][1+dy][1+dz]=1;
            if ((dx==0)||(dy==0)||(dz==0))
            {
              N18s[1+dx][1+dy][1+dz]=1;
            }
            if ((x+dx>=0)&&(y+dy>=0)&&(z+dz>=0)&&(x+dx<imres->width)&&(y+dy<imres->height)&&(z+dz<imres->depth))
            {
              if (imres->mri[x+dx][y+dy][z+dz]!=0)
              {
        	imX[1+dx][1+dy][1+dz]=1;
              }
              else
              {
        	imXb[1+dx][1+dy][1+dz]=1;
              }
            }
            else
            {
              imXb[1+dx][1+dy][1+dz]=1;
            }
          }
          imA1[1+dx][1+dy][1+dz]=N26s[1+dx][1+dy][1+dz]*imX[1+dx][1+dy][1+dz];
          imA2[1+dx][1+dy][1+dz]=N18s[1+dx][1+dy][1+dz]*imXb[1+dx][1+dy][1+dz];
	}
      }
    }   
 
    return  ((CCObjet2(imA2)==1)&&(CCObjet(imA1)==1));

  }

  /****************************************************************************************************************/

  bool simple2(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres)
  {
    char imX[3][3][3];
    char imXb[3][3][3];
    char N26s[3][3][3];
    char N18s[3][3][3];
    char imA1[3][3][3];
    char imA2[3][3][3];
    
    for(int x=0;x<3;x++)
    {
      for(int y=0;y<3;y++)
      {
	for(int z=0;z<3;z++)
	{
	  imX[x][y][z]=0;
	  imXb[x][y][z]=0;
	  N26s[x][y][z]=0;
	  N18s[x][y][z]=0;
	  imA1[x][y][z]=0;
	  imA2[x][y][z]=0;
	}
      }
    }
    
    for(int dx=-1;dx<=1;dx++)
    {
      for(int dy=-1;dy<=1;dy++)
      {
	for(int dz=-1;dz<=1;dz++)
	{
          if ((dx==0)&&(dy==0)&&(dz==0))
          {
            imX[1+dx][1+dy][1+dz]=1;
          }
          else
          {
            N26s[1+dx][1+dy][1+dz]=1;
            if ((dx==0)||(dy==0)||(dz==0))
            {
              N18s[1+dx][1+dy][1+dz]=1;
            }
            if ((x+dx>=0)&&(y+dy>=0)&&(z+dz>=0)&&(x+dx<imres->width)&&(y+dy<imres->height)&&(z+dz<imres->depth))
            {
              if ((imres->mri[x+dx][y+dy][z+dz]==1)||(imres->mri[x+dx][y+dy][z+dz]==-1))
              {
        	imX[1+dx][1+dy][1+dz]=1;
              }
              else
              {
        	imXb[1+dx][1+dy][1+dz]=1;
              }
            }
            else
            {
              imXb[1+dx][1+dy][1+dz]=1;
            }
          }
          imA1[1+dx][1+dy][1+dz]=N26s[1+dx][1+dy][1+dz]*imX[1+dx][1+dy][1+dz];
          imA2[1+dx][1+dy][1+dz]=N18s[1+dx][1+dy][1+dz]*imXb[1+dx][1+dy][1+dz];
	}
      }
    }   
 
    return  ((CCObjet2(imA2)==1)&&(CCObjet(imA1)==1));

  }

  /****************************************************************************************************************/

  void reduction_homotopique(grphic3d *imin,short***imth,short***imtb,grphic3d *imres,int x,int y,int z){
    
    if (imres->mri[x][y][z]==1)
    {
      if ((imin->mri[x][y][z]>imth[x][y][z])||(imin->mri[x][y][z]<imtb[x][y][z]))
      {
	if (simple(x,y,z,imres))
	{
	  imres->mri[x][y][z]=0;
          for(int dx=((x==0)?0:-1);dx<=((x==255)?0:1);dx++)
	  {
            for(int dy=((y==0)?0:-1);dy<=((y==255)?0:1);dy++)
            {
              for(int dz=((z==0)?0:-1); dz<=((z==255)?0:1);dz++)
              {
		reduction_homotopique(imin,imth,imtb,imres,x+dx,y+dy,z+dz);
	      }
	    }
	  }
	}
      }
    }
    return;
  }

  /****************************************************************************************************************/

  void generation_tunnels(grphic3d *imin,short***imth,short***imtb,grphic3d *imres)
  {
    /*
    imth : image seuil haut
    imtb : image seuil bas
    imin : image initiale
    imres : segmentation (binaire !) et image resultat

    Dans imres :
    0 = non segment�
    1 = segment�
    */
    for(int x=0;x<256;x++)
    {
      for(int y=0;y<256;y++)
      {
        for(int z=0;z<256;z++)
        {
	  if (imres->mri[x][y][z]==1)
	  {
	    if ((imin->mri[x][y][z]>imth[x][y][z])||(imin->mri[x][y][z]<imtb[x][y][z]))
	    {
              if ((calcT1(x,y,z,imres)==1)&&(calcT2(x,y,z,imres)==2))
	      {
		imres->mri[x][y][z]=0;
		for(int dx=((x==0)?0:-1);dx<=((x==255)?0:1);dx++)
		{
        	  for(int dy=((y==0)?0:-1);dy<=((y==255)?0:1);dy++)
        	  {
        	    for(int dz=((z==0)?0:-1); dz<=((z==255)?0:1);dz++)
        	    {
		      reduction_homotopique(imin,imth,imtb,imres,x+dx,y+dy,z+dz);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    return;
  }

  /****************************************************************************************************************/

  void separation_composantes(grphic3d *imin,short***imth,short***imtb,grphic3d *imres)
  {
    /*
    imth : image seuil haut
    imtb : image seuil bas
    imin : image initiale
    imres : segmentation (binaire !) et image resultat

    Dans imres :
    0 = non segment�
    1 = segment�
    */
    for(int x=0;x<256;x++)
    {
      for(int y=0;y<256;y++)
      {
        for(int z=0;z<256;z++)
        {
	  if (imres->mri[x][y][z]==1)
	  {
	    if ((imin->mri[x][y][z]>imth[x][y][z])||(imin->mri[x][y][z]<imtb[x][y][z]))
	    {
              if ((calcT1(x,y,z,imres)==2)&&(calcT2(x,y,z,imres)==1))
	      {
		imres->mri[x][y][z]=0;
		for(int dx=((x==0)?0:-1);dx<=((x==255)?0:1);dx++)
		{
        	  for(int dy=((y==0)?0:-1);dy<=((y==255)?0:1);dy++)
        	  {
        	    for(int dz=((z==0)?0:-1); dz<=((z==255)?0:1);dz++)
        	    {
		      reduction_homotopique(imin,imth,imtb,imres,x+dx,y+dy,z+dz);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    return;
  }
  














  /****************************************************************************************************************/

  void segmentationNicolas(grphic3d*imin,short***imth,short***imtb,grphic3d*imres)
  {
    /*
    En entree
    imin : image a segmenter
    imth : image seuil haut
    imtb : image seuil bas
    imres : image resultat
    
    Dans imres :
    0 = non segment�
    1 = segment�
    2 = segment� au bord de l'objet
    3 = non segment� au bord de l'objet
    */

    pt** F; 
    int * N;
    int * B;
    int * B2;

    imres->mri[0][0][0]=1;
    int Ro =9;

    for(int le=2;le<=256;le*=2)
    {
      printf("Reduction (resolution %d) : retrait de points [1,1] (preservation de la topologie)\n", --Ro);
     
      for(int x=(le/2)-1;x>=0;x--)
      {
	for(int y=(le/2)-1;y>=0;y--)
	{
          for(int z=(le/2)-1;z>=0;z--)
          {
	    imres->mri[2*x][2*y][2*z]=imres->mri[x][y][z];
	    imres->mri[2*x][2*y][2*z+1]=imres->mri[x][y][z];
	    imres->mri[2*x][2*y+1][2*z]=imres->mri[x][y][z];
	    imres->mri[2*x][2*y+1][2*z+1]=imres->mri[x][y][z];
	    imres->mri[2*x+1][2*y][2*z]=imres->mri[x][y][z];
	    imres->mri[2*x+1][2*y][2*z+1]=imres->mri[x][y][z];
	    imres->mri[2*x+1][2*y+1][2*z]=imres->mri[x][y][z];
	    imres->mri[2*x+1][2*y+1][2*z+1]=imres->mri[x][y][z];
	  }
	}
      }

      for(int x=0;x<le;x++)
      {
	for(int y=0;y<le;y++)
	{
          for(int z=0;z<le;z++)
          {
	    if (imres->mri[x][y][z]!=0)
	    {
	      short mv=32767;
	      for(int dx=0;dx<=256/le-1;dx++)
	      {
		for(int dy=0;dy<=256/le-1;dy++)
		{
        	  for(int dz=0;dz<=256/le-1;dz++)
        	  {
		    if (imtb[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] > imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz])
		    {
		      if (imtb[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz]-imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz]<=32765)
		      {
			if (2 + imtb[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] - imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] < mv)
			{
		          mv = 2 + imtb[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] - imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz];
			}
		      }
		    }
		    else if (imth[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] < imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz])
		    {
		      if (imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] - imth[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz]<=32765)
		      {
			if (2 + imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] - imth[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] < mv)
			{
		          mv = 2 + imin->mri[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz] - imth[x*(256/le)+dx][y*(256/le)+dy][z*(256/le)+dz];
			}
		      }
		    }
		    else
		    {
		      if (mv > 1)
		      {
		        mv = 1;
		      }
		    }
		  }
		}
	      }	  
	      imres->mri[x][y][z]=mv;
	    }
	  }
	}
      }
      
      short max=1;
      
      for(int x=0;x<le;x++)
      {
	for(int y=0;y<le;y++)
	{
          for(int z=0;z<le;z++)
          {
	    if (imres->mri[x][y][z] > max)
	    {
              max = imres->mri[x][y][z];	  
	    }
          }
	}
      }

      short idx;

      N = (int *)malloc((max+1) * sizeof(int));
      B = (int *)malloc((max+1) * sizeof(int));
      B2 = (int *)malloc((max+1) * sizeof(int));
      for(int i=0;i<=max;i++)
      {
	N[i]=0;
	B[i]=0;
	B2[i]=0;
      }

      for(int x=0;x<le;x++)
      {
	for(int y=0;y<le;y++)
	{
          for(int z=0;z<le;z++)
          {
	    N[imres->mri[x][y][z]]++;
          }
	}
      }

      F = (pt**)malloc((max+1) * sizeof(pt*));
      for(int i=0;i<=max;i++)
      {
	F[i]=(pt*)malloc(N[i]*sizeof(pt));
      }


      for(int x=0;x<le;x++)
      {
	for(int y=0;y<le;y++)
	{
          for(int z=0;z<le;z++)
          {
	    if (imres->mri[x][y][z] != 0)
	    {
	      if ((x==0)||(y==0)||(z==0)||(x==le-1)||(y==le-1)||(z==le-1)){
		F[imres->mri[x][y][z]][(B[imres->mri[x][y][z]]+B2[imres->mri[x][y][z]])%(N[imres->mri[x][y][z]])].x = x;
		F[imres->mri[x][y][z]][(B[imres->mri[x][y][z]]+B2[imres->mri[x][y][z]])%(N[imres->mri[x][y][z]])].y = y;
		F[imres->mri[x][y][z]][(B[imres->mri[x][y][z]]+B2[imres->mri[x][y][z]])%(N[imres->mri[x][y][z]])].z = z;
		B2[imres->mri[x][y][z]]++;
		imres->mri[x][y][z]*=-1; 
	      }
	      else
	      {
		bool bord=false;
		for(int dx=-1;dx<=1;dx++)
		{
		  for(int dy=-1;dy<=1;dy++)
		  {
        	    for(int dz=-1;dz<=1;dz++)
        	    {
		      bord=(bord || (imres->mri[x+dx][y+dy][z+dz]==0));
		    }
        	  }
		}
		if (bord)
		{
		  F[imres->mri[x][y][z]][(B[imres->mri[x][y][z]]+B2[imres->mri[x][y][z]])%(N[imres->mri[x][y][z]])].x = x;
		  F[imres->mri[x][y][z]][(B[imres->mri[x][y][z]]+B2[imres->mri[x][y][z]])%(N[imres->mri[x][y][z]])].y = y;
		  F[imres->mri[x][y][z]][(B[imres->mri[x][y][z]]+B2[imres->mri[x][y][z]])%(N[imres->mri[x][y][z]])].z = z;
		  B2[imres->mri[x][y][z]]++;
		  imres->mri[x][y][z]*=-1; 
		}
	      }
	    }
          }
	}
      }

      idx = max;  

      while (idx>1){
	while ((B2[idx]==0)&&(idx!=1)){
          idx--;
	}
	if (idx!=1)
	{
	  unsigned int x = F[idx][B[idx]].x;
	  unsigned int y = F[idx][B[idx]].y;
	  unsigned int z = F[idx][B[idx]].z;
	  B2[idx]--;
	  B[idx]=(B[idx]+1)%(N[idx]);
          if (simple(x,y,z,imres))
	  {
	    imres->mri[x][y][z]=0;
	    for(int dx=((x==0)?0:-1);dx<=((x==le-1)?0:1);dx++)
	    {
              for(int dy=((y==0)?0:-1);dy<=((y==le-1)?0:1);dy++)
              {
        	for(int dz=((z==0)?0:-1); dz<=((z==le-1)?0:1);dz++)
        	{
		  if (imres->mri[x+dx][y+dy][z+dz]>0)
		  {
		    if (imres->mri[x+dx][y+dy][z+dz]>idx)
		    {
	              idx=imres->mri[x+dx][y+dy][z+dz];
		    }
		    F[imres->mri[x+dx][y+dy][z+dz]][(B[imres->mri[x+dx][y+dy][z+dz]]+B2[imres->mri[x+dx][y+dy][z+dz]])%(N[imres->mri[x+dx][y+dy][z+dz]])].x=x+dx;
		    F[imres->mri[x+dx][y+dy][z+dz]][(B[imres->mri[x+dx][y+dy][z+dz]]+B2[imres->mri[x+dx][y+dy][z+dz]])%(N[imres->mri[x+dx][y+dy][z+dz]])].y=y+dy;
		    F[imres->mri[x+dx][y+dy][z+dz]][(B[imres->mri[x+dx][y+dy][z+dz]]+B2[imres->mri[x+dx][y+dy][z+dz]])%(N[imres->mri[x+dx][y+dy][z+dz]])].z=z+dz;
                    B2[imres->mri[x+dx][y+dy][z+dz]]++;
		    imres->mri[x+dx][y+dy][z+dz]*=-1;
		  }
        	}
              }
	    }

	  }
	  else
	  {
	    imres->mri[x][y][z]*=-1;
	  }
	}
      }

      for(unsigned int x=0;x<le;x++)
      {
	for(unsigned int y=0;y<le;y++)
	{
          for(unsigned int z=0;z<le;z++)
          {
	    if (imres->mri[x][y][z]!=0)
	    {
              imres->mri[x][y][z]=1;
	    }
	  }
	}
      }

      for(int i=0;i<=max;i++)
      {
	free((void*)F[i]);
      }
      free((void*)F);
      free((void*)B);
      free((void*)B2);
      free((void*)N);
    }

    /* bloc a commenter pour ne pas modifier la topologie initiale */
    /*
    printf("Reduction (resolution 1) : retrait de points [2,1] puis [1,1] (separation de composantes)\n");
    separation_composantes(imin,imth,imtb,imres);
    printf("Reduction (resolution 1) : retrait de points [1,2] puis [1,1] (generation de tunnels)\n");
    generation_tunnels(imin,imth,imtb,imres);
    */

    printf("Etiquetage des points\n");
    for(unsigned int x=0;x<256;x++)
    {
      for(unsigned int y=0;y<256;y++)
      {
        for(unsigned int z=0;z<256;z++)
        {
	  if (imres->mri[x][y][z]!=0)
	  {
            if ((x==0)||(y==0)||(z==0)||(x==255)||(y==255)||(z==255))
	    {
              imres->mri[x][y][z]=2;
	    }
	    else
	    {
	      bool bord=false;
	      for(int dx=-1;dx<=1;dx++)
	      {
		for(int dy=-1;dy<=1;dy++)
		{
        	  for(int dz=-1;dz<=1;dz++)
        	  {
		    bord=(bord || (imres->mri[x+dx][y+dy][z+dz]==0));
		  }
        	}
	      }
	      if (bord)
	      {
                imres->mri[x][y][z]=2;
	      }
	      else
	      {
                imres->mri[x][y][z]=1;
	      }
	    }
	  }
	}
      }
    }

    for(unsigned int x=0;x<256;x++)
    {
      for(unsigned int y=0;y<256;y++)
      {
        for(unsigned int z=0;z<256;z++)
        {
	  if (imres->mri[x][y][z]==2)
	  {
	    for(int dx=((x==0)?0:-1);dx<=((x==255)?0:1);dx++)
	    {
              for(int dy=((y==0)?0:-1);dy<=((y==255)?0:1);dy++)
              {
        	for(int dz=((z==0)?0:-1); dz<=((z==255)?0:1);dz++)
        	{
		  if (imres->mri[x+dx][y+dy][z+dz]==0)
		  {
                    imres->mri[x+dx][y+dy][z+dz]=3;
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    return;
  }

  
 //

short int ***alloc_matrix_3dSHORT(int nb_lig, int nb_col, int nb_dep)
{
  short int ***m;
  short int *mt;
  int i, j;

  m = (short int ***)malloc(sizeof(m[0])*nb_lig);

  if (! m)
  {
    return NULL;
  }

  m[0] = (short int **)malloc(sizeof(m[0][0])*nb_lig*nb_col);

  if (! m[0])
  {
    free(m);
    return NULL;
  }

  for(i=1; i<nb_lig; i++)
     m[i] = m[i-1]+nb_col;

  mt = (short int *)malloc(sizeof(m[0][0][0])*nb_lig*nb_col*nb_dep);

  if (! mt)
  {
    free(m[0]);
    free(m);
    return NULL;
  }

  for (i=0; i<nb_lig; i++)
    for (j=0; j<nb_col; j++)
    {
      m[i][j] = mt;
      mt += nb_dep;
    }

  return m;
}

/*
  free_matrix_3d()

  liberation memoire d'une matrice float 3D allouee par alloc_matrix_3d()
*/
int free_matrix_3dSHORT(short int ***m)
{
  free(m[0][0]);
  free(m[0]);
  free(m);

  return(0);
}
void segmentationNicolas2(short***imth,short***imtb,grphic3d*imin,grphic3d*impr,grphic3d*imres,
            Points_3dInt*pointBG,Points_3dInt*pointHD,t_ptrlist*add,t_ptrlist*remove)
  {
    /*
    imth : image seuil haut
    imtb : image seuil haut
    imin : image a segmenter
    impr : segmentation de reference
    imres : image resultat
    pointBG : point de coord. minimales ou les seuils ont ete modifies
    pointHD : point de coord. maximales ou les seuils ont ete modifies
    Add : liste des coord. des points ajoutes
    Add : liste des coord. des points supprimes

    Dans imres et impr :
    0 = non segment�
    1 = segment�
    2 = segment� au bord de l'objet
    3 = non segment� au bord de l'objet
    */

    pt** F; 
    int * N;
    int * B;
    int * B2;

    short ***dist;

    dist=(short***)malloc(256*sizeof(short**));
    for(int i=0;i<256;i++)
    {
      dist[i]=(short**)malloc(256*sizeof(short*));
      for(int j=0;j<256;j++)
      {
        dist[i][j]=(short*)malloc(256*sizeof(short));
      }
    }

    for(int x=0;x<256;x++)
    {
      for(int y=0;y<256;y++)
      {
        for(int z=0;z<256;z++)
        {
	  short mv;
	  
	  if ((impr->mri[x][y][z]==1)||(impr->mri[x][y][z]==2))
	  {
	    if (imtb[x][y][z] > imin->mri[x][y][z])
	    {
	      if (imtb[x][y][z]-imin->mri[x][y][z]<=32765)
	      {
		mv = 2 + imtb[x][y][z] - imin->mri[x][y][z];
	      }
	    }
	    else if (imth[x][y][z] < imin->mri[x][y][z])
	    {
	      if (imin->mri[x][y][z] - imth[x][y][z]<=32765)
	      {
                mv = 2 + imin->mri[x][y][z] - imth[x][y][z];
	      }
	    }
	    else
	    {
	      mv=1;
	    }
	  }
	  else
	  {
	    if ((imtb[x][y][z] <= imin->mri[x][y][z])&&(imth[x][y][z] >= imin->mri[x][y][z]))
	    {
	      if ((imth[x][y][z]-imin->mri[x][y][z]<=32765)&&(imin->mri[x][y][z]-imtb[x][y][z]<=32765))
	      {
	        if (imth[x][y][z]-imin->mri[x][y][z] > imin->mri[x][y][z]-imtb[x][y][z]<=32765)
		{
		  mv = 2 + imth[x][y][z] - imin->mri[x][y][z];
		}
		else
		{
		  mv = 2 + imin->mri[x][y][z] - imtb[x][y][z];
		}
	      }
	    }
 	    else
	    {
	      mv=1;
	    }
	  }
	  dist[x][y][z]=mv;
	}
      }
    }

    short max=1;

    for(int x=0;x<256;x++)
    {
      for(int y=0;y<256;y++)
      {
        for(int z=0;z<256;z++)
        {
	  if (dist[x][y][z] > max)
	  {
            max = dist[x][y][z];	  
	  }
        }
      }
    }

    for(int x=0;x<256;x++)
    {
      for(int y=0;y<256;y++)
      {
        for(int z=0;z<256;z++)
        {
	  if ((impr->mri[x][y][z]==0)||(impr->mri[x][y][z]==3))
	  {
	    imres->mri[x][y][z]=2;
	  }
	  else
	  {
	    imres->mri[x][y][z]=1;
	  }
        }
      }
    }

    N = (int *)malloc((max+1) * sizeof(int));
    B = (int *)malloc((max+1) * sizeof(int));
    B2 = (int *)malloc((max+1) * sizeof(int));

    for(int i=0;i<=max;i++)
    {
      N[i]=0;
      B[i]=0;
      B2[i]=0;
    }
    for(int x=0;x<256;x++)
    {
      for(int y=0;y<256;y++)
      {
        for(int z=0;z<256;z++)
        {
	  N[dist[x][y][z]]++;
        }
      }
    }



    F = (pt**)malloc((max+1) * sizeof(pt*));
    for(int i=0;i<=max;i++)
    {
      F[i]=(pt*)malloc(N[i]*sizeof(pt));
    }

    for(int x=pointBG->x;x<=pointHD->x;x++)
    {
      for(int y=pointBG->y;y<=pointHD->y;y++)
      {
        for(int z=pointBG->z;z<=pointHD->z;z++)
        {
	  if ((impr->mri[x][y][z] == 2)||(impr->mri[x][y][z] == 3))
	  {
	    F[dist[x][y][z]][(B[dist[x][y][z]]+B2[dist[x][y][z]])%(N[dist[x][y][z]])].x = x;
	    F[dist[x][y][z]][(B[dist[x][y][z]]+B2[dist[x][y][z]])%(N[dist[x][y][z]])].y = y;
	    F[dist[x][y][z]][(B[dist[x][y][z]]+B2[dist[x][y][z]])%(N[dist[x][y][z]])].z = z;
	    B2[dist[x][y][z]]++;
	    imres->mri[x][y][z]*=-1; 
	  }
        }
      }
    }

    short idx = max;  
    while (idx>1){
      while ((B2[idx]==0)&&(idx!=1)){
        idx--;
      }
      if (idx!=1)
      {
	unsigned int x = F[idx][B[idx]].x;
	unsigned int y = F[idx][B[idx]].y;
	unsigned int z = F[idx][B[idx]].z;
	B2[idx]--;
	B[idx]=(B[idx]+1)%(N[idx]);
        if (simple2(x,y,z,imres))
	{
	  if (imres->mri[x][y][z]==-2)
	  {
            imres->mri[x][y][z]=-1;
//	    Points_3dInt * point = (Points_3dInt *)malloc(sizeof(Points_3dInt));
//	    point->x = x;
//            point->y = y;
//            point->z = z;
//	    list_push(add,point); 
	  }
	  else
	  {
	    imres->mri[x][y][z]=-2;
//	    Points_3dInt * point = (Points_3dInt *)malloc(sizeof(Points_3dInt));
//	    point->x = x;
////            point->y = y;
//          point->z = z;
//	    list_push(remove,point); 
	  }

	  for(int dx=((x==0)?0:-1);dx<=((x==255)?0:1);dx++)
	  {
            for(int dy=((y==0)?0:-1);dy<=((y==255)?0:1);dy++)
            {
              for(int dz=((z==0)?0:-1); dz<=((z==255)?0:1);dz++)
              {
		if (imres->mri[x+dx][y+dy][z+dz]>0)
		{
                  if (dist[x+dx][y+dy][z+dz]>1)
                  {
                    if (dist[x+dx][y+dy][z+dz]>idx)
                    {
                      idx=dist[x+dx][y+dy][z+dz];
                    }
                    F[dist[x+dx][y+dy][z+dz]][(B[dist[x+dx][y+dy][z+dz]]+B2[dist[x+dx][y+dy][z+dz]])%(N[dist[x+dx][y+dy][z+dz]])].x=x+dx;
                    F[dist[x+dx][y+dy][z+dz]][(B[dist[x+dx][y+dy][z+dz]]+B2[dist[x+dx][y+dy][z+dz]])%(N[dist[x+dx][y+dy][z+dz]])].y=y+dy;
                    F[dist[x+dx][y+dy][z+dz]][(B[dist[x+dx][y+dy][z+dz]]+B2[dist[x+dx][y+dy][z+dz]])%(N[dist[x+dx][y+dy][z+dz]])].z=z+dz;
                    B2[dist[x+dx][y+dy][z+dz]]++;
                    imres->mri[x+dx][y+dy][z+dz]*=-1;
                  }
		}
              }
            }
	  }
          imres->mri[x][y][z]*=-1;
          dist[x][y][z]=1;
	}
	else
	{
	  imres->mri[x][y][z]*=-1;
	}
      }
    }

    printf("Etiquetage des points\n");
    for(unsigned int x=0;x<256;x++)
    {
      for(unsigned int y=0;y<256;y++)
      {
        for(unsigned int z=0;z<256;z++)
        {
          if ((imres->mri[x][y][z]==1)||(imres->mri[x][y][z]==-1))
          {
            imres->mri[x][y][z]=1;
          }
          else
          {
            imres->mri[x][y][z]=0;
          }
        }
      }
    }
    for(unsigned int x=0;x<256;x++)
    {
      for(unsigned int y=0;y<256;y++)
      {
        for(unsigned int z=0;z<256;z++)
        {
	  if (imres->mri[x][y][z]!=0)
	  {
            if ((x==0)||(y==0)||(z==0)||(x==255)||(y==255)||(z==255))
	    {
              imres->mri[x][y][z]=2;
	    }
	    else
	    {
	      bool bord=false;
	      for(int dx=-1;dx<=1;dx++)
	      {
		for(int dy=-1;dy<=1;dy++)
		{
        	  for(int dz=-1;dz<=1;dz++)
        	  {
		    bord=(bord || (imres->mri[x+dx][y+dy][z+dz]==0));
		  }
        	}
	      }
	      if (bord)
	      {
                imres->mri[x][y][z]=2;
	      }
	      else
	      {
                imres->mri[x][y][z]=1;
	      }
	    }
	  }
	}
      }
    }

    for(unsigned int x=0;x<256;x++)
    {
      for(unsigned int y=0;y<256;y++)
      {
        for(unsigned int z=0;z<256;z++)
        {
	  if (imres->mri[x][y][z]==2)
	  {
	    for(int dx=((x==0)?0:-1);dx<=((x==255)?0:1);dx++)
	    {
              for(int dy=((y==0)?0:-1);dy<=((y==255)?0:1);dy++)
              {
        	for(int dz=((z==0)?0:-1); dz<=((z==255)?0:1);dz++)
        	{
		  if (imres->mri[x+dx][y+dy][z+dz]==0)
		  {
                    imres->mri[x+dx][y+dy][z+dz]=3;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    for(int i=0;i<256;i++)
    {
      for(int j=0;j<256;j++)
      {
        free((void*)dist[i][j]);
      }
      free((void*)dist[i]);
    }

    free((void*)dist);
free(N);
free(B);
free(B2);
for(int i=0;i<=max;i++)
    {
      free(F[i]);
    }
free(F);
    return;
  }


//





double FunctionCoutZZ(double *moyenne	              ,
		      double *param                   ,
		      int nb_param		      ,
		      grphic3d * ImageTemplateOsDilate,
		      float *** ImageSeuilBas	      ,
		      double seuilBas                 ,		    
		      int ***masque_param	      )
{
//Determination de la carte de seuil
TranslateBase3dToWholeImage(ImageSeuilBas ,
			    param         ,
			    1		  ,	 
			    nb_param      ,
			    seuilBas     );
//Calcul de la fonction de cout

int resol = BASE3D.resol;
int topD  = (int)pow(2.0,1.0*resol)-1;
double res = 0;
for (int topi=0; topi<topD; topi++)
for (int topj=0; topj<topD; topj++)
for (int topk=0; topk<topD; topk++)
	{
	if (masque_param[topi][topj][topk])		
		{
		int *x0,*x1,*y0,*y1,*z0,*z1,l;
		Points_3dInt pointBG       ; 
		Points_3dInt pointHD       ;
				   
		x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
		l=(int)TOP_conv_ind(topi,topj,topk,nb_param)/3;
		pointBG.x=x0[l];pointHD.x=x1[l];
		pointBG.y=y0[l];pointHD.y=y1[l];
		pointBG.z=z0[l];pointHD.z=z1[l];	
		
		int li = (int)TOP_conv_ind(topi,topj,topk,nb_param);			
		
		for (int i=pointBG.x;i<pointHD.x;i++)
		for (int j=pointBG.y;j<pointHD.y;j++)
		for (int k=pointBG.z;k<pointHD.z;k++)
			{
			if (ImageTemplateOsDilate->mri[i][j][k])
				res += ( (ImageSeuilBas[i][j][k] - moyenne[li] )
					*(ImageSeuilBas[i][j][k]- moyenne[li]));			       
			}																									
		}			   
	}
return res;
}			


double UpdateFunctionCout(int topi,int topj,int topk,
		    double *param	            ,
		    double *moyenne,
		    int nb_param		    ,
		    grphic3d * ImageTemplateOsDilate,
		    float *** ImageSeuilBas	    ,
		    double seuilBas                 ,		    
		    int ***masque_param	            )
{

TranslateBase3dToImage(ImageSeuilBas      ,
		       param              ,
		       1	  	  ,	 
		       nb_param           ,			    
		       topi               ,
		       topj               ,
		       topk               ,
		       seuilBas	 	 );
		       
//Calcul de la fonction de cout

int resol = BASE3D.resol;
int topD  = (int)pow(2.0,1.0*resol)-1;
double res = 0.;



for (int ii=MAXI(topi-1,0);ii<=MINI(topi+1,topD-1);ii++)
for (int jj=MAXI(topj-1,0);jj<=MINI(topj+1,topD-1);jj++)
for (int kk=MAXI(topk-1,0);kk<=MINI(topk+1,topD-1);kk++)
	{
	int *x0,*x1,*y0,*y1,*z0,*z1,l;
	Points_3dInt pointBG       ; 
	Points_3dInt pointHD       ;
				   
	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
	l=(int)TOP_conv_ind(ii,jj,kk,nb_param)/3;
	pointBG.x=x0[l];pointHD.x=x1[l];
	pointBG.y=y0[l];pointHD.y=y1[l];
	pointBG.z=z0[l];pointHD.z=z1[l];	
		
	int li = (int)TOP_conv_ind(ii,jj,kk,nb_param);			
		
	for (int i=pointBG.x;i<pointHD.x;i++)
	for (int j=pointBG.y;j<pointHD.y;j++)
	for (int k=pointBG.z;k<pointHD.z;k++)
		{
		if (ImageTemplateOsDilate->mri[i][j][k])
			res += ( (ImageSeuilBas[i][j][k]- moyenne[li])
			*(ImageSeuilBas[i][j][k]- moyenne[li]));			       
		}
	}																									
return res;
}



int OptimisationLigne(int topi	     	         ,
			 int topj	     	         ,
			 int topk	     	         ,
			 double *param	            	 ,
		    	 double *moyenne		 ,
			 int nb_param		         , 
		    	 grphic3d * ImageTemplateOsDilate,
		    	 float *** ImageSeuilBas	 ,
		    	 double seuilBas                 ,		    
		    	 int ***masque_param	         ,
			 double pas                      ) 
{
double E,Eopt;
int i,j,jMin;
double paramCourant;

int l = TOP_conv_ind(topi,topj,topk,nb_param);
paramCourant = param[l];

//Recherche sur l intervalle parametreCourant + infty
jMin=0;
Eopt =  UpdateFunctionCout(topi,topj,topk,param,moyenne,nb_param,ImageTemplateOsDilate,ImageSeuilBas,seuilBas,masque_param);
	
for (j=-1;j<=1;j++)
	{
	if (j!=0)
		{
		param[l] = paramCourant + (double)j * pas;;
		E =  UpdateFunctionCout(topi,topj,topk,param,moyenne,nb_param,ImageTemplateOsDilate,ImageSeuilBas,seuilBas,masque_param);
		if (E<Eopt)
			{
			Eopt = E;
			jMin = j;
			}
		}	
	}

param[l] = paramCourant + (double)jMin * pas;
E =  UpdateFunctionCout(topi,topj,topk,param,moyenne,nb_param,ImageTemplateOsDilate,ImageSeuilBas,seuilBas,masque_param);

return jMin;
}

void DetermineParameter(double *param	             ,
			int nb_param,
		  	grphic3d * ImageTemplateOsDilate,
			float *** ImageSeuilBas      ,
			double seuilBas               ,
		  	int ***masque_param	      )
{
double functionCoutAncien = 0,E=0;
double *moyenne;
int pari,parj,park,stop,i,j,k,topi,topj,topk,lance,iteration=0;
double Eanc = 0;
double pas = 5.;	
int resol = BASE3D.resol;
int topD  = (int)pow(2.0,1.0*resol)-1;
int tutu = 0;


if((moyenne = (double*)calloc(nb_param, sizeof(double))) == NULL) 
	{
	PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); 
	return; 
	}
	
for (int i=0;i<nb_param;i++)
	moyenne[i] = param[i] +seuilBas ;

Eanc = FunctionCoutZZ(moyenne,param,nb_param,ImageTemplateOsDilate,ImageSeuilBas,seuilBas,masque_param);
printf("iteration : %d %f \n",iteration,Eanc);


int ***masque_param2=(int ***)alloc_imatrix_3d(topD,topD,topD);

for (int i=0;i<topD;i++)
for (int j=0;j<topD;j++)
for (int k=0;k<topD;k++)
	masque_param2[i][j][k] = masque_param[i][j][k];

do
	{	
	stop=0;
	int entre = 0;
	for(pari=MINI(2,resol)-1;pari>=0;pari--) //permet de faire du coloriage
	for(parj=MINI(2,resol)-1;parj>=0;parj--)
	for(park=MINI(2,resol)-1;park>=0;park--)
		{
		for (topi=0; topi<topD; topi++)
		for (topj=0; topj<topD; topj++)
		for (topk=0; topk<topD; topk++)
	   		if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park)&&(masque_param2[topi][topj][topk]))
	  			{
				entre = 1;
				//printf("debut optimisation ligne\n");
				tutu = OptimisationLigne(topi,topj,topk,param,moyenne,nb_param,ImageTemplateOsDilate,
						  ImageSeuilBas,seuilBas,masque_param,pas);
				if (tutu==0)
					masque_param2[topi][topj][topk] = 0;
				else
					{
					for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
 					for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
 					for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
 						masque_param2[i][j][k]=masque_param[i][j][k];
					}						
				}
		}
					   
	iteration ++;	
	
	if (entre == 0)
		{
		pas = pas / 2.;
		for (int i=0;i<topD;i++)
		for (int j=0;j<topD;j++)
		for (int k=0;k<topD;k++)
			masque_param2[i][j][k] = masque_param[i][j][k];		
		}
	entre = 0;		
	if ((iteration % 10) == 0)
		{
		Eanc = FunctionCoutZZ(moyenne,param,nb_param,ImageTemplateOsDilate,ImageSeuilBas,seuilBas,masque_param);
		printf("iteration : %d %f \n",iteration,Eanc);
		}
	if  ( (iteration > 50)|| (pas <0.2) )	
		stop = 1;
	}
	while (stop==0);

//Eanc = FunctionCoutZZ(moyenne,param,nb_param,ImageTemplateOsDilate,ImageSeuilBas,seuilBas,masque_param);
free_imatrix_3d(masque_param2);
free(moyenne);
}		
