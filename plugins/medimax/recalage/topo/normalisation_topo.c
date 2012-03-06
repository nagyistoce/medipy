/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"

#ifndef COMPILE_FOR_MEDIPY
#include "noyau/gui/imx_plot.h"
#endif //fin COMPILE_FOR_MEDIPY


#include "noyau/mani_3d.h"

#ifdef __GTK_GUI
#include "noyau/gui/imx_3dmo.h"
#endif /* __GTK_GUI */

#include "math/imx_matrix.h"
#include "math/oper_3d.h"

#include "outils/imx_sort.h"

#include "traitement/norm_3d.h"

#ifdef  HAS_ITK
#include "interface_itk/itk_algorithms.hpp"
#endif

#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/normalisation_topo.h"
#include "noyau/io/imx_head.h"
#include "noyau/io/imx_file.h"
#include "noyau/zone_3d.h"
#include "outils/imx_misc.h"
#include "traitement/trai_3d.h"





int double_compare_function2(const void *a, const void *b)
{
	return (*(double*)a > *(double*)b);
}


int int_compare_function2(const void *a, const void *b)
{
	return (*(int*)a > *(int*)b);
}



/*******************************************************************************
********************************************************************************
********************************************************************************
******  Fonctions calculant l'histogramme joint par differentes methodes  ******
********************************************************************************
********************************************************************************
********************************************************************************/


/************************** histo_joint_3d() *********************************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void histo_joint_3d(void)
{
int im_1,im_2,im_plot;
grphic3d *im1,*im2;
grphic *implot;

im_1 =GET_PLACE3D(TEXT0224);
im_2 =GET_PLACE3D(TEXT0212);
im_plot=GET_PLACE(TEXT0213);

im1=ptr_img_3d(im_1);
im2=ptr_img_3d(im_2);
implot=ptr_img(im_plot);

implot->width=MAX_WIDTH;
implot->height=MAX_HEIGHT;

histo_joint_3d_p(im1,im2,implot);

imx_iniparaimg(im_plot);
show_picture(im_plot);

}



/**************************** histo_joint_3d_p() *****************************/
/*                                                                           */                         
/*****************************************************************************/

void histo_joint_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i, j, k;
int width, height, depth;
int plot_x,plot_y;
long min_x, max_x,min_y,max_y;

width=im1->width;
height=im1->height;
depth=im1->depth;

imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);

min_x=0;
max_x=im1->max_pixel;
min_y=0;
max_y=im2->max_pixel;

for(i=0;i<implot->width;i++)
    for(j=0;j<implot->height;j++)
    implot->mri[i][j]=0;
		
for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(im1->mri[i][j][k]>0 && im2->mri[i][j][k]>0)
		{
		plot_x=floor(1.0*(im1->mri[i][j][k]/**im1->rcoeff*/-min_x)/(max_x-min_x)*(implot->width-1.0));
		plot_y=floor(1.0*(im2->mri[i][j][k]/**im2->rcoeff*/-min_y)/(max_y-min_y)*(implot->height-1.0));
		implot->mri[plot_x][plot_y]++;
		}			
		
}



/************************** histo_joint_linear_3d() **************************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void histo_joint_linear_3d(void)
{
int im_1,im_2,im_plot;
grphic3d *im1,*im2;
grphic *implot;

im_1= GET_PLACE3D(TEXT0224);
im_2= GET_PLACE3D(TEXT0212);
im_plot=GET_PLACE(TEXT0213);

im1=ptr_img_3d(im_1);
im2=ptr_img_3d(im_2);
implot=ptr_img(im_plot);

implot->width=MAX_WIDTH;
implot->height=MAX_HEIGHT;

histo_joint_linear_3d_p(im1,im2,implot);

imx_iniparaimg(im_plot);
show_picture(im_plot);

}


/**************************** histo_joint_linear_3d_p() **********************/
/*                                                                           */                         
/*****************************************************************************/

void histo_joint_linear_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i, j, k;
int width, height, depth;
int plot_x,plot_y;
long min_x, max_x,min_y,max_y;
double x,y;
double **histo;
grphic3d *cer1,*cer2;


cer1=cr_grphic3d(im1);
cer2=cr_grphic3d(im2);

imx_copie_3d_p(im1,cer1);
imx_copie_3d_p(im2,cer2);


/*for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	if (im1->mask!=NULL)
	if (im1->mask->mri[i][j][k] == 0)
		cer1->mri[i][j][k]=0;
		
for (i=0;i<im2->width;i++)
for (j=0;j<im2->height;j++)
for (k=0;k<im2->depth;k++)
	if (im2->mask!=NULL)
	if (im2->mask->mri[i][j][k] == 0)
		cer2->mri[i][j][k]=0;
	*/	



width=im1->width;
height=im1->height;
depth=im1->depth;

imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);

min_x=0;
max_x=im1->max_pixel;
min_y=0;
max_y=im2->max_pixel;

histo=alloc_dmatrix(implot->width,implot->height); 

for(i=0;i<implot->width;i++)
    for(j=0;j<implot->height;j++)
    histo[i][j]=0;
		
for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(cer1->mri[i][j][k]>0 && cer2->mri[i][j][k]>0) 
		{
		x=1.0*(im1->mri[i][j][k]-min_x)/(max_x-min_x)*(implot->width-1.0);
		y=1.0*(im2->mri[i][j][k]-min_y)/(max_y-min_y)*(implot->height-1.0);
		plot_x=floor(x);
		plot_y=floor(y);
		/*histo[plot_x][plot_y]=histo[plot_x][plot_y]+1.0*(plot_x+1-x)*(plot_y+1-y);
		
		if ((plot_x+1<implot->width)&&(plot_y+1<implot->height))
		histo[plot_x+1][plot_y+1]=histo[plot_x+1][plot_y+1]+1.0*(x-plot_x)*(y-plot_y);
		
		if (plot_x+1<implot->width)
		histo[plot_x+1][plot_y]=histo[plot_x+1][plot_y]+1.0*(x-plot_x)*(plot_y+1-y);
		
		if (plot_y+1<implot->height)
		histo[plot_x][plot_y+1]=histo[plot_x][plot_y+1]+1.0*(plot_x+1-x)*(y-plot_y);*/
		
		if (((plot_x+1)<implot->width)&&((plot_y+1)<implot->height))
			{histo[plot_x][plot_y]=histo[plot_x][plot_y]+1.0*(plot_x+1-x)*(plot_y+1-y);
			 histo[plot_x+1][plot_y+1]=histo[plot_x+1][plot_y+1]+1.0*(x-plot_x)*(y-plot_y);
			 histo[plot_x+1][plot_y]=histo[plot_x+1][plot_y]+1.0*(x-plot_x)*(plot_y+1-y);
			 histo[plot_x][plot_y+1]=histo[plot_x][plot_y+1]+1.0*(plot_x+1-x)*(y-plot_y);
			}
		else
			{histo[plot_x][plot_y]=histo[plot_x][plot_y]+1.0;}
		}			

for(i=0;i<implot->width;i++)
    for(j=0;j<implot->height;j++)
    if (histo[i][j]>0)
		implot->mri[i][j]=(int)floor(histo[i][j]+1);
		else
		implot->mri[i][j]=0;

		
free_dmatrix(histo,implot->width,implot->height);		
free_grphic3d(cer1);
free_grphic3d(cer2);

}

/*void histo_joint_linear_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i, j, k,Ntot,l=0;;
int width, height, depth;
double *data1,*data2;
FILE *fic;

width=im1->width;
height=im1->height;
depth=im1->depth;

Ntot=0;
for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(im1->mri[i][j][k]>0 && im2->mri[i][j][k]>0) 
			Ntot++;
			
data1=(double *)malloc(Ntot*sizeof(double));			
data2=(double *)malloc(Ntot*sizeof(double));			

for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(im1->mri[i][j][k]>0 && im2->mri[i][j][k]>0) 
		{
		data1[l]=im1->mri[i][j][k];
		data2[l]=im2->mri[i][j][k];
		l++;
		}
for (i=0;i<Ntot;i++)
	{
	data1[i]=data1[i]+0.000000001*data2[i];
	}
qsort(data1,Ntot,sizeof(double),double_compare_function2);

for (i=0;i<Ntot;i++)
	{
	l = floor(data1[i]);
	data2[i]=(data1[i]-l)*10000000000;
	data1[i]=l;
	}


// Ouverture du fichier 
fic = fopen("/home/noblet/testhistojoint.txt", "w+");

// Ecriture des donn�s dans le fichier 
for (i=0;i<Ntot;i++)
	fprintf(fic, "%f \t %f \n", data1[i], data2[i]);

// Fermeture du fichier 
fclose(fic);		

free(data1);free(data2);
}*/


/************************** histo_joint_linear_norm_3d() *********************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void histo_joint_linear_norm_3d(void)
{
int im_1,im_2,im_plot,e;
grphic3d *im1,*im2;
grphic *implot;

im_1= GET_PLACE3D(TEXT0224);
im_2= GET_PLACE3D(TEXT0212);
im_plot=GET_PLACE(TEXT0213);


im1=ptr_img_3d(im_1);
im2=ptr_img_3d(im_2);
implot=ptr_img(im_plot);

TOPO_SEUIL_SIGNIFICATIVITE=GET_DOUBLE("seuil de significativite ?",0.05,&e);
histo_joint_linear_norm_3d_p(im1,im2, implot,NULL);
//trimmed_robust_histo (implot,implot,1.96);

imx_iniparaimg(im_plot);
show_picture(im_plot);

}


/**************************** factoriel	**********************/
/*                                                           */                         
/*************************************************************/

int factoriel (int a)
	{
	int i,res;
	res=1;
	for (i=1;i<=a;i++)
		res=res*i;
	return (res);	
	}
	
	
	
/**************************** histo_joint_linear_norm_3d_p() **********************/
/*                                                                           */                         
/*****************************************************************************/

void histo_joint_linear_norm_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot, double **confiance)
{
int i,j,wdth,hght,l;
int *nbx,*nby,Ntot;
double tmp,max,min/*,sigma=1.65*/,p,seuil;


seuil=1-TOPO_SEUIL_SIGNIFICATIVITE;

//implot->width=MAX_WIDTH;
//implot->height=MAX_HEIGHT;

histo_joint_linear_3d_p(im1,im2,implot);
//histo_joint_3d_p(im1,im2,implot);

wdth=implot->width;
hght=implot->height;

// allocation 
nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

// initialisation 
for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}


for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+implot->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+implot->mri[i][j];	
		}

Ntot=0;
for (l=0; l<wdth; l++)
Ntot=Ntot+nbx[l];

/* On garde pour l'histogramme de chaque image les points entre 5% �95% */
/*tmp=0;
i=0;
while (tmp<Ntot*0.05)
{
tmp+=nbx[i];
for (j=0; j<hght; j++)
	implot->mri[i][j]=0;

i++;
}

tmp=0;
i=0;
while (tmp<Ntot*0.05)
{
tmp+=nby[i];
for (j=0; j<wdth; j++)
	implot->mri[j][i]=0;

i++;
}

tmp=0;
i=wdth-1;
while (tmp<Ntot*0.05)
{
tmp+=nbx[i];
for (j=0; j<hght; j++)
	implot->mri[i][j]=0;

i--;
}

tmp=0;
i=hght-1;
while (tmp<Ntot*0.05)
{
tmp+=nby[i];
for (j=0; j<wdth; j++)
	implot->mri[j][i]=0;

i--;
}

*/



for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if ((nbx[i]>0)&&(nby[j]>0)&&(implot->mri[i][j]>0))
		{
		
			
		//if (1.0*nbx[i]*nby[j]/Ntot>10)
		//	{
		//	max=1.0*nbx[i]*nby[j]/Ntot+sigma*sqrt(nbx[i]*nby[j]/Ntot*(1-nbx[i]*nby[j]/(Ntot*Ntot)));
		//	min=1.0*nbx[i]*nby[j]/Ntot-sigma*sqrt(nbx[i]*nby[j]/Ntot*(1-nbx[i]*nby[j]/(Ntot*Ntot)));
		//	}
		//else
			{	

				
				p=1.0*(double)nbx[i]*(double)nby[j]/((double)Ntot*(double)Ntot);
			tmp=eval_binomiale(Ntot,p,implot->mri[i][j]);	
			
		if ((finite(tmp)==0)||(tmp>1.0001))
			printf("y a un bug dans l'evaluation binomial : k : %d   p : %f   proba_tmp %f \n",(int)implot->mri[i][j],p,tmp);
		
		if (tmp<seuil)	
			implot->mri[i][j]=0;
		
		max=min=0;
			
			if (confiance != NULL)
				confiance[i][j]=tmp;
			
			}
			
		min=0;
		if ((implot->mri[i][j]<max)&&(implot->mri[i][j]>min))
			implot->mri[i][j]=0;
	
			
		}


show_picture(implot->pos);
free(nbx);free(nby);		
}


/*********************** eval_binomiale	**********************/
/*                                                           */                         
/*************************************************************/

double eval_binomiale(int N,double p,int k)
{
int i;
double q,res,somme;

q=1.0-p;

/* initialisation de la r�urrence pour i=1 */
res=N*pow(q,N-1)*p;
somme=res;

/*recurrence*/
for (i=2; i<= k; i++)
	{
	res=res*1.0*((double)N-(double)i+1)/(double)i*p/q;
	somme+=res;
	}

return(somme);
}
/*
void histo_joint_linear_norm_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i,j,wdth,hght,nb;
double tmp,IM;
grphic *histo;

implot->width=MAX_WIDTH;
implot->height=MAX_HEIGHT;

histo_joint_linear_3d_p(im1,im2,implot);

wdth=implot->width;
hght=implot->height;
histo=cr_grphic_modif(wdth,hght,0,1,0);

nb=1;
while(nb>0)
{
nb=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	histo->mri[i][j]=implot->mri[i][j];

IM=calcul_IM_histo(histo);

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	if (implot->mri[i][j]>0)
	{
	histo->mri[i][j]=0;
	tmp=calcul_IM_histo(histo);
	
	histo->mri[i][j]=implot->mri[i][j];
	
	printf("%d %d \r",i,j);
	
	if (tmp>IM)
		{implot->mri[i][j]=0;nb++;}
	
	}
	printf("nb : %d\n",nb);

show_picture(implot->pos);
}
free_grphic(histo);

}
*/

double calcul_IM_histo(grphic *histo)
{
int i,j,wdth,hght,l;
int *nbx,*nby,Ntot;
double IM;


wdth=histo->width;
hght=histo->height;

// allocation 
nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

// initialisation 
for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}


for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}

Ntot=0;
for (l=0; l<wdth; l++)
Ntot=Ntot+nbx[l];


IM=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		if (histo->mri[i][j]>0)
		{
		IM+=1.0*histo->mri[i][j]/Ntot*log(1.0*histo->mri[i][j]*Ntot/(1.0*nbx[i]*nby[j]));
		}


free(nbx);free(nby);		
return(IM);
}



/********************* histo_joint_maxentropie_3d_p() ************************/
/*                                                                           */                         
/*****************************************************************************/

void histo_joint_maxentropie_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i, j, k,Ntot,l;
int width, height, depth;
int plot_x,plot_y;
long min_x, max_x,min_y,max_y;
double x,y;
double *nbx,*nby,tmp;

implot->width=MAX_WIDTH;
implot->height=MAX_HEIGHT;

width=im1->width;
height=im1->height;
depth=im1->depth;

/* allocation */
nbx=(double *)malloc(implot->width*sizeof(double));
nby=(double *)malloc(implot->height*sizeof(double));

imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);

min_x=0;
max_x=im1->max_pixel;
min_y=0;
max_y=im2->max_pixel;

for (l=0; l<implot->width; l++)
	nbx[l]=0;

for (l=0; l<implot->height; l++)
	nby[l]=0;
		
for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(im1->mri[i][j][k]>0) 
		{
		x=1.0*(im1->mri[i][j][k]/**im1->rcoeff*/-min_x)/(max_x-min_x)*(implot->width-1.0);
		plot_x=floor(x);
		
		nbx[plot_x]=nbx[plot_x]+1.0*(plot_x+1-x);
		
		if (plot_x+1<implot->width)
		nbx[plot_x]=nbx[plot_x]+1.0*(x-plot_x);	
		}			

for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(im2->mri[i][j][k]>0) 
		{
		y=1.0*(im2->mri[i][j][k]/**im2->rcoeff*/-min_y)/(max_y-min_y)*(implot->height-1.0);
		plot_y=floor(y);
		nby[plot_y]=nby[plot_y]+1.0*(plot_y+1-y);
		
		if (plot_y+1<implot->height)
		nby[plot_y+1]=nby[plot_y+1]+1.0*(y-plot_y);
		}			




Ntot=0;
for (l=0; l<implot->width; l++)
Ntot=Ntot+nbx[l];


for (i=0; i<implot->width; i++)
for (j=0; j<implot->height; j++)
if ((nbx[i]>0)&&(nby[j]>0))
		{
		tmp=(1.0*nbx[i]*nby[j]/Ntot);
		implot->mri[i][j]=(int)(tmp);
		}
else
	implot->mri[i][j]=0;

free(nbx);free(nby);		
}







/*******************************************************************************
********************************************************************************
********************************************************************************
*  Normalisation par m�ange de gaussiennes sans calcul de l'histogramme joint * 
********************************************************************************
********************************************************************************
********************************************************************************/

/*******************************************************************************
**        normalisation_gaussienne2d
**
*******************************************************************************/

void  normalisation_gaussienne2d(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;

	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	nb_classe= GET_INT("nb de classes", 4, &e);

	im1=ptr_img_3d(im_1);
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	imx_inimaxminpixel_3d_p(im1);
	imx_inimaxminpixel_3d_p(im2);
	normalisation_gaussienne2d_p(im1,im2,imres,nb_classe);

	show_picture_3d(im_res);

}

/*******************************************************************************
**         normalisation_gaussienne2d_p
**
*******************************************************************************/

void  normalisation_gaussienne2d_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
Param_Gaussienne_2d *param;

param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));

fit_EM_gaussienne2d(im1,im2,param, nb_classe);

apply_normalisation_gaussienne2d(im1,im2,imres,nb_classe,param);

free(param);
}

/************************** ini_kmeans_histo_joint_3d_p() ********************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void ini_kmeans_histo_joint_3d_p(grphic3d *im1,grphic3d *im2, Param_Gaussienne_2d *param, int nb_classe)
{
int wdth,hght,dpth,i,j,k,l,cl=0,condition,*nb,nb_iter,mini_x,mini_y,maxi_x,maxi_y;
double *mx,*my,*mxold,*myold,tmp,min;
grphic3d *imclasse;
int *pclasse,N;
imclasse=cr_grphic3d(im1);
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

/*allocation memoire*/
mx=(double *)malloc(nb_classe*sizeof(double));
my=(double *)malloc(nb_classe*sizeof(double));
mxold=(double *)malloc(nb_classe*sizeof(double));
myold=(double *)malloc(nb_classe*sizeof(double));
nb=(int *)malloc(nb_classe*sizeof(int));


mini_x=10000;
mini_y=10000;
maxi_x=0;
maxi_y=0;

/* recherche du min et du max*/
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
		{
		mini_x=MINI(mini_x,im1->mri[i][j][k]);
		mini_y=MINI(mini_y,im2->mri[i][j][k]);
		maxi_x=MAXI(maxi_x,im1->mri[i][j][k]);
		maxi_y=MAXI(maxi_y,im2->mri[i][j][k]);
		}


/*initialisation des moyennes*/
for (i=0; i<nb_classe; i++)
	{
	mx[i]=mini_x+1.0*(i+0.5)*(maxi_x-mini_x)/nb_classe;
	my[i]=mini_y+1.0*(i+0.5)*(maxi_y-mini_y)/nb_classe;;
	}

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
	{min=HUGE_VAL;
	for (l=0;l<nb_classe;l++)
		{tmp=(im1->mri[i][j][k]-mx[l])*(im1->mri[i][j][k]-mx[l])+(im2->mri[i][j][k]-my[l])*(im2->mri[i][j][k]-my[l]);
		if (tmp<min) {min=tmp; cl=l;}
		}
		imclasse->mri[i][j][k]=cl;
	}

condition=1;
nb_iter=0;
while (condition>0)
{
/* recopie dans mxold et myold*/
for (l=0;l<nb_classe;l++)
	{
	mxold[l]=mx[l];
	myold[l]=my[l];
	}

/* mise a jour des moyennes */
for (l=0;l<nb_classe;l++)
	{mx[l]=my[l]=nb[l]=0;}

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
			{
			mx[imclasse->mri[i][j][k]]=mx[imclasse->mri[i][j][k]]+im1->mri[i][j][k];
			my[imclasse->mri[i][j][k]]=my[imclasse->mri[i][j][k]]+im2->mri[i][j][k];
			nb[imclasse->mri[i][j][k]]=nb[imclasse->mri[i][j][k]]+1;
			}

for (l=0;l<nb_classe;l++)
	{
	mx[l]=1.0*mx[l]/nb[l];
	my[l]=1.0*my[l]/nb[l];
	}

/* mise a jour de imclasse */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
	{min=HUGE_VAL;
	for (l=0;l<nb_classe;l++)
		{tmp=(im1->mri[i][j][k]-mx[l])*(im1->mri[i][j][k]-mx[l])+(im2->mri[i][j][k]-my[l])*(im2->mri[i][j][k]-my[l]);
		if (tmp<min) {min=tmp; cl=l;}
		}
		imclasse->mri[i][j][k]=cl;
	}


nb_iter++;

condition=0;
for (l=0;l<nb_classe;l++)
	condition=condition+(mx[l]-mxold[l])*(mx[l]-mxold[l])+(my[l]-myold[l])*(my[l]-myold[l]);

if (nb_iter>300)
	condition=0;
}
printf("nombre d'iteration jusqu'a convergence : %d \n",nb_iter);

/* mise a jour de param */
/* allocation */
pclasse=(int *) malloc (nb_classe*sizeof(int));

/* initialisation */
for (l=0;l<nb_classe;l++)
		pclasse[l]=0;

for (i=0;i<nb_classe;i++)
	{
	param[i].mx=0.0;
	param[i].my=0.0;
	param[i].sx=0.0;
	param[i].sy=0.0;
	param[i].sxy=0.0;
	param[i].l=0.0;
	}


/* calcul de la moyenne */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
		{
	param[imclasse->mri[i][j][k]].mx=param[imclasse->mri[i][j][k]].mx+im1->mri[i][j][k];
	param[imclasse->mri[i][j][k]].my=param[imclasse->mri[i][j][k]].my+im2->mri[i][j][k];
	pclasse[imclasse->mri[i][j][k]]=pclasse[imclasse->mri[i][j][k]]+1;
	}

for (l=0;l<nb_classe;l++)
	{
	param[l].mx=1.0*param[l].mx/pclasse[l];
	param[l].my=1.0*param[l].my/pclasse[l];
	}

/* calcul de l'ecart type */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
		{
	param[imclasse->mri[i][j][k]].sx=param[imclasse->mri[i][j][k]].sx+im1->mri[i][j][k]*im1->mri[i][j][k];
	param[imclasse->mri[i][j][k]].sy=param[imclasse->mri[i][j][k]].sy+im2->mri[i][j][k]*im2->mri[i][j][k];
	param[imclasse->mri[i][j][k]].sxy=param[imclasse->mri[i][j][k]].sxy+im1->mri[i][j][k]*im2->mri[i][j][k];
	}

for (l=0;l<nb_classe;l++)
	{
	param[l].sx=1.0*(param[l].sx-param[l].mx*param[l].mx)/pclasse[l];
	param[l].sy=1.0*(param[l].sy-param[l].my*param[l].my)/pclasse[l];
	param[l].sxy=1.0*(param[l].sxy-param[l].mx*param[l].my)/pclasse[l];
	}

N=0;
for (l=0;l<nb_classe;l++)
	N=N+pclasse[l];

for (l=0;l<nb_classe;l++)
	param[l].l=1.0*pclasse[l]/N;

free(pclasse);


free(mx);free(my);free(mxold);free(myold);free(nb);free_grphic3d(imclasse);
}

/*********************************************************************/
/*              Moyenne Ecart-type classe	                           */
/*********************************************************************/
void mean_std_classe(grphic *im,grphic *imcl,int nb_classe, Param_Gaussienne_2d *param)
{
int i,j,l;
int wdth,hght,N;
int *pclasse;

wdth=im->width;
hght=im->height;


/* allocation */
pclasse=(int *) malloc (nb_classe*sizeof(int));

/* initialisation */
for (l=0;l<nb_classe;l++)
		pclasse[l]=0;

for (i=0;i<nb_classe;i++)
	{
	param[i].mx=0.0;
	param[i].my=0.0;
	param[i].sx=0.0;
	param[i].sy=0.0;
	param[i].sxy=0.0;
	param[i].l=0.0;
	}

	
/* calcul de la moyenne */	
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im->mri[i][j]>0)
	{
	param[imcl->mri[i][j]].mx=param[imcl->mri[i][j]].mx+im->mri[i][j]*i;
	param[imcl->mri[i][j]].my=param[imcl->mri[i][j]].my+im->mri[i][j]*j;
	pclasse[imcl->mri[i][j]]=pclasse[imcl->mri[i][j]]+im->mri[i][j];
	}

for (l=0;l<nb_classe;l++)
	{
	param[l].mx=1.0*param[l].mx/pclasse[l];
	param[l].my=1.0*param[l].my/pclasse[l];
	}
	
/* calcul de l'ecart type */
/*for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im->mri[i][j]>0)
	{
	param[imcl->mri[i][j]].sx=param[imcl->mri[i][j]].sx+im->mri[i][j]*i*i;
	param[imcl->mri[i][j]].sy=param[imcl->mri[i][j]].sy+im->mri[i][j]*j*j;
	param[imcl->mri[i][j]].sxy=param[imcl->mri[i][j]].sxy+im->mri[i][j]*i*j;
	}

for (l=0;l<nb_classe;l++)
	{
	param[l].sx=1.0*(param[l].sx-param[l].mx*param[l].mx)/pclasse[l];
	param[l].sy=1.0*(param[l].sy-param[l].my*param[l].my)/pclasse[l];
	param[l].sxy=1.0*(param[l].sxy-param[l].mx*param[l].my)/pclasse[l];
	}	
*/

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im->mri[i][j]>0)
	{
	param[imcl->mri[i][j]].sx=param[imcl->mri[i][j]].sx+im->mri[i][j]*pow((i-param[imcl->mri[i][j]].mx),2.0);
	param[imcl->mri[i][j]].sy=param[imcl->mri[i][j]].sy+im->mri[i][j]*pow((j-param[imcl->mri[i][j]].my),2.0);
	param[imcl->mri[i][j]].sxy=param[imcl->mri[i][j]].sxy+im->mri[i][j]*(i-param[imcl->mri[i][j]].mx)*(j-param[imcl->mri[i][j]].my);
	}

for (l=0;l<nb_classe;l++)
	{
	param[l].sx=1.0*param[l].sx/pclasse[l];
	param[l].sy=1.0*param[l].sy/pclasse[l];
	param[l].sxy=1.0*param[l].sxy/pclasse[l];
	}
		
N=0;
for (l=0;l<nb_classe;l++)
	N=N+pclasse[l];

for (l=0;l<nb_classe;l++)
	param[l].l=1.0*pclasse[l]/N;
		
free(pclasse);	
}

/*******************************************************************************
**      					  fit_EM_gaussienne2d
**
*******************************************************************************/
void fit_EM_gaussienne2d(grphic3d *im1,grphic3d *im2, Param_Gaussienne_2d *param, int nb_classe)
{
int stop,l,nb;
Param_Gaussienne_2d *new_param;

/* allocation memoire */
new_param= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));

/* initialisation */
ini_kmeans_histo_joint_3d_p(im1,im2,param,nb_classe);

nb=0;
stop = 0 ;
printf("%f  \n",vraisemblance_gaussienne_2d(im1,im2, param, nb_classe));
do
{
estime_EM_gaussienne2d(im1,im2,param,new_param,nb_classe);
stop = 1;
for (l=0;l<nb_classe;l++)
	{
	if ((1-1.0*new_param[l].mx/param[l].mx)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].my/param[l].my)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].sx/param[l].sx)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].sy/param[l].sy)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].sxy/param[l].sxy)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].l/param[l].l)>0.05)
		stop=0;
	}

/*for (l=0;l<nb_classe;l++)
	{
	printf("mx 	: %f  %f \n",new_param[l].mx,param[l].mx);
	printf("my 	: %f  %f \n",new_param[l].my,param[l].my);
	printf("sx 	: %f  %f \n",new_param[l].sx,param[l].sx);
	printf("sy 	: %f  %f \n",new_param[l].sy,param[l].sy);
	printf("sxy : %f  %f \n",new_param[l].sxy,param[l].sxy);
	printf("l 	: %f  %f \n",new_param[l].l,param[l].l);
	printf("\n");
	}*/

for (l=0;l<nb_classe;l++)
	{
	param[l].mx=new_param[l].mx;
	param[l].my=new_param[l].my;
	param[l].sx=new_param[l].sx;
	param[l].sy=new_param[l].sy;
	param[l].sxy=new_param[l].sxy;
	param[l].l=new_param[l].l;
	}


printf("%f  \n",vraisemblance_gaussienne_2d(im1,im2, param, nb_classe));
nb++;
if (nb>50)
stop=0;

}
while (stop==0);

free(new_param);
}


/*******************************************************************************
**      					  estime_gaussienne2d_EM
**
*******************************************************************************/
void estime_EM_gaussienne2d(grphic3d *im1,grphic3d *im2,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe)
{
int wdth,hght,dpth,i,j,k,l,Ntot;
double tmp;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;



Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
		Ntot=Ntot+1;



/* mise a jour des moyennes et du facteur de normalisation  */
for (l=0;l<nb_classe;l++)
		{
		new_param[l].mx=0.0;
		new_param[l].my=0.0;
		new_param[l].l=0.0;

		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
			{
			tmp=eval_proba_apriori(im1->mri[i][j][k], im2->mri[i][j][k], param, l,nb_classe);
			new_param[l].l=new_param[l].l+tmp;
			new_param[l].mx=new_param[l].mx+1.0*im1->mri[i][j][k]*tmp;
			new_param[l].my=new_param[l].my+1.0*im2->mri[i][j][k]*tmp;
			}
		}

for (l=0;l<nb_classe;l++)
		{
		new_param[l].l=1.0*new_param[l].l/Ntot;
		}


for (l=0;l<nb_classe;l++)
		{
		new_param[l].mx=1.0*new_param[l].mx/(Ntot*new_param[l].l);
		new_param[l].my=1.0*new_param[l].my/(Ntot*new_param[l].l);
		}


/* mise a jour des variances */
for (l=0;l<nb_classe;l++)
		{
		new_param[l].sx=0.0;
		new_param[l].sy=0.0;
		new_param[l].sxy=0.0;

		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
			{
			tmp=eval_proba_apriori(im1->mri[i][j][k], im2->mri[i][j][k], param, l,nb_classe);

			new_param[l].sx=new_param[l].sx+1.0*(im1->mri[i][j][k]-new_param[l].mx)*(im1->mri[i][j][k]-new_param[l].mx)*tmp;
			new_param[l].sy=new_param[l].sy+1.0*(im2->mri[i][j][k]-new_param[l].my)*(im2->mri[i][j][k]-new_param[l].my)*tmp;
			new_param[l].sxy=new_param[l].sxy+1.0*(im1->mri[i][j][k]-new_param[l].mx)*(im2->mri[i][j][k]-new_param[l].my)*tmp;
			}
		new_param[l].sx=1.0*new_param[l].sx/(Ntot*new_param[l].l);
		new_param[l].sy=1.0*new_param[l].sy/(Ntot*new_param[l].l);
		new_param[l].sxy=1.0*new_param[l].sxy/(Ntot*new_param[l].l);
		}


}


/*******************************************************************************
**      					  vraisemblance_gaussienne2d
**
*******************************************************************************/
double vraisemblance_gaussienne_2d(grphic3d *im1,grphic3d *im2, Param_Gaussienne_2d *param, int nb_classe)
{
double res,tmp;
int i,j,k,l,wdth,hght,dpth;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

res=0;
	for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
			{tmp=0;
			for (l=0;l<nb_classe;l++)
				tmp=tmp+eval_gaussienne_2d(im1->mri[i][j][k],im2->mri[i][j][k],&param[l]);
				res=res+1.0*log(tmp);
			}
return(res);
}

/*******************************************************************************
**      					  eval_proba_apriori
**
*******************************************************************************/
double eval_proba_apriori(double x, double y, Param_Gaussienne_2d *param,int classe, int nb_classe)
{
int l;
double res=0.,tot;
tot=0;
for (l=0;l<nb_classe;l++)
	if (l==classe)
		{
		res = eval_gaussienne_2d(x, y,&param[l]);
		tot=tot + res;
		}
	else
		tot=tot+eval_gaussienne_2d(x, y,&param[l]);

res=1.0*res/tot;

return(res);
}


/*********************************************************************/
/*         apply_normalisation_gaussienne2d						               */
/*********************************************************************/
void apply_normalisation_gaussienne2d(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe, Param_Gaussienne_2d *param)
{
int i,j,k,l,xi;
double max_reca,max_ref,rmax_reca,rmax_ref;
int wdth,hght,dpth;
double x,y,proba_tot,proba_tmp,tmp;
double *var_gaussienne;
double *ft;
int size_ft=TOPO_SIZE_NORM_HISTO;
wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

var_gaussienne = (double *)malloc(nb_classe*sizeof(double));
ft = (double *)malloc(size_ft*sizeof(double));

for (i=0;i<nb_classe;i++)
	var_gaussienne[i]=param[i].sx*param[i].sy-param[i].sxy*param[i].sxy;

/* recherche des max des 2 images */
max_reca=imreca->max_pixel;
max_ref=imref->max_pixel;

rmax_reca=max_reca*imreca->rcoeff;
rmax_ref=max_ref*imref->rcoeff;

if (rmax_reca>rmax_ref)
	imx_copie_param_3d_p(imreca,imres);
else
	imx_copie_param_3d_p(imref,imres);



for (i=0;i<size_ft; i++)
	{	proba_tot=0;
		tmp=0;
		x=1.0*(i+0.5)*max_reca/size_ft;
		for (j=0;j<size_ft;j++)
			{
				y=1.0*(j+0.5)*max_ref/size_ft;

			for (l=0;l<nb_classe;l++)
				{	proba_tmp=eval_gaussienne_2d(x, y, &param[l]);
					if (proba_tmp)
					{proba_tot=proba_tot+proba_tmp;
					tmp=tmp+proba_tmp*y;
					}
				}
			}
	ft[i]=1.0*((tmp*imref->rcoeff/imres->rcoeff)/proba_tot);
	}



for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if (imreca->mri[i][j][k]>0)
	{
	x=imreca->mri[i][j][k]/max_reca*size_ft;
	xi=floor(x);
	imres->mri[i][j][k]=(int)((1.0*(x-xi)*ft[xi]+1.0*(xi+1-x)*ft[xi+1]));
	}
else imres->mri[i][j][k]=0;


free(var_gaussienne);free(ft);

}






/*******************************************************************************
********************************************************************************
********************************************************************************
*  Normalisation par m�ange de gaussiennes avec calcul de l'histogramme joint * 
********************************************************************************
********************************************************************************
********************************************************************************/


/*******************************************************************************
**        normalisation_gaussienne2d_histo_joint                                                
**                                                                                                   
*******************************************************************************/

void  normalisation_gaussienne2d_histo_joint(void)
{
	int im_1,im_2,im_res,im_histo,e;
	int nb_classe,choix;
	grphic3d *im1,*im2,*imres;
	grphic *histo;
/*	int i,j,k,n,nb_gris;*/
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	im_histo=GET_PLACE("histo joint");

	nb_classe= GET_INT("nb de classes", 10, &e);
	
	choix = GET_INT("histo standart(0) histo norme(1) ou histo reclasse(2)", 1, &e);
	
	TOPO_SEUIL_SIGNIFICATIVITE=GET_DOUBLE("seuil de significativite ?",0.05,&e);

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	histo=ptr_img(im_histo);

/*n=0;
for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	if ((im1->mri[i][j][k]!=0)&&(im2->mri[i][j][k]!=0))
		n++;

nb_gris=sqrt(0.01*n);

printf ("Nombre de bins de l'histo : %d\n",nb_gris);

histo->width=histo->height=nb_gris;		
*/
	switch (choix)
{
case 0 : 	normalisation_gaussienne2d_histo_joint_standart_p(im1,im2,imres,histo,nb_classe);break;
case 1 : 	normalisation_gaussienne2d_histo_joint_norm_p(im1,im2,imres,histo,nb_classe); break;
case 2 :	normalisation_gaussienne2d_histo_joint_reclasse_p(im1,im2,imres,histo,nb_classe);;break;
default : normalisation_gaussienne2d_histo_joint_norm_p(im1,im2,imres,histo,nb_classe);;
}
	
show_picture_3d(im_res);

}


void  normalisation_gaussienne2d_histo_joint_standart(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	nb_classe= GET_INT("nb de classes", 10, &e);
	

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

	imx_normalisation_gaussienne2d_histo_joint_standart_p(im1,im2,imres,nb_classe);
	
	show_picture_3d(im_res);

}


void  normalisation_gaussienne_histo(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	nb_classe= GET_INT("nb de classes", 5, &e);
	

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

	#ifdef HAS_ATLAS_PERFUSION
	CPP_normalisation_GM(im1,im2,imres, nb_classe, 0.05);
	#else
	printf("Attention Vous n'avez pas compile avec HAS_ATLAS_PERFUSION ! La 	normalisation Moyenne Ecart Type est utilisee par defaut ...\n");
	imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
	#endif

	
	show_picture_3d(im_res);

}


void  normalisation_quantile_histo(void)
{
#ifdef HAS_ITK

	int im_1,im_2,im_res;
	int nb_classe,e;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	nb_classe= GET_INT("nb de quantile", 7, &e);
	

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

	itk_HistogramMatchingImageFilter_3d_p(im1,im2,imres, 1024, nb_classe);
	
	show_picture_3d(im_res);
#endif
}


void  normalisation_gaussienne2D_histojoint_norm(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	nb_classe= GET_INT("nb de gaussiennes", 10, &e);
	TOPO_SEUIL_SIGNIFICATIVITE=GET_DOUBLE("seuil de significativite ?",0.05,&e);


	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

	imx_normalisation_gaussienne2d_histo_joint_norm_p(im1,im2,imres, nb_classe);
	
	show_picture_3d(im_res);

} 

void  normalisation_moyenne_ecart_type_robust(void)
{
	int im_1,im_2,im_res;
		grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);

	imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
	
	show_picture_3d(im_res);

}


/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_standart_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_standart_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
grphic *histo;
histo=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
histo->pos=0;
normalisation_gaussienne2d_histo_joint_standart_p(im1,im2,imres,histo, nb_classe);
free_grphic(histo);
}




/*******************************************************************************
**         normalisation_gaussienne2d_histo_joint_standart_p                                                
**                                                                                                   
*******************************************************************************/

void  normalisation_gaussienne2d_histo_joint_standart_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe)
{

/*int i,j,k;
double tmp;*/
Param_Gaussienne_2d *param;

param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));

histo_joint_linear_3d_p(im1,im2,histo);
/*tmp = GET_FLOAT("sigma", 1, &j);*/
/*tmp=1;
trimmed_robust_histo (histo,histo,tmp);
*/
imx_inimaxminpixel_p(histo);

/*//fit_gaussienne2d_EM(histo, param, nb_classe);


param[0].mx=100;
param[0].my=150;
param[0].sx=100;
param[0].sy=200;
param[0].sxy=100;
param[0].l=1;

param[1].mx=200;
param[1].my=50;
param[1].sx=400;
param[1].sy=400;
param[1].sxy=300;
param[1].l=0.4;


for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	{tmp=0;
	for (k=0;k<nb_classe;k++)
			{tmp+=	1000000.0*eval_gaussienne_2d(i, j,&param[k]);
			}
	histo->mri[i][j]=(int)tmp;		
	}

*/

fit_gaussienne2d_EM(histo, param, nb_classe);
	
apply_normalisation_gaussienne_histo_joint(im1,im2,imres,histo,nb_classe,param,-1.0,-1.0,1.0);

imx_inimaxminpixel_3d_p(imres);

free(param);
}


/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_norm_symetrique                                                
**                                                                                                   
*******************************************************************************/
void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique(void)
{
	int im_1,im_2,im_res1,im_res2,nb_classe,e;
	grphic3d *im1,*im2,*imres1,*imres2;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res1=GET_PLACE3D(TEXT0006);
	im_res2=GET_PLACE3D(TEXT0006);

	

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres1=ptr_img_3d(im_res1);
	imres2=ptr_img_3d(im_res2);

	nb_classe= GET_INT("nb de classes", 10, &e);
	
	 imx_copie_param_3d_p(im1,imres1);
	 imx_copie_param_3d_p(im2,imres2);

	//imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_p(im1,im2, im1, im2, imres1, imres2, nb_classe);
	imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_iteratif_p(im1,im2, im1, im2, imres1, imres2, nb_classe);

	show_picture_3d(im_res1);
	show_picture_3d(im_res2);
	

} 
/*******************************************************************************
**         imx_norm_seuil_meanecty_symetrique_3d_p                                                
**                                                                                                   
*******************************************************************************/

void imx_norm_seuil_meanecty_symetrique_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *im1d, grphic3d *im2d, grphic3d *imres1, grphic3d *imres2)
{
int i,j,k,l,m;
int height,width,depth;
double moy1,moy2,var1,var2,sygma1,sygma2,moy,sygma;
double rcoeff1,rcoeff2, rcoeff;
double alpha1,alpha2, beta1, beta2;
float max1,min1,max2,min2, max, min;
int err;



m=0;
l=0;
moy1=0;
moy2=0;
var1=0;
var2=0;
height=im1d->height;
width=im1d->width;
depth=im1d->depth;
rcoeff1=im1d->rcoeff;
rcoeff2=im2d->rcoeff;

/*cacul des moyennes des images*/	
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im1d->mri[i][j][k])
      {
        l++;
        moy1 += im1d->mri[i][j][k];
      }
     }		
height=im2d->height;
width=im2d->width;
depth=im2d->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im2d->mri[i][j][k])
      {
        m++;
        moy2 += im2d->mri[i][j][k];
      }
     }

moy1= rcoeff1*moy1/l;
moy2= rcoeff2*moy2/m;

moy=0.5*(moy1+moy2);

/*calcul des variances des images*/
height=im1d->height;
width=im1d->width;
depth=im1d->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im1d->mri[i][j][k])
          var1 +=(rcoeff1*im1d->mri[i][j][k]-moy1)*(rcoeff1*im1d->mri[i][j][k]-moy1);
      }

height=im2d->height;
width=im2d->width;
depth=im2d->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im2d->mri[i][j][k])
          var2 +=(rcoeff2*im2d->mri[i][j][k]-moy2)*(rcoeff2*im2d->mri[i][j][k]-moy2);
      }

var1= var1/l;
var2=var2/m;
sygma1=sqrt(var1);
sygma2=sqrt(var2);

sygma=sqrt(0.5*(var1+var2));

alpha1=rcoeff1*sygma/sygma1;
alpha2=rcoeff2*sygma/sygma2;

beta1=moy-moy1*sygma/sygma1;
beta2=moy-moy2*sygma/sygma2;


imx_copie_param_3d_p(im1,imres1);
imx_copie_param_3d_p(im2,imres2);


max1=(float)(alpha1*(im1->max_pixel)+beta1);
min1=(float)(alpha1*(im1->min_pixel)+beta1);

max2=(float)(alpha2*(im2->max_pixel)+beta2);
min2=(float)(alpha2*(im2->min_pixel)+beta2);

max = MAXI(max1,max2);
min = MINI (min1, min2);

/*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres1);
err=imx_brukermax_3d(max,min,imres2);

rcoeff=imres1->rcoeff;

  
height=im1->height;
width=im1->width;
depth=im1->depth;
for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for(k=0;k<depth;k++)
       if (im1->mri[i][j][k]!=0)
        {imres1->mri[i][j][k]=(TYPEMRI3D)floor((alpha1*im1->mri[i][j][k]+beta1)/rcoeff);}
       else
        {imres1->mri[i][j][k]=0;}
	
/* Modif JP 2/2001 imres->min_pixel=0;*/
imres1->width=width;
imres1->height=height;
imres1->depth=depth;

height=im2->height;
width=im2->width;
depth=im2->depth;
for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for(k=0;k<depth;k++)
       if (im2->mri[i][j][k]!=0)
        {imres2->mri[i][j][k]=(TYPEMRI3D)floor((alpha2*im2->mri[i][j][k]+beta2)/rcoeff);}
       else
        {imres2->mri[i][j][k]=0;}
	
/* Modif JP 2/2001 imres->min_pixel=0;*/
imres2->width=width;
imres2->height=height;
imres2->depth=depth;

/* Mettre le meme rcoeff pour les deux images */
//imx_norm_rcoeff_3d_p(im2,imres);
return;	
}

/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_p(grphic3d *im1,grphic3d *im2, grphic3d *im1d,grphic3d *im2d, grphic3d *im1res,grphic3d *im2res,int nb_classe)
{
grphic *histo1,*histo2;
int nb_gris=TOPO_SIZE_NORM_HISTO;
Param_Gaussienne_2d *param1,*param2;
double max1, max2, N1,N2;
grphic3d *im1dtemp,*im2dtemp,*im1temp,*im2temp;


N1=TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
N2=1-TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
 
printf("Appel de la fonction de normalisation d'intensite symetrique\n");

/* allocation memoire des structures de donnees necessaires */
param1=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param2=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
im1dtemp=cr_grphic3d(im1);
im2dtemp=cr_grphic3d(im2);
im1temp=cr_grphic3d(im1);
im2temp=cr_grphic3d(im2);


imx_norm_seuil_meanecty_symetrique_3d_p(im1, im2, im1, im2, im1temp, im2temp);
imx_norm_seuil_meanecty_symetrique_3d_p(im1d, im2d, im1d, im2d, im1dtemp, im2dtemp);


Troncature_intensite_nsigma_grphic3d_p(im1temp, im1temp, 3);
Troncature_intensite_nsigma_grphic3d_p(im2temp, im2temp, 3);
Troncature_intensite_nsigma_grphic3d_p(im1dtemp, im1dtemp, 3);
Troncature_intensite_nsigma_grphic3d_p(im2dtemp, im2dtemp, 3);




histo1=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo1->pos=0;
histo2=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo2->pos=0;


imx_inimaxminpixel_3d_p(im1dtemp);
imx_inimaxminpixel_3d_p(im2dtemp);
max1=im1dtemp->max_pixel*im1dtemp->rcoeff;
max2=im2dtemp->max_pixel*im2dtemp->rcoeff;


/* calcul de l'histo joint norme entre im1d et im2d */
histo_joint_linear_norm_3d_p(im1dtemp,im2dtemp,histo1,NULL);
trimmed_robust_histo (histo1,histo1,1.96);
imx_inimaxminpixel_p(histo1);

histo_joint_linear_norm_3d_p(im2dtemp,im1dtemp,histo2,NULL);
trimmed_robust_histo (histo2,histo2,1.96);
imx_inimaxminpixel_p(histo2);


/* Estimation du melange de gaussiennes */
fit_gaussienne2d_EM(histo1, param1, nb_classe);
fit_gaussienne2d_EM(histo2, param2, nb_classe);
 
/* mise a jour des images */
apply_normalisation_gaussienne_histo_joint(im1temp,im2temp,im1res,histo1,nb_classe,param1, max1, -1.0,1.0);
apply_normalisation_gaussienne_histo_joint(im2temp,im1temp,im2res,histo2,nb_classe,param2, max2, -1.0,1.0);

imx_mul_coe_3d_p(im1temp,N1,im1dtemp);
imx_mul_coe_3d_p(im2temp,N2,im2dtemp);

imx_mul_coe_3d_p(im1res,N2,im1res);
imx_mul_coe_3d_p(im2res,N1,im2res);


imx_add_3d_p(im1dtemp,im1res,im1res);
imx_add_3d_p(im2dtemp,im2res,im2res);


imx_norm_rcoeff_3d_p(im1res,im2res);

/*save_mri_ipb_3d_p("/home/noblet/im1.ipb", im1) ;	
save_mri_ipb_3d_p("/home/noblet/im2.ipb", im2) ;	
save_mri_ipb_3d_p("/home/noblet/im1res.ipb", im1res) ;	
save_mri_ipb_3d_p("/home/noblet/im2res.ipb", im2res) ;	*/
	
	
	
	
imx_inimaxminpixel_3d_p(im1res);imx_inimaxminpixel_3d_p(im2res);



/* liberation memoire des structures de donnees necessaires */
free_grphic3d(im1dtemp);free_grphic3d(im2dtemp);
free_grphic3d(im1temp);free_grphic3d(im2temp);
free(param1);free(param2);
free_grphic(histo1);free_grphic(histo2);
}



/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_nonsymetrise_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_nonsymetrise_p(grphic3d *im1,grphic3d *im2, grphic3d *im1d,grphic3d *im2d, grphic3d *im1res,grphic3d *im2res,int nb_classe)
{
grphic *histo1,*histo2;
int nb_gris=TOPO_SIZE_NORM_HISTO;
Param_Gaussienne_2d *param1,*param2;
double max1, max2, N1,N2;
grphic3d *im1dtemp,*im2dtemp,*im1temp,*im2temp;


N1=TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
N2=1-TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
 
printf("Appel de la fonction de normalisation d'intensite symetrique\n");

/* allocation memoire des structures de donnees necessaires */
param1=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param2=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
im1dtemp=cr_grphic3d(im1);
im2dtemp=cr_grphic3d(im2);
im1temp=cr_grphic3d(im1);
im2temp=cr_grphic3d(im2);


//imx_norm_seuil_meanecty_symetrique_3d_p(im1, im2, im1, im2, im1temp, im2temp);
//imx_norm_seuil_meanecty_symetrique_3d_p(im1d, im2d, im1d, im2d, im1dtemp, im2dtemp);


imx_norm_seuil_meanecty_3d_p(im1, im2,  im1temp);
imx_norm_seuil_meanecty_3d_p(im1d, im2d,  im1dtemp);
imx_copie_3d_p(im2, im2temp);
imx_copie_3d_p(im2d, im2dtemp);

//Troncature_intensite_nsigma_grphic3d_p(im1temp, im1temp, 3);
//Troncature_intensite_nsigma_grphic3d_p(im2temp, im2temp, 3);
//Troncature_intensite_nsigma_grphic3d_p(im1dtemp, im1dtemp, 3);
//Troncature_intensite_nsigma_grphic3d_p(im2dtemp, im2dtemp, 3);




histo1=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo1->pos=0;
histo2=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo2->pos=0;


imx_inimaxminpixel_3d_p(im1dtemp);
imx_inimaxminpixel_3d_p(im2dtemp);
max1=im1dtemp->max_pixel*im1dtemp->rcoeff;
max2=im2dtemp->max_pixel*im2dtemp->rcoeff;


/* calcul de l'histo joint norme entre im1d et im2d */
histo_joint_linear_norm_3d_p(im1dtemp,im2dtemp,histo1,NULL);
trimmed_robust_histo (histo1,histo1,1.96);
imx_inimaxminpixel_p(histo1);



/* Estimation du melange de gaussiennes */
fit_gaussienne2d_EM(histo1, param1, nb_classe);
 
/* mise a jour des images */
apply_normalisation_gaussienne_histo_joint(im1temp,im2temp,im1res,histo1,nb_classe,param1, max1, -1.0,1.0);


imx_copie_3d_p(im2temp, im2res);



imx_norm_rcoeff_3d_p(im1res,im2res);

/*save_mri_ipb_3d_p("/home/miv/noblet/validation_recalage_sym/im1.ipb", im1) ;	
save_mri_ipb_3d_p("/home/miv/noblet/validation_recalage_sym/im2.ipb", im2) ;	
save_mri_ipb_3d_p("/home/miv/noblet/validation_recalage_sym/im1res.ipb", im1res) ;	
save_mri_ipb_3d_p("/home/miv/noblet/validation_recalage_sym/im2res.ipb", im2res) ;	
save_mri_ipb_3d_p("/home/miv/noblet/validation_recalage_sym/im1d.ipb", im1d) ;	
save_mri_ipb_3d_p("/home/miv/noblet/validation_recalage_sym/im2d.ipb", im2d) ;	*/

	
	
	
	
imx_inimaxminpixel_3d_p(im1res);imx_inimaxminpixel_3d_p(im2res);



/* liberation memoire des structures de donnees necessaires */
free_grphic3d(im1dtemp);free_grphic3d(im2dtemp);
free_grphic3d(im1temp);free_grphic3d(im2temp);
free(param1);free(param2);
free_grphic(histo1);free_grphic(histo2);
}





/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_iteratif_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_norm_symetrique_iteratif_p(grphic3d *im1,grphic3d *im2, grphic3d *im1d,grphic3d *im2d, grphic3d *im1res,grphic3d *im2res,int nb_classe)
{
grphic *histo1,*histo2;
int nb_gris=TOPO_SIZE_NORM_HISTO,l;
Param_Gaussienne_2d *param1,*param2;
double max1, max2, N1,N2;
grphic3d *im1dtemp,*im2dtemp,*im1restemp,*im2restemp,*im1temp,*im2temp;


N1=TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
N2=1-TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
 
printf("Appel de la fonction de normalisation d'intensite symetrique\n");

/* allocation memoire des structures de donnees necessaires */
param1=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param2=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
im1dtemp=cr_grphic3d(im1);
im2dtemp=cr_grphic3d(im2);
im1restemp=cr_grphic3d(im1);
im2restemp=cr_grphic3d(im2);
im1temp=cr_grphic3d(im1);
im2temp=cr_grphic3d(im2);

histo1=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo1->pos=0;
histo2=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo2->pos=0;



imx_copie_3d_p(im1d,im1dtemp);
imx_copie_3d_p(im2d,im2dtemp);

imx_copie_3d_p(im1,im1temp);
imx_copie_3d_p(im2,im2temp);

for (l=0;l<10;l++) // boucle sur le nombre d'iteration
{
imx_inimaxminpixel_3d_p(im1dtemp);
imx_inimaxminpixel_3d_p(im2dtemp);
max1=im1dtemp->max_pixel*im1dtemp->rcoeff;
max2=im2dtemp->max_pixel*im2dtemp->rcoeff;


/* calcul de l'histo joint norme entre im1d et im2d */
histo_joint_linear_norm_3d_p(im1dtemp,im2dtemp,histo1,NULL);
trimmed_robust_histo (histo1,histo1,1.96);
imx_inimaxminpixel_p(histo1);

histo_joint_linear_norm_3d_p(im2dtemp,im1dtemp,histo2,NULL);
trimmed_robust_histo (histo2,histo2,1.96);
imx_inimaxminpixel_p(histo2);


/* Estimation du melange de gaussiennes */
fit_gaussienne2d_EM(histo1, param1, nb_classe);
fit_gaussienne2d_EM(histo2, param2, nb_classe);
 
/* mise a jour des images originale */
apply_normalisation_gaussienne_histo_joint(im1temp,im2temp,im1res,histo1,nb_classe,param1, max1, -1.0,1.0);
apply_normalisation_gaussienne_histo_joint(im2temp,im1temp,im2res,histo2,nb_classe,param2, max2, -1.0,1.0);

imx_mul_coe_3d_p(im1temp,N1,im1restemp);
imx_mul_coe_3d_p(im2temp,N2,im2restemp);

imx_mul_coe_3d_p(im1res,N2,im1res);
imx_mul_coe_3d_p(im2res,N1,im2res);


imx_add_3d_p(im1restemp,im1res,im1res);
imx_add_3d_p(im2restemp,im2res,im2res);


imx_norm_rcoeff_3d_p(im1res,im2res);

imx_copie_3d_p(im1res,im1temp);
imx_copie_3d_p(im2res,im2temp);


/* mise a jour des images deformees */
apply_normalisation_gaussienne_histo_joint(im1dtemp,im2dtemp,im1res,histo1,nb_classe,param1, max1, -1.0,1.0);
apply_normalisation_gaussienne_histo_joint(im2dtemp,im1dtemp,im2res,histo2,nb_classe,param2, max2, -1.0,1.0);

imx_mul_coe_3d_p(im1dtemp,N1,im1restemp);
imx_mul_coe_3d_p(im2dtemp,N2,im2restemp);

imx_mul_coe_3d_p(im1res,N2,im1res);
imx_mul_coe_3d_p(im2res,N1,im2res);


imx_add_3d_p(im1restemp,im1res,im1res);
imx_add_3d_p(im2restemp,im2res,im2res);


imx_norm_rcoeff_3d_p(im1res,im2res);

imx_copie_3d_p(im1res,im1dtemp);
imx_copie_3d_p(im2res,im2dtemp);
	
	
	
	
imx_inimaxminpixel_3d_p(im1res);imx_inimaxminpixel_3d_p(im2res);


}
imx_copie_3d_p(im1temp,im1res);
imx_copie_3d_p(im2temp,im2res);


/* liberation memoire des structures de donnees necessaires */
free_grphic3d(im1dtemp);free_grphic3d(im2dtemp);
free_grphic3d(im1temp);free_grphic3d(im2temp);
free_grphic3d(im1restemp);free_grphic3d(im2restemp);
free(param1);free(param2);
free_grphic(histo1);free_grphic(histo2);


}



/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_norm_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_norm_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
grphic *histo;
int nb_gris=TOPO_SIZE_NORM_HISTO;
/*int i,j,k,n;
n=0;
for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	if ((im1->mri[i][j][k]!=0)&&(im2->mri[i][j][k]!=0))
		n++;

nb_gris=sqrt(0.01*n);

printf ("Nombre de bins de l'histo : %d\n",nb_gris);
*/		
histo=cr_grphic_modif(nb_gris,nb_gris,0,1,0);
histo->pos=0;
normalisation_gaussienne2d_histo_joint_norm_p(im1,im2,imres,histo, nb_classe);
free_grphic(histo);
}

/*******************************************************************************
**         normalisation_gaussienne2d_histo_joint_norm_p                                                
**                                                                                                   
*******************************************************************************/

void  normalisation_gaussienne2d_histo_joint_norm_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe)
{ 

grphic3d *im1temp,*im2temp;
im1temp=cr_grphic3d(im1);
im2temp=cr_grphic3d(im2);

Troncature_intensite_nsigma_grphic3d_p(im1, im1temp, 3);
Troncature_intensite_nsigma_grphic3d_p(im2, im2temp, 3);

imx_norm_seuil_meanecty_3d_p(im1temp,im2temp,im1temp);
	

//save_mri_ipb_3d_p("/home/noblet/im1.ipb", im1temp) ;	
//save_mri_ipb_3d_p("/home/noblet/im2.ipb", im2temp) ;	


//imx_copie_3d_p(im1,im1temp);
//imx_copie_3d_p(im2,im2temp);

Param_Gaussienne_2d *param;

param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));

histo_joint_linear_norm_3d_p(im1temp,im2temp,histo,NULL);

trimmed_robust_histo (histo,histo,1.96);

imx_inimaxminpixel_p(histo);

fit_gaussienne2d_EM(histo, param, nb_classe);

apply_normalisation_gaussienne_histo_joint(im1temp,im2temp,imres,histo,nb_classe,param,-1.0,-1.0,1.0);


imx_inimaxminpixel_3d_p(imres);

free_grphic3d(im1temp);free_grphic3d(im2temp);

free(param);
}

/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_recalsse_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_reclasse_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
grphic *histo;
histo=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
//histo=cr_grphic_modif(64,64,0,1,0);
histo->pos=0;
normalisation_gaussienne2d_histo_joint_reclasse_p(im1,im2,imres,histo, nb_classe);
free_grphic(histo);
}

/*******************************************************************************
**         normalisation_gaussienne2d_histo_joint_reclasse_p                                                
**                                                                                                   
*******************************************************************************/

void  normalisation_gaussienne2d_histo_joint_reclasse_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe)
{
Param_Gaussienne_2d *param;

param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));

histo_joint_maxinfo_3d_p(im1,im2,histo);

imx_inimaxminpixel_p(histo);

fit_gaussienne2d_EM(histo, param, nb_classe);

apply_normalisation_gaussienne_histo_joint(im1,im2,imres,histo,nb_classe,param,-1.0,-1.0,1.0);

imx_inimaxminpixel_3d_p(imres);

free(param);
}

/*******************************************************************************
**      					  fit_gaussienne2d_EM                                                
**                                                                                                   
*******************************************************************************/
void fit_gaussienne2d_EM(grphic *histo, Param_Gaussienne_2d *param, int nb_classe)
{
int stop,l,i,j,k,kopt=0,nb,Ntot;
grphic *histo_cl,*histo_tmp/*,*trimmedhisto*/;
Param_Gaussienne_2d *new_param,*param_opt;
double J,Jmin;
/* allocation memoire */
histo_cl=cr_grphic_modif(histo->width,histo->height,0,1,0);
histo_tmp=cr_grphic_modif(histo->width,histo->height,histo->icomp,histo->rcoeff,histo->max_pixel);
//trimmedhisto=cr_grphic_modif(histo->width,histo->height,histo->icomp,histo->rcoeff,histo->max_pixel);

new_param= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param_opt= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));

/*for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo_tmp->mri[i][j]=histo->mri[i][j];


// initialisation 
kmeans_histo_joint_3d_p(histo,histo_cl,nb_classe);
mean_std_classe(histo,histo_cl,nb_classe, param);

stop = 0 ;
//printf("%f  \n",vraisemblance_gaussienne2d(histo, param, nb_classe));

do
{

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];

for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
		plot_iso_gaussienne_2d(histo,&param[l],1, histo->max_pixel);

if (histo->pos!=0)
show_picture(histo->pos);

//trimmed_histo_gaussienne_2d(histo_tmp,trimmedhisto,param,nb_classe,0.6);

estime_gaussienne2d_EM(histo_tmp,param,new_param,nb_classe);
//estime_gaussienne2d_isosigma_EM(histo_tmp,param,new_param,nb_classe);
//estime_gaussienne2d_isosxsy_EM(histo_tmp,param,new_param,nb_classe);

stop = 1;
for (l=0;l<nb_classe;l++)
	{
	if ((1-1.0*new_param[l].mx/param[l].mx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].my/param[l].my)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sx/param[l].sx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sy/param[l].sy)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sxy/param[l].sxy)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].l/param[l].l)>0.05)
		stop=0;	
	}

for (l=0;l<nb_classe;l++)
	{	
	param[l].mx=new_param[l].mx;
	param[l].my=new_param[l].my;
	param[l].sx=new_param[l].sx;
	param[l].sy=new_param[l].sy;
	param[l].sxy=new_param[l].sxy;
	param[l].l=new_param[l].l;
	}
//printf("%f  \n",vraisemblance_gaussienne2d(histo, param, nb_classe));

}
while (stop==0);
*/
Ntot=0;
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	Ntot+=histo->mri[i][j];


Jmin=HUGE_VAL;
for (k=nb_classe;k<=nb_classe;k++)  // dans le cas monomodal, on ne choisit pas automatiquement le meilleur nombre de classe voir fonction fit_gaussienne2d_EM_multi
{
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo_tmp->mri[i][j]=histo->mri[i][j];


// initialisation 
kmeans_histo_joint_3d_p(histo,histo_cl,k);
//kmeans_L1_histo_joint_3d_p(histo,histo_cl,k);
mean_std_classe(histo,histo_cl,k, param);

stop = 0 ;
nb=0;
do
{

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];

for (l=0;l<k;l++)
	if (param[l].l>0)
		plot_iso_gaussienne_2d(histo,&param[l],1, histo->max_pixel);

if (histo->pos!=0)
show_picture(histo->pos);


estime_gaussienne2d_EM(histo_tmp,param,new_param,k);
//robust_estime_gaussienne2d_EM(histo_tmp,param,new_param,k);


stop = 1;
for (l=0;l<k;l++)
	{
	if ((1-1.0*new_param[l].mx/param[l].mx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].my/param[l].my)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sx/param[l].sx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sy/param[l].sy)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sxy/param[l].sxy)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].l/param[l].l)>0.05)
		stop=0;	
	}

for (l=0;l<k;l++)
	{	
	param[l].mx=new_param[l].mx;
	param[l].my=new_param[l].my;
	param[l].sx=new_param[l].sx;
	param[l].sy=new_param[l].sy;
	param[l].sxy=new_param[l].sxy;
	param[l].l=new_param[l].l;
	}
//printf("%f  \n",vraisemblance_gaussienne2d(histo_tmp, param, k));

nb++;
if (nb>100)
	stop=1;
}
while (stop==0);



J=0;
/*for (l=0;l<k;l++)
	if (param[l].l>0)
	J=J+0.5*param[l].l*log(param[l].sx*param[l].sy-param[l].sxy*param[l].sxy)-param[l].l*log(param[l].l);*/

/*J=-1.0*vraisemblance_gaussienne2d(histo_tmp, param, k)+0.5*(6*k-1)*log(Ntot)*100;

printf ("k= %d  J = %f  Vraisemblance %f     Penalisation : %f\n",k,J,vraisemblance_gaussienne2d(histo_tmp, param, k),0.5*(6*k-1)*log(Ntot)); */ 
if (J<Jmin)
	{Jmin=J; kopt=k;
	for (l=0;l<k;l++)
		{	
		param_opt[l].mx=param[l].mx;
		param_opt[l].my=param[l].my;
		param_opt[l].sx=param[l].sx;
		param_opt[l].sy=param[l].sy;
		param_opt[l].sxy=param[l].sxy;
		param_opt[l].l=param[l].l;
		}
	}
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];
	
}

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];


for (l=0;l<kopt;l++)
	{	
	param[l].mx=param_opt[l].mx;
	param[l].my=param_opt[l].my;
	param[l].sx=param_opt[l].sx;
	param[l].sy=param_opt[l].sy;
	param[l].sxy=param_opt[l].sxy;
	param[l].l=param_opt[l].l;
	}

for (l=kopt;l<nb_classe;l++)
	{	
	param[l].mx=0;
	param[l].my=0;
	param[l].sx=0;
	param[l].sy=0;
	param[l].sxy=0;
	param[l].l=0;
	}

printf("Nombre de gaussienne optimal : %d \n",kopt);	

for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
	plot_iso_gaussienne_2d(histo,&param[l],1, histo->max_pixel);

if (histo->pos!=0)
show_picture(histo->pos);

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];


free_grphic(histo_cl);
free_grphic(histo_tmp);
//free_grphic(trimmedhisto);
free(new_param);
free(param_opt);
}


/*******************************************************************************
**      					  fit_gaussienne2d_EM_multi                                                
**                                                                                                   
*******************************************************************************/
void fit_gaussienne2d_EM_multi(grphic *histo, Param_Gaussienne_2d *param, int nb_classe)
{
int stop,l,i,j,k,kopt=0,Ntot;
grphic *histo_cl,*histo_tmp/*,*trimmedhisto*/;
Param_Gaussienne_2d *new_param,*param_opt;
double J,Jmin,ML;
/* allocation memoire */
histo_cl=cr_grphic_modif(histo->width,histo->height,0,1,0);
histo_tmp=cr_grphic_modif(histo->width,histo->height,histo->icomp,histo->rcoeff,histo->max_pixel);
//trimmedhisto=cr_grphic_modif(histo->width,histo->height,histo->icomp,histo->rcoeff,histo->max_pixel);

new_param= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param_opt= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));

/*for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo_tmp->mri[i][j]=histo->mri[i][j];


// initialisation 
kmeans_histo_joint_3d_p(histo,histo_cl,nb_classe);
mean_std_classe(histo,histo_cl,nb_classe, param);

stop = 0 ;
//printf("%f  \n",vraisemblance_gaussienne2d(histo, param, nb_classe));

do
{

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];

for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
		plot_iso_gaussienne_2d(histo,&param[l],1, histo->max_pixel);

if (histo->pos!=0)
show_picture(histo->pos);

//trimmed_histo_gaussienne_2d(histo_tmp,trimmedhisto,param,nb_classe,0.6);

estime_gaussienne2d_EM(histo_tmp,param,new_param,nb_classe);
//estime_gaussienne2d_isosigma_EM(histo_tmp,param,new_param,nb_classe);
//estime_gaussienne2d_isosxsy_EM(histo_tmp,param,new_param,nb_classe);

stop = 1;
for (l=0;l<nb_classe;l++)
	{
	if ((1-1.0*new_param[l].mx/param[l].mx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].my/param[l].my)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sx/param[l].sx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sy/param[l].sy)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sxy/param[l].sxy)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].l/param[l].l)>0.05)
		stop=0;	
	}

for (l=0;l<nb_classe;l++)
	{	
	param[l].mx=new_param[l].mx;
	param[l].my=new_param[l].my;
	param[l].sx=new_param[l].sx;
	param[l].sy=new_param[l].sy;
	param[l].sxy=new_param[l].sxy;
	param[l].l=new_param[l].l;
	}
//printf("%f  \n",vraisemblance_gaussienne2d(histo, param, nb_classe));

}
while (stop==0);
*/

Ntot=0;
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	Ntot+=histo->mri[i][j];



Jmin=HUGE_VAL;
for (k=1;k<=nb_classe;k++)
{
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo_tmp->mri[i][j]=histo->mri[i][j];


// initialisation 
kmeans_histo_joint_3d_p(histo,histo_cl,k);
mean_std_classe(histo,histo_cl,k, param);

stop = 0 ;
//printf("%f  \n",vraisemblance_gaussienne2d(histo, param, k));

do
{

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];

for (l=0;l<k;l++)
	if (param[l].l>0)
		plot_iso_gaussienne_2d(histo,&param[l],1, histo->max_pixel);

if (histo->pos!=0)
show_picture(histo->pos);


estime_gaussienne2d_EM(histo_tmp,param,new_param,k);

stop = 1;
for (l=0;l<k;l++)
	{
	if ((1-1.0*new_param[l].mx/param[l].mx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].my/param[l].my)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sx/param[l].sx)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sy/param[l].sy)>0.01)
		stop=0;
	if ((1-1.0*new_param[l].sxy/param[l].sxy)>0.05)
		stop=0;
	if ((1-1.0*new_param[l].l/param[l].l)>0.05)
		stop=0;	
	}

for (l=0;l<k;l++)
	{	
	param[l].mx=new_param[l].mx;
	param[l].my=new_param[l].my;
	param[l].sx=new_param[l].sx;
	param[l].sy=new_param[l].sy;
	param[l].sxy=new_param[l].sxy;
	param[l].l=new_param[l].l;
	}
//printf("%f  \n",vraisemblance_gaussienne2d(histo, param, k));

}
while (stop==0);
J=0;

/*for (l=0;l<k;l++)
	if (param[l].l>0)
			J=J+0.5*param[l].l*log(param[l].sx*param[l].sy-param[l].sxy*param[l].sxy)-param[l].l*log(param[l].l);
*/

ML=vraisemblance_gaussienne2d(histo, param, k);
J=0.5*log(Ntot)*(3*k+2)-ML; /* MDL : Minimun description length */

printf ("k= %d  J = %f\n",k,J);  
if (J<Jmin)
	{Jmin=J; kopt=k;
	for (l=0;l<k;l++)
		{	
		param_opt[l].mx=param[l].mx;
		param_opt[l].my=param[l].my;
		param_opt[l].sx=param[l].sx;
		param_opt[l].sy=param[l].sy;
		param_opt[l].sxy=param[l].sxy;
		param_opt[l].l=param[l].l;
		}
	}
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];
	
}

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];


for (l=0;l<kopt;l++)
	{	
	param[l].mx=param_opt[l].mx;
	param[l].my=param_opt[l].my;
	param[l].sx=param_opt[l].sx;
	param[l].sy=param_opt[l].sy;
	param[l].sxy=param_opt[l].sxy;
	param[l].l=param_opt[l].l;
	}

for (l=kopt;l<nb_classe;l++)
	{	
	param[l].mx=0;
	param[l].my=0;
	param[l].sx=0;
	param[l].sy=0;
	param[l].sxy=0;
	param[l].l=0;
	}

printf("Nombre de gaussienne optimal : %d \n",kopt);	

for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
	plot_iso_gaussienne_2d(histo,&param[l],1, histo->max_pixel);

if (histo->pos!=0)
show_picture(histo->pos);

free_grphic(histo_cl);
free_grphic(histo_tmp);
//free_grphic(trimmedhisto);
free(new_param);
free(param_opt);
}

/*******************************************************************************
**      					  estime_gaussienne2d_EM                                                
**                                                                                                   
*******************************************************************************/
void estime_gaussienne2d_EM(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe)
{
double ***proba_apriori;
int wdth,hght,i,j,k,Ntot;
double tmp;
wdth=histo->width;
hght=histo->height;

/* allocation memoire */
proba_apriori=alloc_dmatrix_3d(wdth,hght,nb_classe);

/* remplissage de proba_apriori */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		for (k=0;k<nb_classe;k++)
			{proba_apriori[i][j][k]=	eval_gaussienne_2d(i, j,&param[k]);
			 if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
			}

Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		Ntot=Ntot+histo->mri[i][j];

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (histo->mri[i][j]>0)
	{tmp=0;
	for (k=0;k<nb_classe;k++)
	tmp=tmp+proba_apriori[i][j][k];
	
	for (k=0;k<nb_classe;k++)
		{
		proba_apriori[i][j][k]=1.0*proba_apriori[i][j][k]/tmp;
		if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
		}
	
	}



/* mise a jour du facteur de normalisation */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].l=0.0;
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].l=new_param[k].l+1.0*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].l=1.0*new_param[k].l/Ntot;
		}

/* mise a jour des moyennes */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].mx=0.0;
		new_param[k].my=0.0;
		
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].mx=new_param[k].mx+1.0*i*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[k].my=new_param[k].my+1.0*j*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].mx=1.0*new_param[k].mx/(Ntot*new_param[k].l);
		new_param[k].my=1.0*new_param[k].my/(Ntot*new_param[k].l);		
		}

/* mise a jour des variances */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].sx=0.0;
		new_param[k].sy=0.0;
		new_param[k].sxy=0.0;
		
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].sx=new_param[k].sx+1.0*(i-new_param[k].mx)*(i-new_param[k].mx)*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[k].sy=new_param[k].sy+1.0*(j-new_param[k].my)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[k].sxy=new_param[k].sxy+1.0*(i-new_param[k].mx)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].sx=1.0*new_param[k].sx/(Ntot*new_param[k].l);
		new_param[k].sy=1.0*new_param[k].sy/(Ntot*new_param[k].l);		
		new_param[k].sxy=1.0*new_param[k].sxy/(Ntot*new_param[k].l);		
		}

		
free_dmatrix_3d(proba_apriori);
}


/*******************************************************************************
**      					  vraisemblance_gaussienne2d                                                
**                                                                                                   
*******************************************************************************/
double vraisemblance_gaussienne2d(grphic *histo, Param_Gaussienne_2d *param, int nb_classe)
{
double res,tmp;
int i,j,l,wdth,hght;
wdth=histo->width;
hght=histo->height;

res=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
	{tmp=0;
	for (l=0;l<nb_classe;l++)
		tmp=tmp+eval_gaussienne_2d(i,j,&param[l]);
	res=res+1.0*histo->mri[i][j]*log(tmp);
	}
return(res);
}

/*******************************************************************************
**      					  eval_gaussienne_2d                                                
**                                                                                                   
*******************************************************************************/
/*double eval_gaussienne_2d(double x, double y, Param_Gaussienne_2d *param)  

///  Attention la gfeneralisation pour que ca marche en 1D aussi entraine un bug 
{
double sigma,res;

if ((param->sx>0)&&(param->sy>0))
	{
	sigma=param->sx*param->sy-param->sxy*param->sxy;
	res=1.0*param->l/sqrt(sigma)*exp(-1.0*((x-param->mx)*(x-param->mx)*param->sy+(y-param->my)*(y-param->my)*param->sx
			-2.0*(x-param->mx)*(y-param->my)*param->sxy)
			/(2.0*sigma));
	}
else
	{
	if (param->sx>0)
		{
		res=1.0*param->l/sqrt(param->sx)*exp(-1.0*((x-param->mx)*(x-param->mx))/(2.0*param->sx));
		}
	if (param->sy>0)
		{
		res=1.0*param->l/sqrt(param->sy)*exp(-1.0*((y-param->my)*(y-param->my))/(2.0*param->sy));
		}
	
	if ((param->sx==0)&&(param->sy==0))
	res=0;
	
	}
return(res);

}*/

double eval_gaussienne_2d(double x, double y, Param_Gaussienne_2d *param)
{
double sigma,res;

sigma=param->sx*param->sy-param->sxy*param->sxy;
	if (sigma>0)
	res=1.0*param->l/sqrt(sigma)*exp(-1.0*((x-param->mx)*(x-param->mx)*param->sy+(y-param->my)*(y-param->my)*param->sx
			-2.0*(x-param->mx)*(y-param->my)*param->sxy)
			/(2.0*sigma));
	else
	res=0;
	
return(res);
}

double eval_mahalanobis_2d(double x, double y, Param_Gaussienne_2d *param)
{
double sigma,res;

sigma=param->sx*param->sy-param->sxy*param->sxy;
	if (sigma>0)
	res=1.0*((x-param->mx)*(x-param->mx)*param->sy+(y-param->my)*(y-param->my)*param->sx
			-2.0*(x-param->mx)*(y-param->my)*param->sxy)
			/(2.0*sigma);
	else
	res=HUGE_VAL;
	
return(res);
}

/*********************************************************************/
/*         apply_normalisation_gaussienne_histo_joint                */
/*********************************************************************/
void apply_normalisation_gaussienne_histo_joint(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param,double maxreca, double maxref, double alpha)
{
int i,j,k,l,xi;
double max_reca,max_ref,rmax_reca,rmax_ref,max;
int wdth,hght,dpth;
double x,y,proba_tot,proba_tmp,tmp;
double *ft,ftmax,ftmin;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

ft = (double *)malloc((histo->width+1)*sizeof(double));

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

if (maxreca==-1.0)
rmax_reca=max_reca*imreca->rcoeff;
else
rmax_reca=maxreca;

if (maxref==-1.0)
rmax_ref=max_ref*imref->rcoeff;
else
rmax_ref=maxref;


imx_copie_param_3d_p(imref,imres);


	
/* On applique la valeur m�iane */ 
ftmax=rmax_ref;
ftmin=HUGE_VAL;

for (xi=0;xi<=histo->width; xi++)
	{	max=0;tmp=0;proba_tot=0;
		
		for (y=0;y<histo->height;y++)
			for (l=0;l<nb_classe;l++)
			if (param[l].l>0)
				proba_tot+=eval_gaussienne_2d(xi, y, &param[l]);
				proba_tmp=0;
				proba_tot=proba_tot*0.5;
				
			y=-1;
			while (proba_tmp<proba_tot)
				{y++;
				for (l=0;l<nb_classe;l++)
				if (param[l].l>0)
					{
					proba_tmp+=eval_gaussienne_2d(xi, y, &param[l]);
					}
				}
	ft[xi]=1.0*(((1.0*y*rmax_ref)/(histo->height-1)));
	
	ftmin=MINI(ftmin,ft[xi]);
	ftmax=MAXI(ftmax,ft[xi]);
	
	
	}


imx_brukermax_3d(ftmax,ftmin,imres);

for (xi=0;xi<=histo->width; xi++)
	ft[xi]=ft[xi]/imres->rcoeff;

/* correction de la FT pour les valeurs sans �hantillons */
ft[0]=0;
ft[histo->width]=1.0*rmax_ref/imres->rcoeff;
ft[histo->width-1]=1.0*rmax_ref/imres->rcoeff;

xi=1;
while (xi<histo->width-1)
{
l=0;
for (j=0;j<histo->height;j++)
	l+=histo->mri[xi][j];

if (l>0)
	{xi++;}
else
	{
	ft[xi]=-1;
	xi++;
	/*i=xi;
	//printf("xi : %d\n",xi);
	while ((l==0)&&(i<histo->width-1))
		{i++;
		for (j=0;j<histo->height;j++)
				l+=histo->mri[i][j];
		}
	//printf("i : %d\n",i);
	for (j=xi; j<i; j++)
		ft[j]=ft[xi-1]+1.0*(ft[i]-ft[xi-1])*(j-xi+1)/(i-xi+1);

	xi=i;
	*/
	}
}



	
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if (imreca->mri[i][j][k]>0)
	{
	x=imreca->mri[i][j][k]*imreca->rcoeff/rmax_reca*(histo->width-1);
	x=MINI(histo->width-1,x);

	xi=floor(x);
	
	if (ft[xi]<0)
		{imres->mri[i][j][k]=(int)((imreca->mri[i][j][k]*imreca->rcoeff)/imres->rcoeff);
		}
	else
		{tmp=(1.0*(x-xi)*ft[xi]+1.0*(xi+1-x)*ft[xi+1])*imres->rcoeff-imreca->mri[i][j][k]*imreca->rcoeff;
		imres->mri[i][j][k]=(int)((imreca->mri[i][j][k]*imreca->rcoeff+alpha*tmp)/imres->rcoeff);}
	
	}
else imres->mri[i][j][k]=0;

imx_inimaxminpixel_3d_p(imres);

free(ft);

}

/*********************************************************************/
/*         apply_normalisation_gaussienne_histo_joint_multimodal     */
/*********************************************************************/
void apply_normalisation_gaussienne_histo_joint_multimodal(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param)
{
int i,j,k,l,lmax;
double max_reca,max_ref,rmax_reca,rmax_ref,max;
int wdth,hght,dpth,u,v;
double x,y,proba_tot,ft,tmp;
double *stdx, *stdy;
double ***proba_apriori;
double **proba_marg;



wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

stdx=(double *)malloc(nb_classe*sizeof(double));
stdy=(double *)malloc(nb_classe*sizeof(double));

for (l=0; l<nb_classe; l++)
			{
			stdx[l]=sqrt(param[l].sx);
			stdy[l]=sqrt(param[l].sy);
			}	
/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

rmax_reca=max_reca*imreca->rcoeff;
rmax_ref=max_ref*imref->rcoeff;

imx_copie_param_3d_p(imref,imres);

/* allocation memoire */
proba_apriori=alloc_dmatrix_3d(histo->width,histo->height,nb_classe);
proba_marg=alloc_dmatrix(histo->width,histo->height);

/* remplissage de proba_apriori */
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	if (histo->mri[i][j]>0)
		for (k=0;k<nb_classe;k++)
			{proba_apriori[i][j][k]=	eval_gaussienne_2d(i, j,&param[k]);
			 if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
			}


for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
if (histo->mri[i][j]>0)
	{tmp=0;
	for (k=0;k<nb_classe;k++)
	tmp=tmp+proba_apriori[i][j][k];
	
	for (k=0;k<nb_classe;k++)
		{
		proba_apriori[i][j][k]=1.0*proba_apriori[i][j][k]/tmp;
		if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
		}
	
	}


for (i=0;i<histo->width;i++)
{
proba_tot=0;
for (k=0;k<nb_classe;k++)
	{tmp=0;
	for (j=0;j<histo->height;j++)
		tmp+=proba_apriori[i][j][k];
		
		proba_marg[i][k]=tmp;
		if(isnan(proba_marg[i][k])) 
		 		proba_marg[i][k]=0;
	
	proba_tot+=proba_marg[i][k];
	}

for (k=0;k<nb_classe;k++)
	{proba_marg[i][k]=proba_marg[i][k]/proba_tot;

	if(isnan(proba_marg[i][k])) 
		 		proba_marg[i][k]=0;}

}



	



for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if ((imreca->mri[i][j][k]>0)&&(imref->mri[i][j][k]>0))
		{
		printf("%d %d %d\r",i,j,k);
		x=imreca->mri[i][j][k]/max_reca*(histo->width-1);
		u=floor(x);
 		y=imref->mri[i][j][k]/max_ref*(histo->height-1);
		v=floor(y);
 		ft=0;
		proba_tot=0;
		if (histo->mri[u][v]>0)
			{
			max=0; lmax=0;
			for (l=0; l<nb_classe; l++)
				if (param[l].l>0)
				{
				if (proba_apriori[u][v][l]>max)
					{max=proba_apriori[u][v][l]; lmax=l;}
				//ft+=1.0*((u-param[l].mx)/stdx[l]*stdy[l]+param[l].my)*proba_apriori[u][v][l];
				}
			ft+=1.0*((u-param[lmax].mx)/stdx[lmax]*stdy[lmax]+param[lmax].my);
			}
		else
			{
			max=0; lmax=0;
		
			for (l=0; l<nb_classe; l++)
				if (param[l].l>0)
				{
				if (proba_marg[u][l]>max)
					{max=proba_marg[u][l]; lmax=l;}
				//ft+=1.0*((u-param[l].mx)/stdx[l]*stdy[l]+param[l].my)*proba_marg[u][l];
				}
			ft+=1.0*((u-param[lmax].mx)/stdx[lmax]*stdy[lmax]+param[lmax].my);
		
			} 
		imres->mri[i][j][k]=(int)(((1.0*ft*rmax_ref/imres->rcoeff)/(histo->height-1)));
	/*	if (imres->mri[i][j][k]==0)
			printf("y a un bug\n");*/
		}
	else 
	imres->mri[i][j][k]=0;

free(stdx);
free(stdy);
free_dmatrix(proba_marg,histo->width,histo->height);
free_dmatrix_3d(proba_apriori);
}

/*void apply_normalisation_gaussienne_histo_joint(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param)
{
int i,j,k,l,xi;
double max_reca,max_ref,rmax_reca,rmax_ref;
int wdth,hght,dpth;
double x,y,proba_tot,proba_tmp,tmp;
//double *var_gaussienne;
double *ft;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

//var_gaussienne = (double *)malloc(nb_classe*sizeof(double));
ft = (double *)malloc((histo->width+1)*sizeof(double));

//for (i=0;i<nb_classe;i++)
//	var_gaussienne[i]=param[i].sx*param[i].sy-param[i].sxy*param[i].sxy;

max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

rmax_reca=max_reca*imreca->rcoeff;
rmax_ref=max_ref*imref->rcoeff;

imx_copie_param_3d_p(imref,imres);

 

for (xi=0;xi<=histo->width; xi++)
	{	proba_tot=0;
		tmp=0;
		for (y=0;y<histo->height;y++)
			for (l=0;l<nb_classe;l++)
			if (param[l].l>0)
		{		proba_tmp=eval_gaussienne_2d(xi, y, &param[l]);
				proba_tot=proba_tot+proba_tmp;
			tmp=tmp+proba_tmp*y;
		}
	ft[xi]=1.0*(((1.0*tmp*rmax_ref/imres->rcoeff)/(histo->height-1))/proba_tot);
	}
	
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if (imreca->mri[i][j][k]>0)
	{
	x=imreca->mri[i][j][k]/max_reca*(histo->width-1);
	xi=floor(x);
	imres->mri[i][j][k]=(int)((1.0*(x-xi)*ft[xi]+1.0*(xi+1-x)*ft[xi+1]));
	}
else imres->mri[i][j][k]=0;


	

//free(var_gaussienne);
free(ft);

}*/






/*******************************************************************************
********************************************************************************
********************************************************************************
**********  Normalisation par m�ange de gaussiennes (cas multimodal)  ********* 
********************************************************************************
********************************************************************************
********************************************************************************/

/*******************************************************************************
**        normalisation_gaussienne2d_histo_joint_multimodal                                                
**                                                                                                   
*******************************************************************************/

void  normalisation_gaussienne2d_histo_joint_multimodal(void)
{
	int im_1,im_2,im_res,im_histo,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	grphic *histo;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	im_histo=GET_PLACE("histo joint");

	nb_classe= GET_INT("nb de classes", 4, &e);
	
	TOPO_SEUIL_SIGNIFICATIVITE=GET_DOUBLE("seuil de significativite ?",0.05,&e);

	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	histo=ptr_img(im_histo);

	normalisation_gaussienne2d_histo_joint_multimodal_p(im1,im2,imres,histo,nb_classe);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
**         imx_normalisation_gaussienne2d_histo_joint_multimodal_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_gaussienne2d_histo_joint_multimodal_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
grphic *histo;

//histo=cr_grphic_modif(im1->width,im1->height,0,1,0);
histo=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
histo->pos=0;
 normalisation_gaussienne2d_histo_joint_multimodal_p(im1,im2,imres,histo, nb_classe);
 free_grphic(histo);
}

/*******************************************************************************
**         normalisation_gaussienne2d_histo_joint_multimodal_p                                                
**                                                                                                   
*******************************************************************************/
void  normalisation_gaussienne2d_histo_joint_multimodal_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe)
{
Param_Gaussienne_2d *param;
grphic3d *im1_cl,*im2_cl,*mask1,*mask2;
double *moyenne1, *std1, *moyenne2, *std2,*covariance,**confiance;
int i,j,k,*nb1,*nb2/*,nb_label,nb_label_prev*/;

/* allocations m�oires */
param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
im1_cl=cr_grphic3d(im1);
im2_cl=cr_grphic3d(im2);
moyenne1= (double *) malloc (nb_classe*sizeof(double));
std1= (double *) malloc (nb_classe*sizeof(double));
moyenne2= (double *) malloc (nb_classe*sizeof(double));
std2= (double *) malloc (nb_classe*sizeof(double));
covariance= (double *) malloc (nb_classe*sizeof(double));
nb1 = (int *) malloc (nb_classe*sizeof(int));
nb2 = (int *) malloc (nb_classe*sizeof(int));

confiance=alloc_dmatrix(histo->width, histo->height);


//histo_tmp=cr_grphic_modif(histo->width,histo->height,0,1,0);

/* calcul de l'histo joint */
histo_joint_linear_norm_3d_p(im1,im2,histo,confiance);

imx_inimaxminpixel_p(histo);
imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);


			
			


/* estimation du melange de gaussiennes */
fit_gaussienne2d_EM(histo, param, nb_classe);




/* segmentation initial a partir du melange de gaussiennes*/
init_segment_gaussienne2d(im1,im2,im1_cl, im2_cl, histo,nb_classe, param	);


im1_cl->mask=cr_grphic3d(im1_cl);imx_copie_3d_p(im1_cl,im1_cl->mask);
im2_cl->mask=cr_grphic3d(im2_cl);imx_copie_3d_p(im2_cl,im2_cl->mask);

for (i=0;i<im1_cl->width;i++)
for (j=0;j<im1_cl->height;j++)
for (k=0;k<im1_cl->depth;k++)
	{
	im1_cl->mask->mri[i][j][k]=im1_cl->mri[i][j][k];
	}

for (i=0;i<im2_cl->width;i++)
for (j=0;j<im2_cl->height;j++)
for (k=0;k<im2_cl->depth;k++)
	{
	im2_cl->mask->mri[i][j][k]=im2_cl->mri[i][j][k];
	}
				


segment_mean_std_classe(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
segment_mean_std_classe(im2,im2_cl,moyenne2,std2,nb2,nb_classe);


raffine_segment_croissance_region2(im1,im1_cl,moyenne1,std1,nb_classe);
raffine_segment_croissance_region2(im2,im2_cl,moyenne2,std2,nb_classe);
 

/* segmentation des deux images*/
//segment_gaussienne2d_recherche_voisinage(im1,im2,im1_cl, im2_cl, histo,nb_classe, param	,confiance);




/* Calcul des param�res pour chacune des classes */
segment_mean_std_classe(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
segment_mean_std_classe(im2,im2_cl,moyenne2,std2,nb2,nb_classe);
segment_covariance_classe(im1,im1_cl,im2,im2_cl,covariance,moyenne1,moyenne2,std1,std2,nb_classe);
/* Application de la normalisation d'intensit�moyenne/ecart-type par classe */
apply_norm_mean_std_classe_reg(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);

//segment_mean_std_classe_robust(im1,im1_cl,moyenne1,std1,nb_classe);
//segment_mean_std_classe_robust(im2,im2_cl,moyenne2,std2,nb_classe);

//segment_signe_covariance_classe_robuste(im1,im1_cl,im2,im2_cl,covariance,moyenne1,moyenne2,std1,std2,nb_classe,histo);
//apply_norm_mean_std_classe(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);
//apply_norm_mean_std_classe_optcov(im1,im1_cl,im2,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);
//apply_norm_regression_lineaire(im1,im1_cl,im2,imres,histo, nb_classe);


imx_copie_param_3d_p(im2,imres);
imx_inimaxminpixel_3d_p(imres);


if (imres->mask!=NULL)
	remplissage_carte_de_probabilite(im1,im1_cl,im2,im2_cl,imres,histo,moyenne1,moyenne2,std1,std2,nb1,nb2,nb_classe);

if ((im1->pos>0)&&(im1->pos<5)&&(im2->pos>0)&&(im2->pos<5))
{
 mask1 = ptr_mask_3d(im1->pos);
 mask2 = ptr_mask_3d(im2->pos);
 
mask1->max_pixel = nb_classe;
mask1->min_pixel = 0;
mask1->cutoff_max = nb_classe;
mask1->cutoff_min = 0;

#ifndef COMPILE_FOR_MEDIPY
ptr_mask_activate_3d(im1->pos, TRUE);
#endif

mask2->max_pixel = nb_classe;
mask2->min_pixel = 0;
mask2->cutoff_max = nb_classe;
mask2->cutoff_min = 0;

#ifndef COMPILE_FOR_MEDIPY
ptr_mask_activate_3d(im2->pos, TRUE);
#endif


for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	{
	mask1->mri[i][j][k]=im1_cl->mri[i][j][k];
	mask2->mri[i][j][k]=im2_cl->mri[i][j][k];
	}

 Refresh_3d();
} 

/* liberation memoire */
free_dmatrix(confiance,histo->width, histo->height);
free_grphic3d(im1_cl);
free_grphic3d(im2_cl);
free(param);
free(moyenne1);
free(moyenne2);
free(std1);
free(std2);
free(covariance);
free(nb1);
free(nb2);
}

/*******************************************************************************
**         normalisation_gaussienne2d_histo_joint_multimodal2_p                                                
**                                                    (back-up de la version qui marche)                                               
*******************************************************************************/


void  normalisation_gaussienne2d_histo_joint_multimodal2_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe)
{
Param_Gaussienne_2d *param;
grphic3d *im1_cl,*im2_cl,*mask1,*mask2;
double *moyenne1, *std1, *moyenne2, *std2,*covariance;
int i,j,k,*nb1,*nb2/*,nb_label,nb_label_prev*/;
//double sigma_gaussien=10;
//int wnd_gaussien=9;
//grphic *histo_tmp;

/* allocations m�oires */
param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));
im1_cl=cr_grphic3d(im1);
im2_cl=cr_grphic3d(im2);
moyenne1= (double *) malloc (nb_classe*sizeof(double));
std1= (double *) malloc (nb_classe*sizeof(double));
moyenne2= (double *) malloc (nb_classe*sizeof(double));
std2= (double *) malloc (nb_classe*sizeof(double));
covariance= (double *) malloc (nb_classe*sizeof(double));
nb1 = (int *) malloc (nb_classe*sizeof(int));
nb2 = (int *) malloc (nb_classe*sizeof(int));

//histo_tmp=cr_grphic_modif(histo->width,histo->height,0,1,0);

/* calcul de l'histo joint */
histo_joint_linear_norm_3d_p(im1,im2,histo,NULL);

/*for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo_tmp->mri[i][j]=histo->mri[i][j];*/

/* filtrage de l'histojoint */
//histo->cmapinfo.noColormap  = 1;
//imx_gaussian_filter_p(histo,histo, sigma_gaussien, wnd_gaussien);

imx_inimaxminpixel_p(histo);

imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);


			
			


/* estimation du melange de gaussiennes */
//fit_gaussienne2d_EM_multi(histo, param, nb_classe);
fit_gaussienne2d_EM(histo, param, nb_classe);

/*for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	histo->mri[i][j]=histo_tmp->mri[i][j];*/

/* segmentation initial a partir du melange de gaussiennes*/
//histo_joint_linear_3d_p(im1,im2,histo);
init_segment_gaussienne2d(im1,im2,im1_cl, im2_cl, histo,nb_classe, param	);


im1_cl->mask=cr_grphic3d(im1_cl);imx_copie_3d_p(im1_cl,im1_cl->mask);
im2_cl->mask=cr_grphic3d(im2_cl);imx_copie_3d_p(im2_cl,im2_cl->mask);

for (i=0;i<im1_cl->width;i++)
for (j=0;j<im1_cl->height;j++)
for (k=0;k<im1_cl->depth;k++)
	{
	im1_cl->mask->mri[i][j][k]=im1_cl->mri[i][j][k];
	}

for (i=0;i<im2_cl->width;i++)
for (j=0;j<im2_cl->height;j++)
for (k=0;k<im2_cl->depth;k++)
	{
	im2_cl->mask->mri[i][j][k]=im2_cl->mri[i][j][k];
	}
				

/* Calcul des param�res pour chacune des classes */
segment_mean_std_classe(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
segment_mean_std_classe(im2,im2_cl,moyenne2,std2,nb2,nb_classe);


/*nb_label_prev=0;
nb_label=nb_classe;

while (nb_label!=nb_label_prev) 
{
nb_label_prev=nb_label;
nb_label=merge_classe(im1,im1_cl,im2,im2_cl,moyenne1,moyenne2,std1,std2,nb1,nb2,nb_classe);
}
*/

/*  Raffinement de la segmentation par croissance de r�ion  */
//raffine_segment_MAP(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
//raffine_segment_MAP(im2,im2_cl,moyenne2,std2,nb2,nb_classe);


raffine_segment_croissance_region2(im1,im1_cl,moyenne1,std1,nb_classe);
raffine_segment_croissance_region2(im2,im2_cl,moyenne2,std2,nb_classe);
 

//raffine_segment_potts_phi_variable(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
//raffine_segment_potts_phi_variable(im2,im2_cl,moyenne2,std2,nb2,nb_classe);


//raffine_segment_ML(im1,im1_cl,moyenne1,std1,nb_classe);
//raffine_segment_ML(im2,im2_cl,moyenne2,std2,nb_classe);

//affine_segment_croissance_de_region(im1,im1_cl,nb_classe);
//raffine_segment_croissance_de_region(im2,im2_cl,nb_classe);


/* Calcul des param�res pour chacune des classes */
segment_mean_std_classe(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
segment_mean_std_classe(im2,im2_cl,moyenne2,std2,nb2,nb_classe);


segment_covariance_classe(im1,im1_cl,im2,im2_cl,covariance,moyenne1,moyenne2,std1,std2,nb_classe);


/* Application de la normalisation d'intensit�moyenne/ecart-type par classe */
//apply_norm_mean_std_classe(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);
//apply_norm_mean_std_classe_vp(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);
apply_norm_mean_std_classe_reg(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);
//apply_norm_quantile_classe(im1,im1_cl,im2,im2_cl,imres,covariance, nb_classe, nb_quantile);
//apply_norm_mean_std_classe_proba_reg(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param,nb_classe);

/*
for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	{
	imres->mri[i][j][k]=imres->mri[i][j][k]-im2->mri[i][j][k];
	}
*/
imx_copie_param_3d_p(im2,imres);

//segment_mean_std_classe(imres,im2_cl,moyenne2,std2,nb_classe);


imx_inimaxminpixel_3d_p(imres);

//apply_normalisation_gaussienne_histo_joint2(im1,im2,imres,histo,nb_classe,param);




if (imres->mask!=NULL)
	remplissage_carte_de_probabilite(im1,im1_cl,im2,im2_cl,imres,histo,moyenne1,moyenne2,std1,std2,nb1,nb2,nb_classe);







if ((im1->pos>0)&&(im1->pos<5)&&(im2->pos>0)&&(im2->pos<5))
{
 mask1 = ptr_mask_3d(im1->pos);
 mask2 = ptr_mask_3d(im2->pos);
 
mask1->max_pixel = nb_classe;
mask1->min_pixel = 0;
mask1->cutoff_max = nb_classe;
mask1->cutoff_min = 0;

#ifndef COMPILE_FOR_MEDIPY
ptr_mask_activate_3d(im1->pos, TRUE);
#endif

mask2->max_pixel = nb_classe;
mask2->min_pixel = 0;
mask2->cutoff_max = nb_classe;
mask2->cutoff_min = 0;

#ifndef COMPILE_FOR_MEDIPY
ptr_mask_activate_3d(im2->pos, TRUE);
#endif



for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	{
	mask1->mri[i][j][k]=im1_cl->mri[i][j][k];
	mask2->mri[i][j][k]=im2_cl->mri[i][j][k];
	}

 Refresh_3d();
} 


 
/* liberation memoire */
free_grphic3d(im1_cl);
free_grphic3d(im2_cl);
free(param);
free(moyenne1);
free(moyenne2);
free(std1);
free(std2);
free(covariance);
free(nb1);
free(nb2);

 //free_grphic(histo_tmp);
}


/*********************************************************************/
/*      				calcul_fonction_partition_potts				  		         */
/*********************************************************************/
/*void calcul_fonction_partition_potts (void)
{
int im_1,nb_classe,err;
double phi;
grphic3d *im1;

nb_classe= GET_INT("Nombre de classe", 10, &err);
phi = GET_DOUBLE("Phi",2,&err);

im_1 = GET_PLACE3D(TEXT0231);
im1=ptr_img_3d(im_1);

swendsen_wang(im1,phi,nb_classe);
 Refresh_3d();
show_picture_3d(im_1);

}
*/
/*********************************************************************/
/*      								calcul_Uz				  									         */
/*********************************************************************/
int calcul_Uz (grphic3d *im,grphic3d *mask)
{
int i,j,k,wdth,hght,dpth,Uz;

wdth=im->width-1;
hght=im->height-1;
dpth=im->depth-1;

Uz=0;
for (i=0;i<im->width;i++)
for (j=0;j<im->height;j++)
for (k=0;k<im->depth;k++)
	if (mask->mri[i][j][k]>-1)
		{
		if (i<wdth)
		if (im->mri[i+1][j][k]==im->mri[i][j][k]) Uz++;
		
		if (j<hght)
		if (im->mri[i][j+1][k]==im->mri[i][j][k]) Uz++;
		
		if (k<dpth)
		if (im->mri[i][j][k+1]==im->mri[i][j][k]) Uz++;
		}

return(Uz);

}
/*********************************************************************/
/*      				calcul_fonction_partition_potts				  		         */
/*********************************************************************/
void calcul_fonction_partition_potts(void)
{
double phi_max, incr,*Z,reel;
int nb_classe,nb_elt,entier,err;
char     nom[256];
FILE     *fichier;
HEADER header;

/*ecriture du fichier*/
/*carateristiques du fichier*/   
header.exists=0;
header.maxlines = 2000;
header.numlines = 0;
header.lines = (char **)malloc(2000 * sizeof(char *));;
  


phi_max = GET_DOUBLE("Phi_max",2.0,&err);
incr = GET_DOUBLE("Increment",0.1,&err);
nb_classe= GET_INT("Nombre de classe", 10, &err);


nb_elt=floor(1.0*phi_max/incr)+1;
Z=(double *) malloc(nb_elt*sizeof(double));

calcul_fonction_partition_potts_p(Z, incr, phi_max ,nb_classe);



strcpy(nom,"fonction_repartition_potts.header");
fichier=fopen(nom,"r");
if (fichier!=NULL) {fclose(fichier); remove(nom);}


reel = phi_max;	put_header(&header, "phi_max=", REEL, &reel, (int)NULL);
reel = incr;	put_header(&header, "incr=", REEL, &reel, (int)NULL);
entier =nb_elt;put_header(&header, "nb_elt=", ENTIER, &entier,(int)NULL);
entier =nb_classe;put_header(&header, "nb_classe=", ENTIER, &entier,(int)NULL);

save_header(&header/*???*/,nom);

strcpy(nom,"fonction_repartition_potts.param");
fichier=fopen(nom,"wb");
imx_fwrite(Z,nb_elt,sizeof(double),fichier);
fclose(fichier);
  

free(header.lines);			
free(Z);
}

/*********************************************************************/
/*      				calcul_fonction_partition_potts				  		         */
/*********************************************************************/
void calcul_fonction_partition_potts_p(double *Z, double incr, double phi_max ,int nb_classe)
{
int i,j,k,l,u,v,w,iter,nb_chaine=4,Nmax=100000,size=128,**Uz,stop,nb,nb_per_loop=10,nb_elt,debut,fin,total,nb_clique,size1;
double tirage,phi, *EUz, mean, precision=0.001,*Ztmp,mean_prec;
grphic3d **chaine;
grphic3d *im1,*im2,*im3,*im4;

size1=size-1;

im1=ptr_img_3d(1);
im2=ptr_img_3d(2);
im3=ptr_img_3d(3);
im4=ptr_img_3d(4);



chaine=cr_ptr_grphic3d_array(nb_chaine, size, size, size);
Uz=alloc_imatrix(Nmax, nb_chaine);
nb_elt=floor(1.0*phi_max/incr)+1;
EUz= (double *) malloc(nb_chaine*sizeof(double));
Ztmp= (double *) malloc(nb_elt*sizeof(double));

nb_clique=0;
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	if (i<size1)
		nb_clique++;
	if (j<size1)
		nb_clique++;
	if (k<size1)
		nb_clique++;	
	}

/*  Initialisation des chaines  */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	tirage=drand48();
	chaine[0]->mri[i][j][k]=(int)floor(tirage*nb_classe);
	
	tirage=drand48();
	chaine[1]->mri[i][j][k]=(int)floor(tirage*nb_classe);
	
	tirage=drand48();
	chaine[2]->mri[i][j][k]=(int)floor(tirage*nb_classe);
	
	tirage=drand48();
	chaine[3]->mri[i][j][k]=(int)floor(tirage*nb_classe);
	}

debut=(int)(0.25*size);
fin=(int)(0.75*size);

printf("debut %d  fin %d \n",debut,fin);

for (i=debut;i<fin;i++)
for (j=debut;j<fin;j++)
for (k=debut;k<fin;k++)
	{
	chaine[1]->mri[i][j][k]=0;
	}



for (i=debut;i<fin;i++)
for (j=debut;j<fin;j++)
for (k=debut;k<fin;k++)
	{
	chaine[2]->mri[i][j][k]=(int)floor(i*nb_classe/size);
	}

for (i=debut;i<fin;i++)
for (j=debut;j<fin;j++)
for (k=debut;k<fin;k++)
	{
	chaine[3]->mri[i][j][k]=(int)(k%nb_classe);
	}

for (i=0;i<Nmax;i++)
	for (j=0;j<nb_chaine;j++)
		Uz[i][j]=0;


/* copie dans les visu */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
		{
		im1->mri[i][j][k]=chaine[0]->mri[i][j][k];
		im2->mri[i][j][k]=chaine[1]->mri[i][j][k];
		im3->mri[i][j][k]=chaine[2]->mri[i][j][k];
		im4->mri[i][j][k]=chaine[3]->mri[i][j][k];
		}
Refresh_3d();

for (iter=0; iter<nb_elt;iter++)
{
phi=1.0*iter*incr;


mean_prec=0;		
stop=1;
nb=0;
while (stop>0)
{

/* update de la chaine */
for (l=0;l<nb_per_loop;l++)
	{j=nb*nb_per_loop+l;
	for (i=0;i<nb_chaine;i++)
		{
		swendsen_wang(chaine[i],phi, nb_classe);
		Uz[j][i]=  calcul_Uz (chaine[i],chaine[i]);
		}
		
		
		/* copie dans les visu */
		for (u=0;u<size;u++)
		for (v=0;v<size;v++)
		for (w=0;w<size;w++)
			{
			im1->mri[u][v][w]=chaine[0]->mri[u][v][w];
			im2->mri[u][v][w]=chaine[1]->mri[u][v][w];
			im3->mri[u][v][w]=chaine[2]->mri[u][v][w];
			im4->mri[u][v][w]=chaine[3]->mri[u][v][w];
			}	
			Refresh_3d();
				
	printf("Phi %f Nb iter %d  \t chaine1 %d \t chaine2 %d \t chaine3 %d \t chaine4 %d \n",phi,j,Uz[j][0],Uz[j][1],Uz[j][2],Uz[j][3]);
	}
nb++;

/* calcul de l'esperance sur 80% des elements*/
	for (i=0;i<nb_chaine;i++)
		EUz[i]=0.0;
	
	debut=(int)floor(nb*nb_per_loop*0.2);
	fin=nb*nb_per_loop;
	total=0;
	
	for (i=debut;i<fin;i++)
		{
		total++;
		for (j=0;j<nb_chaine;j++)
			EUz[j]+=Uz[i][j];
		}
	
		for (j=0;j<nb_chaine;j++)
			EUz[j]=EUz[j]/total;
	
		mean=0;
		for (j=0;j<nb_chaine;j++)
			mean+=EUz[j];
		
		mean=mean/nb_chaine;
		
		stop=0;
		for (j=0;j<nb_chaine;j++)
			if  (fabs((EUz[j]-mean)/mean)<precision)
			stop++;
		
		printf("Phi : %f Esperance \t chaine1 %f \t chaine2 %f \t chaine3 %f \t chaine4 %f \n",phi,EUz[0],EUz[1],EUz[2],EUz[3]);

		printf("nombre de chaine verifiant la condition %d / %d  \n",stop,nb_chaine);
		
	if ((stop==nb_chaine)&&(fabs((mean_prec-mean)/mean)<precision))
			stop=0;
		else {stop=1;	mean_prec=mean;}	

		
}

Ztmp[iter]=mean;
}


/* integration */
Z[0]=Ztmp[0];
for (l=1;l<nb_elt;l++)
	{
	Z[l]=Z[l-1]+0.5*(Ztmp[l]+Ztmp[l-1])*incr;
	}


/* normalisation par le nombre de clique */
for (l=0;l<nb_elt;l++)
	Z[l]=Z[l]/nb_clique;
	


for (l=0;l<nb_elt;l++)
{phi=1.0*l*incr;
printf("Phi  %f \t Z: %f  Nombre de voisin identique moyen : %f\n",phi,Z[l],3*Ztmp[l]/nb_clique);
}


free(EUz);free(Ztmp);
free_imatrix(Uz, Nmax, nb_chaine);
free_ptr_grphic3d_array(chaine,nb_chaine);
}


/*********************************************************************/
/*      				calcul_fonction_partition_potts				  		         */
/*********************************************************************/
void calcul_fonction_partition_potts2 (void)
{
double phi_max, incr,*EU,*Z,tirage,phi,*proba,proba_tot,moyenne,reel;
int nb_classe,nb_elt,Nmoyennage=10,entier,nb_clique,err;
int i,j,k,l,t,size=128,stop,*voisinage,size1,nb,newEU1,newEU_prec1,newEU2,newEU_prec2,im_1,im_2;
grphic3d *im1,*im2;
char     nom[256];
FILE     *fichier;
HEADER header;

 /*ecriture du fichier*/
  /*carateristiques du fichier*/   
   header.exists=0;
   header.maxlines = 2000;
   header.numlines = 0;
   header.lines = (char **)malloc(2000 * sizeof(char *));;
  

size1=size-1;

phi_max = GET_DOUBLE("Phi_max",2,&err);
incr = GET_DOUBLE("Increment",0.1,&err);
nb_classe= GET_INT("Nombre de classe", 10, &err);


strcpy(nom,"fonction_repartition_potts.header");
//fichier=fopen(nom,"r");
//if (fichier!=NULL) {fclose(fichier); remove(nom);}



nb_elt=floor(1.0*phi_max/incr)+1;
EU=(double *) malloc(nb_elt*sizeof(double));
Z=(double *) malloc(nb_elt*sizeof(double));

for (i=0; i<nb_elt; i++)
	{EU[i]=0.0;Z[i]=0.0;}


voisinage=(int*) malloc(nb_classe*sizeof(int));
proba=(double*) malloc(nb_classe*sizeof(double));
im_1 = GET_PLACE3D(TEXT0231);
im_2 = GET_PLACE3D(TEXT0031);
	
im1=ptr_img_3d(im_1);	
im2=ptr_img_3d(im_2);

//im1 = cr_grphic3d_modif(size, size, size, 0.0, 1, nb_classe-1);
//im2 = cr_grphic3d_modif(size, size, size, 0.0, 1, nb_classe-1);

im1->width=size;
im1->height=size;
im1->depth=size;
im1->icomp=0;
im1->rcoeff=1;
im1->max_pixel=nb_classe;
im1->min_pixel=0;
im1->x3dc=(int)size/2;
im1->y3dc=(int)size/2;
im1->z3dc=(int)size/2;


im2->width=size;
im2->height=size;
im2->depth=size;
im2->icomp=0;
im2->rcoeff=1;
im2->max_pixel=nb_classe;
im2->min_pixel=0;
im2->x3dc=(int)size/2;
im2->y3dc=(int)size/2;
im2->z3dc=(int)size/2;


//raffine_segment_potts_phi_variable(im1,im2,NULL,NULL,nb_classe);

EU[0]=log(nb_classe);
EU[0]=0;

/* initialisation */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	im2->mri[i][j][k]=0;
	tirage=drand48()*nb_classe;
	im1->mri[i][j][k]=floor(tirage);
	}

Refresh_3d();

for (l=1;l<nb_elt;l++)
{
phi=1.0*l*incr;

/* initialisation */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	im2->mri[i][j][k]=0;
	tirage=drand48()*nb_classe;
	im1->mri[i][j][k]=floor(tirage);
	}

stop=0;
nb=0; 

nb_clique=0;
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	if (i<size1)
		nb_clique++;
	if (j<size1)
		nb_clique++;
	if (k<size1)
		nb_clique++;	
	}
printf ("Nombre de clique : %d \n",nb_clique);

newEU1=0;newEU2=0;
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	
	if (i<size1)
	if (im1->mri[i][j][k]==im1->mri[i+1][j][k])
		newEU1++;
	
	if (j<size1)
	if (im1->mri[i][j][k]==im1->mri[i][j+1][k])
		newEU1++;
	
	if (k<size1)
	if (im1->mri[i][j][k]==im1->mri[i][j][k+1])
		newEU1++;

	if (i<size1)
	if (im2->mri[i][j][k]==im2->mri[i+1][j][k])
		newEU2++;

	if (j<size1)
	if (im2->mri[i][j][k]==im2->mri[i][j+1][k])
		newEU2++;

	if (k<size1)
	if (im2->mri[i][j][k]==im2->mri[i][j][k+1])
		newEU2++;
	}
	
printf("Phi  %f \t Nombre iteration : %d \t   Esperance moyenne EU1:  %f \t   Esperance moyenne EU2:  %f \n",phi,nb,1.0*newEU1/(size1*size1*size1),1.0*newEU2/(size1*size1*size1));

newEU_prec1=0;
newEU_prec2=0;
while (stop==0)
{

/* tirage d'une nouvelle carte de segmentation1 */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	for (t=0; t<nb_classe; t++)
		voisinage[t]=0;
	
	if (i>0) voisinage[im1->mri[i-1][j][k]]++;
	if (j>0) voisinage[im1->mri[i][j-1][k]]++;
	if (k>0) voisinage[im1->mri[i][j][k-1]]++;
	if (i<size1) voisinage[im1->mri[i+1][j][k]]++;
	if (j<size1) voisinage[im1->mri[i][j+1][k]]++;
	if (k<size1) voisinage[im1->mri[i][j][k+1]]++;
	
	for (t=0; t<nb_classe; t++)
		proba[t]=exp(phi*voisinage[t]);

			proba_tot=0;
			for (t=0;t<nb_classe;t++)
				proba_tot+=proba[t];
							
			for (t=0;t<nb_classe;t++)
				proba[t]=proba[t]/proba_tot;
								
			for (t=1;t<nb_classe;t++)
				proba[t]=proba[t]+proba[t-1];
								
			tirage=drand48();
						
			t=0;
			while((tirage>proba[t])&&(t<(nb_classe-1)))
							t++;
							
			im1->mri[i][j][k]=t;
							
	}

/* tirage d'une nouvelle carte de segmentation2 */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	for (t=0; t<nb_classe; t++)
		voisinage[t]=0;
	
	if (i>0) voisinage[im2->mri[i-1][j][k]]++;
	if (j>0) voisinage[im2->mri[i][j-1][k]]++;
	if (k>0) voisinage[im2->mri[i][j][k-1]]++;
	if (i<size1) voisinage[im2->mri[i+1][j][k]]++;
	if (j<size1) voisinage[im2->mri[i][j+1][k]]++;
	if (k<size1) voisinage[im2->mri[i][j][k+1]]++;
	
	for (t=0; t<nb_classe; t++)
		proba[t]=exp(phi*voisinage[t]);

			proba_tot=0;
			for (t=0;t<nb_classe;t++)
				proba_tot+=proba[t];
							
			for (t=0;t<nb_classe;t++)
				proba[t]=proba[t]/proba_tot;
								
			for (t=1;t<nb_classe;t++)
				proba[t]=proba[t]+proba[t-1];
								
			tirage=drand48();
						
			t=0;
			while((tirage>proba[t])&&(t<(nb_classe-1)))
							t++;
							
			im2->mri[i][j][k]=t;
							
	}
/* calcul de l'esp�ance de U */
newEU1=0;
newEU2=0;
for (i=0;i<size1;i++)
for (j=0;j<size1;j++)
for (k=0;k<size1;k++)
	{
	if (im1->mri[i][j][k]==im1->mri[i+1][j][k])
		newEU1++;
	if (im1->mri[i][j][k]==im1->mri[i][j+1][k])
		newEU1++;
	if (im1->mri[i][j][k]==im1->mri[i][j][k+1])
		newEU1++;
	
	if (im2->mri[i][j][k]==im2->mri[i+1][j][k])
		newEU2++;
	if (im2->mri[i][j][k]==im2->mri[i][j+1][k])
		newEU2++;
	if (im2->mri[i][j][k]==im2->mri[i][j][k+1])
		newEU2++;
	}
nb++;
Refresh_3d();

printf("Phi  %f \t Nombre iteration : %d \t   Esperance moyenne EU1:  %f \t   Esperance moyenne EU2:  %f \n",phi,nb,1.0*newEU1/(size1*size1*size1),1.0*newEU2/(size1*size1*size1));

if ((fabs(1.0*(newEU_prec1-newEU1)/newEU1)<0.001)&&(fabs(1.0*(newEU_prec2-newEU2)/newEU2)<0.001)&&(fabs(1.0*(newEU1-newEU2)/(newEU1+newEU2))<0.001))
	{stop=1;
	printf("condition 1 : %f \t  condition 2 : %f \t  condition 3 : %f	\n",fabs(1.0*(newEU_prec1-newEU1)/newEU1),fabs(1.0*(newEU_prec2-newEU2)/newEU2),fabs(1.0*(newEU1-newEU2)/(newEU1+newEU2)));	
	}
	
	
newEU_prec1=newEU1;
newEU_prec2=newEU2;
	
}

moyenne=0;
for (stop=0; stop<Nmoyennage; stop++)
{

/* tirage d'une nouvelle carte de segmentation1 */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	for (t=0; t<nb_classe; t++)
		voisinage[t]=0;
	
	if (i>0) voisinage[im1->mri[i-1][j][k]]++;
	if (j>0) voisinage[im1->mri[i][j-1][k]]++;
	if (k>0) voisinage[im1->mri[i][j][k-1]]++;
	if (i<size1) voisinage[im1->mri[i+1][j][k]]++;
	if (j<size1) voisinage[im1->mri[i][j+1][k]]++;
	if (k<size1) voisinage[im1->mri[i][j][k+1]]++;
	
	for (t=0; t<nb_classe; t++)
		proba[t]=exp(phi*voisinage[t]);

			proba_tot=0;
			for (t=0;t<nb_classe;t++)
				proba_tot+=proba[t];
							
			for (t=0;t<nb_classe;t++)
				proba[t]=proba[t]/proba_tot;
								
			for (t=1;t<nb_classe;t++)
				proba[t]=proba[t]+proba[t-1];
								
			tirage=drand48();
						
			t=0;
			while((tirage>proba[t])&&(t<(nb_classe-1)))
							t++;
							
			im1->mri[i][j][k]=t;
							
	}

/* tirage d'une nouvelle carte de segmentation2 */
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	for (t=0; t<nb_classe; t++)
		voisinage[t]=0;
	
	if (i>0) voisinage[im2->mri[i-1][j][k]]++;
	if (j>0) voisinage[im2->mri[i][j-1][k]]++;
	if (k>0) voisinage[im2->mri[i][j][k-1]]++;
	if (i<size1) voisinage[im2->mri[i+1][j][k]]++;
	if (j<size1) voisinage[im2->mri[i][j+1][k]]++;
	if (k<size1) voisinage[im2->mri[i][j][k+1]]++;
	
	for (t=0; t<nb_classe; t++)
		proba[t]=exp(phi*voisinage[t]);

			proba_tot=0;
			for (t=0;t<nb_classe;t++)
				proba_tot+=proba[t];
							
			for (t=0;t<nb_classe;t++)
				proba[t]=proba[t]/proba_tot;
								
			for (t=1;t<nb_classe;t++)
				proba[t]=proba[t]+proba[t-1];
								
			tirage=drand48();
						
			t=0;
			while((tirage>proba[t])&&(t<(nb_classe-1)))
							t++;
							
			im2->mri[i][j][k]=t;
							
	}
/* calcul de l'esp�ance de U */
newEU1=0;
newEU2=0;
for (i=0;i<size1;i++)
for (j=0;j<size1;j++)
for (k=0;k<size1;k++)
	{
	if (im1->mri[i][j][k]==im1->mri[i+1][j][k])
		newEU1++;
	if (im1->mri[i][j][k]==im1->mri[i][j+1][k])
		newEU1++;
	if (im1->mri[i][j][k]==im1->mri[i][j][k+1])
		newEU1++;
	
	if (im2->mri[i][j][k]==im2->mri[i+1][j][k])
		newEU2++;
	if (im2->mri[i][j][k]==im2->mri[i][j+1][k])
		newEU2++;
	if (im2->mri[i][j][k]==im2->mri[i][j][k+1])
		newEU2++;
	}

moyenne+=newEU1+newEU2;
printf("Phi  %f \t Calcul de la moyenne sur les echantillons: %d / %d \t   EU1:  %f \t   EU2:  %f \n",phi,stop,Nmoyennage,1.0*newEU1/(size1*size1*size1),1.0*newEU2/(size1*size1*size1));


	
newEU_prec1=newEU1;
newEU_prec2=newEU2;
	
}

moyenne=0.5*moyenne/Nmoyennage;
EU[l]=moyenne/nb_clique;
printf("Phi %f   \t  esperance moyenne de U(z) %f\n",phi,EU[l]);
Refresh_3d();

}

/* on int�re */
Z[0]=EU[0];
for (l=1;l<nb_elt;l++)
	{
	Z[l]=Z[l-1]+0.5*(EU[l]+EU[l-1])*incr;
	}

for (l=0;l<nb_elt;l++)
{phi=1.0*l*incr;
printf("Phi  %f \t Z: %f  \n",phi,Z[l]);

}

reel = phi_max;	put_header(&header, "phi_max=", REEL, &reel, (int)NULL);
reel = incr;	put_header(&header, "incr=", REEL, &reel, (int)NULL);
entier =nb_elt;put_header(&header, "nb_elt=", ENTIER, &entier,(int)NULL);
entier =nb_classe;put_header(&header, "nb_classe=", ENTIER, &entier,(int)NULL);

save_header(&header/*???*/,nom);

strcpy(nom,"fonction_repartition_potts.param");
fichier=fopen(nom,"wb");
imx_fwrite(Z,nb_elt,sizeof(double),fichier);
fclose(fichier);
  
//save_mri_ipb_3d_p("/home/noblet/test_phi_08",im1);

free(header.lines);			
free_grphic3d(im1);
free_grphic3d(im2);
free(EU);free(Z);
free(voisinage);free(proba);
}

/**********************************************************************/
/*      					  raffine_segment_potts_phi_variable							  */
/**********************************************************************/

void raffine_segment_potts_phi_variable(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int *nb_tot,int nb_classe)
{
int i,j,k,l,lopt,wdth,hght,dpth,wdth1,hght1,dpth1,nb_clique;
int *voisinage,nb,nb_elt,nb_classe2;
double phi=1,*proba,proba_max,tmp,potts,vraisemblance,tirage,proba_tot,phi_max,incr,*Z,*proba_phi,Uz;
HEADER header;
char nom[256];
FILE     *fichier;

/* recuperation des parametres des fichiers fonction_repartition_potts.header et fonction_repartition_potts.param  */

strcpy(nom,"fonction_repartition_potts.header");
fichier=fopen(nom,"r");

if (fichier==NULL)
  {  printf("Error in opening fonction_repartition_potts.header in load_field\n");
    exit(1);
  }

load_header(nom,&header/*???*/);
phi_max=atof(getheader_interfile(nom,"phi_max=",REEL,(int)NULL));
incr=atof(getheader_interfile(nom,"incr=",REEL,(int)NULL));
nb_elt=atoi(getheader_interfile(nom,"nb_elt=",ENTIER,(int)NULL));
nb_classe2=atoi(getheader_interfile(nom,"nb_classe=",ENTIER,(int)NULL));

if (nb_classe!=nb_classe2)
	{printf("Attention, le nombre de classes dans le fichier fonction_repartition_potts.header ne correspond au nombre de classes de la segmentation courante !!!\n");
		exit(1);	
	}
fclose(fichier);


Z=(double*) malloc(nb_elt*sizeof(double));
proba_phi=(double*) malloc(nb_elt*sizeof(double));

strcpy(nom,"fonction_repartition_potts.param");
fichier=fopen(nom,"rb");

if (fichier==NULL)
  {  printf("Error in opening fonction_repartition_potts.param in load_field\n");
   exit(1);
  }

imx_fread(Z,nb_elt,sizeof(double),fichier);	

fclose(fichier);



wdth=im->width;
hght=im->height;
dpth=im->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;


voisinage=(int *)malloc(nb_classe*sizeof(int));
proba=(double *)malloc(nb_classe*sizeof(double));



nb_clique=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im_cl->mask->mri[i][j][k]<0)
	{
	if (i<wdth1)
		nb_clique++;
	if (j<hght1)
		nb_clique++;
	if (k<dpth1)
		nb_clique++;	
	}

for (l=0;l<nb_elt;l++)
	Z[l]=Z[l]*nb_clique; 

nb=0;
do
{


Uz=0;
/* Calcul de Uz */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im_cl->mask->mri[i][j][k]<0)
		{
		if (i<wdth1)
		if (im_cl->mri[i+1][j][k]==im_cl->mri[i][j][k]) Uz++;
		
		if (j<hght1)
		if (im_cl->mri[i][j+1][k]==im_cl->mri[i][j][k]) Uz++;
		
		if (k<dpth1)
		if (im_cl->mri[i][j][k+1]==im_cl->mri[i][j][k]) Uz++;
		}

/* Estimation de phi  */
for (l=0;l<nb_elt;l++)
		proba_phi[l]=Uz*l*incr-Z[l];

proba_max=proba_phi[0];
lopt=0;
for (l=1;l<nb_elt;l++)
		if (proba_phi[l]>proba_max)
			{proba_max=proba_phi[l];lopt=l;}
		
		phi=lopt*incr;
			
		for (i=1;i<wdth1;i++)
			for (j=1;j<hght1;j++)
			for (k=1;k<dpth1;k++)
				if (im_cl->mask->mri[i][j][k]<0)
					{
						
						for (l=0;l<nb_classe;l++)
							voisinage[l]=0;
							
						if (im_cl->mri[i-1][j][k]>0) voisinage[im_cl->mri[i-1][j][k]-1]++;
						if (im_cl->mri[i][j-1][k]>0) voisinage[im_cl->mri[i][j-1][k]-1]++;
						if (im_cl->mri[i][j][k-1]>0) voisinage[im_cl->mri[i][j][k-1]-1]++;
						if (im_cl->mri[i+1][j][k]>0) voisinage[im_cl->mri[i+1][j][k]-1]++;
						if (im_cl->mri[i][j+1][k]>0) voisinage[im_cl->mri[i][j+1][k]-1]++;
						if (im_cl->mri[i][j][k+1]>0) voisinage[im_cl->mri[i][j][k+1]-1]++;
	
						for (l=0;l<nb_classe;l++)
							{
							potts=phi*voisinage[l];
							vraisemblance=-0.5*(im->mri[i][j][k]-moyenne[l])*(im->mri[i][j][k]-moyenne[l])/(std[l]*std[l])-log(std[l]);
							tmp=exp(potts+vraisemblance);
							if (isnan(tmp))
								proba[l]=0;
							else
							proba[l]=tmp;
							}
						
									
						proba_tot=0;
						for (l=0;l<nb_classe;l++)
							proba_tot+=proba[l];
							
						for (l=0;l<nb_classe;l++)
								proba[l]=proba[l]/proba_tot;
								
						for (l=1;l<nb_classe;l++)
								proba[l]=proba[l]+proba[l-1];
								
						tirage=drand48();
						
						l=0;
						while((tirage>proba[l])&&(l<(nb_classe-1)))
							l++;
						
						im_cl->mri[i][j][k]=l+1;
						
					
					}
nb++;
segment_mean_std_classe(im,im_cl,moyenne,std,nb_tot,nb_classe);					
printf("Phi %f    Nombre d'iteration %d \r",phi,nb);
}while (nb<40);




free(voisinage);
free(proba);
free(Z);free(proba_phi);
}




/*********************************************************************/
/*      					  raffine_segment_potts										         */
/*********************************************************************/

void raffine_segment_potts(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int nb_classe)
{
int i,j,k,l,lopt=0,wdth,hght,dpth,wdth1,hght1,dpth1;
int *voisinage,*nb_tot,nb,nb_voisin;
double psi=2,*proba,tmp,potts,vraisemblance,tirage,proba_tot;
grphic3d **label;



wdth=im->width;
hght=im->height;
dpth=im->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;

label=cr_ptr_grphic3d_array(nb_classe, wdth, hght, dpth);

		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
		for (l=0;l<nb_classe;l++)
			label[l]->mri[i][j][k]=0;

voisinage=(int *)malloc(nb_classe*sizeof(int));
nb_tot=(int *)malloc(nb_classe*sizeof(int));
proba=(double *)malloc(nb_classe*sizeof(double));

nb=0;
do
{
		for (i=1;i<wdth1;i++)
			for (j=1;j<hght1;j++)
			for (k=1;k<dpth1;k++)
				if (im_cl->mask->mri[i][j][k]!=0)
					{
					nb_voisin=0;
					if ((im_cl->mri[i-1][j][k]>0)&&(im_cl->mri[i-1][j][k]!=im_cl->mri[i][j][k])) nb_voisin++;
					if ((im_cl->mri[i][j-1][k]>0)&&(im_cl->mri[i][j-1][k]!=im_cl->mri[i][j][k])) nb_voisin++;
					if ((im_cl->mri[i][j][k-1]>0)&&(im_cl->mri[i][j][k-1]!=im_cl->mri[i][j][k])) nb_voisin++;
					if ((im_cl->mri[i+1][j][k]>0)&&(im_cl->mri[i+1][j][k]!=im_cl->mri[i][j][k])) nb_voisin++;
					if ((im_cl->mri[i][j+1][k]>0)&&(im_cl->mri[i][j+1][k]!=im_cl->mri[i][j][k])) nb_voisin++;
					if ((im_cl->mri[i][j][k+1]>0)&&(im_cl->mri[i][j][k+1]!=im_cl->mri[i][j][k])) nb_voisin++;
	
					if (nb_voisin>0)
						{
						
						for (l=0;l<nb_classe;l++)
							voisinage[l]=0;
							
						if (im_cl->mri[i-1][j][k]>0) voisinage[im_cl->mri[i-1][j][k]-1]++;
						if (im_cl->mri[i][j-1][k]>0) voisinage[im_cl->mri[i][j-1][k]-1]++;
						if (im_cl->mri[i][j][k-1]>0) voisinage[im_cl->mri[i][j][k-1]-1]++;
						if (im_cl->mri[i+1][j][k]>0) voisinage[im_cl->mri[i+1][j][k]-1]++;
						if (im_cl->mri[i][j+1][k]>0) voisinage[im_cl->mri[i][j+1][k]-1]++;
						if (im_cl->mri[i][j][k+1]>0) voisinage[im_cl->mri[i][j][k+1]-1]++;
	
						for (l=0;l<nb_classe;l++)
							{
							potts=psi*voisinage[l];
							vraisemblance=-0.5*(im->mri[i][j][k]-moyenne[l])*(im->mri[i][j][k]-moyenne[l])/(std[l]*std[l])-log(std[l]);
							tmp=exp(potts+vraisemblance);
							if (isnan(tmp))
								proba[l]=0;
							else
							proba[l]=tmp;
							}
						
						/*proba_max=0;lopt=-1;
						for (l=0;l<nb_classe;l++)
							if (proba[l]>proba_max)
								{proba_max=proba[l]; lopt=l;}
						
						if (im_cl->mri[i][j][k]!=(lopt+1))
							{nb++;im_cl->mri[i][j][k]=lopt+1;}*/
							
							proba_tot=0;
						for (l=0;l<nb_classe;l++)
							proba_tot+=proba[l];
							
						for (l=0;l<nb_classe;l++)
								proba[l]=proba[l]/proba_tot;
								
						for (l=1;l<nb_classe;l++)
								proba[l]=proba[l]+proba[l-1];
								
						tirage=drand48();
						
						l=0;
						while((tirage>proba[l])&&(l<(nb_classe-1)))
							l++;
						
						im_cl->mri[i][j][k]=l+1;
						label[l]->mri[i][j][k]++;	
							
						}
					
					}
nb++;
segment_mean_std_classe(im,im_cl,moyenne,std,nb_tot,nb_classe);					
printf("Psi %f    Nombre d'iteration %d \n",psi,nb);
}while (nb<40);


/* choix de la segmentation finale */
for (i=1;i<wdth1;i++)
			for (j=1;j<hght1;j++)
			for (k=1;k<dpth1;k++)
				if (im_cl->mask->mri[i][j][k]!=0)
					{
					nb=0;
					for (l=0;l<nb_classe;l++)
						if (label[l]->mri[i][j][k]>nb)
							{nb=label[l]->mri[i][j][k];lopt=l;}
					
					im_cl->mri[i][j][k]=lopt+1;
					}


free_ptr_grphic3d_array(label,nb_classe);
free(voisinage);
free(proba);
free(nb_tot);
}


/*********************************************************************/
/*        remplissage_carte_de_probabilite						  		         */
/*********************************************************************/
void remplissage_carte_de_probabilite(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,grphic3d *imres,grphic *histo, double *moyenne1,double *moyenne2,double *std1,double *std2,int *nb1,int *nb2,int nb_classe)
{
int i,j,k,wdth,hght,dpth;
long min_x, max_x,min_y,max_y;
int plot_x,plot_y;

wdth=imres->width;
hght=imres->height;
dpth=imres->depth;


min_x=0;
max_x=im1->max_pixel;
min_y=0;
max_y=im2->max_pixel;

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
for (k=0; k<dpth; k++)
	{
		plot_x=floor(1.0*(im1->mri[i][j][k]-min_x)/(max_x-min_x)*(histo->width-1.0));
		plot_y=floor(1.0*(im2->mri[i][j][k]-min_y)/(max_y-min_y)*(histo->height-1.0));

	/*if ((im1->mri[i][j][k]>0)&&(histo->mri[plot_x][plot_y]==0))
	{
	
	ptot=0;
	for (l=0;l<nb_classe; l++)
		{
		tmp=nb1[l]/std1[l]*exp(-0.5*(im1->mri[i][j][k]-moyenne1[l])*(im1->mri[i][j][k]-moyenne1[l])/(std1[l]*std1[l]));
		ptot+=tmp;
		if (l==(im1_cl->mri[i][j][k]-1))
			p1=tmp;
		}
		p1=p1/ptot;
	
		
	ptot=0;
	for (l=0;l<nb_classe; l++)
		{
		tmp=nb2[l]/std2[l]*exp(-0.5*(im2->mri[i][j][k]-moyenne2[l])*(im2->mri[i][j][k]-moyenne2[l])/(std2[l]*std2[l]));
		ptot+=tmp;
		if (l==(im2_cl->mri[i][j][k]-1))
			p2=tmp;
		}
		p2=p2/ptot;
	
	imres->mask->mri[i][j][k]=(int)(10000.0*p1*p2);	
	
		
	imres->mask->mri[i][j][k]=500;
	}
	else
	imres->mask->mri[i][j][k]=1000;		*/
	
	if (im1_cl->mri[i][j][k]<=0)
		imres->mask->mri[i][j][k]=0;
	else
		imres->mask->mri[i][j][k]=1;
	}

}


/*********************************************************************/
/*        							 merge_classe									  		         */
/*********************************************************************/
int merge_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,double *moyenne1,double *moyenne2,double *std1,double *std2,int *nb1,int *nb2,int nb_classe)
{
int i,j,k,l,label,new_label,nb_label;
int **independance,nb_sigma=1;
double t1,t2,seuil=1;


independance=alloc_imatrix(nb_classe,nb_classe);

for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	{
	
	t1=(MINI(moyenne1[i]+nb_sigma*std1[i],moyenne1[j]+nb_sigma*std1[j])-MAXI(moyenne1[i]-nb_sigma*std1[i],moyenne1[j]-nb_sigma*std1[j]))
		/(MINI(nb_sigma*std1[i],nb_sigma*std1[j])); 
	t2=(MINI(moyenne2[i]+nb_sigma*std2[i],moyenne2[j]+nb_sigma*std2[j])-MAXI(moyenne2[i]-nb_sigma*std2[i],moyenne2[j]-nb_sigma*std2[j]))
		/(MINI(nb_sigma*std2[i],nb_sigma*std2[j])); 
	
	
	if ((t1>seuil)&&(t2>seuil))
		independance[i][j]=1;
	else
		independance[i][j]=0;
	}

for (i=0;i<im1_cl->width;i++)
for (j=0;j<im1_cl->height;j++)
for (k=0;k<im1_cl->depth;k++)
	if (im1_cl->mri[i][j][k]>0)
	{
	label=im1_cl->mri[i][j][k]-1;
	new_label=label;
	for (l=label; l<nb_classe;l++)
		if (independance[label][l]==1)
		new_label=MAXI(l,new_label);
		
	im1_cl->mri[i][j][k]=new_label-1;
	}

for (i=0;i<im2_cl->width;i++)
for (j=0;j<im2_cl->height;j++)
for (k=0;k<im2_cl->depth;k++)
	if (im2_cl->mri[i][j][k]>0)
	{
	label=im2_cl->mri[i][j][k]-1;
	new_label=label;
	for (l=label; l<nb_classe;l++)
		if (independance[label][l]==1)
		new_label=MAXI(l,new_label);
		
	im2_cl->mri[i][j][k]=new_label-1;
	}


segment_mean_std_classe(im1,im1_cl,moyenne1,std1,nb1,nb_classe);
segment_mean_std_classe(im2,im2_cl,moyenne2,std2,nb2,nb_classe);

nb_label=0;
for (l=0; l<nb_classe;l++)
	if (nb1[l]>0)
	nb_label++; 

free_imatrix(independance,nb_classe,nb_classe);
printf("Nombre de classe restante apres merge : %d sur %d \n",nb_label,nb_classe);
return(nb_label);
}


/*********************************************************************/
/*        					 apply_norm_mean_std_classe						           */
/*********************************************************************/
void apply_norm_mean_std_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int
nb_classe)
{
int i,j,k,l,wdth,hght,dpth;
double tmp;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im1->mri[i][j][k]>0)
		{
		l=im1_cl->mri[i][j][k]-1;
		
		if (covariance[l]>=0)
			tmp=1.0*(im1->mri[i][j][k]-moyenne1[l])/std1[l]*std2[l]+moyenne2[l];
			else
			tmp=1.0*(moyenne1[l]-im1->mri[i][j][k])/std1[l]*std2[l]+moyenne2[l];
		
			
		imres->mri[i][j][k]=(int)tmp;
		}
		else
		imres->mri[i][j][k]=0;

}

/*********************************************************************/
/*        					 apply_norm_mean_std_classe_optcov						           */
/*********************************************************************/
void apply_norm_mean_std_classe_optcov(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2, grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int
nb_classe)
{
int i,j,k,l,t,wdth,hght,dpth;
double tmp,Emin,Etmp;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

/* on cherche �estimer le signe de la covariance */

for (l=0; l<nb_classe;l++)
	covariance[l]=1;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im1->mri[i][j][k]>0)
		{
		l=im1_cl->mri[i][j][k]-1;
		
		if (covariance[l]>=0)
			tmp=1.0*(im1->mri[i][j][k]-moyenne1[l])/std1[l]*std2[l]+moyenne2[l];
			else
			tmp=1.0*(moyenne1[l]-im1->mri[i][j][k])/std1[l]*std2[l]+moyenne2[l];
		
			
		imres->mri[i][j][k]=(int)tmp;
		}
		else
		imres->mri[i][j][k]=0;

Emin=imx_energie_inter_images_3d_p(im2,imres);


for (t=0; t<nb_classe;t++)
{
printf("label %d sur %d \r",t,nb_classe);
covariance[t]=-1;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im1->mri[i][j][k]>0)
		{
		l=im1_cl->mri[i][j][k]-1;
		
		if (covariance[l]>=0)
			tmp=1.0*(im1->mri[i][j][k]-moyenne1[l])/std1[l]*std2[l]+moyenne2[l];
			else
			tmp=1.0*(moyenne1[l]-im1->mri[i][j][k])/std1[l]*std2[l]+moyenne2[l];
		
			
		imres->mri[i][j][k]=(int)tmp;
		}
		else
		imres->mri[i][j][k]=0;

Etmp=imx_energie_inter_images_3d_p(im2,imres);

if (Etmp<Emin)
	Emin=Etmp;
else
	covariance[t]=1;

}


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im1->mri[i][j][k]>0)
		{
		l=im1_cl->mri[i][j][k]-1;
		
		if (covariance[l]>=0)
			tmp=1.0*(im1->mri[i][j][k]-moyenne1[l])/std1[l]*std2[l]+moyenne2[l];
			else
			tmp=1.0*(moyenne1[l]-im1->mri[i][j][k])/std1[l]*std2[l]+moyenne2[l];
		
			
		imres->mri[i][j][k]=(int)tmp;
		}
		else
		imres->mri[i][j][k]=0;

}


/*********************************************************************/
/*        					 apply_norm_regression_lineaire						           */
/*********************************************************************/
void apply_norm_regression_lineaire(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *imres,grphic * histo, int nb_classe)
{
int i,j,k,l,wdth,hght,dpth,*nb,u,v;
double *xk, *yk, *xk2, *xkyk,*a,*b,tmp;
double max_reca, max_ref,x,y;

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (im1->mri[i][j][k]>max_reca) max_reca=im1->mri[i][j][k];
	if (im2->mri[i][j][k]>max_ref) max_ref=im2->mri[i][j][k];
	}
	

xk=(double *)malloc(nb_classe * sizeof(double));
yk=(double *)malloc(nb_classe * sizeof(double));
xk2=(double *)malloc(nb_classe * sizeof(double));
xkyk=(double *)malloc(nb_classe * sizeof(double));
a=(double *)malloc(nb_classe * sizeof(double));
b=(double *)malloc(nb_classe * sizeof(double));
nb=(int *)malloc(nb_classe * sizeof(int));

for (l=0;l<nb_classe;l++)
	a[l]=b[l]=xk[l]=yk[l]=xk2[l]=xkyk[l]=0.0;

for (l=0;l<nb_classe;l++)
		nb[l]=0;
		
		
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im1->mri[i][j][k]>0)&&(im1_cl->mri[i][j][k]>0))
			{
			x=1.0*im1->mri[i][j][k]/max_reca*(histo->width-1.0);
			y=1.0*im2->mri[i][j][k]/max_ref*(histo->height-1.0);
			u=floor(x);
			v=floor(y);
			
			if ((u>=0)&&(v>=0))
				if (histo->mri[u][v]>0)
				{
				l=(int)(im1_cl->mri[i][j][k]-1);
				nb[l]++;
				xk[l]=xk[l]+(double)im1->mri[i][j][k];
				yk[l]=yk[l]+(double)im2->mri[i][j][k];
				xk2[l]=xk2[l]+(double)(im1->mri[i][j][k]*im1->mri[i][j][k]);
				xkyk[l]=xkyk[l]+(double)(im1->mri[i][j][k]*im2->mri[i][j][k]);
				}
			}	

/* calcul des coeffients de la droite de regression lin�ire */
for (l=0;l<nb_classe;l++)
	if (nb[l]>0)
	{
	a[l]=((double)nb[l]*xkyk[l]-xk[l]*yk[l])/((double)nb[l]*xk2[l]-xk[l]*xk[l]);
	b[l]=(xk2[l]*yk[l]-xk[l]*xkyk[l])/((double)nb[l]*xk2[l]-xk[l]*xk[l]);
	
	printf("regression lineaire ax+b :  a = %f  \t  b = %f\n",a[l],b[l]);
	}
	else
	{a[l]=b[l]=0;}


/* on applique la transformation lin�ire sur les intensit� */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im1->mri[i][j][k]>0)&&(im1_cl->mri[i][j][k]>0))
			{
			l=im1_cl->mri[i][j][k]-1;
			tmp=im1->mri[i][j][k]*a[l]+b[l];
			imres->mri[i][j][k]=(int)(tmp);
			}
			else
			imres->mri[i][j][k]=0;



free(xk);
free(yk);
free(xk2);
free(xkyk);
free(nb);
free(a);free(b);
}


/*********************************************************************/
/*        					 apply_norm_quantile_classe						           */
/*********************************************************************/
void apply_norm_quantile_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,grphic3d *imres,double *covariance,int nb_classe,int nb_quantile)
{
int i,j,k,l,q,wdth,hght,dpth,nb;
double tmp,alpha;
double** quantile1, **quantile2;
double *data;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

quantile1 = alloc_dmatrix(nb_classe,nb_quantile);
quantile2 = alloc_dmatrix(nb_classe,nb_quantile);
data=(double *) malloc(wdth*hght*dpth*sizeof(double));

/* remplissage de quantile */
for (l=1; l<=nb_classe; l++)
{
nb=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im1->mri[i][j][k]>0)&&(im1_cl->mri[i][j][k]==l))
		{
		data[nb]=im1->mri[i][j][k];
		nb++;
		}

qsort(data,nb,sizeof(double),double_compare_function2);

nb --;

for (i=0;i<nb_quantile;i++)
	{
	quantile1[l-1][i]=data[(int)(i*nb/(nb_quantile-1))];
	}	
}



for (l=1; l<=nb_classe; l++)
{
nb=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im2->mri[i][j][k]>0)&&(im2_cl->mri[i][j][k]==l))
		{
		data[nb]=im2->mri[i][j][k];
		nb++;
		}

qsort(data,nb,sizeof(double),double_compare_function2);

nb --;

for (i=0;i<nb_quantile;i++)
	{
	quantile2[l-1][i]=data[(int)(i*nb/(nb_quantile-1))];
	}	
}

/* Inversion dans quantile2 des lignes pour lesquelles la covariance est n�ative */
for (l=0; l<nb_classe; l++)
if (covariance[l]<0)
	{
	for (i=0;i<floor(nb_quantile/2);i++)
		{
		tmp=quantile2[l][i];
		quantile2[l][i]=quantile2[l][nb_quantile-1-i];
		quantile2[l][nb_quantile-1-i]=tmp;
		}
	}


/* Application de la normalisation */

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im1->mri[i][j][k]>0)&&(im1_cl->mri[i][j][k]>0))
		{
		l=im1_cl->mri[i][j][k]-1;
		
		q=0;
		while (quantile1[l][q]<im1->mri[i][j][k])
			q++;
		
		if (q>0)
			{alpha=(im1->mri[i][j][k]-quantile1[l][q-1])/(quantile1[l][q]-quantile1[l][q-1]);
			tmp = quantile2[l][q-1]+alpha*(quantile2[l][q]-quantile2[l][q-1]);}
		else
			tmp=quantile2[l][0];
			
		imres->mri[i][j][k]=(int)tmp;
		}
		else
		imres->mri[i][j][k]=0;


free_dmatrix(quantile1,nb_classe,nb_quantile);
free_dmatrix(quantile2,nb_classe,nb_quantile);
free(data);
}


/*********************************************************************/
/*        					 apply_norm_mean_std_classe_reg						           */
/*********************************************************************/
void apply_norm_mean_std_classe_reg(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int
nb_classe)
{
int i,j,k,l,wdth,hght,dpth;
double tmp;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im1->mri[i][j][k]>0)&&(im1_cl->mri[i][j][k]>0))
		{
		l=im1_cl->mri[i][j][k]-1;
		
		tmp=1.0*(im1->mri[i][j][k]-moyenne1[l])*covariance[l]/(std1[l]*std1[l])+moyenne2[l];
		
			
		imres->mri[i][j][k]=(int)tmp;
		}
		else
		imres->mri[i][j][k]=0;

}

/*********************************************************************/
/*        		 apply_norm_mean_std_classe_proba_reg				           */
/*********************************************************************/
void apply_norm_mean_std_classe_proba_reg(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int nb_classe)
{
int i,j,k,l,wdth,hght,dpth,*label,nb_cluster,*nb_nonrecale,*nb_total;
double tirage;
grphic3d *im_label;

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;


im_label = cr_grphic3d_modif(wdth,hght,dpth,0.0,1,1);


//imx_labelcroix_3d_p(im1_cl, im_label);
imx_labelcube_3d_p(im1_cl, im_label);
nb_cluster=im_label->max_pixel;

label=(int *) malloc (nb_cluster*sizeof(int));
nb_nonrecale=(int *) malloc (nb_cluster*sizeof(int));
nb_total=(int *) malloc (nb_cluster*sizeof(int));

for (l=0; l<nb_cluster;l++)
	nb_nonrecale[l]=0;


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im_label->mri[i][j][k]>0)
		{
		l=im_label->mri[i][j][k]-1;
		
		if (im1_cl->mask->mri[i][j][k]<0)
			nb_nonrecale[l]++;
		
		nb_total[l]++;
		
		label[l]=	im1_cl->mri[i][j][k];
		}


for (l=0; l<nb_cluster;l++)
	{
	tirage=drand48();
	
	if (tirage<1.0*nb_nonrecale[l]/nb_total[l])
		label[l]=floor(drand48()*nb_classe+1);
	}

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im_label->mri[i][j][k]>0)
		{
		l=im_label->mri[i][j][k]-1;
		
		im1_cl->mri[i][j][k]=label[l];
		}

		
apply_norm_mean_std_classe_reg(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param, nb_classe);

free_grphic3d(im_label);
free(label); free(nb_nonrecale);free(nb_total);
}


/*********************************************************************/
/*        					 apply_norm_mean_std_classe_vp				           */
/*********************************************************************/
void apply_norm_mean_std_classe_vp(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int
nb_classe)
{
int i,j,k,l,wdth,hght,dpth,t,wdth1,hght1,dpth1,v;
int label[7];
double tmp,proba_tot,proba,u;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if (im1->mri[i][j][k]>0)
		{
		
		
		label[0]=im1_cl->mri[i][j][k]-1;
		l=1;
		
		if (i>0)  {v=0; for (t=0;t<l;t++) if ((im1_cl->mri[i-1][j][k]!=(label[t]+1))&&(im1_cl->mri[i-1][j][k]>0)) v++;
									if (v==l) {label[l]=im1_cl->mri[i-1][j][k]-1;l++;}}
	
		if (j>0) {v=0;for (t=0;t<l;t++) if ((im1_cl->mri[i][j-1][k]!=(label[t]+1))&&(im1_cl->mri[i][j-1][k]>0)) v++;
									if (v==l) {label[l]=im1_cl->mri[i][j-1][k]-1;l++;}}
		
		if (k>0) {v=0;for (t=0;t<l;t++) if ((im1_cl->mri[i][j][k-1]!=(label[t]+1))&&(im1_cl->mri[i][j][k-1]>0)) v++;
									if (v==l) {label[l]=im1_cl->mri[i][j][k-1]-1;l++;}}
		
		if (i<wdth1) {v=0;for (t=0;t<l;t++) if ((im1_cl->mri[i+1][j][k]!=(label[t]+1))&&(im1_cl->mri[i+1][j][k]>0)) v++;
									if (v==l) {label[l]=im1_cl->mri[i+1][j][k]-1;l++;}}
		
		if (j<hght1) {v=0;for (t=0;t<l;t++) if ((im1_cl->mri[i][j+1][k]!=(label[t]+1))&&(im1_cl->mri[i][j+1][k]>0)) v++;
									if (v==l) {label[l]=im1_cl->mri[i][j+1][k]-1;l++;}}
		
		if (k<dpth1) {v=0;for (t=0;t<l;t++) if ((im1_cl->mri[i][j][k+1]!=(label[t]+1))&&(im1_cl->mri[i][j][k+1]>0)) v++;
									if (v==l) {label[l]=im1_cl->mri[i][j][k+1]-1;l++;}}
		
		proba_tot=0;
		tmp=0;
		for (t=0; t<l; t++)
			{
		u=((double)im1->mri[i][j][k]-moyenne1[label[t]]);	
		proba=exp(-u*u/(2*std1[label[t]]*std1[label[t]]))/std1[label[t]];
		proba_tot+=proba;
		
		if (covariance[label[t]]>=0)
			tmp+=(1.0*u/std1[label[t]]*std2[label[t]]+moyenne2[label[t]])*proba;
			else
			tmp+=(-1.0*u/std1[label[t]]*std2[label[t]]+moyenne2[label[t]])*proba;
			}
			
		imres->mri[i][j][k]=(int)(tmp/proba_tot);
		}
		else
		imres->mri[i][j][k]=0;

}
/*********************************************************************/
/*        		segment_covariance_classe   								           */
/*********************************************************************/
void segment_covariance_classe(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,double *covariance,double *moyenne1, double *moyenne2, double *std1 , double
*std2, int nb_classe)
{
int i,j,k,l,label,wdth,hght,dpth,nb;
double cov;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

for (l=1; l<=nb_classe; l++)
	{
	label=l-1;
	cov=0;
	nb=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if ((im1_cl->mri[i][j][k]==l)&&(im2_cl->mri[i][j][k]==l))
			{
			cov+=(im1->mri[i][j][k]-moyenne1[label])*(im2->mri[i][j][k]-moyenne2[label]);
			nb++;
			}
	if (nb>0)		
	cov=cov/nb;
	else
	cov=0;
	
	covariance[label]=cov;
	}
}


/*********************************************************************/
/*        		segment_covariance_classe   								           */
/*********************************************************************/
void segment_signe_covariance_classe_robuste(grphic3d *im1,grphic3d *im1_cl,grphic3d *im2,grphic3d *im2_cl,double *covariance,double *moyenne1, double *moyenne2, double *std1 , double
*std2, int nb_classe, grphic *histo)
{
int i,j,k,l,label,wdth,hght,dpth,nb,u,v;
double cov,tmp, max_ref,max_reca,x,y;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (im1->mri[i][j][k]>max_reca) max_reca=im1->mri[i][j][k];
	if (im2->mri[i][j][k]>max_ref) max_ref=im2->mri[i][j][k];
	}
	

for (l=1; l<=nb_classe; l++)
	{
	label=l-1;
	cov=0;
	nb=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		{
			x=1.0*im1->mri[i][j][k]/max_reca*(histo->width-1.0);
			y=1.0*im2->mri[i][j][k]/max_ref*(histo->height-1.0);
			u=floor(x);
			v=floor(y);
		if ((u>=0)&&(v>=0))	
		if ((im1_cl->mri[i][j][k]==l)&&(im2_cl->mri[i][j][k]==l)&&(histo->mri[u][v]>0))
			{
			tmp=(im1->mri[i][j][k]-moyenne1[label])*(im2->mri[i][j][k]-moyenne2[label]);
			
			if (tmp>=0)
				cov++;
			
			if (tmp<=0)
			 cov--; 
			}
		}
	covariance[label]=cov;
	}
}


/*********************************************************************/
/*        					 segment_mean_std_classe   										           */
/*********************************************************************/
void segment_mean_std_classe(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int *nbtot,int nb_classe)
{
int i,j,k,l,wdth,hght,dpth,nb;
double tmp_mean, tmp_var;
wdth=im->width;
hght=im->height;
dpth=im->depth;

for (l=1; l<=nb_classe; l++)
	{
	tmp_mean=0;
	nb=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if (im_cl->mri[i][j][k]==l)
			{
			tmp_mean+=im->mri[i][j][k];
			nb++;
			}
	if (nb>0)		
	tmp_mean=tmp_mean/nb;
	else
	tmp_mean=0;
	
	tmp_var=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if (im_cl->mri[i][j][k]==l)
			tmp_var+=(im->mri[i][j][k]-tmp_mean)*(im->mri[i][j][k]-tmp_mean);
			
	
	if (nb>1)		
	tmp_var=tmp_var/(nb-1);
	else
	tmp_var=0;
	 
	moyenne[l-1]=tmp_mean;
	std[l-1]=sqrt(tmp_var);
	nbtot[l-1]=nb;
	
	}

}


/*********************************************************************/
/*        					 segment_mean_std_classe_robust   										           */
/*********************************************************************/
void segment_mean_std_classe_robust(grphic3d *im,grphic3d *im_cl,double *moyenne,double *std,int nb_classe)
{
int i,j,k,l,wdth,hght,dpth,nb;
double tmp_mean, tmp_var,*data;
wdth=im->width;
hght=im->height;
dpth=im->depth;

data=(double *)malloc(wdth*hght*dpth*sizeof(double));

for (l=1; l<=nb_classe; l++)
	{
	nb=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if (im_cl->mri[i][j][k]==l)
			{
			data[nb]=im->mri[i][j][k];
			nb++;
			}
	if (nb>0)		
		{
		qsort(data,nb,sizeof(double),double_compare_function2);
		tmp_mean=data[(int)(nb/2)];			nb++;

		}
	else
	tmp_mean=0;
	
	nb=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if (im_cl->mri[i][j][k]==l)
			{
			data[nb]=fabs(im->mri[i][j][k]-tmp_mean);
			nb++;
			}
			
	if (nb>0)		
		{
		qsort(data,nb,sizeof(double),double_compare_function2);
		tmp_var=1.4826*data[(int)(nb/2)];
		}
	else
	tmp_var=0;
	
	moyenne[l-1]=tmp_mean;
	std[l-1]=tmp_var;
	
	}

free(data);
}



/*********************************************************************/
/*         raffine_segment_croissance_region         							   */
/*********************************************************************/
void raffine_segment_croissance_region(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int nb_classe )
{
int i,j,k,wdth,hght,dpth,voisinage,seuil,label,nb,nb_sigma=1/*,Ntot*/;
double  min,sigma;
wdth=im->width;
hght=im->height;
dpth=im->depth;

sigma=HUGE_VAL;

for (i=0;i<nb_classe; i++)
if (std[i]>0)
sigma= MINI(std[i],sigma);

sigma=0.5*sigma;
//Ntot=1;
do
{
/*Ntot=0;
for (i=1;i<wdth-1;i++)
for (j=1;j<hght-1;j++)
for (k=1;k<dpth-1;k++)
	if (im_cl->mri[i][j][k]==-1)
		{
		Ntot++;
		if (nb_sigma==20)
			printf("y a un bug\n");
		}*/
seuil=6;

	do
	{
		do
		{
			nb=0;
			
			for (i=1;i<wdth-1;i++)
			for (j=1;j<hght-1;j++)
			for (k=1;k<dpth-1;k++)
				if (im_cl->mri[i][j][k]==-1)
					{
					voisinage=0;
					if (im_cl->mri[i-1][j][k]>0) voisinage++;
					if (im_cl->mri[i][j-1][k]>0) voisinage++;
					if (im_cl->mri[i][j][k-1]>0) voisinage++;
					if (im_cl->mri[i+1][j][k]>0) voisinage++;
					if (im_cl->mri[i][j+1][k]>0) voisinage++;
					if (im_cl->mri[i][j][k+1]>0) voisinage++;
		
					if (voisinage>=seuil)
						{
						min=HUGE_VAL;
						label=-1;
						if ((fabs(im->mri[i-1][j][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i-1][j][k]>0)) {min=fabs(im->mri[i-1][j][k]-im->mri[i][j][k]); label=im_cl->mri[i-1][j][k];}
						if ((fabs(im->mri[i][j-1][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j-1][k]>0)) {min=fabs(im->mri[i][j-1][k]-im->mri[i][j][k]); label=im_cl->mri[i][j-1][k];}
						if ((fabs(im->mri[i][j][k-1]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j][k-1]>0)) {min=fabs(im->mri[i][j][k-1]-im->mri[i][j][k]); label=im_cl->mri[i][j][k-1];}
						if ((fabs(im->mri[i+1][j][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i+1][j][k]>0)) {min=fabs(im->mri[i+1][j][k]-im->mri[i][j][k]); label=im_cl->mri[i+1][j][k];}
						if ((fabs(im->mri[i][j+1][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j+1][k]>0)) {min=fabs(im->mri[i][j+1][k]-im->mri[i][j][k]); label=im_cl->mri[i][j+1][k];}
						if ((fabs(im->mri[i][j][k+1]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j][k+1]>0)) {min=fabs(im->mri[i][j][k+1]-im->mri[i][j][k]); label=im_cl->mri[i][j][k+1];}
			
						/* validation de l'affectation */
						if (fabs(im->mri[i][j][k]-moyenne[label-1])<nb_sigma*std[label-1])
						//if (fabs(im->mri[i][j][k]-moyenne[label-1])<nb_sigma*sigma)
						//if (min<nb_sigma*sigma)
							{
							im_cl->mri[i][j][k]=label;
							nb++;
							}
						}
			
					}
			//printf ("nb_sigma %d   seuil  %d   nb  %d \n",nb_sigma,seuil, nb);
			}while (nb>0);
		seuil--;
	}
	while (seuil>0);
nb_sigma++;
}
while (nb_sigma<20);
}

/*********************************************************************/
/*         raffine_segment_croissance_region         							   */
/*********************************************************************/
void raffine_segment_croissance_region2(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int nb_classe )
{
int i,j,k,wdth,hght,dpth,voisinage,label,nb,nb_sigma=1,max_nbsigma,Ntot;
int wdth1,hght1,dpth1;
double  min,sigma,incr;
int *nbtot;

incr=0.25;
//incr=0.01;

nbtot=(int *) malloc(nb_classe*sizeof(int));

wdth=im->width;
hght=im->height;
dpth=im->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;

sigma=HUGE_VAL;

for (i=0;i<nb_classe; i++)
if (std[i]>0)
sigma= MINI(std[i],sigma);


max_nbsigma=im->max_pixel/sigma;

printf("sigma : %f   max_sigma %d \n",sigma,max_nbsigma);

do
{
	do
		{
			nb=0;
			for (i=0;i<wdth;i++)
			for (j=0;j<hght;j++)
			for (k=0;k<dpth;k++)
				if (im_cl->mri[i][j][k]==-1)
					{
					voisinage=0;
					if (i>0) if (im_cl->mri[i-1][j][k]>0) voisinage++;
					if (j>0) if (im_cl->mri[i][j-1][k]>0) voisinage++;
					if (k>0) if (im_cl->mri[i][j][k-1]>0) voisinage++;
					if (i<wdth1) if (im_cl->mri[i+1][j][k]>0) voisinage++;
					if (j<hght1) if (im_cl->mri[i][j+1][k]>0) voisinage++;
					if (k<dpth1) if (im_cl->mri[i][j][k+1]>0) voisinage++;
		
					if (voisinage>=1)
						{
						min=HUGE_VAL;
						label=-1;
						if (i>0) if ((fabs(im->mri[i-1][j][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i-1][j][k]>0)) {min=fabs(im->mri[i-1][j][k]-im->mri[i][j][k]); label=im_cl->mri[i-1][j][k];}
						if (j>0) if ((fabs(im->mri[i][j-1][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j-1][k]>0)) {min=fabs(im->mri[i][j-1][k]-im->mri[i][j][k]); label=im_cl->mri[i][j-1][k];}
						if (k>0) if ((fabs(im->mri[i][j][k-1]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j][k-1]>0)) {min=fabs(im->mri[i][j][k-1]-im->mri[i][j][k]); label=im_cl->mri[i][j][k-1];}
						if (i<wdth1) if ((fabs(im->mri[i+1][j][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i+1][j][k]>0)) {min=fabs(im->mri[i+1][j][k]-im->mri[i][j][k]); label=im_cl->mri[i+1][j][k];}
						if (j<hght1) if ((fabs(im->mri[i][j+1][k]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j+1][k]>0)) {min=fabs(im->mri[i][j+1][k]-im->mri[i][j][k]); label=im_cl->mri[i][j+1][k];}
						if (k<dpth1) if ((fabs(im->mri[i][j][k+1]-im->mri[i][j][k])<min)&&(im_cl->mri[i][j][k+1]>0)) {min=fabs(im->mri[i][j][k+1]-im->mri[i][j][k]); label=im_cl->mri[i][j][k+1];}
			
						/* validation de l'affectation */
						//if (fabs(im->mri[i][j][k]-moyenne[label-1])<nb_sigma*sigma)
						if (fabs(im->mri[i][j][k]-moyenne[label-1])<incr*nb_sigma*std[label-1])
								{
							im_cl->mri[i][j][k]=label;
							nb++;
							}
						}
			
					}
			}while (nb>0);

Ntot=0;
	for (i=0;i<wdth;i++)
			for (j=0;j<hght;j++)
			for (k=0;k<dpth;k++)
				if (im_cl->mri[i][j][k]==-1)
						Ntot++;
nb_sigma++;
segment_mean_std_classe(im,im_cl,moyenne,std,nbtot,nb_classe);

//printf("nb_sigma : %d   Not %d \n",nb_sigma,Ntot);
}
while ((Ntot>0)&&(nb_sigma<100));
//while (nb_sigma<max_nbsigma);
//while (nb_sigma<16);

free(nbtot);
}


/*********************************************************************/
/*         raffine_segment_croissance_region         							   */
/*********************************************************************/
void raffine_segment_croissance_region3(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int nb_classe )
{
int i,j,k,wdth,hght,dpth,voisinage,label,nb,nb_sigma=1,voisinage_max;
int wdth1,hght1,dpth1;
double  min,tmp;
wdth=im->width;
hght=im->height;
dpth=im->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;
voisinage_max=6;

do
{
nb_sigma=1;
do
{
	do
		{
			nb=0;
			for (i=0;i<wdth;i++)
			for (j=0;j<hght;j++)
			for (k=0;k<dpth;k++)
				if (im_cl->mri[i][j][k]==-1)
					{
					voisinage=0;
					if (i>0) if (im_cl->mri[i-1][j][k]>0) voisinage++;
					if (j>0) if (im_cl->mri[i][j-1][k]>0) voisinage++;
					if (k>0) if (im_cl->mri[i][j][k-1]>0) voisinage++;
					if (i<wdth1) if (im_cl->mri[i+1][j][k]>0) voisinage++;
					if (j<hght1) if (im_cl->mri[i][j+1][k]>0) voisinage++;
					if (k<dpth1) if (im_cl->mri[i][j][k+1]>0) voisinage++;
		
					if (voisinage>=5)
						{
						min=HUGE_VAL;
						label=-1;
						if (i>0) if (im_cl->mri[i-1][j][k]>0) 
									{tmp=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i-1][j][k]-1])/std[im_cl->mri[i-1][j][k]-1];
									if (tmp<min){min=tmp; label=im_cl->mri[i-1][j][k]-1;}}
					
						if (j>0) if (im_cl->mri[i][j-1][k]>0) 
									{tmp=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j-1][k]-1])/std[im_cl->mri[i][j-1][k]-1];
									if (tmp<min){min=tmp; label=im_cl->mri[i][j-1][k]-1;}}
										
						if (k>0) if (im_cl->mri[i][j][k-1]>0) 
									{tmp=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j][k-1]-1])/std[im_cl->mri[i][j][k-1]-1];
									if (tmp<min){min=tmp; label=im_cl->mri[i][j][k-1]-1;}}
					
						if (i<wdth1) if (im_cl->mri[i+1][j][k]>0) 
									{tmp=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i+1][j][k]-1])/std[im_cl->mri[i+1][j][k]-1];
									if (tmp<min){min=tmp; label=im_cl->mri[i+1][j][k]-1;}}
					
						if (j<hght1) if (im_cl->mri[i][j+1][k]>0) 
									{tmp=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j+1][k]-1])/std[im_cl->mri[i][j+1][k]-1];
									if (tmp<min){min=tmp; label=im_cl->mri[i][j+1][k]-1;}}
										
						if (k<dpth1) if (im_cl->mri[i][j][k+1]>0) 
									{tmp=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j][k+1]-1])/std[im_cl->mri[i][j][k+1]-1];
									if (tmp<min){min=tmp; label=im_cl->mri[i][j][k+1]-1;}}

					
						/* validation de l'affectation */
						//if (fabs(im->mri[i][j][k]-moyenne[label-1])<nb_sigma*sigma)
							if (fabs(im->mri[i][j][k]-moyenne[label])<nb_sigma*std[label])
								{
							im_cl->mri[i][j][k]=label+1;
							nb++;
							}
						}
			
					}
					printf("voisinage %d  nb_sigma %d  nb_segemente %d \n",voisinage_max,nb_sigma,nb);
			}while (nb>0);
nb_sigma++;
//segment_mean_std_classe(im,im_cl,moyenne,std,nb_classe);
}
while (nb_sigma<5);
voisinage_max--;

}
while(voisinage_max>=1);
}
/*********************************************************************/
/*         raffine_segment_ML         							   */
/*********************************************************************/
void raffine_segment_ML(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int nb_classe )
{
int i,j,k,l,wdth,hght,dpth,label=0;
double  min,tmp;
wdth=im->width;
hght=im->height;
dpth=im->depth;


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if ((im->mri[i][j][k]>0)&&(im_cl->mri[i][j][k]==-1))
		{
			min=HUGE_VAL;
			for (l=0;l<nb_classe;l++)
				{
				tmp=fabs(im->mri[i][j][k]-moyenne[l])/std[l];
				if (tmp<min)
					{min=tmp; label=l;}
				
				}	
				im_cl->mri[i][j][k]=label+1;	
		}

}

/*********************************************************************/
/*         raffine_segment_MAP         							   */
/*********************************************************************/
void raffine_segment_MAP(grphic3d *im,grphic3d *im_cl, double *moyenne, double *std , int *nb, int nb_classe )
{
int i,j,k,l,wdth,hght,dpth,label=0;
double  max,tmp;
wdth=im->width;
hght=im->height;
dpth=im->depth;


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if ((im->mri[i][j][k]>0)&&(im_cl->mri[i][j][k]==-1))
		{
			max=0;
			for (l=0;l<nb_classe;l++)
				{
				tmp=fabs(im->mri[i][j][k]-moyenne[l])/std[l];
				tmp=exp(-0.5*tmp*tmp)/std[l]*nb[l];
				if (tmp>max)
					{max=tmp; label=l;}
				
				}	
				im_cl->mri[i][j][k]=label+1;	
		}

}



/*********************************************************************/
/*         segment_gaussienne2d_recherche_voisinage               */
/*********************************************************************/
void segment_gaussienne2d_recherche_voisinage(grphic3d *imreca,grphic3d *imref,grphic3d *imreca_cl,grphic3d *imref_cl, 
grphic *histo, int nb_classe, Param_Gaussienne_2d *param,double ** confiance)
{
int i,j,k,l,wdth,hght,dpth,u,v,classe=0,nb,Ntot,voisin_max=10,label;
double x,y,max_reca, max_ref,proba_max,proba_tmp,proba_tot,max_confiance,min,distance_min;
grphic *histo_cl;
int itmp,jtmp,ktmp,ideb,ifin,jdeb,jfin,kdeb,kfin,uopt,vopt,utmp,vtmp;
double xtmp,ytmp,seuil=0.95;
double **proba_ij;

proba_ij=alloc_dmatrix(histo->width, histo->height);


if (BASE3D.x1!=NULL)
voisin_max=(int)((BASE3D.x1[0]-BASE3D.x0[0])/2);

voisin_max=MINI(voisin_max, 4);

printf("taille du voisinage dans la segment_gaussienne2d_recherche_voisinage : %d \n",voisin_max);


min=1;
for (l=0;l<nb_classe;l++)
	{
	if (param[l].l<min)
		{min=param[l].l; label=l;}
	}


label=label+1;
label=-1;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

histo_cl=cr_grphic_modif(histo->width,histo->height,0,1,0);

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

Ntot=nb=0;


/* classification de l'histogramme joint */
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
	{
	Ntot++;
	proba_max=0;
	classe=0;
	for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
		{
		proba_tmp=eval_gaussienne_2d(i, j, &param[l]);
		if (proba_tmp>proba_max){proba_max=proba_tmp; classe=l;}
		}
		histo_cl->mri[i][j]=classe+1;nb++;
	}


/* on remplit proba_ij avec  la proba p(i,j) obtenue par l'estimation du m�ange de gaussiennes */

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
		{
		proba_tot=0;
		for (k=0;k<nb_classe;k++)
			{proba_tmp=	eval_gaussienne_2d(i, j,&param[k]);
			 if(isnan(proba_tmp)) 
			 		proba_tmp=0;
				proba_tot+=proba_tmp;
			}
		proba_ij[i][j]=proba_tot;
		}


/*  Segmentation des deux images */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	printf("voxels segment� %d %d %d \r",i,j,k);
	if ((imref->mri[i][j][k]<=0)&&(imreca->mri[i][j][k]<=0))
		{ imref_cl->mri[i][j][k]=0;imreca_cl->mri[i][j][k]=0;}
	else
		{
			x=1.0*imreca->mri[i][j][k]/max_reca*(histo->width-1.0);
			y=1.0*imref->mri[i][j][k]/max_ref*(histo->height-1.0);
			u=floor(x);
			v=floor(y);
			
			if (x<=0)
					imreca_cl->mri[i][j][k]=0;
					
			if (y<=0)
					imref_cl->mri[i][j][k]=0;
					
			if ((x>0)&&(y>0))
				{
				max_confiance=confiance[u][v];
				vopt=v;
				uopt=u;
				imreca_cl->mri[i][j][k]=histo_cl->mri[u][v];
				imref_cl->mri[i][j][k]=histo_cl->mri[u][v];
				}
				else
				{max_confiance=0;uopt=vopt=0;}
			
		
	

if (max_confiance<seuil)
			{
				if (x>0)
				{
					ideb=MAXI(0,i-voisin_max);
					ifin=MINI(wdth,i+voisin_max)+1;
					jdeb=MAXI(0,j-voisin_max);
					jfin=MINI(hght,j+voisin_max)+1;
					kdeb=MAXI(0,k-voisin_max);
					kfin=MINI(dpth,k+voisin_max)+1;
			
					distance_min=HUGE_VAL;
				
					for (itmp=ideb;itmp<ifin;itmp++)
					for (jtmp=jdeb;jtmp<jfin;jtmp++)
					for (ktmp=kdeb;ktmp<kfin;ktmp++)
						if (imref->mri[itmp][jtmp][ktmp]>0)
						{
						ytmp=1.0*imref->mri[itmp][jtmp][ktmp]/max_ref*(histo->height-1.0);
						vtmp=floor(ytmp);
						
						for (l=0;l<nb_classe;l++)
							{
							proba_tmp=eval_mahalanobis_2d(x, ytmp, &param[l]);
							//proba_tmp=-1.0*eval_gaussienne_2d(x, ytmp, &param[l]);
							
							if (proba_tmp<distance_min)
								{classe=l; distance_min=proba_tmp;}
							
							}
						}
					imreca_cl->mri[i][j][k]=classe;	
					}
				
									
				if (y>0)
				{ 
					ideb=MAXI(0,i-voisin_max);
					ifin=MINI(wdth,i+voisin_max)+1;
					jdeb=MAXI(0,j-voisin_max);
					jfin=MINI(hght,j+voisin_max)+1;
					kdeb=MAXI(0,k-voisin_max);
					kfin=MINI(dpth,k+voisin_max)+1;
				
					distance_min=HUGE_VAL;
				
					for (itmp=ideb;itmp<ifin;itmp++)
					for (jtmp=jdeb;jtmp<jfin;jtmp++)
					for (ktmp=kdeb;ktmp<kfin;ktmp++)
						if (imreca->mri[itmp][jtmp][ktmp]>0)
						{
						xtmp=1.0*imreca->mri[itmp][jtmp][ktmp]/max_reca*(histo->width-1.0);
						utmp=floor(xtmp);
					
						for (l=0;l<nb_classe;l++)
							{
							proba_tmp=eval_mahalanobis_2d(xtmp, y, &param[l]);
							//proba_tmp=-1.0*eval_gaussienne_2d(xtmp, y, &param[l]);
							
							if (proba_tmp<distance_min)
								{classe=l; distance_min=proba_tmp;}
							
							}
						}
							imref_cl->mri[i][j][k]=classe;
					}
				}
	/*	if (max_confiance<seuil)
			{
				if (x>0)
				{
				voisin=1;
				while ((max_confiance<seuil)&&(voisin<voisin_max))
					{
					incr=2*voisin;
					ideb=MAXI(0,i-voisin);
					ifin=MINI(wdth,i+voisin);
					jdeb=MAXI(0,j-voisin);
					jfin=MINI(hght,j+voisin);
					kdeb=MAXI(0,k-voisin);
					kfin=MINI(dpth,k+voisin);
			
					for (itmp=ideb;itmp<ifin;itmp+=incr)
					for (jtmp=jdeb;jtmp<jfin;jtmp+=incr)
					for (ktmp=kdeb;ktmp<kfin;ktmp+=incr)
						if (imref->mri[itmp][jtmp][ktmp]>0)
						{
						ytmp=1.0*imref->mri[itmp][jtmp][ktmp]/max_ref*(histo->height-1.0);
						vtmp=floor(ytmp);
						
						if (confiance[u][vtmp]>max_confiance)
							{vopt=vtmp; max_confiance=confiance[u][vtmp];}
						}
			
					voisin++;
					}
				
				if (max_confiance>seuil)	
				imreca_cl->mri[i][j][k]=histo_cl->mri[u][vopt];
				else
					{
					ideb=MAXI(0,i-voisin_max);
					ifin=MINI(wdth,i+voisin_max);
					jdeb=MAXI(0,j-voisin_max);
					jfin=MINI(hght,j+voisin_max);
					kdeb=MAXI(0,k-voisin_max);
					kfin=MINI(dpth,k+voisin_max);
					
					for (itmp=ideb;itmp<ifin;itmp++)
					for (jtmp=jdeb;jtmp<jfin;jtmp++)
					for (ktmp=kdeb;ktmp<kfin;ktmp++)
						if (imref->mri[itmp][jtmp][ktmp]>0)
						{
						ytmp=1.0*imref->mri[itmp][jtmp][ktmp]/max_ref*(histo->height-1.0);
						vtmp=floor(ytmp);
						
						if (proba_ij[u][vtmp]>proba_max)
							{vopt=vtmp; proba_max=proba_ij[u][vtmp];}
						}
					
						imreca_cl->mri[i][j][k]=histo_cl->mri[u][vopt];
			
					}
					
				
				}
		
				if (y>0)
				{ 
				voisin=1;
	 			max_confiance=confiance_tmp;
				proba_max=proba_tmp;
				while ((max_confiance<seuil)&&(voisin<voisin_max))
					{
					incr=2*voisin;
					ideb=MAXI(0,i-voisin);
					ifin=MINI(wdth,i+voisin);
					jdeb=MAXI(0,j-voisin);
					jfin=MINI(hght,j+voisin);
					kdeb=MAXI(0,k-voisin);
					kfin=MINI(dpth,k+voisin);
				
					for (itmp=ideb;itmp<ifin;itmp+=incr)
					for (jtmp=jdeb;jtmp<jfin;jtmp+=incr)
					for (ktmp=kdeb;ktmp<kfin;ktmp+=incr)
						if (imreca->mri[itmp][jtmp][ktmp]>0)
						{
						xtmp=1.0*imreca->mri[itmp][jtmp][ktmp]/max_reca*(histo->width-1.0);
						utmp=floor(xtmp);
					
						if (confiance[utmp][v]>max_confiance)
							{uopt=utmp; max_confiance=confiance[utmp][v];}
						}
					voisin++;	
					}
				
				if (max_confiance>seuil)
				imref_cl->mri[i][j][k]=histo_cl->mri[uopt][v];
				else
					{
					ideb=MAXI(0,i-voisin_max);
					ifin=MINI(wdth,i+voisin_max);
					jdeb=MAXI(0,j-voisin_max);
					jfin=MINI(hght,j+voisin_max);
					kdeb=MAXI(0,k-voisin_max);
					kfin=MINI(dpth,k+voisin_max);
					
					for (itmp=ideb;itmp<ifin;itmp++)
					for (jtmp=jdeb;jtmp<jfin;jtmp++)
					for (ktmp=kdeb;ktmp<kfin;ktmp++)
						if (imreca->mri[itmp][jtmp][ktmp]>0)
						{
						xtmp=1.0*imreca->mri[itmp][jtmp][ktmp]/max_reca*(histo->height-1.0);
						utmp=floor(xtmp);
						
						if (proba_ij[utmp][v]>proba_max)
							{uopt=utmp; proba_max=proba_ij[utmp][v];}
						}
					
						imref_cl->mri[i][j][k]=histo_cl->mri[uopt][v];
			
					}
				
				}
			}*/
			
		
		
			
			
		}
	}
free_grphic(histo_cl);
free_dmatrix(proba_ij, histo->width, histo->height);
}











/*********************************************************************/
/*         init_segment_gaussienne2d               */
/*********************************************************************/
void init_segment_gaussienne2d(grphic3d *imreca,grphic3d *imref,grphic3d *imreca_cl,grphic3d *imref_cl, grphic *histo, int nb_classe, Param_Gaussienne_2d *param)
{
int i,j,k,l,wdth,hght,dpth,u,v,classe,nb,Ntot;
double x,y,max_reca, max_ref,proba_max,proba_tmp,PROBA_SEUIL=HUGE_VAL;
grphic *histo_cl;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

histo_cl=cr_grphic_modif(histo->width,histo->height,0,1,0);

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

Ntot=nb=0;
/* classification de l'histogramme joint */
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
if (histo->mri[i][j]>0)
	{
	Ntot++;
	proba_max=0;
	classe=0;
	for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
		{
		proba_tmp=eval_gaussienne_2d(i, j, &param[l]);
		if (proba_tmp>proba_max){proba_max=proba_tmp; classe=l;}
		}
	proba_max=eval_mahalanobis_2d(i, j, &param[classe]);
	if (proba_max<PROBA_SEUIL)
		{histo_cl->mri[i][j]=classe+1;nb++;}
	else
	histo_cl->mri[i][j]=-1;
	
	}
else
	histo_cl->mri[i][j]=-1;


printf("Nombre de site segment� : %d sur un total de %d \n",nb,Ntot);

/*  Segmentation des deux images */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if ((imref->mri[i][j][k]<=0)&&(imreca->mri[i][j][k]<=0))
				{ imref_cl->mri[i][j][k]=0;imreca_cl->mri[i][j][k]=0;}
			
	if ((imreca->mri[i][j][k]>0)&&(imref->mri[i][j][k]<=0))
				{ imref_cl->mri[i][j][k]=0;imreca_cl->mri[i][j][k]=-1;}
		
	if ((imreca->mri[i][j][k]<=0)&&(imref->mri[i][j][k]>0))
				{ imref_cl->mri[i][j][k]=-1;imreca_cl->mri[i][j][k]=0;}
	
	
	if ((imreca->mri[i][j][k]>0)&&(imref->mri[i][j][k]>0))
		{
			x=1.0*imreca->mri[i][j][k]/max_reca*(histo->width-1.0);
			y=1.0*imref->mri[i][j][k]/max_ref*(histo->height-1.0);
			u=floor(x);
			v=floor(y);
			imref_cl->mri[i][j][k]=histo_cl->mri[u][v];
			imreca_cl->mri[i][j][k]=histo_cl->mri[u][v];
		}
	}

free_grphic(histo_cl);

}



/*********************************************************************/
/*         apply_normalisation_gaussienne_histo_joint2               */
/*********************************************************************/
void apply_normalisation_gaussienne_histo_joint2(grphic3d *imreca,grphic3d *imref,grphic3d *imres,grphic *histo, int nb_classe, Param_Gaussienne_2d *param)
{
int i,j,k,l,xi,yi,classe/*,l1,l2*/;
double max_reca,max_ref,rmax_reca,rmax_ref;
int wdth,hght,dpth;
double x,y,proba_tot,proba_tmp,tmp,proba_max;
double **ft,**proba_ft;
int **is_classe,**nb_per_classe,max_classe;
grphic *histo_cl;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

ft = alloc_dmatrix(histo->width,nb_classe);
proba_ft = alloc_dmatrix(histo->width,nb_classe);
is_classe = alloc_imatrix(histo->width,nb_classe);
nb_per_classe = alloc_imatrix(histo->width,nb_classe);

	
histo_cl=cr_grphic_modif(histo->width,histo->height,0,1,0);


/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

rmax_reca=max_reca*imreca->rcoeff;
rmax_ref=max_ref*imref->rcoeff;

imx_copie_param_3d_p(imref,imres);



for (xi=0;xi<histo->width; xi++)
	{	
		for (l=0;l<nb_classe;l++)
		if (param[l].l>0)
			{
			proba_tot=0;
			tmp=0;
			for (y=0;y<histo->height;y++)
				{
				proba_tmp=eval_gaussienne_2d(xi, y, &param[l]);
				proba_tot=proba_tot+proba_tmp;
				tmp=tmp+proba_tmp*y;
				}
			if (proba_tot>0)
			ft[xi][l]=1.0*(((tmp*rmax_ref/imres->rcoeff)/(histo->height-1))/proba_tot);
			else
			ft[xi][l]=0;
			
			proba_ft[xi][l]=proba_tot;
			}
	}

/* classification de l'histo */
for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
{
proba_max=0;
classe=0;
for (l=0;l<nb_classe;l++)
if (param[l].l>0)
	{
	proba_tmp=eval_gaussienne_2d(i, j, &param[l]);
	if (proba_tmp>proba_max){proba_max=proba_tmp; classe=l;}
	}
histo_cl->mri[i][j]=classe;
}


/* initialisation de is_classe*/
for (i=0;i<histo->width;i++)
for (j=0;j<nb_classe;j++)
is_classe[i][j]=0;

/* initialisation de nb_per_classe*/
for (i=0;i<histo->width;i++)
for (j=0;j<nb_classe;j++)
nb_per_classe[i][j]=0;

for (i=0;i<histo->width;i++)
for (j=0;j<histo->height;j++)
nb_per_classe[i][histo_cl->mri[i][j]]=nb_per_classe[i][histo_cl->mri[i][j]]+histo->mri[i][j];	



for (i=0;i<histo->width;i++)
{	max_classe=0;
	tmp=0;
	for (j=0;j<nb_classe;j++)
		if (proba_ft[i][j]>tmp) {tmp=proba_ft[i][j]; max_classe=j; }
		
	is_classe[i][max_classe]=1;
}

for (i=0;i<histo->width;i++)
{	max_classe=0;
	tmp=0;
	for (j=0;j<nb_classe;j++)
		if (is_classe[i][j]<1)
			if (proba_ft[i][j]>tmp) 
				{tmp=proba_ft[i][j]; max_classe=j; }
		
	is_classe[i][max_classe]=2;
	}
/*for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if (imreca->mri[i][j][k]>0)
	{
	x=imreca->mri[i][j][k]/max_reca*(histo->width-1);
	y=imref->mri[i][j][k]/max_ref*(histo->height-1);
	xi=floor(x);
	tmp=0;
	proba_tot=0;
	proba_max=0;
	for (l=0;l<nb_classe;l++)
	if (param[l].l>0)
		{
		if (xi<histo->width-1)
		tmp=tmp+((1.0*(x-xi)*ft[xi][l]+1.0*(xi+1-x)*ft[xi+1][l]))*proba_ft[xi][l];
		else
		tmp=tmp+(1.0*(x-xi)*ft[xi][l])*proba_ft[xi][l];
		
		proba_tot=proba_tot+proba_ft[xi][l];
		}
		
	imres->mri[i][j][k]=(int)(1.0*tmp/proba_tot);
	}
else imres->mri[i][j][k]=0;
*/
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if (imreca->mri[i][j][k]>0)
	{
	x=imreca->mri[i][j][k]/max_reca*(histo->width-1);
	y=imref->mri[i][j][k]/max_ref*(histo->height-1);
	xi=floor(x);
	yi=floor(y);
	 
/*	proba_tmp=0;
	for (l=0;l<nb_classe;l++)
		{tmp=eval_gaussienne_2d(x, y, &param[l]);
		if (tmp>proba_tmp) {proba_tmp=tmp; max_classe=l;}
		}
	if (xi<histo->width-1)	
	imres->mri[i][j][k]=(int)((1.0*(x-xi)*ft[xi][max_classe]+1.0*(xi+1-x)*ft[xi+1][max_classe]));
	else
	imres->mri[i][j][k]=(int)(1.0*(x-xi)*ft[xi][max_classe]);*/

/* Test de l'appartenance a au moins une classe */
/*tmp=0;
for (l=0;l<nb_classe;l++)
		if (param[l].l>0)
			{
				proba_tmp=distance_mahalanobis_2d(x, y, &param[l]);
				if (proba_tmp<2)
					tmp=1;
			}*/
			
	
	if  (imref->mri[i][j][k]>0)//((histo->mri[xi][yi]>0)||(tmp>0)) //(imref->mri[i][j][k]>0)
		{
			proba_tot=0;
			tmp=0;
			for (l=0;l<nb_classe;l++)
			if (param[l].l>0)
			if (is_classe[xi][l]>0)
			{
				proba_tmp=eval_gaussienne_2d(x, y, &param[l]);
				proba_tot=proba_tot+proba_tmp;
		
				if (xi<histo->width-1)	
					tmp=tmp+proba_tmp*(1.0*(x-xi)*ft[xi][l]+1.0*(xi+1-x)*ft[xi+1][l]);
				else
					tmp=tmp+proba_tmp*(1.0*(x-xi)*ft[xi][l]);
				}
			tmp=1.0*tmp/proba_tot;
			imres->mri[i][j][k]=(int)tmp;
			}
		else
			/*{
			if (imref->mri[i][j][k]>0)
			{
			tmp=0;
			proba_tot=0;
			for (l=0;l<nb_classe;l++)
			if (param[l].l>0)
			if (is_classe[xi][l]>0)
				{
				proba_tmp=proba_ft[xi][l];
				proba_tot=proba_tot+proba_tmp;
		
				if (xi<histo->width-1)	
					tmp=tmp+proba_tmp*(1.0*(x-xi)*ft[xi][l]+1.0*(xi+1-x)*ft[xi+1][l]);
				else
					tmp=tmp+proba_tmp*(1.0*(x-xi)*ft[xi][l]);
				}
				tmp=1.0*tmp/proba_tot;
				tmp=0;
				imres->mri[i][j][k]=(int)tmp;
				
			}
			else*/	
			{
			tmp=0;
			proba_tot=0;
			for (l=0;l<nb_classe;l++)
			if (param[l].l>0)
			if (is_classe[xi][l]>0)
				{
				proba_tmp=proba_ft[xi][l];
				proba_tot=proba_tot+proba_tmp;
		
				if (xi<histo->width-1)	
					tmp=tmp+proba_tmp*(1.0*(x-xi)*ft[xi][l]+1.0*(xi+1-x)*ft[xi+1][l]);
				else
					tmp=tmp+proba_tmp*(1.0*(x-xi)*ft[xi][l]);
				}
				tmp=1.0*tmp/proba_tot;
				imres->mri[i][j][k]=(int)tmp;
			}
			//}			
	}
else imres->mri[i][j][k]=0;


/* segmentation des deux images */
/*for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
if ((imreca->mri[i][j][k]>0)&&(imref->mri[i][j][k]>0))
	{
	x=imreca->mri[i][j][k]/max_reca*(histo->width-1);
	y=imref->mri[i][j][k]/max_ref*(histo->height-1);
	xi=floor(x);
	yi=floor(y);
	if (histo->mri[xi][yi]>0)
	imref->mri[i][j][k]=imreca->mri[i][j][k]=histo_cl->mri[xi][yi]+1;
	else
	imref->mri[i][j][k]=imreca->mri[i][j][k]=0;
	
	}
else
	imref->mri[i][j][k]=imreca->mri[i][j][k]=-1;

imref->rcoeff=1;
imreca->rcoeff=1;

imx_iniparaimg(imref->pos);
show_picture(imref->pos);	
imx_iniparaimg(imreca->pos);
show_picture(imreca->pos);	
*/


free_grphic(histo_cl);

//free(var_gaussienne);
free_dmatrix(ft,histo->width,nb_classe);
free_dmatrix(proba_ft,histo->width,nb_classe);
free_imatrix(is_classe,histo->width,nb_classe);
free_imatrix(nb_per_classe,histo->width,nb_classe);
}





/*******************************************************************************
********************************************************************************
********************************************************************************
***********   Quelques variantes dans l'estimation des gaussiennes  ************ 
********************************************************************************
********************************************************************************
********************************************************************************/

/*******************************************************************************
**      					  estime_gaussienne2d_isosigma_EM                                                
**                                                                                                   
*******************************************************************************/
void estime_gaussienne2d_isosigma_EM(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe)
{
double ***proba_apriori;
int wdth,hght,i,j,k,Ntot;
double tmp;
wdth=histo->width;
hght=histo->height;

/* allocation memoire */
proba_apriori=alloc_dmatrix_3d(wdth,hght,nb_classe);

/* remplissage de proba_apriori */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		for (k=0;k<nb_classe;k++)
			{proba_apriori[i][j][k]=	eval_gaussienne_2d(i, j,&param[k]);
			 if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
			}

Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		Ntot=Ntot+histo->mri[i][j];

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (histo->mri[i][j]>0)
	{tmp=0;
	for (k=0;k<nb_classe;k++)
	tmp=tmp+proba_apriori[i][j][k];
	
	for (k=0;k<nb_classe;k++)
	proba_apriori[i][j][k]=1.0*proba_apriori[i][j][k]/tmp;
	
	}



/* mise a jour du facteur de normalisation */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].l=0.0;
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].l=new_param[k].l+1.0*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].l=1.0*new_param[k].l/Ntot;
		}

/* mise a jour des moyennes */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].mx=0.0;
		new_param[k].my=0.0;
		
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].mx=new_param[k].mx+1.0*i*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[k].my=new_param[k].my+1.0*j*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].mx=1.0*new_param[k].mx/(Ntot*new_param[k].l);
		new_param[k].my=1.0*new_param[k].my/(Ntot*new_param[k].l);		
		}

/* mise a jour des variances */
		new_param[0].sx=0.0;
		new_param[0].sy=0.0;
		new_param[0].sxy=0.0;
		
	for (k=0;k<nb_classe;k++)
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[0].sx=new_param[0].sx+1.0*(i-new_param[k].mx)*(i-new_param[k].mx)*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[0].sy=new_param[0].sy+1.0*(j-new_param[k].my)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[0].sxy=new_param[0].sxy+1.0*(i-new_param[k].mx)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[0].sx=1.0*new_param[0].sx/Ntot;
		new_param[0].sy=1.0*new_param[0].sy/Ntot;		
		new_param[0].sxy=1.0*new_param[0].sxy/Ntot;		
	
	for (k=1;k<nb_classe;k++)
		{
		new_param[k].sx=1.0*new_param[0].sx;
		new_param[k].sy=1.0*new_param[0].sy;		
		new_param[k].sxy=1.0*new_param[0].sxy;		
		}
		
free_dmatrix_3d(proba_apriori);
}

/*******************************************************************************
**      					  estime_gaussienne2d_isosxsy_EM                                                
**                                                                                                   
*******************************************************************************/
void estime_gaussienne2d_isosxsy_EM(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe)
{
double ***proba_apriori;
int wdth,hght,i,j,k,Ntot;
double tmp;
wdth=histo->width;
hght=histo->height;

/* allocation memoire */
proba_apriori=alloc_dmatrix_3d(wdth,hght,nb_classe);

/* remplissage de proba_apriori */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		for (k=0;k<nb_classe;k++)
			{proba_apriori[i][j][k]=	eval_gaussienne_2d(i, j,&param[k]);
			 if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
			}

Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		Ntot=Ntot+histo->mri[i][j];

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (histo->mri[i][j]>0)
	{tmp=0;
	for (k=0;k<nb_classe;k++)
	tmp=tmp+proba_apriori[i][j][k];
	
	for (k=0;k<nb_classe;k++)
	proba_apriori[i][j][k]=1.0*proba_apriori[i][j][k]/tmp;
	
	}



/* mise a jour du facteur de normalisation */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].l=0.0;
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].l=new_param[k].l+1.0*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].l=1.0*new_param[k].l/Ntot;
		}

/* mise a jour des moyennes */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].mx=0.0;
		new_param[k].my=0.0;
		
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].mx=new_param[k].mx+1.0*i*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[k].my=new_param[k].my+1.0*j*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[k].mx=1.0*new_param[k].mx/(Ntot*new_param[k].l);
		new_param[k].my=1.0*new_param[k].my/(Ntot*new_param[k].l);		
		}

/* mise a jour des variances */
		new_param[0].sx=0.0;
		new_param[0].sy=0.0;
		
	for (k=0;k<nb_classe;k++)
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[0].sx=new_param[0].sx+1.0*(i-new_param[k].mx)*(i-new_param[k].mx)*histo->mri[i][j]*proba_apriori[i][j][k];
			new_param[0].sy=new_param[0].sy+1.0*(j-new_param[k].my)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		new_param[0].sx=1.0*new_param[0].sx/Ntot;
		new_param[0].sy=1.0*new_param[0].sy/Ntot;		
	
	for (k=1;k<nb_classe;k++)
		{
		new_param[k].sx=1.0*new_param[0].sx;
		new_param[k].sy=1.0*new_param[0].sy;		
		}

/* mise a jour des variances */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].sxy=0.0;
		
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].sxy=new_param[k].sxy+1.0*(i-new_param[k].mx)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			}
			new_param[k].sxy=1.0*new_param[k].sxy/(Ntot*new_param[k].l);		
		}

		
free_dmatrix_3d(proba_apriori);
}

/*******************************************************************************
**        				 trimmed_histo_gaussienne_2d                                                
**                                                                                                   
*******************************************************************************/

void  trimmed_histo_gaussienne_2d(grphic *histo,grphic *trimmedhisto, Param_Gaussienne_2d *param,int nb_classe, double percentage)
{
int i,j,l,wdth,hght,cl=0,*classe,k;
grphic *histo_cl;
double max_prob,tmp,*proba_classe,**histo_proba,*seuil;

wdth=histo->width;
hght=histo->height;

/* allocation */
histo_cl=cr_grphic_modif(wdth,hght,0,1,0);
classe=(int *)malloc(nb_classe*sizeof(int));
seuil=(double *)malloc(nb_classe*sizeof(double));
histo_proba=alloc_dmatrix(wdth,hght);

for (l=0; l<nb_classe;l++)
		classe[l]=0;

/* classification au sens du maximum de vraisemblance */
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	{
	max_prob=0.0;
	for (l=0; l<nb_classe;l++)
		{
		tmp=eval_gaussienne_2d(i, j,&param[l]);
		if (tmp>max_prob) {max_prob=tmp; cl=l;}
		}
	histo_cl->mri[i][j]=cl;
	histo_proba[i][j]=max_prob;
	classe[cl]++;
	}

for (l=0; l<nb_classe;l++)
	{
	proba_classe=(double* )malloc(classe[l]*sizeof(double));
	k=0;
		for (i=0; i<wdth; i++)
		for (j=0; j<hght; j++)
			if (histo_cl->mri[i][j]==l)
				{
				proba_classe[k]=histo_proba[i][j];
				k++;
				}
	tri_rapide_double(proba_classe,0,classe[l]-1);
	seuil[l]=proba_classe[(int)floor(classe[l]*(1.0-percentage))];		
	free(proba_classe);
	}

/* seuillage de l'histo */
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		if (histo_proba[i][j]>=seuil[histo_cl->mri[i][j]])
		trimmedhisto->mri[i][j]=histo->mri[i][j];
		else
		trimmedhisto->mri[i][j]=0;


free_grphic(histo_cl);free(classe);free(seuil);
free_dmatrix(histo_proba,wdth,hght);
}

/*******************************************************************************
**        				 trimmed_histo_global_gaussienne_2d                                                
**                                                                                                   
*******************************************************************************/

void  trimmed_histo_global_gaussienne_2d(grphic *histo,grphic *trimmedhisto, Param_Gaussienne_2d *param,int nb_classe, double percentage)
{
int i,j,l,wdth,hght,Ntot,k,size;
double tmp,**histo_proba,**histo_erreur,*liste_erreur,seuil;

wdth=histo->width;
hght=histo->height;

/* allocation */
histo_proba=alloc_dmatrix(wdth,hght);
histo_erreur=alloc_dmatrix(wdth,hght);
liste_erreur=(double *) malloc(wdth*hght*sizeof(double));

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	{tmp=0;
	for (l=0; l<nb_classe; l++)
		tmp=tmp+eval_gaussienne_2d(i, j,&param[l]);
	
	histo_proba[i][j]=tmp;
	}

tmp=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
		tmp=tmp+histo_proba[i][j];

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
		histo_proba[i][j]=histo_proba[i][j]/tmp;

Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
		Ntot=Ntot+histo->mri[i][j];

size=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		{histo_erreur[i][j]=fabs(histo_proba[i][j]/histo->mri[i][j]*Ntot-1);size++;}
	else	
		histo_erreur[i][j]=-1;

k=0;	
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo_erreur[i][j]>-1)
		{liste_erreur[k]=histo_erreur[i][j];
		k++;}
			
tri_rapide_double(liste_erreur,0,size-1);

seuil=liste_erreur[(int)floor(size*percentage)];

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if ((histo->mri[i][j]>0)&&(histo_erreur[i][j]<=seuil))
	trimmedhisto->mri[i][j]=histo->mri[i][j];
	else
	trimmedhisto->mri[i][j]=0;	

free_dmatrix(histo_proba,wdth,hght);
free_dmatrix(histo_erreur,wdth,hght);
free(liste_erreur);
}

/*******************************************************************************
**        				 trimmed_histo_gaussienne_2d                                                
**                                                                                                   
*******************************************************************************/

void  trimmed_histo_local_gaussienne_2d(grphic *histo,grphic *trimmedhisto, Param_Gaussienne_2d *param,int nb_classe, double percentage)
{
int i,j,l,wdth,hght,cl=0,*classe,k,Ntot;
grphic *histo_cl;
double max_prob,tmp,*proba_classe,**histo_proba,*seuil,ptot;

wdth=histo->width;
hght=histo->height;

/* allocation */
histo_cl=cr_grphic_modif(wdth,hght,0,1,0);
classe=(int *)malloc(nb_classe*sizeof(int));
seuil=(double *)malloc(nb_classe*sizeof(double));
histo_proba=alloc_dmatrix(wdth,hght);

for (l=0; l<nb_classe;l++)
		classe[l]=0;

/* classification au sens du maximum de vraisemblance */
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	{
	max_prob=0.0;
	for (l=0; l<nb_classe;l++)
		{
		tmp=eval_gaussienne_2d(i, j,&param[l]);
		if (tmp>max_prob) {max_prob=tmp; cl=l;}
		}
	histo_cl->mri[i][j]=cl;
	histo_proba[i][j]=max_prob;
	classe[cl]++;
	}

for (l=0; l<nb_classe;l++)
	{
	proba_classe=(double* )malloc(classe[l]*sizeof(double));
	k=0;
	Ntot=0;
	for (i=0; i<wdth; i++)
		for (j=0; j<hght; j++)
			if ((histo_cl->mri[i][j]==l)&&(histo->mri[i][j]>0))
			Ntot++;
	ptot=0;
	for (i=0; i<wdth; i++)
		for (j=0; j<hght; j++)
			if (histo_cl->mri[i][j]==l)
				ptot=ptot+histo_proba[i][j];
			
	for (i=0; i<wdth; i++)
		for (j=0; j<hght; j++)
			if ((histo_cl->mri[i][j]==l)&&(histo->mri[i][j]>0))
				{
				proba_classe[k]=fabs(histo_proba[i][j]/histo->mri[i][j]*Ntot/ptot-1);
				k++;
				}
	tri_rapide_double(proba_classe,0,Ntot-1);
	seuil[l]=proba_classe[(int)floor(Ntot*percentage)];	
	for (i=0; i<wdth; i++)
		for (j=0; j<hght; j++)
		if ((histo_cl->mri[i][j]==l)&&(histo->mri[i][j]>0))
			{	
				if (fabs(histo_proba[i][j]/histo->mri[i][j]*Ntot/ptot-1)<=seuil[l])
				trimmedhisto->mri[i][j]=histo->mri[i][j];
				else
				trimmedhisto->mri[i][j]=0;
			}
			else
				trimmedhisto->mri[i][j]=0;	
	free(proba_classe);
	}

/* seuillage de l'histo */


free_grphic(histo_cl);free(classe);free(seuil);
free_dmatrix(histo_proba,wdth,hght);
}


/*******************************************************************************
**      					  robust_estime_gaussienne2d_EM                                                
**                                                                                                   
*******************************************************************************/
void robust_estime_gaussienne2d_EM(grphic *histo,Param_Gaussienne_2d *param,Param_Gaussienne_2d *new_param,int nb_classe)
{
double ***proba_apriori;
int wdth,hght,i,j,k,Ntot,l=0,m;
double tmp,max,*datax,*datay;
double ui,vi,Dem1,Dem2,Num;
wdth=histo->width;
hght=histo->height;

/* allocation memoire */
proba_apriori=alloc_dmatrix_3d(wdth,hght,nb_classe);

/* remplissage de proba_apriori */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		for (k=0;k<nb_classe;k++)
			{proba_apriori[i][j][k]=	eval_gaussienne_2d(i, j,&param[k]);
			 if(isnan(proba_apriori[i][j][k])) 
			 		proba_apriori[i][j][k]=0;
			}

Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		Ntot=Ntot+histo->mri[i][j];

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (histo->mri[i][j]>0)
	{tmp=0;
	for (k=0;k<nb_classe;k++)
	tmp=tmp+proba_apriori[i][j][k];
	
	for (k=0;k<nb_classe;k++)
	proba_apriori[i][j][k]=1.0*proba_apriori[i][j][k]/tmp;
	
	}



/* binarisation de la proba a posteriori */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
	if (histo->mri[i][j]>0)
		{
		max=0;
		for (k=0;k<nb_classe;k++)
			if (proba_apriori[i][j][k]>max)
				{	max=proba_apriori[i][j][k]; l=k;}
		
		for (k=0;k<nb_classe;k++)
			if (k==l)
				proba_apriori[i][j][k]=1;
			else
				proba_apriori[i][j][k]=0;		
		}
		
		
		
/* mise a jour du facteur de normalisation */
for (k=0;k<nb_classe;k++)
		{
		new_param[k].l=0.0;
		for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].l=new_param[k].l+1.0*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		}

/* mise a jour des moyennes */
for (k=0;k<nb_classe;k++)
		{
		
	datax=(double *)malloc(new_param[k].l*sizeof(double));
	datay=(double *)malloc(new_param[k].l*sizeof(double));
	
	l=0;	
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		if (proba_apriori[i][j][k]>0)
		for (m=0;m<histo->mri[i][j];m++)
			{
			datax[l]=i;
			datay[l]=j;
			l++;
			}
	qsort(datax,new_param[k].l,sizeof(double),double_compare_function2);
	qsort(datay,new_param[k].l,sizeof(double),double_compare_function2);
	
	new_param[k].mx=datax[(int)new_param[k].l/2];
	new_param[k].my=datay[(int)new_param[k].l/2];
	
	for (l=0;l<new_param[k].l;l++)
		{
		datax[l]=fabs(datax[l]-new_param[k].mx);
		datay[l]=fabs(datay[l]-new_param[k].my);
		}	
		
	qsort(datax,new_param[k].l,sizeof(double),double_compare_function2);
	qsort(datay,new_param[k].l,sizeof(double),double_compare_function2);
	
	//new_param[k].sx=pow(1.4826*datax[(int)new_param[k].l/2],2.0);
	//new_param[k].sy=pow(1.4826*datay[(int)new_param[k].l/2],2.0);
	
	new_param[k].sx=datax[(int)new_param[k].l/2];
	new_param[k].sy=datay[(int)new_param[k].l/2];
	


		new_param[k].sxy=0.0;
		
		/*for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		if (histo->mri[i][j]>0)
			{
			new_param[k].sxy=new_param[k].sxy+1.0*(i-new_param[k].mx)*(j-new_param[k].my)*histo->mri[i][j]*proba_apriori[i][j][k];
			}
		
		new_param[k].sxy=1.0*new_param[k].sxy/(new_param[k].l);		
		*/
		Num=Dem1=Dem2=0;
		for (i=0;i<wdth;i++)
			for (j=0;j<hght;j++)
				if (proba_apriori[i][j][k]>0)
				{
				ui=(1.0*i-new_param[k].mx)/(9.0*new_param[k].sx);
				vi=(1.0*j-new_param[k].my)/(9.0*new_param[k].sy);
				
				if (fabs(ui)<1)
					Dem1+= (1.0-ui*ui)*(1.0-5.0*ui*ui)*histo->mri[i][j];
					
				if (fabs(vi)<1)
					Dem2+= (1.0-vi*vi)*(1.0-5.0*vi*vi)*histo->mri[i][j];
					
				if ((fabs(ui)<1)&&(fabs(vi)<1))
					Num+= (1.0*i-new_param[k].mx)*pow((1.0-ui*ui),2.0)*(1.0*j-new_param[k].my)*pow((1.0-vi*vi),2.0)*histo->mri[i][j];
				}
			new_param[k].sxy=new_param[k].l*Num/(Dem1*Dem2);

		
		new_param[k].sx=pow(1.4826*new_param[k].sx,2.0);
		new_param[k].sy=pow(1.4826*new_param[k].sy,2.0);
		free(datax);free(datay);	
		}


for (k=0;k<nb_classe;k++)
		{
		new_param[k].l=1.0*new_param[k].l/Ntot;
		}
		
free_dmatrix_3d(proba_apriori);
}







/*******************************************************************************
********************************************************************************
********************************************************************************
*****  Suppression de points non significatifs dans l'histogramme joint  ******* 
********************************************************************************
********************************************************************************
********************************************************************************/

/*******************************************************************************
**        				 trimmed_histo                                                
**                                                                                                   
*******************************************************************************/
void  trimmed_histo(grphic *histo,grphic *trimmedhisto,double percentage)
{
int i,j,k,wdth,hght,Ntot;
long int *liste,seuil;

wdth=histo->width;
hght=histo->height;

Ntot=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	if (histo->mri[i][j]>0)
	Ntot++;
	
liste=(long int *)malloc(Ntot*sizeof(long int));

k=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (histo->mri[i][j]>0)
	{
	liste[k]=histo->mri[i][j];
	k++;
	}

tri_rapide(liste,0,Ntot-1);
seuil=liste[(int)floor(Ntot*(1.0-percentage))];

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (histo->mri[i][j]>=seuil)
	trimmedhisto->mri[i][j]=histo->mri[i][j];
else
	trimmedhisto->mri[i][j]=0;

}

/*******************************************************************************
**        				 trimmed_info_histo                                                
**                                                                                                   
*******************************************************************************/
void  trimmed_info_histo(grphic *histo,grphic *trimmedhisto,double nb_sigma)
{
int i,j,wdth,hght,l;
int *nbx,*nby;
double *seuilx,*seuily;

wdth=histo->width;
hght=histo->height;

/* allocation */
nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));
seuilx=(double *)malloc(wdth*sizeof(double));
seuily=(double *)malloc(hght*sizeof(double));

/* initialisation */
for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}


for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}


/* calcul des seuils */
for (l=0; l<wdth; l++)
	{seuilx[l]=1.0*nbx[l]/hght+1.0*nb_sigma*sqrt(nbx[l]*(hght-1)/hght/hght);}
for (l=0; l<hght; l++)
	{seuily[l]=1.0*nby[l]/wdth+1.0*nb_sigma*sqrt(nby[l]*(wdth-1)/wdth/wdth);}


for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (histo->mri[i][j]>=MAXI(seuilx[i],seuily[j]))
	trimmedhisto->mri[i][j]=histo->mri[i][j];
else
	trimmedhisto->mri[i][j]=0;


free(nbx);free(nby);
free(seuilx);free(seuily);
}

/*******************************************************************************
**        				 trimmed_info_histo                                                
**                                                                                                   
*******************************************************************************/
void  trimmed_info_xy_histo(grphic *histo,grphic *trimmedhisto,double nb_sigma)
{
int i,j,wdth,hght,l;
int *nbx,*nby,Ntot;
double moyenne,seuil;

wdth=histo->width;
hght=histo->height;

/* allocation */
nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

/* initialisation */
for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}


for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}


Ntot=0;
for (l=0; l<wdth; l++)
Ntot=Ntot+nbx[l];




for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		//moyenne=1.0*nbx[i]*nby[j]/(1.0*wdth*hght*Ntot);
		moyenne=1.0*nbx[i]*nby[j]/(1.0*Ntot);
		
		seuil=moyenne+1.0*nb_sigma*sqrt(moyenne*(1.0-1.0*moyenne/Ntot));
		
			if (histo->mri[i][j]>seuil)
			trimmedhisto->mri[i][j]=histo->mri[i][j]-seuil;	
			else
			trimmedhisto->mri[i][j]=0;

		}

free(nbx);free(nby);
}






/*******************************************************************************
********************************************************************************
********************************************************************************
*******  Reclassement de points non significatifs de l'histogramme joint  ****** 
********************************************************************************
********************************************************************************
********************************************************************************/

/**************************** histo_joint_maxinfo_3d_p() *****************************/
/*                                                                           */                         
/*****************************************************************************/

void histo_joint_maxinfo_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i,j,l,wdth,hght,m,n,i0,j0,continu,nb,Ntot/*,Npart*/;
double tmp,**produit,toto/*,**proba_histo*/;
double /*pij,pin,pmj,pmn,pijt,pint,pmjt,pmnt,*//*,moyenne,sigma*/seuil,min,max;
double *mIj,*mJi,delta;
int *nby,*nbx;
//grphic *histo_nonrecale;

//implot->width=MAX_WIDTH;
//implot->height=MAX_HEIGHT;
//implot->width=64;
//implot->height=64;

histo_joint_linear_3d_p(im1,im2,implot);
//histo_joint_linear_norm_3d_p(im1,im2,implot);
//histo_joint_maxentropie_3d_p(im1,im2,implot);


if (implot->pos!=0)
show_picture(implot->pos);

wdth=implot->width;
hght=implot->height;

//histo_nonrecale=cr_grphic_modif(wdth,hght,0,1,0);
//histo_joint_maxentropie_3d_p(im1,im2,histo_nonrecale);

produit=alloc_dmatrix(wdth,hght);
//proba_histo=alloc_dmatrix(wdth,hght);

/*for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (histo_nonrecale->mri[i][j]>0)
			proba_histo[i][j]=1.0*implot->mri[i][j]/histo_nonrecale->mri[i][j];
else
			proba_histo[i][j]=1.0*implot->mri[i][j];  
*/			
			
mIj=(double *)malloc(hght*sizeof(double));
mJi=(double *)malloc(wdth*sizeof(double));
nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

for (l=0; l<hght; l++)
	{nby[l]=0;}

for (l=0; l<wdth; l++)
	{nbx[l]=0;}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+implot->mri[i][j];	
		}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+implot->mri[i][j];	
		}
			

Ntot=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
 Ntot=Ntot+implot->mri[i][j];



max=0;min=HUGE_VAL;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
{
if (implot->mri[i][j]>max)
	{max=implot->mri[i][j];}

if ((implot->mri[i][j]<min)&&(implot->mri[i][j]>0))
	{min=implot->mri[i][j];}

}


/*max=0;min=HUGE_VAL;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
{
if (proba_histo[i][j]>max)
	{max=proba_histo[i][j];}

if (proba_histo[i][j]<min)
	{min=proba_histo[i][j];}

}
*/

/* calcul de mIj */
for (j=0; j<hght; j++)
	{
	mIj[j]=0;
	for (i=0; i<wdth; i++)
	mIj[j]=mIj[j]+i*implot->mri[i][j];
	
	mIj[j]=mIj[j]/nby[j];
	
	}

for (i=0; i<wdth; i++)
	{
	mJi[i]=0;
	for (j=0; j<hght; j++)
	mJi[i]=mJi[i]+j*implot->mri[i][j];
	
	mJi[i]=mJi[i]/nbx[i];
	
	}


seuil=min;
//seuil=max;

continu=0;
nb=0;
while((continu==0)/*&&(nb<500)*/)
{
/*Npart=0;
moyenne=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (implot->mri[i][j]>0)
{moyenne=moyenne+implot->mri[i][j];
Npart++;
}
moyenne=1.0*moyenne/Npart;

seuil=moyenne;
*/
/*seuil=0;
for (i=0; i<wdth; i++)
{
min=HUGE_VAL;
for (j=0; j<hght; j++)
if (implot->mri[i][j]>0)
if (implot->mri[i][j]<min)
	{min=implot->mri[i][j];}

seuil=MINI(seuil,min);
}

seuil=(0.9*seuil+0.1*moyenne);
*/
continu=1;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
{
		
if ((implot->mri[i][j]>0)&&(implot->mri[i][j]<=seuil))
		{
		printf("%d  %d \r",i,j);
		for (m=0; m<wdth; m++)
			for (n=0; n<hght; n++)
			/*for (m=MAXI(0,i-voisinage); m<MINI(wdth,i+voisinage); m++)
			for (n=MAXI(0,j-voisinage); n<MINI(hght,j+voisinage); n++)*/
			if ((m!=i)&&(n!=j)&&(implot->mri[m][n]>0)&&(implot->mri[i][n]>0)&&(implot->mri[m][j]>0)&&(implot->mri[i][n]>seuil)&&(implot->mri[m][j]>seuil))
			{
			tmp=1.0*MINI(implot->mri[i][j],implot->mri[m][n]);
			delta=1.0*tmp/Ntot;
			/*pij=1.0*implot->mri[i][j]/Ntot;
			pijt=1.0*(implot->mri[i][j]-tmp)/Ntot;
			if (pijt==0) pijt=1;
			
			pmn=1.0*implot->mri[m][n]/Ntot;
			pmnt=1.0*(implot->mri[m][n]-tmp)/Ntot;
			if (pmnt==0) pmnt=1;
			
			pin=1.0*implot->mri[i][n]/Ntot;
			pint=1.0*(implot->mri[i][n]+tmp)/Ntot;
			pmj=1.0*implot->mri[m][j]/Ntot;
			pmjt=1.0*(implot->mri[m][j]+tmp)/Ntot;
			
			toto=(pijt)*log(pijt)-(pij)*log(pij)
					+(pmnt)*log(pmnt)-(pmn)*log(pmn)
					+(pint)*log(pint)-(pin)*log(pin)
					+(pmjt)*log(pmjt)-(pmj)*log(pmj);*/
			
			toto=2.0*delta*(i-m)*(mIj[j]-mIj[n])-1.0*(i-m)*(i-m)*delta*delta*Ntot*(1.0/nby[j]+1.0/nby[n]);
			toto=toto+2.0*delta*(j-n)*(mJi[i]-mJi[m])-1.0*(j-n)*(j-n)*delta*delta*Ntot*(1.0/nbx[i]+1.0/nbx[m]);
			produit[m][n]=toto;		
			
			}
			else
			produit[m][n]=0;
		
		/* recherche du min */
		i0=i;
		j0=j;
		//tmp=implot->mri[i][j]*implot->mri[j][j];
		tmp=0;
		//tmp=2.0*implot->mri[i][j];
		for (m=0; m<wdth; m++)
			for (n=0; n<hght; n++)
				if (produit[m][n]<tmp)
					{tmp=produit[m][n]; i0=m;j0=n;}
		
		if ((i0!=i)&&(j0!=j))
			{
			//printf("diminution de %f \n",tmp);
			tmp=MINI(implot->mri[i][j],implot->mri[i0][j0]);
			delta=1.0*tmp;
			implot->mri[i][j]=implot->mri[i][j]-tmp;
			implot->mri[i0][j0]=implot->mri[i0][j0]-tmp;
			implot->mri[i][j0]=implot->mri[i][j0]+tmp;
			implot->mri[i0][j]=implot->mri[i0][j]+tmp;
			mIj[j]=mIj[j]+1.0*(i0-i)*delta/nby[j];
			mIj[j0]=mIj[j0]+1.0*(i-i0)*delta/nby[j0];
			mJi[i]=mJi[i]+1.0*(j0-j)*delta/nbx[i];
			mJi[i0]=mJi[i0]+1.0*(j-j0)*delta/nbx[i0];
			
			continu=0;
			}			
		}
	}
if ((continu==1)||(nb>10))
	{seuil=seuil+(max-min)/100;continu=0;nb=0;}
					
//printf("seuil %f \n",seuil);

if (seuil>max)
	continu=1;
	
if (implot->pos!=0)
show_picture(implot->pos);

nb++;
}

/*for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]-implot->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]-implot->mri[i][j];	
		}

for (i=0; i<wdth; i++)
	if (nbx[i]!=0)
		printf("Y a un bug !!!!!!!!!!!\n");
		
for (i=0; i<hght; i++)
	if (nby[i]!=0)
		printf("Y a un bug !!!!!!!!!!!\n");
*/		
	
//free_grphic(histo_nonrecale);	
free_dmatrix(produit,wdth,hght);
//free_dmatrix(proba_histo,wdth,hght);
free(mIj);free(nby);
free(mJi);free(nbx);
}


/**************************** histo_joint_maxinfo_3d_p() *****************************/
/*                                                                           */                         
/*****************************************************************************/

void histo_joint_maxinfo_multimodal_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i,j,l,wdth,hght,m,n,i0,j0,continu,nb,nb_classe=4,Ntot,V=15,x,y=0,indx[4],indy[4],signe;
double tmp,**produit,toto,tutu;
double min,max,seuil;
Param_Gaussienne_2d *param,*param_tmp;

grphic *histo_cl;

implot->width=MAX_WIDTH;
implot->height=MAX_HEIGHT;
implot->width=64;
implot->height=64;

histo_joint_linear_3d_p(im1,im2,implot);
implot->zoom=(int)1.0*MAX_WIDTH/implot->width;

if (implot->pos!=0)
show_picture(implot->pos);

wdth=implot->width;
hght=implot->height;

histo_cl=cr_grphic_modif(wdth,hght,0,1,0);

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	histo_cl->mri[i][j]=0;
	  
produit=alloc_dmatrix(wdth,hght);
			
param=(Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param_tmp=(Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));

//mean_std_classe(implot,histo_cl,nb_classe,param);

for (l=0;l<nb_classe;l++)
	{	
	param_tmp[l].mx=param[l].mx;
	param_tmp[l].my=param[l].my;
	param_tmp[l].sx=param[l].sx;
	param_tmp[l].sy=param[l].sy;
	param_tmp[l].sxy=param[l].sxy;
	param_tmp[l].l=param[l].l;
	}

Ntot=0;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
	Ntot=Ntot+implot->mri[i][j];
	
max=0;min=HUGE_VAL;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
{
if (implot->mri[i][j]>max)
	{max=implot->mri[i][j];}

if ((implot->mri[i][j]<min)&&(implot->mri[i][j]>0))
	{min=implot->mri[i][j];}

}





seuil=min;

continu=0;
nb=0;
while (continu==0)
{
/*kmeans_histo_joint_3d_p(implot,histo_cl,nb_classe);*/
			
continu=1;
for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
{
indx[0]=i;
indy[0]=j;
indx[2]=i;
indy[3]=j;		

if ((implot->mri[i][j]>0)&&(implot->mri[i][j]<=seuil))
		{
		printf("%d  %d \r",i,j);
		for (m=0; m<wdth; m++)
			for (n=0; n<hght; n++)
			if ((m!=i)&&(n!=j)&&(implot->mri[m][n]>0)&&(implot->mri[i][n]>0)&&(implot->mri[m][j]>0))
			{
			tmp=1.0*MINI(implot->mri[i][j],implot->mri[m][n]);
			
			/* calcul du critere */
			indx[1]=m;
			indy[1]=n;
			indy[2]=n;
			indx[3]=m;		
		toto=tutu=0;
		for (l=0;l<4;l++)
			{
			
			if (l<2)
			signe=-1;
			else
			signe=1;
			
			param[l].mx=0;
			param[l].my=0;
			param[l].sx=0;
			param[l].sy=0;
			param[l].sxy=0;
			param[l].l=0;
			
			for (x=MAXI(indx[l]-V,0);x<MINI(indx[l]+V,wdth);x++)
			for (y=MAXI(indy[l]-V,0);y<MINI(indy[l]+V,hght);y++)
				{
				param[l].l=param[l].l+implot->mri[x][y];
				param[l].mx=param[l].mx+implot->mri[x][y]*x;
				param[l].my=param[l].my+implot->mri[x][y]*y;
				param[l].sx=param[l].sx+implot->mri[x][y]*x*x;
				param[l].sy=param[l].sy+implot->mri[x][y]*y*y;
				param[l].sxy=param[l].sxy+implot->mri[x][y]*x*y;
				}
				
			param[l].mx=1.0*param[l].mx/param[l].l;
			param[l].my=1.0*param[l].my/param[l].l;
			param[l].sx=1.0*param[l].sx/param[l].l-param[l].mx*param[l].mx;
			param[l].sy=1.0*param[l].sy/param[l].l-param[l].my*param[l].my;
			param[l].sxy=1.0*param[l].sxy/param[l].l-param[l].mx*param[l].my;
			
			
			param_tmp[l].l=param[l].l+signe*tmp;
			param_tmp[l].mx=1.0*(param[l].l*param[l].mx+signe*tmp*x)/param_tmp[l].l;
			param_tmp[l].my=1.0*(param[l].l*param[l].my+signe*tmp*y)/param_tmp[l].l;
			param_tmp[l].sx=1.0*(param[l].l*param[l].sx+signe*tmp*x*x)/param_tmp[l].l-param_tmp[l].mx*param_tmp[l].mx;
			param_tmp[l].sy=1.0*(param[l].l*param[l].sy+signe*tmp*y*y)/param_tmp[l].l-param_tmp[l].my*param_tmp[l].my;
			param_tmp[l].sxy=1.0*(param[l].l*param[l].sxy+signe*tmp*x*y)/param_tmp[l].l-param_tmp[l].mx*param_tmp[l].my;
			
			toto=toto+param_tmp[l].l*(param_tmp[l].sx*param_tmp[l].sy-param_tmp[l].sxy*param_tmp[l].sxy);
			tutu=tutu+param[l].l*(param[l].sx*param[l].sy-param[l].sxy*param[l].sxy);
			}
				
			
			/*************************/
			if (toto<tutu)
			produit[m][n]=toto;		
			else
			produit[m][n]=HUGE_VAL;
			}
			else
			produit[m][n]=HUGE_VAL;
		
		/* recherche du min */
		i0=i;
		j0=j;
		tmp=HUGE_VAL;
		for (m=0; m<wdth; m++)
			for (n=0; n<hght; n++)
				if (produit[m][n]<tmp)
					{tmp=produit[m][n]; i0=m;j0=n;}
		
		if ((i0!=i)&&(j0!=j))
			{
			//printf("diminution de %f \n",tmp);
			tmp=MINI(implot->mri[i][j],implot->mri[i0][j0]);
			implot->mri[i][j]=implot->mri[i][j]-tmp;
			implot->mri[i0][j0]=implot->mri[i0][j0]-tmp;
			implot->mri[i][j0]=implot->mri[i][j0]+tmp;
			implot->mri[i0][j]=implot->mri[i0][j]+tmp;
			
			/* mise a jour de param */
		
			
			
			/*************************/
			continu=0;
			}			
		}
	}
if ((continu==1)||(nb>10))
	{seuil=seuil+(max-min)/100;continu=0;nb=0;}
					
//printf("seuil %f \n",seuil);

if (seuil>max)
	continu=1;
	
if (implot->pos!=0)
show_picture(implot->pos);

nb++;
}

	
free_grphic(histo_cl);	
free_dmatrix(produit,wdth,hght);
free(param);free(param_tmp);
}






/*******************************************************************************
********************************************************************************
********************************************************************************
*******************   Divers fonctions relatives a l'EM  *********************** 
********************************************************************************
********************************************************************************
********************************************************************************/

/*******************************************************************************
**        				 plot_iso_gaussienne_2d                                                
**                                                                                                   
*******************************************************************************/

void  plot_iso_gaussienne_2d(grphic *histo,Param_Gaussienne_2d *param,double n_iso, int value)
{
double sigma,delta,y1,y2;
int i;
sigma=param->sx*param->sy-param->sxy*param->sxy;

for (i=0;i<histo->width;i++)
	{
	delta=(2.0*(i-param->mx)*param->sxy)*(2.0*(i-param->mx)*param->sxy)-4.0*param->sx*(param->sy*(i-param->mx)*(i-param->mx)-1.0*n_iso*sigma);
	
	if (delta>=0)
		{
		y1=1.0*(2.0*(i-param->mx)*param->sxy-sqrt(delta))/(2.0*param->sx)+param->my;
		y2=1.0*(2.0*(i-param->mx)*param->sxy+sqrt(delta))/(2.0*param->sx)+param->my;
		
		if ((y1>0)&&(y1<histo->height))
			histo->mri[i][(int)floor(y1)]=value;
		
		if ((y2>0)&&(y2<histo->height))
			histo->mri[i][(int)floor(y2)]=value;
				
		}
	
	}
}

/*******************************************************************************
**      					  distance_mahalanobis_2d                                                
**                                                                                                   
*******************************************************************************/
double distance_mahalanobis_2d(double x, double y, Param_Gaussienne_2d *param)
{
double sigma,res;
sigma=param->sx*param->sy-param->sxy*param->sxy;
res=sqrt(1.0*((x-param->mx)*(x-param->mx)*param->sy+(y-param->my)*(y-param->my)*param->sx
			-2.0*(x-param->mx)*(y-param->my)*param->sxy)
			/(1.0*sigma));

return(res);
}

/*******************************************************************************
**      					 verif_EM                                               
**                                                                                                   
*******************************************************************************/

void verif_EM(void)
{
grphic *histo;
Param_Gaussienne_2d *param,*est_param;
int nb_classe=3,l,i,j;
double tmp;
/* allocation memoire */
histo=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
param= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
est_param= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));

est_param[0].mx=50;
est_param[0].my=110;
est_param[0].sx=15;
est_param[0].sy=20;
est_param[0].sxy=5;
est_param[0].l=0.5;

est_param[1].mx=100;
est_param[1].my=40;
est_param[1].sx=20;
est_param[1].sy=40;
est_param[1].sxy=10;
est_param[1].l=0.3;
 
est_param[2].mx=110;
est_param[2].my=90;
est_param[2].sx=16;
est_param[2].sy=40;
est_param[2].sxy=10;
est_param[2].l=0.2;


/*param[0].mx=50;
param[0].my=110;
param[0].sx=20;
param[0].sy=20;
param[0].sxy=0;
param[0].l=0.7;

param[1].mx=100;
param[1].my=40;
param[1].sx=20;
param[1].sy=40;
param[1].sxy=0;
param[1].l=0.3;*/
/*
est_param[2].mx=50;
est_param[2].my=200;
est_param[2].sx=1600;
est_param[2].sy=1600;
est_param[2].sxy=1000;
est_param[2].l=0.2;
*/
for (i=0;i<TOPO_SIZE_NORM_HISTO;i++)
for (j=0;j<TOPO_SIZE_NORM_HISTO;j++)
	{tmp=0;
	for (l=0;l<nb_classe;l++)
	tmp=tmp+1000000*eval_gaussienne_2d(i,j,&est_param[l]);
	
	histo->mri[i][j]=(int)tmp;
	}



fit_gaussienne2d_EM(histo, param, nb_classe);


for (l=0;l<nb_classe;l++)
	{
	printf("mx 	: %f  %f \n",est_param[l].mx,param[l].mx);
	printf("my 	: %f  %f \n",est_param[l].my,param[l].my);
	printf("sx 	: %f  %f \n",est_param[l].sx,param[l].sx);
	printf("sy 	: %f  %f \n",est_param[l].sy,param[l].sy);
	printf("sxy : %f  %f \n",est_param[l].sxy,param[l].sxy);
	printf("l 	: %f  %f \n",est_param[l].l,param[l].l);
	printf("\n");
	
	}
	
free_grphic(histo);
free(param);free(est_param);

}





/*******************************************************************************
********************************************************************************
********************************************************************************
****   Fonctions interfacant la normalisation avec le recalage non-rigide  ***** 
********************************************************************************
********************************************************************************
********************************************************************************/

void TOP_norm_mean (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
imx_norm_seuil_mean_3d_p(im1,im2,imres);
}

void TOP_norm_meanecty (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
}

void TOP_norm_rl (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
imx_norm_seuil_rl_3d_p(im1,im2,imres);
}

void TOP_norm_mean_ratio (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
imx_norm_seuil_3d_p(im1,im2,imres);
}

void TOP_norm_minmax (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
imx_norm_seuil_minmax_3d_p(im1,im2,imres);
}

void TOP_norm_spline (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
#ifdef HAS_IMLIB3D
CPP_normalisation_histo(im1,im2,imres);
#else
printf("Attention Vous n'avez pas compile avec HAS_IMLIB3D ! La normalisation Moyenne Ecart Type est utilisee par defaut ...\n");
imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
#endif
}

void TOP_norm_GM (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
#ifdef HAS_ATLAS_PERFUSION

CPP_normalisation_GM(im1,im2,imres, nb_classe, 0.05);

#else
printf("Attention Vous n'avez pas compile avec HAS_ATLAS_PERFUSION ! La normalisation Moyenne Ecart Type est utilisee par defaut ...\n");
imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
#endif
}

void TOP_HistogramMatchingImageFilter_3d_p (grphic3d *im1, grphic3d *im2, grphic3d *imres, int nb_classe)
{
#ifdef  HAS_ITK
itk_HistogramMatchingImageFilter_3d_p(im1,im2,imres, 1024, nb_classe);
#else
printf("Attention Vous n'avez pas compile avec ITK ! La normalisation Moyenne Ecart Type est utilisee par defaut ...\n");
imx_norm_seuil_meanecty_3d_p(im1,im2,imres);
#endif
}
 

void norm_histojoint_spline_Bosc(void)
{
	int im_1,im_2,im_res;

	im_1=GET_PLACE3D(TEXT0231);
	im_2=GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	#ifdef HAS_IMLIB3D
	CPP_normalisation_histo(ptr_img_3d(im_1),ptr_img_3d(im_2),ptr_img_3d(im_res));
	#else
	printf("Attention Vous n'avez pas compile avec HAS_IMLIB3D ! La normalisation Moyenne Ecart Type est utilisee par defaut ...\n");
	imx_norm_seuil_meanecty_3d_p(ptr_img_3d(im_1),ptr_img_3d(im_2),ptr_img_3d(im_res));
	#endif

	show_picture_3d(im_res);
}


/*******************************************************************************
********************************************************************************
********************************************************************************
*******************   Differentes methodes de normalisation   ****************** 
********************************************************************************
********************************************************************************
********************************************************************************/

/*******************************************************************************
**        normalisation_histo_joint                                                
**                                                                                                   
*******************************************************************************/

void normalisation_histo_joint(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	nb_classe= GET_INT("nb de classes", 0, &e);
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	normalise_histo_joint(im1,im2,imres,nb_classe);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
** 	normalise_histo_joint(imreca,imref,imres,nb_classe)                                   
**                                                                   
**	Cette fonction permet dnormaliser l'histogramme conjoint entre 2 images       
*******************************************************************************/
void normalise_histo_joint(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe)
{
int **histojoint,**histotraite;
double *ft;
double max_reca,max_ref;
int i,j,k,wdth,hght,dpth;
int plot_x,plot_y,nbi,nbtmp,nbmid,itmp=0;
wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

if ((wdth!=imref->width)||(hght!=imref->height)||(dpth!=imref->depth))
	{printf("Les images ne sont pas de meme tailles !!!\n");
	return;
	}

if (nb_classe==0)
	{nb_classe=floor(2.0*imreca->max_pixel/1000);
	//printf("Nb de classe par defaut  : %d\n",nb_classe);
	}
	
/* Allocation memoire */
histojoint=alloc_imatrix(nb_classe,nb_classe);
histotraite=alloc_imatrix(nb_classe,nb_classe);
ft=malloc(nb_classe*sizeof(double));

/*initialisation*/
for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	histojoint[i][j]=0;
	
for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	histotraite[i][j]=0;
	
for (i=0;i<nb_classe;i++)
	ft[i]=0.0;

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

max_reca=max_reca*imreca->rcoeff;
max_ref=max_ref*imref->rcoeff;

if (max_reca>max_ref)
	imx_copie_param_3d_p(imreca,imres);
else
	imx_copie_param_3d_p(imref,imres);

/*remplissage de l'histogramme joint*/
for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
	if(imreca->mri[i][j][k]>0 && imref->mri[i][j][k]>0)
		{
		plot_x=floor(1.0*(imreca->mri[i][j][k]*imreca->rcoeff)/(max_reca)*(nb_classe-1.0));
		plot_y=floor(1.0*(imref->mri[i][j][k]*imref->rcoeff)/(max_ref)*(nb_classe-1.0));
		histojoint[plot_x][plot_y]++;
		}	

/*remplissage de histo traite */
for (i=0;i<nb_classe;i++)
	{
	nbi=0;
	for (j=0;j<nb_classe;j++)
	nbi=nbi+histojoint[i][j];
	
	nbmid=floor(1.0*nbi/2+0.5);
	nbtmp=histojoint[i][0];
	j=0;
	while(nbtmp<nbmid)
		{j++;
		nbtmp=nbtmp+histojoint[i][j];
		}
	histotraite[i][j]=histotraite[i][j]+nbi;
	}

for (i=0;i<nb_classe;i++)
	{
	nbi=0;
	for (j=0;j<nb_classe;j++)
	nbi=nbi+histojoint[j][i];
	
	nbmid=floor(1.0*nbi/2+0.5);
	nbtmp=histojoint[0][i];
	j=0;
	while(nbtmp<nbmid)
		{j++;
		nbtmp=nbtmp+histojoint[j][i];
		}
	histotraite[j][i]=histotraite[j][i]+nbi;
	}	
	
/*remplissage de la fonction de transfert */
for (i=0;i<nb_classe;i++)
	{
	nbi=0;
	for (j=0;j<nb_classe;j++)
		{
		ft[i]=ft[i]+j*histotraite[i][j];
		nbi=nbi+histotraite[i][j];
		}
		if (nbi!=0)
			ft[i]=1.0*ft[i]/nbi;
		ft[i]=ft[i]+1;
	}

/*On veut que la fonction de transfert soit croissante*/
for (i=0;i<nb_classe-1;i++)
	{
	if (ft[i]>=ft[i+1])
		{
		itmp=i+1;
		while((ft[itmp]<=ft[i])&&(itmp<nb_classe))
			itmp++;
		
		if (itmp==nb_classe)
			{ft[itmp-1]=nb_classe-1;itmp=itmp-1;}
	
		for (j=i+1;j<itmp;j++)
			{ft[j]=ft[i]+1.0*(j-i)*(ft[itmp]-ft[i])/(itmp-i);
			}
		}
	}

/* On applique la fonction de transfert a imreca */
for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
	if(imreca->mri[i][j][k]>0 /*&& imref->mri[i][j][k]>0*/)
		{
		plot_x=floor(1.0*(imreca->mri[i][j][k]*imreca->rcoeff)/(max_reca)*(nb_classe-1.0));
		
		if (plot_x==0)
		imres->mri[i][j][k]=(int)(1.0*imreca->rcoeff/imres->rcoeff*max_ref/max_reca*(ft[0]*(1.0*imreca->mri[i][j][k])));
		else
		imres->mri[i][j][k]=(int)(1.0/imres->rcoeff*max_ref/(nb_classe-1.0)*(ft[plot_x-1]+(ft[plot_x]-ft[plot_x-1])*(1.0*imreca->mri[i][j][k]*imreca->rcoeff*(nb_classe-1.0)/max_reca-1.0*plot_x)));
		}
	else
		imres->mri[i][j][k]=(int)(1.0*imreca->mri[i][j][k]*imreca->rcoeff/imres->rcoeff); 

/*liberation memoire*/
free_imatrix(histojoint,nb_classe,nb_classe);
free_imatrix(histotraite,nb_classe,nb_classe);
free(ft);
}


/*******************************************************************************
**        normalisation_kmeans                                                
**                                                                                                   
*******************************************************************************/

void normalisation_kmeans(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	nb_classe= GET_INT("nb de classes", 10, &e);
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	normalise_kmeans_p(im1,im2,imres,nb_classe);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
** 	normalise_kmeans_p                                  
**                                                                   
*******************************************************************************/
void normalise_kmeans_p(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe)
{
int **histojoint,*p_ref,*p_reca;
int i,j,k,wdth,hght,dpth,l,nb,nb_reca,nb_ref;
double *mean_reca,*mean_ref,*std_reca,*std_ref,max_reca,max_ref,tmp,/*p_tot_ref,p_tot_reca,*/p_tot,tutu;
grphic3d *imreca_cl,*imref_cl;
wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

if ((wdth!=imref->width)||(hght!=imref->height)||(dpth!=imref->depth))
	{printf("Les images ne sont pas de meme tailles !!!\n");
	return;
	}

	
/* Allocation memoire */
histojoint=alloc_imatrix(nb_classe,nb_classe);
p_ref=(int *)malloc(nb_classe*sizeof(int));
p_reca=(int *)malloc(nb_classe*sizeof(int));

mean_reca=(double *)malloc(nb_classe*sizeof(double));
mean_ref=(double *)malloc(nb_classe*sizeof(double));
std_reca=(double *)malloc(nb_classe*sizeof(double));
std_ref=(double *)malloc(nb_classe*sizeof(double));

imref_cl=cr_grphic3d(imref);
imreca_cl=cr_grphic3d(imreca);

/*initialisation*/
for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	histojoint[i][j]=0;

for (i=0;i<nb_classe;i++)
	p_ref[i]=p_reca[i]=0;

for (i=0;i<nb_classe;i++)
	mean_ref[i]=mean_reca[i]=0;

for (i=0;i<nb_classe;i++)
	std_ref[i]=std_reca[i]=0;

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

max_reca=max_reca*imreca->rcoeff;
max_ref=max_ref*imref->rcoeff;

if (max_reca>max_ref)
	imx_copie_param_3d_p(imreca,imres);
else
	imx_copie_param_3d_p(imref,imres);

/* kmeans sur les deux images */
#ifndef COMPILE_FOR_MEDIPY
imx_kmean_p(imref, imref_cl,nb_classe, 250);
imx_kmean_p(imreca, imreca_cl,nb_classe, 250);
#else
printf("imx_kmean_p not defined in normalise_kmeans_p (normalisation_topo.c) due to medipy compatibility!\n");
return;
#endif

nb=0;
/*remplissage de l'histogramme joint*/
for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
	if((imreca->mri[i][j][k]>0) || (imref->mri[i][j][k]>0))
		{
		histojoint[imreca_cl->mri[i][j][k]][imref_cl->mri[i][j][k]]++;
		nb++;
		}	

/* proba marginal*/
for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	p_reca[i]=p_reca[i]+histojoint[i][j];

for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	p_ref[i]=p_ref[i]+histojoint[j][i];

nb_reca=nb_ref=0;	
/* moyenne des classes */
 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imreca->mri[i][j][k]>0)
			{
			mean_reca[imreca_cl->mri[i][j][k]]=mean_reca[imreca_cl->mri[i][j][k]]+imreca->mri[i][j][k];
			if (imreca_cl->mri[i][j][k]==0) nb_reca++;
			}

 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imref->mri[i][j][k]>0)
			{
			mean_ref[imref_cl->mri[i][j][k]]=mean_ref[imref_cl->mri[i][j][k]]+imref->mri[i][j][k];
			if (imref_cl->mri[i][j][k]==0) nb_ref++;
			}


if (nb_reca>0) mean_reca[0]=1.0*mean_reca[0]*imreca->rcoeff/nb_reca;
if (nb_ref>0) mean_ref[0]=1.0*mean_ref[0]*imref->rcoeff/nb_ref;

for (i=1;i<nb_classe;i++)
	{
	if (p_reca[i]!=0) mean_reca[i]=1.0*mean_reca[i]*imreca->rcoeff/p_reca[i];
	if (p_ref[i]!=0) mean_ref[i]=1.0*mean_ref[i]*imref->rcoeff/p_ref[i];
	}
	

/*ecart-type */
 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imreca->mri[i][j][k]>0)
			{
			std_reca[imreca_cl->mri[i][j][k]]=std_reca[imreca_cl->mri[i][j][k]]+(imreca->mri[i][j][k]-mean_reca[imreca_cl->mri[i][j][k]])*(imreca->mri[i][j][k]-mean_reca[imreca_cl->mri[i][j][k]]);
			}

for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imref->mri[i][j][k]>0)
			{
			std_ref[imref_cl->mri[i][j][k]]=std_ref[imref_cl->mri[i][j][k]]+(imref->mri[i][j][k]-mean_ref[imref_cl->mri[i][j][k]])*(imref->mri[i][j][k]-mean_ref[imref_cl->mri[i][j][k]]);
		}

if (nb_reca>0) std_reca[0]=sqrt(1.0*std_reca[0]*imreca->rcoeff*imreca->rcoeff/nb_reca);
if (nb_ref>0) std_ref[0]=sqrt(1.0*std_ref[0]*imref->rcoeff*imref->rcoeff/nb_ref);

for (i=1;i<nb_classe;i++)
	{
	if (p_reca[i]!=0) std_reca[i]=sqrt(1.0*std_reca[i]*imreca->rcoeff*imreca->rcoeff/p_reca[i]);
	if (p_ref[i]!=0) std_ref[i]=sqrt(1.0*std_ref[i]*imref->rcoeff*imref->rcoeff/p_ref[i]);
	}

/*remplissage de imres*/
/* for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
			if(imreca->mri[i][j][k]>0)
				{tmp=0;
				for (l=0;l<nb_classe;l++)
					{
					tmp=tmp+1.0*
					(
						(
							(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[imreca_cl->mri[i][j][k]])/std_reca[imreca_cl->mri[i][j][k]]*std_ref[l]
								+mean_ref[l])*histojoint[imreca_cl->mri[i][j][k]][l]/p_reca[imreca_cl->mri[i][j][k]]);
					}
				imres->mri[i][j][k]=(int)(1.0*tmp/imres->rcoeff);
				} 
				else imres->mri[i][j][k]=0;
*/
 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
			if(imreca->mri[i][j][k]>0)
				{tmp=0;
				p_tot=0;
				/*p_tot_ref=0;
				p_tot_reca=0;
				for (l=0;l<nb_classe;l++)
					{
					p_tot_ref=p_tot_ref+1.0/(sqrt(2.0*PI)*std_ref[l])
							*exp(-0.5*(imref->mri[i][j][k]*imref->rcoeff-mean_ref[l])*(imref->mri[i][j][k]*imref->rcoeff-mean_ref[l])/(std_ref[l]*std_ref[l]));
					}
				for (l=0;l<nb_classe;l++)
					{
					p_tot_reca=p_tot_reca+1.0/(sqrt(2.0*PI)*std_reca[l])
							*exp(-0.5*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])/(std_reca[l]*std_reca[l]));
					}
				*/
				/*for (l=0;l<nb_classe;l++)
					for (t=0;t<nb_classe;t++)
					if ((std_ref[l]!=0)&&(std_reca[t]!=0))
					{
					
							tutu=1.0*histojoint[t][l]/p_reca[t]
							*1.0/(sqrt(2.0*PI)*std_reca[t])
							*exp(-0.5*((imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[t])*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[t])/(std_reca[t]*std_reca[t])));
					tmp=tmp+1.0*
					(
						(
							(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[t])/std_reca[t]*std_ref[l]
								+mean_ref[l])*tutu);
					p_tot=p_tot+tutu;		
					}
				imres->mri[i][j][k]=(int)(1.0*tmp/(p_tot*imres->rcoeff));
				} 
				else imres->mri[i][j][k]=0;
				*/for (l=0;l<nb_classe;l++)
					if ((std_ref[l]!=0)&&(std_reca[l]!=0))
					{
					
							tutu=1.0/(sqrt(2.0*PI)*std_reca[l])
							*exp(-0.5*((imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])/(std_reca[l]*std_reca[l])));
					tmp=tmp+1.0*
					(
						(
							(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])/std_reca[l]*std_ref[l]
								+mean_ref[l])*tutu);
					p_tot=p_tot+tutu;		
					}
				imres->mri[i][j][k]=(int)(1.0*tmp/(p_tot*imres->rcoeff));
				} 
				else imres->mri[i][j][k]=0;
/*liberation memoire*/
free_imatrix(histojoint,nb_classe,nb_classe);
free(p_ref);free(p_reca);
free(mean_reca);free(mean_ref);free(std_reca);free(std_ref);
free_grphic3d(imref_cl);
free_grphic3d(imreca_cl);
}

/*******************************************************************************
**        normalisation_markov                                                
**                                                                                                   
*******************************************************************************/

void normalisation_markov(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	nb_classe= GET_INT("nb de classes", 10, &e);
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	 normalise_markov_p(im1,im2,imres,nb_classe);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
** 	normalise_markov_p                                  
**                                                                   
*******************************************************************************/
void normalise_markov_p(grphic3d *imreca,grphic3d *imref,grphic3d *imres,int nb_classe)
{
int **histojoint,*p_ref,*p_reca;
int i,j,k,wdth,hght,dpth,l,nb,nb_reca,nb_ref;
double *mean_reca,*mean_ref,*std_reca,*std_ref,max_reca,max_ref,tmp,/*p_tot_ref,p_tot_reca,*/p_tot,tutu;
grphic3d *imreca_cl,*imref_cl;
wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

if ((wdth!=imref->width)||(hght!=imref->height)||(dpth!=imref->depth))
	{printf("Les images ne sont pas de meme tailles !!!\n");
	return;
	}

	
/* Allocation memoire */
histojoint=alloc_imatrix(nb_classe,nb_classe);
p_ref=(int *)malloc(nb_classe*sizeof(int));
p_reca=(int *)malloc(nb_classe*sizeof(int));

mean_reca=(double *)malloc(nb_classe*sizeof(double));
mean_ref=(double *)malloc(nb_classe*sizeof(double));
std_reca=(double *)malloc(nb_classe*sizeof(double));
std_ref=(double *)malloc(nb_classe*sizeof(double));

imref_cl=cr_grphic3d(imref);
imreca_cl=cr_grphic3d(imreca);

/*initialisation*/
for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	histojoint[i][j]=0;

for (i=0;i<nb_classe;i++)
	p_ref[i]=p_reca[i]=0;

for (i=0;i<nb_classe;i++)
	mean_ref[i]=mean_reca[i]=0;

for (i=0;i<nb_classe;i++)
	std_ref[i]=std_reca[i]=0;

/* recherche des max des 2 images */
max_reca=max_ref=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (imreca->mri[i][j][k]>max_reca) max_reca=imreca->mri[i][j][k];
	if (imref->mri[i][j][k]>max_ref) max_ref=imref->mri[i][j][k];
	}

max_reca=max_reca*imreca->rcoeff;
max_ref=max_ref*imref->rcoeff;

if (max_reca>max_ref)
	imx_copie_param_3d_p(imreca,imres);
else
	imx_copie_param_3d_p(imref,imres);

/* kmeans sur les deux images */
#ifndef COMPILE_FOR_MEDIPY
imx_segm_markov_p(imref, imref_cl,nb_classe,5, 1, 1, 0,1);
imx_segm_markov_p(imreca, imreca_cl,nb_classe,5, 1, 1, 0,1);
#else
printf("imx_segm_markov_p not defined in normalise_markov_p (normalisation_topo.c) due to medipy compatibility!\n");
return;
#endif


for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
	{
	if (imreca_cl->mri[i][j][k]>0)
	imreca_cl->mri[i][j][k]=imreca_cl->mri[i][j][k]-1;
	if (imref_cl->mri[i][j][k]>0)
	imref_cl->mri[i][j][k]=imref_cl->mri[i][j][k]-1;
	}
	
nb=0;
/*remplissage de l'histogramme joint*/
for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
	if((imreca->mri[i][j][k]>0) || (imref->mri[i][j][k]>0))
		{
		histojoint[imreca_cl->mri[i][j][k]][imref_cl->mri[i][j][k]]++;
		nb++;
		}	

/* proba marginal*/
for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	p_reca[i]=p_reca[i]+histojoint[i][j];

for (i=0;i<nb_classe;i++)
for (j=0;j<nb_classe;j++)
	p_ref[i]=p_ref[i]+histojoint[j][i];

nb_reca=nb_ref=0;	
/* moyenne des classes */
 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imreca->mri[i][j][k]>0)
			{
			mean_reca[imreca_cl->mri[i][j][k]]=mean_reca[imreca_cl->mri[i][j][k]]+imreca->mri[i][j][k];
			if (imreca_cl->mri[i][j][k]==0) nb_reca++;
			}

 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imref->mri[i][j][k]>0)
			{
			mean_ref[imref_cl->mri[i][j][k]]=mean_ref[imref_cl->mri[i][j][k]]+imref->mri[i][j][k];
			if (imref_cl->mri[i][j][k]==0) nb_ref++;
			}


if (nb_reca>0) mean_reca[0]=1.0*mean_reca[0]*imreca->rcoeff/nb_reca;
if (nb_ref>0) mean_ref[0]=1.0*mean_ref[0]*imref->rcoeff/nb_ref;

for (i=1;i<nb_classe;i++)
	{
	if (p_reca[i]!=0) mean_reca[i]=1.0*mean_reca[i]*imreca->rcoeff/p_reca[i];
	if (p_ref[i]!=0) mean_ref[i]=1.0*mean_ref[i]*imref->rcoeff/p_ref[i];
	}
	

/*ecart-type */
 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imreca->mri[i][j][k]>0)
			{
			std_reca[imreca_cl->mri[i][j][k]]=std_reca[imreca_cl->mri[i][j][k]]+imreca->mri[i][j][k]*imreca->mri[i][j][k];
			}

for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
		if(imref->mri[i][j][k]>0)
			{
			std_ref[imref_cl->mri[i][j][k]]=std_ref[imref_cl->mri[i][j][k]]+imref->mri[i][j][k]*imref->mri[i][j][k];
			}

if (nb_reca>0) std_reca[0]=sqrt(1.0*std_reca[0]*imreca->rcoeff*imreca->rcoeff/nb_reca-mean_reca[0]*mean_reca[0]);
if (nb_ref>0) std_ref[0]=sqrt(1.0*std_ref[0]*imref->rcoeff*imref->rcoeff/nb_ref-mean_ref[0]*mean_ref[0]);

for (i=1;i<nb_classe;i++)
	{
	if (p_reca[i]!=0) std_reca[i]=sqrt(1.0*std_reca[i]*imreca->rcoeff*imreca->rcoeff/p_reca[i]-mean_reca[i]*mean_reca[i]);
	if (p_ref[i]!=0) std_ref[i]=sqrt(1.0*std_ref[i]*imref->rcoeff*imref->rcoeff/p_ref[i]-mean_ref[i]*mean_ref[i]);
	}

/*remplissage de imres*/
/* for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
			if(imreca->mri[i][j][k]>0)
				{tmp=0;
				for (l=0;l<nb_classe;l++)
					{
					tmp=tmp+1.0*
					(
						(
							(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[imreca_cl->mri[i][j][k]])/std_reca[imreca_cl->mri[i][j][k]]*std_ref[l]
								+mean_ref[l])*histojoint[imreca_cl->mri[i][j][k]][l]/p_reca[imreca_cl->mri[i][j][k]]);
					}
				imres->mri[i][j][k]=(int)(1.0*tmp/imres->rcoeff);
				} 
				else imres->mri[i][j][k]=0;
*/
 for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
		for(k=0;k<dpth;k++)
			if(imreca->mri[i][j][k]>0)
				{tmp=0;
				p_tot=0;
				/*p_tot_ref=0;
				p_tot_reca=0;
				for (l=0;l<nb_classe;l++)
					{
					p_tot_ref=p_tot_ref+1.0/(sqrt(2.0*PI)*std_ref[l])
							*exp(-0.5*(imref->mri[i][j][k]*imref->rcoeff-mean_ref[l])*(imref->mri[i][j][k]*imref->rcoeff-mean_ref[l])/(std_ref[l]*std_ref[l]));
					}
				for (l=0;l<nb_classe;l++)
					{
					p_tot_reca=p_tot_reca+1.0/(sqrt(2.0*PI)*std_reca[l])
							*exp(-0.5*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])/(std_reca[l]*std_reca[l]));
					}
				*/
				/*for (l=0;l<nb_classe;l++)
					for (t=0;t<nb_classe;t++)
					if ((std_ref[l]!=0)&&(std_reca[t]!=0))
					{
					
							tutu=1.0*histojoint[t][l]/p_reca[t]
							*1.0/(sqrt(2.0*PI)*std_reca[t])
							*exp(-0.5*((imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[t])*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[t])/(std_reca[t]*std_reca[t])));
					tmp=tmp+1.0*
					(
						(
							(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[t])/std_reca[t]*std_ref[l]
								+mean_ref[l])*tutu);
					p_tot=p_tot+tutu;		
					}
				imres->mri[i][j][k]=(int)(1.0*tmp/(p_tot*imres->rcoeff));
				} 
				else imres->mri[i][j][k]=0;
				*/for (l=0;l<nb_classe;l++)
					if ((std_ref[l]!=0)&&(std_reca[l]!=0))
					{
					
							tutu=1.0/(sqrt(2.0*PI)*std_reca[l])
							*exp(-0.5*((imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])*(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])/(std_reca[l]*std_reca[l])));
					tmp=tmp+1.0*
					(
						(
							(imreca->mri[i][j][k]*imreca->rcoeff-mean_reca[l])/std_reca[l]*std_ref[l]
								+mean_ref[l])*tutu);
					p_tot=p_tot+tutu;		
					}
				imres->mri[i][j][k]=(int)(1.0*tmp/(p_tot*imres->rcoeff));
				} 
				else imres->mri[i][j][k]=0;
/*liberation memoire*/
free_imatrix(histojoint,nb_classe,nb_classe);
free(p_ref);free(p_reca);
free(mean_reca);free(mean_ref);free(std_reca);free(std_ref);
free_grphic3d(imref_cl);
free_grphic3d(imreca_cl);
}


/************************** kmeans_histo_joint_3d() **************************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void kmeans_histo_joint_3d(void)
{
int im_1,im_res,nb_classe,e;
grphic *im1,*imres;
Param_Gaussienne_2d *param;

im_1 = GET_PLACE("histo joint");
im_res=GET_PLACE("resultat");
nb_classe= GET_INT("nb de classes", 4, &e);

im1=ptr_img(im_1);
imres=ptr_img(im_res);

param=(Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));

kmeans_histo_joint_3d_p(im1,imres,nb_classe);
mean_std_classe(im1,imres,nb_classe,param);

for (e=0;e<nb_classe;e++)
printf ("%d   mx:%f  my:%f    sx:%f  sy:%f  sxy:%f    l%f\n",e,param[e].mx,param[e].my,param[e].sx,param[e].sy,param[e].sxy,param[e].l);


for (e=0;e<nb_classe;e++)
plot_iso_gaussienne_2d(im1,&param[e],1, im1->max_pixel);
	
imx_iniparaimg(im_res);
show_picture(im_res);
show_picture(im_1);
free(param);
}


/************************** kmeans_histo_joint_3d_p() **************************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void kmeans_histo_joint_3d_p(grphic *im1,grphic *imres,int nb_classe)
{
int wdth,hght,i,j,l,cl=0,condition,*nb,nb_iter/*,mini_x,mini_y,maxi_x,maxi_y*/,Ntot,ntmp,ideb,ifin;
double *mx,*my,*mxold,*myold,tmp,min;
wdth=im1->width;
hght=im1->height;

/*allocation memoire*/
mx=(double *)malloc(nb_classe*sizeof(double));
my=(double *)malloc(nb_classe*sizeof(double));
mxold=(double *)malloc(nb_classe*sizeof(double));
myold=(double *)malloc(nb_classe*sizeof(double));
nb=(int *)malloc(nb_classe*sizeof(int));

/*initialisation des moyennes*/

/*mini_x=wdth;
mini_y=hght;
maxi_x=0;
maxi_y=0;

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		if (im1->mri[i][j]>0)
			{
			if (i<mini_x) mini_x=i;
			if (j<mini_y) mini_y=j;
			if (i>maxi_x) maxi_x=i;
			if (j>maxi_y) maxi_y=j;
			
			}
for (i=0; i<nb_classe; i++)
	{
	mx[i]=mini_x+1.0*(i+0.5)*(maxi_x-mini_x)/(nb_classe+1);
	my[i]=mini_y+1.0*(i+0.5)*(maxi_y-mini_y)/(nb_classe+1);
	}
*/

Ntot=0;
for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		Ntot+=im1->mri[i][j];
		

ideb=0;
l=0;


while (l<nb_classe)
{	
	i=ideb;
	ntmp=0;
	while ((ntmp<(Ntot/nb_classe))&&(i<wdth))
	{	
	for (j=0;j<hght;j++)
		ntmp+=im1->mri[i][j];

	i++;
	}

	ifin=i-1;
	mx[l]=0;
	my[l]=0;

for (i=ideb;i<ifin;i++)
for (j=0;j<hght;j++)
if (im1->mri[i][j]>0)
	{
	mx[l]+=i*im1->mri[i][j];
	my[l]+=j*im1->mri[i][j];
	}
mx[l]=mx[l]/ntmp;
my[l]=my[l]/ntmp;

l++;
ideb=ifin;
}


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im1->mri[i][j]>0)
	{min=HUGE_VAL;
	for (l=0;l<nb_classe;l++)
		{tmp=(i-mx[l])*(i-mx[l])+(j-my[l])*(j-my[l]);
		if (tmp<min) {min=tmp; cl=l;}
		}
		imres->mri[i][j]=cl;
	}

condition=1;
nb_iter=0;
while (condition>0)
{
/* recopie dans mxold et myold*/
for (l=0;l<nb_classe;l++)
	{
	mxold[l]=mx[l];
	myold[l]=my[l];
	}
	
/* mise a jour des moyennes */
for (l=0;l<nb_classe;l++)
	{mx[l]=my[l]=nb[l]=0;}

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		if (im1->mri[i][j]>0)
			{
			mx[imres->mri[i][j]]=mx[imres->mri[i][j]]+i*im1->mri[i][j];
			my[imres->mri[i][j]]=my[imres->mri[i][j]]+j*im1->mri[i][j];
			nb[imres->mri[i][j]]=nb[imres->mri[i][j]]+im1->mri[i][j];
			}

for (l=0;l<nb_classe;l++)
	{
	mx[l]=1.0*mx[l]/nb[l];
	my[l]=1.0*my[l]/nb[l];	
	}	

/* mise a jour de imres */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im1->mri[i][j]>0)
	{min=HUGE_VAL;
	for (l=0;l<nb_classe;l++)
		{tmp=(i-mx[l])*(i-mx[l])+(j-my[l])*(j-my[l]);
		if (tmp<min) {min=tmp; cl=l;}
		}
		imres->mri[i][j]=cl;
	}

nb_iter++;

condition=0;
for (l=0;l<nb_classe;l++)
	condition=condition+(mx[l]-mxold[l])*(mx[l]-mxold[l])+(my[l]-myold[l])*(my[l]-myold[l]);

if (nb_iter>300)
	condition=0;	
}
printf("nombre d'iteration jusqu'a convergence : %d \n",nb_iter);
free(mx);free(my);free(mxold);free(myold);free(nb);
}


/************************** kmeans_L1_histo_joint_3d_p() **************************/
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

void kmeans_L1_histo_joint_3d_p(grphic *im1,grphic *imres,int nb_classe)
{
int wdth,hght,i,j,l,cl=0,condition,*nb,nb_iter,mini_x,mini_y,maxi_x,maxi_y;
double *mx,*my,*mxold,*myold,tmp,min;
wdth=im1->width;
hght=im1->height;

/*allocation memoire*/
mx=(double *)malloc(nb_classe*sizeof(double));
my=(double *)malloc(nb_classe*sizeof(double));
mxold=(double *)malloc(nb_classe*sizeof(double));
myold=(double *)malloc(nb_classe*sizeof(double));
nb=(int *)malloc(nb_classe*sizeof(int));


mini_x=wdth;
mini_y=hght;
maxi_x=0;
maxi_y=0;

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		if (im1->mri[i][j]>0)
			{
			if (i<mini_x) mini_x=i;
			if (j<mini_y) mini_y=j;
			if (i>maxi_x) maxi_x=i;
			if (j>maxi_y) maxi_y=j;
			
			}
/*initialisation des moyennes*/
for (i=0; i<nb_classe; i++)
	{
	mx[i]=mini_x+1.0*(i+0.5)*(maxi_x-mini_x)/nb_classe;
	my[i]=mini_y+1.0*(i+0.5)*(maxi_y-mini_y)/nb_classe;;
	}

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im1->mri[i][j]>0)
	{min=HUGE_VAL;
	for (l=0;l<nb_classe;l++)
		{tmp=fabs(i-mx[l])+fabs(j-my[l]);
		if (tmp<min) {min=tmp; cl=l;}
		}
		imres->mri[i][j]=cl;
	}

condition=1;
nb_iter=0;
while (condition>0)
{
/* recopie dans mxold et myold*/
for (l=0;l<nb_classe;l++)
	{
	mxold[l]=mx[l];
	myold[l]=my[l];
	}
	
/* mise a jour des moyennes */
for (l=0;l<nb_classe;l++)
	{mx[l]=my[l]=nb[l]=0;}

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		if (im1->mri[i][j]>0)
			{
			mx[imres->mri[i][j]]=mx[imres->mri[i][j]]+i*im1->mri[i][j];
			my[imres->mri[i][j]]=my[imres->mri[i][j]]+j*im1->mri[i][j];
			nb[imres->mri[i][j]]=nb[imres->mri[i][j]]+im1->mri[i][j];
			}

for (l=0;l<nb_classe;l++)
	{
	mx[l]=1.0*mx[l]/nb[l];
	my[l]=1.0*my[l]/nb[l];	
	}	

/* mise a jour de imres */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
if (im1->mri[i][j]>0)
	{min=HUGE_VAL;
	for (l=0;l<nb_classe;l++)
		{tmp=fabs(i-mx[l])+fabs(j-my[l]);
		if (tmp<min) {min=tmp; cl=l;}
		}
		imres->mri[i][j]=cl;
	}

nb_iter++;

condition=0;
for (l=0;l<nb_classe;l++)
	condition=condition+(mx[l]-mxold[l])*(mx[l]-mxold[l])+(my[l]-myold[l])*(my[l]-myold[l]);

if (nb_iter>300)
	condition=0;	
}
printf("nombre d'iteration jusqu'a convergence : %d \n",nb_iter);
free(mx);free(my);free(mxold);free(myold);free(nb);
}


/*******************************************************************************
**        normalisation_kmeans_histo_joint                                                
**                                                                                                   
*******************************************************************************/

void normalisation_kmeans_histo_joint(void)
{
	int im_1,im_2,im_res,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);

	nb_classe= GET_INT("nb de classes", 10, &e);
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	normalise_kmeans_histo_joint_p(im1,im2,imres,nb_classe);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
**        normalisation_kmeans_histo_joint_p                                                
**                                                                                                   
*******************************************************************************/

void normalise_kmeans_histo_joint_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
grphic *histo,*histo_cl;
Param_Gaussienne_2d *param;

histo=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
histo_cl=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));

histo_joint_3d_p(im1,im2,histo);

kmeans_histo_joint_3d_p(histo,histo_cl,nb_classe);

mean_std_classe(histo,histo_cl,nb_classe, param);

apply_normalisation_gaussienne_histo_joint(im1,im2,imres,histo,nb_classe,param,-1.0,-1.0,1.0);

free_grphic(histo);free_grphic(histo_cl);
free(param);
}


/*******************************************************************************
** 	void normalisation_local_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)

**
**  Normalisation locale par moyenne avec continuite ds intensite
*******************************************************************************/

void normalisation_local_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres)
{
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
	int bx0,by0,bz0,width,height,depth,dx,dy,dz;
	int resol,nb_param,topD,i,j,k,topi,topj,topk,nb,l;
	double ***m1,***m2,***dm;
	double moyenne,rcoeff1,rcoeff2,rcoeff,max,min;
	double a,b,c,d,e,f,g,h;
	double x,y,z;

	nb_param=3*BASE3D.nb_func;
	resol=BASE3D.resol;
	x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	dx=(x1[0]-x0[0])/2;
	dy=(y1[0]-y0[0])/2;
	dz=(z1[0]-z0[0])/2;
		
	width=im1->width;
	height=im1->height;
	depth=im1->depth;
	
	topD=TOP_D(nb_param);
	m1=alloc_dmatrix_3d(topD,topD,topD);
	m2=alloc_dmatrix_3d(topD,topD,topD);
	dm=alloc_dmatrix_3d(topD+1,topD+1,topD+1);
	
	rcoeff1=im1->rcoeff;
	rcoeff2=im2->rcoeff;
	imx_copie_param_3d_p(im1,imres);
	
	/* Calcul de m1 et m2 contenant les moyennes des images sur chaque sous boite i,j,k */
	/* Ce calcul peut etre optimise en divisant l'image en sous sous bloc ...*/
	for (topi=0;topi<topD;topi++)
		for (topj=0;topj<topD;topj++)
			for (topk=0;topk<topD;topk++)
			{	
			l=1.0*TOP_conv_ind(topi,topj,topk,nb_param)/3;		
			moyenne=0;nb=0;
			for (i=x0[l];i<x1[l];i++)
				for (j=y0[l];j<y1[l];j++)
					for (k=z0[l];k<z1[l];k++)
						if (im1->mri[i][j][k]!=0)
							{moyenne=moyenne+im1->mri[i][j][k];
								nb++;
							}
					if (nb>0)
					m1[topi][topj][topk]=1.0*rcoeff1*moyenne/nb;
					else
					m1[topi][topj][topk]=0.0;
			}
	
	for (topi=0;topi<topD;topi++)
		for (topj=0;topj<topD;topj++)
			for (topk=0;topk<topD;topk++)
			{
			l=1.0*TOP_conv_ind(topi,topj,topk,nb_param)/3;	
			moyenne=0;nb=0;
			for (i=x0[l];i<x1[l];i++)
				for (j=y0[l];j<y1[l];j++)
					for (k=z0[l];k<z1[l];k++)
						if (im2->mri[i][j][k]!=0)
							{moyenne=moyenne+im2->mri[i][j][k];
								nb++;
							}
					if (nb>0)
					m2[topi][topj][topk]=1.0*rcoeff2*moyenne/nb;
					else
					m2[topi][topj][topk]=0.0;
			}
			
	for (topi=0;topi<=topD;topi++)
		for (topj=0;topj<=topD;topj++)
			for (topk=0;topk<=topD;topk++)
				dm[topi][topj][topk]=0.0;
	
	for (topi=0;topi<topD;topi++)
		for (topj=0;topj<topD;topj++)
			for (topk=0;topk<topD;topk++)
				dm[topi][topj][topk]=m2[topi][topj][topk]-m1[topi][topj][topk];
				


/*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
min=HUGE_VAL;
max=-HUGE_VAL;

for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
		for (k=0;k<topD;k++)
			{
			if (dm[i][j][k]>max) max=dm[i][j][k];
			if (dm[i][j][k]<min) min=dm[i][j][k];
			}

max=im1->max_pixel*rcoeff1+max;
min=im1->min_pixel*rcoeff1+min;

imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;



	/* Calcul de Imres */
		for (topi=0;topi<=topD;topi++)
			for (topj=0;topj<=topD;topj++)
				for (topk=0;topk<=topD;topk++)
					{
					/* calcul des coefficients du polynome */
					
					/* a */
					if ((topi>0)&&(topj>0)&&(topk>0))
						a=dm[topi-1][topj-1][topk-1];
					else a=0.0;
					
					/* b */
					if ((topi>0)&&(topj>0)&&(topk>0))
						b=2.0*(dm[topi][topj-1][topk-1]-dm[topi-1][topj-1][topk-1]);
					else
						if ((topj>0)&&(topk>0))
							b=2.0*dm[topi][topj-1][topk-1];
						else b=0.0;
						
					
					/* c */
					if ((topi>0)&&(topj>0)&&(topk>0))
						c=2.0*(dm[topi-1][topj][topk-1]-dm[topi-1][topj-1][topk-1]);
					
					else if ((topi>0)&&(topk>0))
							c=2.0*dm[topi-1][topj][topk-1];
					
					else c=0.0;
						
					
					/* d */
					if ((topi>0)&&(topj>0)&&(topk>0))
						d=2.0*(dm[topi-1][topj-1][topk]-dm[topi-1][topj-1][topk-1]);
					
					else if ((topi>0)&&(topj>0))
						d=2.0*dm[topi-1][topj-1][topk];
					
					else d=0.0;


					/* e */
					if ((topi>0)&&(topj>0)&&(topk>0))
						e=4.0*(dm[topi][topj][topk-1]+dm[topi-1][topj-1][topk-1]-dm[topi][topj-1][topk-1]-dm[topi-1][topj][topk-1]);
					
					else if (topk>0)
							{
								e= 4.0*(dm[topi][topj][topk-1]);	
								if (topi>0)
								e=4.0*(dm[topi][topj][topk-1]-dm[topi-1][topj][topk-1]);
								if (topj>0)
								e=4.0*(dm[topi][topj][topk-1]-dm[topi][topj-1][topk-1]);
							}
					else e=0.0;					
					


					/* f */
					if ((topi>0)&&(topj>0)&&(topk>0))
						f=4.0*(dm[topi][topj-1][topk]+dm[topi-1][topj-1][topk-1]-dm[topi][topj-1][topk-1]-dm[topi-1][topj-1][topk]);
					else if (topj>0)
							{
								f= 4.0*(dm[topi][topj-1][topk]);	
								if (topi>0)
								f=4.0*(dm[topi][topj-1][topk]-dm[topi-1][topj-1][topk]);
								if (topk>0)
								f=4.0*(dm[topi][topj-1][topk]-dm[topi][topj-1][topk-1]);
								
							}
					else f=0.0;					
					

					/* g */
					if ((topi>0)&&(topj>0)&&(topk>0))
						g=4.0*(dm[topi-1][topj][topk]+dm[topi-1][topj-1][topk-1]-dm[topi-1][topj][topk-1]-dm[topi-1][topj-1][topk]);
					
					else if (topi>0)
							{
								g=4.0*(dm[topi-1][topj][topk]);	
								if (topj>0)
								g=4.0*(dm[topi-1][topj][topk]-dm[topi-1][topj-1][topk]);
								if (topk>0)
								g=4.0*(dm[topi-1][topj][topk]-dm[topi-1][topj][topk-1]);
								
							}
					else g=0.0;					
					
					
					/* h */
					if ((topi>0)&&(topj>0)&&(topk>0))
						h=8.0*(dm[topi][topj][topk]+dm[topi][topj-1][topk-1]+dm[topi-1][topj][topk-1]+dm[topi-1][topj-1][topk]
								-dm[topi][topj][topk-1]-dm[topi][topj-1][topk]-dm[topi-1][topj][topk]-dm[topi-1][topj-1][topk-1]);
					else if ((topi>0)&&(topj>0))
						h=8.0*(dm[topi][topj][topk]+dm[topi-1][topj-1][topk]-dm[topi][topj-1][topk]-dm[topi-1][topj][topk]);
								
					else if ((topi>0)&&(topk>0))
						h=8.0*(dm[topi][topj][topk]+dm[topi-1][topj][topk-1]-dm[topi][topj][topk-1]-dm[topi-1][topj][topk]);
						
					else if ((topj>0)&&(topk>0))
						h=8.0*(dm[topi][topj][topk]+dm[topi][topj-1][topk-1]-dm[topi][topj][topk-1]-dm[topi][topj-1][topk]);
					
					else if (topi>0)
						h=8.0*(dm[topi][topj][topk]-dm[topi-1][topj][topk]);	
					
					else if (topj>0)
						h=8.0*(dm[topi][topj][topk]-dm[topi][topj-1][topk]);
									
					else if (topk>0)
						h=8.0*(dm[topi][topj][topk]-dm[topi][topj][topk-1]);	
						
					else h=8.0*dm[topi][topj][topk];		
					
					/* remplissage de imres */
					l=1.0*TOP_conv_ind(MINI(topi,topD-1),MINI(topj,topD-1),MINI(topk,topD-1),nb_param)/3;		
					for (i=0;i<dx;i++)
						for (j=0;j<dy;j++)
							for (k=0;k<dz;k++)
								{
								if (topi<topD) bx0=x0[l]+i; 
								else {bx0=x0[l]+dx+i;}
								
								if (topj<topD)
								by0=y0[l]+j;
								else {by0=y0[l]+dy+j;}
								
								if (topk<topD)
								bz0=z0[l]+k;
								else {bz0=z0[l]+dz+k;}
								
								x=1.0*i/(2.0*dx);
								y=1.0*j/(2.0*dy);
								z=1.0*k/(2.0*dz);
								
								if (im1->mri[bx0][by0][bz0]!=0)
								imres->mri[bx0][by0][bz0]=floor(0.5+1.0*(im1->mri[bx0][by0][bz0]*rcoeff1+a+1.0*b*x+1.0*c*y+1.0*d*z
																								+1.0*e*x*y+1.0*f*x*z+1.0*g*y*z+1.0*h*x*y*z)/rcoeff);

								}
										

		}

				
imx_norm_rcoeff_3d_p(im2,imres);
				
	
free_dmatrix_3d(m1);free_dmatrix_3d(m2);free_dmatrix_3d(dm);
}




/*******************************************************************************
********************************************************************************
********************************************************************************
***************   normalisation avec segmentation sous-jacente  **************** 
********************************************************************************
********************************************************************************
********************************************************************************/


/*******************************************************************************
**        normalisation_segmentation
**
*******************************************************************************/

void  normalisation_segmentation(void)
{
	int im_1,im_2,im_res,im_histo,e;
	int nb_classe;
	grphic3d *im1,*im2,*imres;
	grphic *histo;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	im_histo=GET_PLACE("histo joint");

	nb_classe= GET_INT("nb de classes", 4, &e);
	
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	histo=ptr_img(im_histo);
	
			if ((im1->mask==NULL)||(im1->mask==NULL))
			printf ("\n\n Attention, il n'y a pas de masque du cerveau associ�aux images !!! Cela peut d��iorer la qualit�de la normalisation d'intensit�n\n");
		
	
	normalisation_segmentation_p(im1,im2,imres,histo,nb_classe);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
**         imx_normalisation_segmentation_p                                                
**                                                                                                   
*******************************************************************************/

void  imx_normalisation_segmentation_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
grphic *histo;
histo=cr_grphic_modif(TOPO_SIZE_NORM_HISTO,TOPO_SIZE_NORM_HISTO,0,1,0);
histo->pos=0;
normalisation_segmentation_p(im1,im2,imres,histo, nb_classe);
free_grphic(histo);
}

/*******************************************************************************
**         normalisation_segmentation_p
**
*******************************************************************************/

void  normalisation_segmentation_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,grphic *histo,int nb_classe)
{
Param_Gaussienne_2d *param;
grphic3d *imref_cl,*imreca_cl;
grphic3d *cer1,*cer2;

int i,j,k;
param=(Param_Gaussienne_2d *) malloc(nb_classe*sizeof(Param_Gaussienne_2d));

//histo_joint_linear_3d_p(im1,im2,histo);
//imx_inimaxminpixel_p(histo);

//segmentation_histo(histo,param,nb_classe);

imref_cl=cr_grphic3d(im1);
imreca_cl=cr_grphic3d(im2);

cer1=cr_grphic3d(im1);
cer2=cr_grphic3d(im2);

imx_copie_3d_p(im1,cer1);
imx_copie_3d_p(im2,cer2);


for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	if (im1->mask->mri[i][j][k] == 0)
		cer1->mri[i][j][k]=0;
		
for (i=0;i<im2->width;i++)
for (j=0;j<im2->height;j++)
for (k=0;k<im2->depth;k++)
	if (im2->mask->mri[i][j][k] == 0)
		cer2->mri[i][j][k]=0;
		
#ifndef COMPILE_FOR_MEDIPY
imx_segm_markov_p(cer1, imref_cl,nb_classe,5, 1, 1, 0,1);
imx_segm_markov_p(cer2, imreca_cl,nb_classe,5, 1, 1, 0,1);
#else
printf("imx_segm_markov_p not defined in normalisation_segmentation_p (normalisation_topo.c) due to medipy compatibility!\n");
return;
#endif

for (i=0;i<im1->width;i++)
for (j=0;j<im1->height;j++)
for (k=0;k<im1->depth;k++)
	if (imref_cl->mri[i][j][k]==imreca_cl->mri[i][j][k])
		{
		//imref_cl->mri[i][j][k]=im1->mri[i][j][k];
		//imreca_cl->mri[i][j][k]=im2->mri[i][j][k];		
		imref_cl->mri[i][j][k]=cer1->mri[i][j][k];
		imreca_cl->mri[i][j][k]=cer2->mri[i][j][k];		
		
		}
	else
		{
		imref_cl->mri[i][j][k]=0;
		imreca_cl->mri[i][j][k]=0;		
		}


histo_joint_linear_3d_p(imref_cl,imreca_cl,histo);
imx_inimaxminpixel_p(histo);


fit_gaussienne2d_EM(histo, param, nb_classe);

apply_normalisation_gaussienne_histo_joint(im1,im2,imres,histo,nb_classe,param,-1.0,-1.0,1.0);

imx_inimaxminpixel_3d_p(imres);
free_grphic3d(imref_cl);
free_grphic3d(imreca_cl);
free_grphic3d(cer1);
free_grphic3d(cer2);

free(param);
}




/*******************************************************************************
**         normalisation_segmentation_p
**
*******************************************************************************/

void  segmentation_histo(grphic *histo, Param_Gaussienne_2d *param, int nb_classe)
{
int wdth,hght,i,j,l,m,lmax=0,mmax=0,stop,Ntot;
Param_Gaussienne_2d *param_x,*param_y;
grphic *histo_x,*histo_y;
double **proba_cond_x,**proba_cond_y,**proba_classe;
double tot,tmp,max;
float *x,*y;

wdth=histo->width;
hght=histo->height;


/* allocation m�oire */
param_x= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
param_y= (Param_Gaussienne_2d *)malloc(nb_classe*sizeof(Param_Gaussienne_2d));
histo_x=cr_grphic_modif(wdth,1,0,1,0);
histo_y=cr_grphic_modif(1,hght,0,1,0);

proba_cond_x=alloc_dmatrix(wdth,nb_classe);
proba_cond_y=alloc_dmatrix(hght,nb_classe);
proba_classe=alloc_dmatrix(nb_classe,nb_classe);

x=(float *)malloc(wdth*sizeof(float));
y=(float *)malloc(wdth*sizeof(float));

//do
//{
stop=0;
for (i=0;i<wdth;i++)
	histo_x->mri[i][0]=0;

for (j=0;j<hght;j++)
	histo_y->mri[0][j]=0;

Ntot=0;
for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	{
	histo_x->mri[i][0]+=histo->mri[i][j];
	histo_y->mri[0][j]+=histo->mri[i][j];
	Ntot+=histo->mri[i][j];
	}

fit_gaussienne2d_EM(histo_x,param_x, nb_classe);
fit_gaussienne2d_EM(histo_y,param_y, nb_classe);


for (i=0;i<wdth;i++)
	{x[i]=i;
	y[i]=histo_x->mri[i][0];
	}

#ifndef COMPILE_FOR_MEDIPY
  	plot_histo(x,y,wdth,"Histogram of the ROI","ROI");
#endif

for (i=0;i<wdth;i++)
	y[i]=0;

for (i=0;i<wdth;i++)
	for(l=0;l<nb_classe;l++)
	{
	y[i]+=Ntot*eval_gaussienne_2d(i,0,&param_x[l]);
	}

#ifndef COMPILE_FOR_MEDIPY
	plot_histo(x,y,wdth,"Histogram of the ROI","ROI");
#endif


for (l=0;l<nb_classe;l++)
	{
	tot=0;
	for (i=0;i<wdth;i++)
		{
		proba_cond_x[i][l]=eval_gaussienne_2d(i,0,&param_x[l]);
		tot+=proba_cond_x[i][l];
		}
	for (i=0;i<wdth;i++)	
		proba_cond_x[i][l]=1.0*proba_cond_x[i][l]/tot;
	
	tot=0;
	for (i=0;i<hght;i++)
		{
		proba_cond_y[i][l]=eval_gaussienne_2d(0,i,&param_y[l]);
		tot+=proba_cond_y[i][l];
		}
	for (i=0;i<hght;i++)	
		proba_cond_y[i][l]=1.0*proba_cond_y[i][l]/tot;
	}

	
	for (i=0;i<nb_classe;i++)
	for (j=0;j<nb_classe;j++)
		proba_classe[i][j]=0;
		
	
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (l=0;l<nb_classe;l++)
	for (m=0;m<nb_classe;m++)
		{
		proba_classe[l][m]+=histo->mri[i][j]*proba_cond_x[i][l]*proba_cond_y[j][m];
		}	
	
	tot=0;
	for (l=0;l<nb_classe;l++)
	for (m=0;m<nb_classe;m++)
		tot+=proba_classe[l][m];
		
	for (l=0;l<nb_classe;l++)
	for (m=0;m<nb_classe;m++)
		proba_classe[l][m]=proba_classe[l][m]/tot;
		
	
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		{
		max=0;tmp=0;tot=0;
		for (l=0;l<nb_classe;l++)
		for (m=0;m<nb_classe;m++)
			{
			tmp=proba_cond_x[i][l]*proba_cond_y[j][m];
			if (tmp>max)
				{max=tmp; lmax=l; mmax=m;}
			}
			
		histo->mri[i][j]=histo->mri[i][j]*proba_classe[lmax][mmax];
			/*tmp=proba_cond_x[i][l]*proba_cond_y[j][m];
			max+=tmp*proba_classe[l][m];
			tot+=tmp;
			}
			histo->mri[i][j]=histo->mri[i][j]*(10*max/tot);*/
		}

	/*for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		{
		tmp=0;
		for (l=0;l<nb_classe;l++)
			tmp+=Ntot*proba_cond_x[i][l]*proba_cond_y[j][l]*(double)histo->mri[i][j];
		
		histo->mri[i][j]=(int)tmp;
		}*/
			
if (histo->pos!=0)
show_picture(histo->pos);
//}
//while (stop==1);

free_dmatrix(proba_classe,nb_classe,nb_classe);
free_dmatrix(proba_cond_x,wdth,nb_classe);
free_dmatrix(proba_cond_y,hght,nb_classe);
free_grphic(histo_x);
free_grphic(histo_y);
free(param_x);
free(param_y);
free(x); free(y);
}


/*******************************************************************************
**        normalisation_equalisation_histo
**
*******************************************************************************/

void  normalisation_equalisation_histo(void)
{
	int im_1,im_2,im_res,e;
	int nb_bin;
	grphic3d *im1,*im2,*imres;

	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
	nb_bin= GET_INT("nb de bins", 1024, &e);
	
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	normalisation_equalisation_histo_p(im1,im2,imres,nb_bin);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
**         normalisation_equalisation_histo_p
**
*******************************************************************************/

void  normalisation_equalisation_histo_p(grphic3d *im1,grphic3d *im2,grphic3d *imres,int nb_classe)
{
int i,j,k,l,wdth,hght,dpth,Ntot,r;
double *tab1,*tab2;
double *diff1,*diff2;
double *sec1,*sec2;

double *t1,*t2,tmp,dg,dd;
  FILE * fic;

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);
imx_copie_param_3d_p(im2,imres);

Ntot=0;

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
 	Ntot++;

r=floor(Ntot/nb_classe)+1;
Ntot=r*nb_classe;
tab1=(double *) malloc((Ntot+2)*sizeof(double));
tab2=(double *) malloc((Ntot+2)*sizeof(double));
diff1=(double *) malloc((Ntot+2)*sizeof(double));
diff2=(double *) malloc((Ntot+2)*sizeof(double));
sec1=(double *) malloc((Ntot+2)*sizeof(double));
sec2=(double *) malloc((Ntot+2)*sizeof(double));

t1=(double *) malloc((nb_classe+1)*sizeof(double));
t2=(double *) malloc((nb_classe+1)*sizeof(double));


l=1;
tab1[0]=MAXI(im1->min_pixel,0);
tab2[0]=MAXI(im2->min_pixel,0);
t1[0]=MAXI(im1->min_pixel,0);
t2[0]=MAXI(im2->min_pixel,0);


for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if ((im1->mri[i][j][k]>0)&&(im2->mri[i][j][k]>0))
 		{
		tab1[l]=im1->mri[i][j][k];
		tab2[l]=im2->mri[i][j][k];
		l++;
		}

for (i=l;i<=Ntot;i++)
{
tab1[i]=tab1[2*l-i];
tab2[i]=tab2[2*l-i];
}

tab1[Ntot+1]=im1->max_pixel;
tab2[Ntot+1]=im2->max_pixel;
	
/*printf("debut tri 1 \n");
tri_rapide_double(tab1,0,Ntot+2);
printf("debut tri 2 \n");
tri_rapide_double(tab2,0,Ntot+2);
printf("fin tri 2 \n");*/

printf("debut tri 1 \n");
qsort(tab1,Ntot+2,sizeof(double),double_compare_function2);
printf("debut tri 2 \n");
qsort(tab2,Ntot+2,sizeof(double),double_compare_function2);
printf("fin tri 2 \n");


for (i=1;i<nb_classe;i++)
{tmp=0;
for (j=-r;j<r;j++)
	tmp+=tab1[i*r+j];
	
t1[i]=0.5*tmp/r;

tmp=0;
for (j=-r;j<r;j++)
	tmp+=tab2[i*r+j];
	
t2[i]=0.5*tmp/r;

}

t1[nb_classe]=im1->max_pixel;
t2[nb_classe]=im2->max_pixel;

/*
for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if (im1->mri[i][j][k]>0)
 	{
	l=0;
	while (im1->mri[i][j][k]>tab1[l])
		l++;
	
	imres->mri[i][j][k]=(int)(tab2[l-1]*(im1->mri[i][j][k]-tab1[l-1])+tab2[l]*(tab1[l]-im1->mri[i][j][k]))/(tab1[l]-tab1[l-1]);	
	}
	else
	imres->mri[i][j][k]=0;
*/

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if (im1->mri[i][j][k]>0)
 	{
	l=0;
	while (im1->mri[i][j][k]>t1[l])
		l++;
	
	imres->mri[i][j][k]=(int)(t2[l-1]*(im1->mri[i][j][k]-t1[l-1])+t2[l]*(t1[l]-im1->mri[i][j][k]))/(t1[l]-t1[l-1]);	
	}
	else
	imres->mri[i][j][k]=0;

/* calcul de la d�iv�*/
/*for (i=0;i<Ntot+2;i++)
	{
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=tab1[j];
		
	dg=dg/(i-MAXI(0,i-r));
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,Ntot+1);j++)
		dd+=tab1[j];
	
	dd=dd/(MINI(i+r,Ntot+1)-i);	
	
	diff1[j]=dd-dg;
	
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=tab2[j];
		
	dg=dg/(i-MAXI(0,i-r));
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,Ntot+1);j++)
		dd+=tab2[j];
	
	dd=dd/(MINI(i+r,Ntot+1)-i);	
	
	diff2[j]=dd-dg;
	}
*/
r=floor(nb_classe/10);

for (i=0;i<nb_classe;i++)
	{
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=t1[j];
		
	dg=dg/MAXI((i-MAXI(0,i-r)),1);
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,nb_classe-1);j++)
		dd+=t1[j];
	
	dd=dd/MAXI((MINI(i+r,nb_classe-1)-i),1);	
	
	diff1[i]=dd-dg;
	
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=t2[j];
		
	dg=dg/MAXI((i-MAXI(0,i-r)),1);
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,nb_classe-1);j++)
		dd+=t2[j];
	
	dd=dd/MAXI((MINI(i+r,nb_classe-1)-i),1);	
	
	diff2[i]=dd-dg;
	}

/* calcul de la d�iv�seconde */
/*for (i=0;i<Ntot+2;i++)
	{
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=diff1[j];
		
	dg=dg/(i-MAXI(0,i-r));
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,Ntot+1);j++)
		dd+=diff1[j];
	
	dd=dd/(MINI(i+r,Ntot+1)-i);	
	
	sec1[j]=dd-dg;
	
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=diff2[j];
		
	dg=dg/(i-MAXI(0,i-r));
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,Ntot+1);j++)
		dd+=diff2[j];
	
	dd=dd/(MINI(i+r,Ntot+1)-i);	
	
	sec2[j]=dd-dg;
	}
*/
for (i=0;i<nb_classe;i++)
	{
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=diff1[j];
		
	dg=dg/MAXI((i-MAXI(0,i-r)),1);
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,nb_classe-1);j++)
		dd+=diff1[i];
	
	dd=dd/MAXI((MINI(i+r,nb_classe-1)-i),1);	
	
	sec1[j]=dd-dg;
	
	if (isnan(sec1[j]))
		dg=0;
		
	dg=0;
	for (j=MAXI(0,i-r);j<i;j++)
		dg+=diff2[j];
		
	dg=dg/MAXI((i-MAXI(0,i-r)),1);
	
	dd=0;
	for (j=i+1;j<=MINI(i+r,nb_classe-1);j++)
		dd+=diff2[j];
	
	dd=dd/MAXI((MINI(i+r,nb_classe-1)-i),1);	
	
	sec2[i]=dd-dg;
	}


/* Ouverture du fichier */
fic = fopen("test_histo.dat", "w+");

/* Ecriture des donn�s dans le fichier */
for (i=0;i<nb_classe;i++)
	fprintf(fic, "%f %f %f %f %f %f\n", t1[i], t2[i],diff1[i],diff2[i],sec1[i],sec2[i]);

/* Fermeture du fichier */
fclose(fic);

imx_inimaxminpixel_3d_p(imres);
free(tab1);
free(tab2);
free(t1);
free(t2);
free(diff1);
free(diff2);
free(sec1);
free(sec2);


}




/*******************************************************************************
**        imx_norm_seuil_meanecty_robust_3d
**
*******************************************************************************/

void  imx_norm_seuil_meanecty_robust_3d(void)
{
	int im_1,im_2,im_res;
	grphic3d *im1,*im2,*imres;
		
	im_1 = GET_PLACE3D(TEXT0231);
	im_2 = GET_PLACE3D(TEXT0031);
	im_res=GET_PLACE3D(TEXT0006);
	
		
	
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	imx_norm_seuil_meanecty_robust_3d_p(im1,im2,imres,0);
	
	show_picture_3d(im_res);

}
/******************************************************************
*** imx_norm_seuil_meanecty_robust_3d_p 
**/
/*!  NORMALISATION D'UNE IMAGE (im1) PAR RAPPORT A UNE AUTRE (im2) SEUILLEE
***
***	\param im1 : image source
***	\param im2 : image ref
***	\param imres : image resultat (E/S)
**************************************************************/
void imx_norm_seuil_meanecty_robust_3d_p(grphic3d *im1, grphic3d *im2, grphic3d *imres,int nb_classe)
{
int i,j,k,l,m;
int height,width,depth;
double *data1,*data2;
double moy1,moy2,var1,var2,sygma1,sygma2;
double rcoeff1,rcoeff2,rcoeff;
double aux1,aux2;
float max,min;
int err;



m=0;
l=0;
moy1=0;
moy2=0;
var1=0;
var2=0;
height=im1->height;
width=im1->width;
depth=im1->depth;
rcoeff1=im1->rcoeff;
rcoeff2=im2->rcoeff;

/*cacul des moyennes des images en utilisant la m�iane*/	
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im1->mri[i][j][k])
      {
        l++;
      }
     }		
height=im2->height;
width=im2->width;
depth=im2->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im2->mri[i][j][k])
      {
        m++;
      }
     }


/* Vectorisation de l'image */
data1=(double *)malloc(l*sizeof(double));
data2=(double *)malloc(m*sizeof(double));

l=0;
height=im1->height;
width=im1->width;
depth=im1->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im1->mri[i][j][k])
  				 {data1[l]=im1->mri[i][j][k]; l++;
					 }
	    }

m=0;
height=im2->height;
width=im2->width;
depth=im2->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im2->mri[i][j][k])
      {
        data2[m]=im2->mri[i][j][k];
				m++;
      }
     }

qsort(data1,l,sizeof(double),double_compare_function2);
qsort(data2,m,sizeof(double),double_compare_function2);


moy1= data1[(int)l/2];
moy2= data2[(int)m/2];

l=0;
height=im1->height;
width=im1->width;
depth=im1->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im1->mri[i][j][k])
  				 {data1[l]=fabs(im1->mri[i][j][k]-moy1); l++;
					 }
	    }

m=0;
height=im2->height;
width=im2->width;
depth=im2->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im2->mri[i][j][k])
      {
        data2[m]=fabs(im2->mri[i][j][k]-moy2);
				m++;
      }
     }

qsort(data1,l,sizeof(double),double_compare_function2);
qsort(data2,m,sizeof(double),double_compare_function2);


sygma1=1.4826*rcoeff1*data1[(int)l/2];
sygma2=1.4826*rcoeff2*data2[(int)m/2];

moy1=moy1*rcoeff1;
moy2=moy2*rcoeff2;


printf("image1 : %f +/- %f \n",moy1,sygma1);
printf("image2 : %f +/- %f \n",moy2,sygma2);

	
aux1=rcoeff1*sygma2/sygma1;
aux2=moy2-moy1*sygma2/sygma1;

imx_copie_param_3d_p(im1,imres);
max=(float)(aux1*(im1->max_pixel)+aux2);
min=(float)(aux1*(im1->min_pixel)+aux2);

/*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);
rcoeff=imres->rcoeff;

  
height=im1->height;
width=im1->width;
depth=im1->depth;
for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for(k=0;k<depth;k++)
       if (im1->mri[i][j][k]!=0)
        {imres->mri[i][j][k]=(TYPEMRI3D)floor((aux1*im1->mri[i][j][k]+aux2)/rcoeff);}
       else
        {imres->mri[i][j][k]=0;}
	
/* Modif JP 2/2001 imres->min_pixel=0;*/
imres->width=width;
imres->height=height;
imres->depth=depth;

/* Mettre le meme rcoeff pour les deux images */
imx_norm_rcoeff_3d_p(im2,imres);
free(data1);
free(data2);

return;	
}

/*******************************************************************************
**        				 trimmed_robust_histo                                                
**                                                                                                   
*******************************************************************************/

void  trimmed_robust_histo(grphic *histo,grphic *trimmedhisto,double nb_sigma)
{
int i,j,wdth,hght,l,seuil,tmp;
int moy,std,*nbx,*nby;
double seuilhaut,seuilbas;
grphic* mask;

wdth=histo->width;
hght=histo->height;

mask=cr_grphic_modif(wdth,hght,0,1,0);

nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		mask->mri[i][j]=0;

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}


for (i=0; i<wdth; i++)
	{
	moy=0;
	seuil=(int)nbx[i]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<hght-1))
		{
		tmp+=histo->mri[i][moy];
		moy++;
		}
		
	
	std=0;
	tmp=histo->mri[i][moy];
	while (tmp<seuil)
		{
		std++;
		if (moy+std<wdth)
			tmp+=histo->mri[i][moy+std];
		
		if (moy-std>=0)
			tmp+=histo->mri[i][moy-std];		
		}
	
	std=std*1.4826;
	
	seuilbas=moy-1.0*nb_sigma*std;
	seuilhaut=moy+1.0*nb_sigma*std;
	
	
	
	for (j=0; j<hght; j++)
		if ((j>seuilbas)&&(j<seuilhaut))
		mask->mri[i][j]++;
	}


for (j=0; j<hght; j++)
	{
	moy=0;
	seuil=(int)nby[j]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<wdth-1))
		{
		tmp+=histo->mri[moy][j];
		moy++;
		}
		
	
	std=0;
	tmp=histo->mri[moy][j];
	while (tmp<seuil)
		{
		std++;
		if (moy+std<hght)
			tmp+=histo->mri[moy+std][j];
		
		if (moy-std>=0)
			tmp+=histo->mri[moy-std][j];		
		}
	
	std=std*1.4826;
	
	seuilbas=moy-1.0*nb_sigma*std;
	seuilhaut=moy+1.0*nb_sigma*std;
	
	
	
	for (i=0; i<wdth; i++)
		if ((i>seuilbas)&&(i<seuilhaut))
		mask->mri[i][j]++;
	}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (mask->mri[i][j]==2)
	trimmedhisto->mri[i][j]=histo->mri[i][j];
else
	trimmedhisto->mri[i][j]=0;


free(nbx);free(nby);

free_grphic(mask);
}

/************************************************************************************/

/*void  trimmed_robust_histo(grphic *histo,grphic *trimmedhisto,double nb_sigma)
{
int i,j,k,wdth,hght,l,seuil,tmp;
int moy,std,*nbx,*nby,Ntot;
double seuilhaut,seuilbas,*data;
grphic* mask;

wdth=histo->width;
hght=histo->height;

mask=cr_grphic_modif(wdth,hght,0,1,0);

nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		mask->mri[i][j]=0;

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}

Ntot=0;
for (j=0; j<hght; j++)
		{
		Ntot+=nby[j];	
		}


// Vectorisation de l'image 
data=(double *)malloc(Ntot*sizeof(double));


l=0;
for (i=0; i<wdth; i++)
	{
	moy=0;
	seuil=(int)nbx[i]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<hght-1))
		{
		tmp+=histo->mri[i][moy];
		moy++;
		}
		
	

	for (j=0; j<hght; j++)
		for (k=0;k<histo->mri[i][j];k++)
			{data[l]=fabs(moy-j); l++;}
	}

qsort(data,Ntot,sizeof(double),double_compare_function2);
std=data[(int)(Ntot/2)];
std=std*1.4826;

for (i=0; i<wdth; i++)
{moy=0;
	seuil=(int)nbx[i]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<hght-1))
		{
		tmp+=histo->mri[i][moy];
		moy++;
		}
seuilbas=moy-1.0*nb_sigma*std;
seuilhaut=moy+1.0*nb_sigma*std;
for (j=0; j<hght; j++)
	if ((j>seuilbas)&&(j<seuilhaut))
		mask->mri[i][j]++;
}	

l=0;
for (j=0; j<hght; j++)
	{
	moy=0;
	seuil=(int)nby[j]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<wdth-1))
		{
		tmp+=histo->mri[moy][j];
		moy++;
		}
		
	

	for (i=0; i<wdth; i++)
		for (k=0;k<histo->mri[i][j];k++)
			{data[l]=fabs(moy-i); l++;}
	}

qsort(data,Ntot,sizeof(double),double_compare_function2);
std=data[(int)(Ntot/2)];
std=std*1.4826;

for (j=0; j<hght; j++)
	{
	moy=0;
	seuil=(int)nby[j]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<wdth-1))
		{
		tmp+=histo->mri[moy][j];
		moy++;
		}
seuilbas=moy-1.0*nb_sigma*std;
seuilhaut=moy+1.0*nb_sigma*std;
for (i=0; i<wdth; i++)
		if ((i>seuilbas)&&(i<seuilhaut))
		mask->mri[i][j]++;
	}	

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (mask->mri[i][j]==2)
	trimmedhisto->mri[i][j]=histo->mri[i][j];
else
	trimmedhisto->mri[i][j]=0;


free(nbx);free(nby);
free(data);
free_grphic(mask);
}*/

/*
void  trimmed_robust_histo(grphic *histo,grphic *trimmedhisto,double nb_sigma)
{
int i,j,k,wdth,hght,l,seuil,tmp;
int moy,std,*nbx,*nby,Ntot;
double seuilhaut,seuilbas,*sx,*sy;
grphic* mask;

wdth=histo->width;
hght=histo->height;

mask=cr_grphic_modif(wdth,hght,0,1,0);

nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		mask->mri[i][j]=0;

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}

Ntot=0;
for (j=0; j<hght; j++)
		{
		Ntot+=nby[j];	
		}


sx=(double *)malloc(wdth*sizeof(double));
sy=(double *)malloc(hght*sizeof(double));


for (i=0; i<wdth; i++)
	{
	moy=0;
	seuil=(int)nbx[i]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<hght-1))
		{
		tmp+=histo->mri[i][moy];
		moy++;
		}
		
	

	std=0;
	tmp=histo->mri[i][moy];
	while (tmp<seuil)
		{
		std++;
		if (moy+std<wdth)
			tmp+=histo->mri[i][moy+std];
		
		if (moy-std>=0)
			tmp+=histo->mri[i][moy-std];		
		}
	
	std=std*1.4826;
	sx[i]=std*std;
	}

std=0;
for (i=0; i<wdth; i++)
	std+=sx[i]*nbx[i];
	
std=sqrt(1.0*std/Ntot);

for (i=0; i<wdth; i++)
{moy=0;
	seuil=(int)nbx[i]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<hght-1))
		{
		tmp+=histo->mri[i][moy];
		moy++;
		}
seuilbas=moy-1.0*nb_sigma*std;
seuilhaut=moy+1.0*nb_sigma*std;
for (j=0; j<hght; j++)
	if ((j>seuilbas)&&(j<seuilhaut))
		mask->mri[i][j]++;
}	


for (j=0; j<hght; j++)
	{
	moy=0;
	seuil=(int)nby[j]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<wdth-1))
		{
		tmp+=histo->mri[moy][j];
		moy++;
		}
		
	
	std=0;
	tmp=histo->mri[moy][j];
	while (tmp<seuil)
		{
		std++;
		if (moy+std<hght)
			tmp+=histo->mri[moy+std][j];
		
		if (moy-std>=0)
			tmp+=histo->mri[moy-std][j];		
		}
	
	std=std*1.4826;
	sy[j]=std*std;

	}

std=0;
for (j=0; j<hght; j++)
	std+=sy[j]*nby[j];
	
std=sqrt(1.0*std/Ntot);

for (j=0; j<hght; j++)
	{
	moy=0;
	seuil=(int)nby[j]/2;
	tmp=0;
	while ((tmp<seuil)&&(moy<wdth-1))
		{
		tmp+=histo->mri[moy][j];
		moy++;
		}
seuilbas=moy-1.0*nb_sigma*std;
seuilhaut=moy+1.0*nb_sigma*std;
for (i=0; i<wdth; i++)
		if ((i>seuilbas)&&(i<seuilhaut))
		mask->mri[i][j]++;
	}	

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (mask->mri[i][j]==2)
	trimmedhisto->mri[i][j]=histo->mri[i][j];
else
	trimmedhisto->mri[i][j]=0;


free(nbx);free(nby);
free(sx);
free(sy);

free_grphic(mask);
}*/
/*void  trimmed_robust_histo(grphic *histo,grphic *trimmedhisto,double nb_sigma)
{
int i,j,wdth,hght,l;
int *nbx,*nby;
double moy,std,seuilhaut,seuilbas;
grphic* mask;

wdth=histo->width;
hght=histo->height;

mask=cr_grphic_modif(wdth,hght,0,1,0);

nbx=(int *)malloc(wdth*sizeof(int));
nby=(int *)malloc(hght*sizeof(int));

for (l=0; l<wdth; l++)
	{nbx[l]=0;}
for (l=0; l<hght; l++)
	{nby[l]=0;}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		mask->mri[i][j]=0;

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
		{
		nbx[i]=nbx[i]+histo->mri[i][j];	
		}

for (j=0; j<hght; j++)
for (i=0; i<wdth; i++)
		{
		nby[j]=nby[j]+histo->mri[i][j];	
		}


for (i=0; i<wdth; i++)
	{
	moy=0;
	for (j=0; j<hght; j++)
		moy+=j*histo->mri[i][j];
		
	moy=1.0*moy/nbx[i];
	
	std=0;
	for (j=0; j<hght; j++)
		std+=pow((j-moy),2.0)*histo->mri[i][j];
		
	std=sqrt(1.0*std/nbx[i]);
	
	seuilbas=moy-1.0*nb_sigma*std;
	seuilhaut=moy+1.0*nb_sigma*std;
	
	if(isnan(seuilbas))
		seuilbas=hght; 
	
	if(isnan(seuilhaut))
		seuilhaut=0; 
	
	
	for (j=0; j<hght; j++)
		if ((j>seuilbas)&&(j<seuilhaut))
		mask->mri[i][j]++;
	}


for (j=0; j<hght; j++)
	{
	moy=0;
	for (i=0; i<wdth; i++)
		moy+=i*histo->mri[i][j];
		
	moy=1.0*moy/nby[j];
	
	std=0;
	for (i=0; i<wdth; i++)
		std+=pow((i-moy),2.0)*histo->mri[i][j];
		
	std=sqrt(1.0*std/nby[j]);
	
	seuilbas=moy-1.0*nb_sigma*std;
	seuilhaut=moy+1.0*nb_sigma*std;
	
	if(isnan(seuilbas))
		seuilbas=wdth; 
	
	if(isnan(seuilhaut))
		seuilhaut=0; 
	
	
	for (i=0; i<wdth; i++)
		if ((i>seuilbas)&&(i<seuilhaut))
		mask->mri[i][j]++;
	}

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
if (mask->mri[i][j]==2)
	trimmedhisto->mri[i][j]=histo->mri[i][j];
else
	trimmedhisto->mri[i][j]=0;


free(nbx);free(nby);

free_grphic(mask);
}
*/

/*******************************************************************************
**        				 equalisation_histogramme                                                
**                                                                                                   
*******************************************************************************/

void equalisation_histogramme(grphic3d *imsrc,grphic3d *imres,double *ft,int size_ft)
{
int i,j,k,l,wdth,hght,dpth,min,max;
double dt;

wdth=imsrc->width;
hght=imsrc->height;
dpth=imsrc->depth;

imx_inimaxminpixel_3d_p(imsrc);
imx_copie_param_3d_p(imsrc,imres);


l=0;
for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if ((imsrc->mri[i][j][k]>0)&&(l<size_ft))
 		{
		ft[l]=imsrc->mri[i][j][k];
		l++;
		}



printf("debut tri  \n");
qsort(ft,size_ft,sizeof(double),double_compare_function2);
printf("fin tri  \n");

min=ft[0];
max=ft[size_ft-1];
dt=1.0*(max-min)/(size_ft-1.0);

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if (imsrc->mri[i][j][k]>0)
 	{
	l=0;
	while ((imsrc->mri[i][j][k]>ft[l]))
		{l=l+1000;
		if (l>size_ft)
			{l=l-1000;break;}
	
		}
	
	
	while ((imsrc->mri[i][j][k]>ft[l])&&(l<size_ft))
		l++;
		
	
	imres->mri[i][j][k]=(int)(min+1.0*l*dt);	
	}
	else
	imres->mri[i][j][k]=0;
}

/*******************************************************************************
**        				 inv_equalisation_histogramme                                                
**                                                                                                   
*******************************************************************************/

void inv_equalisation_histogramme(grphic3d *imsrc,grphic3d *imres,double *ft,int size_ft)
{
int i,j,k,wdth,hght,dpth,min,max,tmp;
double dt;

wdth=imsrc->width;
hght=imsrc->height;
dpth=imsrc->depth;

imx_inimaxminpixel_3d_p(imsrc);
imx_copie_param_3d_p(imsrc,imres);


min=ft[0];
max=ft[size_ft-1];
dt=1.0*(max-min)/(size_ft-1.0);

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if (imsrc->mri[i][j][k]>0)
 	{
	tmp=floor(1.0*(imsrc->mri[i][j][k]-min)/dt+0.5);
	imres->mri[i][j][k]=ft[tmp];
	}
	else
	imres->mri[i][j][k]=0;
}

/**************************** histo_joint_3d_p() *****************************/
/*                                                                           */                         
/*****************************************************************************/

void Update_image_from_histojointnorm_3d_p(grphic3d *im1, grphic3d *im2, grphic *implot)
{
int i, j, k;
int width, height, depth;
int plot_x,plot_y;
long min_x, max_x,min_y,max_y;

width=im1->width;
height=im1->height;
depth=im1->depth;

imx_inimaxminpixel_3d_p(im1);
imx_inimaxminpixel_3d_p(im2);

min_x=0;
max_x=im1->max_pixel;
min_y=0;
max_y=im2->max_pixel;

		
for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		if(im1->mri[i][j][k]>0 && im2->mri[i][j][k]>0)
		{
		plot_x=floor(1.0*(im1->mri[i][j][k]/**im1->rcoeff*/-min_x)/(max_x-min_x)*(implot->width-1.0));
		plot_y=floor(1.0*(im2->mri[i][j][k]/**im2->rcoeff*/-min_y)/(max_y-min_y)*(implot->height-1.0));
		
		if (implot->mri[plot_x][plot_y]==0)
			{im1->mri[i][j][k]=im2->mri[i][j][k]=0;}
		
		}			
		
}


/**************************** Troncature_intensite_nsigma_grphic3d_p() *****************************/
/*                                                                           */                         
/*****************************************************************************/
void Troncature_intensite_nsigma_grphic3d_p(grphic3d *im, grphic3d *imres, int nb_sigma)
{
int i,j,k,l;
int height,width,depth, seuil_haut, seuil_bas, seuil_haut_res, seuil_bas_res;
double moy,var,sigma;
double rcoeff,rcoeff_res,rapport_rcoeff;
float max,min;
int err;


l=0;
moy=0;
var=0;
height=im->height;
width=im->width;
depth=im->depth;
rcoeff=im->rcoeff;

/*cacul des moyennes des images*/	
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im->mri[i][j][k]>0)
      {
        l++;
        moy += im->mri[i][j][k];
      }
     }	
     	

moy= rcoeff*moy/l;

/*calcul des variances des images*/
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im->mri[i][j][k]>0)
          var +=(rcoeff*im->mri[i][j][k]-moy)*(rcoeff*im->mri[i][j][k]-moy);
      }

var= var/l;
sigma=sqrt(var);


seuil_bas= (int)((moy - nb_sigma*sigma)/rcoeff);
seuil_haut= (int)((moy + nb_sigma*sigma)/rcoeff);


imx_copie_param_3d_p(im,imres);
max=seuil_haut*rcoeff;
min=seuil_bas*rcoeff;


/*   Calcul de imres->maxpixel et imres->icomp et imres->rcoeff a partir de max*/
err=imx_brukermax_3d(max,min,imres);

rcoeff_res=imres->rcoeff;
rapport_rcoeff=rcoeff/rcoeff_res;
seuil_bas_res=(int)(1.0*seuil_bas*rcoeff/rcoeff_res);
seuil_haut_res=(int)(1.0*seuil_haut*rcoeff/rcoeff_res);

for (i=0;i<width;i++)
   for (j=0;j<height;j++)
      for(k=0;k<depth;k++)
       if (im->mri[i][j][k]!=0)
        {
	if (im->mri[i][j][k]>seuil_haut) imres->mri[i][j][k]=seuil_haut_res;
	if (im->mri[i][j][k]<seuil_bas) imres->mri[i][j][k]=seuil_bas_res;
	if ((im->mri[i][j][k]<=seuil_haut)&&(im->mri[i][j][k]>=seuil_bas)) imres->mri[i][j][k]=(int)(1.0*im->mri[i][j][k]*rapport_rcoeff);
	}
       else
        {imres->mri[i][j][k]=0;}


imx_inimaxminpixel_3d_p(imres);

}


/******************************************************************************
 ** -- equalisation_histo_3d ---------------------------------------------------------------
 **
 ******************************************************************************/
void equalisation_histo_3d(void)
{
  grphic3d *im1,*imres;
  int im_1,im_res;

  im_1=GET_PLACE3D( "Histogram equalization from");
  im_res=GET_PLACE3D( "to :");
  im1=ptr_img_3d(im_1);
  imres=ptr_img_3d(im_res);
  equalisation_histo_3d_p(im1,imres);
  show_picture_3d(im_res);
}

/******************************************************************************
 ** -- equalisation_histo_3d_p ------------------------------------------------
 **
 ******************************************************************************/
 void equalisation_histo_3d_p(grphic3d *im1,grphic3d *imres)
 {
int i,j,k,l,wdth,hght,dpth,Ntot,*tab1;
double dt;

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

imx_inimaxminpixel_3d_p(im1);
imx_copie_param_3d_p(im1,imres);
Ntot=0;

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if ((im1->mri[i][j][k]>0))
 	Ntot++;

tab1=(int *) malloc((Ntot+2)*sizeof(int));

l=1;
tab1[0]=MAXI(im1->min_pixel,0);


for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if ((im1->mri[i][j][k]>0))
 	{
	tab1[l]=im1->mri[i][j][k];
	l++;
	}


tab1[Ntot+1]=im1->max_pixel;
	
dt=(tab1[Ntot+1] - tab1[0])/(Ntot+1.0);
	
printf("debut tri 1 \n");
qsort(tab1,Ntot+2,sizeof(int),int_compare_function2);
printf("fin tri 1 \n");


for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
 if (im1->mri[i][j][k]>0)
 	{
	l=0;
	while (im1->mri[i][j][k]>tab1[l])
		l++;
	
	imres->mri[i][j][k]=(int)(tab1[0]+(double)l*dt);	
	}
	else
	imres->mri[i][j][k]=0;


imx_inimaxminpixel_3d_p(imres);
free(tab1);
}


 
 
