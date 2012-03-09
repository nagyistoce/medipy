/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>
#include <time.h>

#include "recalage/chps_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/divers_topo.h"
#include "math/imx_matrix.h"
#include "noyau/gui/imx_picture3d.h"

#ifndef COMPILE_FOR_MEDIPY
#include "noyau/gui/imx_plot.h"
#endif //fin COMPILE_FOR_MEDIPY

#include "math/imx_bspline.h"	/* pour les fonctions d'interpolation */
#include "noyau/mani_3d.h"
#ifdef __GTK_GUI
#include "noyau/gui/imx_3dmo.h"
#endif /* __GTK_GUI */
#include "noyau/io/imx_log.h"
#include "math/oper_3d.h"
#include "math/operation.h"
#include "traitement/trai_3d.h"
#include "math/ana_3d.h"
#include "noyau/io/imx_export_file.h"
#include "outils/imx_misc.h"
#include "recalage/topo/normalisation_topo.h"
#include "outils/imx_sort.h"
#include "morpho/morpho_3d.h"
#include "noyau/io/imx_file.h"
#include "noyau/io/imx_head.h"
#include "noyau/io/file_3d.h"

int NB_BLOC_POINT_FIXE;
int INV_BSPLINE_RESOL=-1;
int INV_BSPLINE_COEFF2L=-1;
double PRECISION_NUMERIQUE_INV_BSPLINE;

/*******************************************************************************
** 	int inversion_mat3d(mat3d *A, mat3d* invA);
**
**  Inversion d'une matrice 3*3
*******************************************************************************/
int inversion_mat3d(mat3d* A, mat3d* invA)
{
double d1,d2,d3,d;
int i,j;

d1=A->H[1][1]*A->H[2][2]-A->H[2][1]*A->H[1][2];
d2=A->H[0][2]*A->H[2][1]-A->H[0][1]*A->H[2][2];
d3=A->H[0][1]*A->H[1][2]-A->H[0][2]*A->H[1][1];
d=A->H[0][0]*d1+A->H[1][0]*d2+A->H[2][0]*d3;

if (d==0) return(-1);

invA->H[0][0]=d1;
invA->H[0][1]=d2;
invA->H[0][2]=d3;

invA->H[1][0]=A->H[2][0]*A->H[1][2]-A->H[1][0]*A->H[2][2];
invA->H[1][1]=A->H[0][0]*A->H[2][2]-A->H[0][2]*A->H[2][0];
invA->H[1][2]=A->H[1][0]*A->H[0][2]-A->H[0][0]*A->H[1][2];

invA->H[2][0]=A->H[1][0]*A->H[2][1]-A->H[2][0]*A->H[1][1];
invA->H[2][1]=A->H[0][1]*A->H[2][0]-A->H[0][0]*A->H[2][1];
invA->H[2][2]=A->H[0][0]*A->H[1][1]-A->H[0][1]*A->H[1][0];



for (i=0;i<3;i++)
	for (j=0;j<3;j++)
		invA->H[i][j]=1.0*invA->H[i][j]/d;

if ((d<0)||(A->H[2][2]<=0)||(d1<=0)) return(0);

return(1);
}

/*******************************************************************************
** 	double determinant_mat3d(mat3d *A);
**
**  Inversion d'une matrice 3*3
*******************************************************************************/
double determinant_mat3d(mat3d *A)
{
double d1,d2,d3,d;

d1=A->H[1][1]*A->H[2][2]-A->H[2][1]*A->H[1][2];
d2=A->H[0][2]*A->H[2][1]-A->H[0][1]*A->H[2][2];
d3=A->H[0][1]*A->H[1][2]-A->H[0][2]*A->H[1][1];
d=A->H[0][0]*d1+A->H[1][0]*d2+A->H[2][0]*d3;


if (isnan(d))
	printf("y a un bug determinant\n");
	

return(d);
}

/*******************************************************************************
**      jacobien_bspline1_3d(transfo,imres)
**	
**      jacobien d'une transformation Bspline d'ordre 1 
*******************************************************************************/

int	jacobien_bspline1_3d(transf3d *transfo, grphic3d *imres)
{
 	double  *p;
 	int     nb_param;
 	int	wdth,hght,dpth;
 	double  ***J;
 	double maxJ=-100000,minJ=100000,rcoeff,tmpJ;
	int singularite=0,i,j,k; 
 	double **toto,Jtot=0;
	int tdebut,tfin;

 toto=alloc_dmatrix(3,3);
 
	tdebut=time(NULL);

 if (transfo->degre!=1) return(0);
 wdth=transfo->width;hght=transfo->height;dpth=transfo->depth;
 J=alloc_dmatrix_3d(wdth,hght,dpth);

 p=CALLOC(transfo->nb_param,double);  
 
 nb_param=transfo->nb_param;
 
for (i=0;i<nb_param/3;i++)
 { 
 p[3*i]=(transfo->param)[3*i]/transfo->dx;
 p[3*i+1]=(transfo->param)[3*i+1]/transfo->dy;
 p[3*i+2]=(transfo->param)[3*i+2]/transfo->dz;
 }
 

 for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
	   {
  	   
			
		//test de la nouvelle fonction d'évaluation du jacobien 
		eval_jacobien_transfo_bspline1_3d(p, nb_param, wdth, hght, dpth, i, j, k, toto);
		tmpJ=toto[0][0]*(toto[1][1]*toto[2][2]-toto[1][2]*toto[2][1])
					-toto[1][0]*(toto[0][1]*toto[2][2]-toto[0][2]*toto[2][1])
					+toto[2][0]*(toto[0][1]*toto[1][2]-toto[0][2]*toto[1][1]);
		
   
   	Jtot=Jtot+tmpJ;
	 J[i][j][k]=tmpJ;
		if (tmpJ>maxJ) maxJ=tmpJ;
		if (tmpJ<minJ) minJ=tmpJ;
   
	 }
  
 printf("Jacobien minimum: %f  \t maximum : %f \n",minJ,maxJ);



 free(p);
 //calcul de l'image
 imx_brukermax_3d(maxJ,minJ,imres);
 rcoeff=imres->rcoeff;
  
    
  for (i=0;i<imres->width;i++)
   for (j=0;j<imres->height;j++)
    for (k=0;k<imres->depth;k++)
    {
     imres->mri[i][j][k]=(int)(J[i][j][k]/rcoeff);
		 
		 
		 if (J[i][j][k]<0)
		 	singularite++;
    }

 //printf("Jacobien minimum: %.2f \n",minJ);
 printf("%d singularite sur %d  soit %f pourcent \n",singularite,imres->width*imres->height*imres->depth,100.0*singularite/(imres->width*imres->height*imres->depth));

Jtot=Jtot/(wdth*hght*dpth);
   printf("Moyenne du Jacobien : %f \n",Jtot);
 
 imx_inimaxminpixel_3d_p(imres);
 free_dmatrix_3d(J);
 free_dmatrix(toto,3,3);
 
		  tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
 
 return(1);
}


/*******************************************************************************
**      integrale_jacobien_bspline1_compute(transfo,imres)
**	
**      integrale jacobien d'une transformation Bspline d'ordre 1 
*******************************************************************************/

int	integrale_jacobien_bspline1_compute(transf3d *transfo, grphic3d *imres)
{

int topD,i,j,k,l;
TSlpqr Slpqr[8];
int topi,topj,topk, nb_param,aux0,aux1;
int *x00,*x11,*y00,*y11,*z00,*z11;
int width,height,depth,resol;
double xm,xM,ym,yM,zm,zM;
scal_func_t scal_func=Bsplined1;
double *param_norm;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double auxdbl[20],auxdbl21a[21], auxdbl21b[21];   
double  ***J;
double maxJ=-100000,minJ=100000,rcoeff,tmpJ,xt,yt,zt,xt1,yt1,zt1,aux2,Jtot=0;
int singularite=0,vol;
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k

nb_param=transfo->nb_param;

param_norm=malloc(nb_param*sizeof(double));
   

nb_param=init_base_3d(transfo->width,transfo->height,transfo->depth,scal_func);

for (l=1;l<transfo->resol;l++)
     nb_param=base_resol_up_3d(param_norm,nb_param);


width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
vol=width*height*depth;
resol=BASE3D.resol; 
x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;

 J=alloc_dmatrix_3d(width,height,depth);
 
for (i=0;i<8;i++)
	{
	Slpqr[i].Jm=-HUGE_VAL;
	Slpqr[i].JM=HUGE_VAL;
	} 
  

for (l=0;l<nb_param/3;l++)
	{
	param_norm[3*l]=transfo->param[3*l]/width/transfo->dx;
	param_norm[3*l+1]=transfo->param[3*l+1]/height/transfo->dy;
	param_norm[3*l+2]=transfo->param[3*l+2]/depth/transfo->dz;
	}


topD=(int)pow(2.0,1.0*transfo->resol)-1;


for (topi=0; topi<topD; topi=topi+2)
for (topj=0; topj<topD; topj=topj+2)
for (topk=0; topk<topD; topk=topk+2)
{
	l = TOP_conv_ind(topi,topj,topk,nb_param);
  	//---------------------------------------------------------------
  	//----- initialisation Slpqr            -------------------------
  	//---------------------------------------------------------------

  	xm = 1.0*x00[l/3]/width; xM = 1.0*x11[l/3]/width; ym = 1.0*y00[l/3]/height; 
  	yM = 1.0*y11[l/3]/height; zm = 1.0*z00[l/3]/depth; zM =1.0*z11[l/3]/depth;

  	Slpqr[0].b.xm = xm; Slpqr[0].b.xM = (xm+xM)/2;
  	Slpqr[1].b.xm = xm; Slpqr[1].b.xM = (xm+xM)/2;
  	Slpqr[2].b.xm = xm; Slpqr[2].b.xM = (xm+xM)/2;
  	Slpqr[3].b.xm = xm; Slpqr[3].b.xM = (xm+xM)/2;
  	Slpqr[4].b.xM = xM; Slpqr[4].b.xm = (xm+xM)/2;
  	Slpqr[5].b.xM = xM; Slpqr[5].b.xm = (xm+xM)/2;
  	Slpqr[6].b.xM = xM; Slpqr[6].b.xm = (xm+xM)/2;
  	Slpqr[7].b.xM = xM; Slpqr[7].b.xm = (xm+xM)/2;

  	Slpqr[0].b.ym = ym; Slpqr[0].b.yM = (ym+yM)/2;
  	Slpqr[1].b.ym = ym; Slpqr[1].b.yM = (ym+yM)/2;
  	Slpqr[4].b.ym = ym; Slpqr[4].b.yM = (ym+yM)/2;
  	Slpqr[5].b.ym = ym; Slpqr[5].b.yM = (ym+yM)/2;
  	Slpqr[2].b.yM = yM; Slpqr[2].b.ym = (ym+yM)/2;
  	Slpqr[3].b.yM = yM; Slpqr[3].b.ym = (ym+yM)/2;
  	Slpqr[6].b.yM = yM; Slpqr[6].b.ym = (ym+yM)/2;
  	Slpqr[7].b.yM = yM; Slpqr[7].b.ym = (ym+yM)/2;

  	Slpqr[0].b.zm = zm; Slpqr[0].b.zM = (zm+zM)/2;
  	Slpqr[2].b.zm = zm; Slpqr[2].b.zM = (zm+zM)/2;
  	Slpqr[4].b.zm = zm; Slpqr[4].b.zM = (zm+zM)/2;
  	Slpqr[6].b.zm = zm; Slpqr[6].b.zM = (zm+zM)/2;
  	Slpqr[1].b.zM = zM; Slpqr[1].b.zm = (zm+zM)/2;
  	Slpqr[3].b.zM = zM; Slpqr[3].b.zm = (zm+zM)/2;
  	Slpqr[5].b.zM = zM; Slpqr[5].b.zm = (zm+zM)/2;
  	Slpqr[7].b.zM = zM; Slpqr[7].b.zm = (zm+zM)/2;



  //-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

  TOP_init_slpqr(Slpqr  ,topi-1,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+1,topi-1,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+2,topi-1,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+3,topi-1,topj  ,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+4,topi  ,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+5,topi  ,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+6,topi  ,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+7,topi  ,topj  ,topk  ,D,nb_param,param_norm);


  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
	for (aux0=0; aux0<8; aux0++)
  	{ //--- initialisation du terme constant ------------------------------
    	TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    	aux2 = pow(2.0,-1.0 * resol);
		poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,auxdbl21b);
    	for (aux1=0; aux1<20; aux1++) Slpqr[aux0].alpha[aux1].ordorg = auxdbl21b[aux1+1];
  
	ideb=Slpqr[aux0].b.xm*width;
	ifin=Slpqr[aux0].b.xM*width;
	jdeb=Slpqr[aux0].b.ym*height;
	jfin=Slpqr[aux0].b.yM*height;
	kdeb=Slpqr[aux0].b.zm*depth;
	kfin=Slpqr[aux0].b.zM*depth;

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	{
	
        
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		xt1=(double)(i+1)/width;yt1=(double)(j+1)/height;zt1=(double)(k+1)/depth;
		eval_integrale_fun(auxdbl21b,  xt, yt, zt,  xt1, yt1, zt1, &tmpJ);
		
		Jtot=Jtot+tmpJ;	
		tmpJ=tmpJ*vol;	
		
	     J[i][j][k]=tmpJ;
		if (tmpJ>maxJ) maxJ=tmpJ;
		if (tmpJ<minJ) minJ=tmpJ;
	    
	    }
	}   



  
}

 //calcul de l'image
 imx_brukermax_3d(maxJ,minJ,imres);
 rcoeff=imres->rcoeff;

for (i=0;i<imres->width;i++)
   for (j=0;j<imres->height;j++)
    for (k=0;k<imres->depth;k++)
    {
     imres->mri[i][j][k]=(int)(J[i][j][k]/rcoeff);
		 
		 
		 if (J[i][j][k]<0)
		 	singularite++;
    }

 printf("Jacobien minimum: %f  \t maximum : %f \n",minJ,maxJ);
  printf("%d singularite sur %d  soit %f pourcent \n",singularite,imres->width*imres->height*imres->depth,100.0*singularite/(imres->width*imres->height*imres->depth));
 printf("Moyenne du jacobien sur l'ensemble de l'image :  %f \n",Jtot);
 imx_inimaxminpixel_3d_p(imres);
 free_dmatrix_3d(J);
 free(param_norm);
 end_base_3d();


 
 return(1);
}

/*******************************************************************************
**      borne_jacobien_bspline1_local_3d()
**	
**      jacobien d'une transformation Bspline d'ordre 1 
*******************************************************************************/
int	borne_jacobien_bspline1_local_3d(int nb_param, double *p, double *param_norm,int topi, int topj, int topk,double *min,double *max)
{
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
 	int	wdth,hght,dpth,l;
 	TSlpqr Slpqr[8];
 	double xm,xM,ym,yM,zm,zM;
 	int resol,topD,topDx,topDy,topDz;
	double *fx,*fy,*fz;
	int i,j,k;
	int aux0,aux1;  								// variables auxiliaires
	double auxdbl[20],aux2;
	double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
	int D ; 			     
	double tmpJ,mymono[21],xt,yt,zt;
	int ideb,ifin,jdeb,jfin,kdeb,kfin;
	int bx0,bx1,by0,by1,bz0,bz1;
	double maxJ=-100000,minJ=100000;
	
 
 topD=D=TOP_D(nb_param);
 
 wdth=BASE3D.width;
 hght=BASE3D.height;
 dpth=BASE3D.depth;
  
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
 resol=BASE3D.resol; 
 l = TOP_conv_ind(topi,topj,topk,nb_param)/3;
  	
			  //---------------------------------------------------------------
  			//----- initialisation Slpqr            -------------------------
  			//---------------------------------------------------------------

  			xm = 1.0*x0[l]/wdth; xM = 1.0*x1[l]/wdth; ym = 1.0*y0[l]/hght; 
  			yM = 1.0*y1[l]/hght; zm = 1.0*z0[l]/dpth; zM =1.0*z1[l]/dpth;

  			Slpqr[0].b.xm = xm; Slpqr[0].b.xM = (xm+xM)/2;
  			Slpqr[1].b.xm = xm; Slpqr[1].b.xM = (xm+xM)/2;
  			Slpqr[2].b.xm = xm; Slpqr[2].b.xM = (xm+xM)/2;
  			Slpqr[3].b.xm = xm; Slpqr[3].b.xM = (xm+xM)/2;
  			Slpqr[4].b.xM = xM; Slpqr[4].b.xm = (xm+xM)/2;
  			Slpqr[5].b.xM = xM; Slpqr[5].b.xm = (xm+xM)/2;
  			Slpqr[6].b.xM = xM; Slpqr[6].b.xm = (xm+xM)/2;
  			Slpqr[7].b.xM = xM; Slpqr[7].b.xm = (xm+xM)/2;

  			Slpqr[0].b.ym = ym; Slpqr[0].b.yM = (ym+yM)/2;
  			Slpqr[1].b.ym = ym; Slpqr[1].b.yM = (ym+yM)/2;
  			Slpqr[4].b.ym = ym; Slpqr[4].b.yM = (ym+yM)/2;
  			Slpqr[5].b.ym = ym; Slpqr[5].b.yM = (ym+yM)/2;
  			Slpqr[2].b.yM = yM; Slpqr[2].b.ym = (ym+yM)/2;
  			Slpqr[3].b.yM = yM; Slpqr[3].b.ym = (ym+yM)/2;
  			Slpqr[6].b.yM = yM; Slpqr[6].b.ym = (ym+yM)/2;
  			Slpqr[7].b.yM = yM; Slpqr[7].b.ym = (ym+yM)/2;

  			Slpqr[0].b.zm = zm; Slpqr[0].b.zM = (zm+zM)/2;
  			Slpqr[2].b.zm = zm; Slpqr[2].b.zM = (zm+zM)/2;
  			Slpqr[4].b.zm = zm; Slpqr[4].b.zM = (zm+zM)/2;
  			Slpqr[6].b.zm = zm; Slpqr[6].b.zM = (zm+zM)/2;
  			Slpqr[1].b.zM = zM; Slpqr[1].b.zm = (zm+zM)/2;
  			Slpqr[3].b.zM = zM; Slpqr[3].b.zm = (zm+zM)/2;
  			Slpqr[5].b.zM = zM; Slpqr[5].b.zm = (zm+zM)/2;
  			Slpqr[7].b.zM = zM; Slpqr[7].b.zm = (zm+zM)/2;

			  


  //-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

  TOP_init_slpqr(Slpqr  ,topi-1,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+1,topi-1,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+2,topi-1,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+3,topi-1,topj  ,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+4,topi  ,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+5,topi  ,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+6,topi  ,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+7,topi  ,topj  ,topk  ,D,nb_param,param_norm);


bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;


// On balaye chacun des pixels de la boite

  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,auxdbl21b);
    for (aux1=0; aux1<20; aux1++) Slpqr[aux0].alpha[aux1].ordorg = auxdbl21b[aux1+1];
  
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	   {
    xt=(double)i/wdth;yt=(double)j/hght;zt=(double)k/dpth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &tmpJ, mymono);
		if (tmpJ>maxJ) maxJ=tmpJ;
		if (tmpJ<minJ) minJ=tmpJ;
      }

   }
   
*max=maxJ;
*min=minJ;

 return(1);
}

/*******************************************************************************
**      test_correlation
**	
**      elimine les blocs ou les residus sont nuls
*******************************************************************************/


void test_correlation(grphic3d *imref, grphic3d *imres,int nb_param,int ***masque_param)
{
int D,wdth,hght,dpth,taille;
int topi,topj,topk,i,j,k,ind,nb,nb_tot,N;
int *x0,*x1,*y0,*y1,*z0,*z1;
double moyenne,sigma,percentage,seuil;
int ***difference,m;
double *residu;


wdth=imref->width;
hght=imref->height;
dpth=imref->depth;

D=TOP_D(nb_param);
x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

nb_tot= nb=0;
printf("debut test_correlation\n");	

	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if ((imref->mri[i][j][k]!=0)&&(imres->mri[i][j][k]!=0))
			nb_tot++;
			
residu=(double *)malloc(nb_tot*sizeof(double));

m=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if ((imref->mri[i][j][k]!=0)&&(imres->mri[i][j][k]!=0))
			{residu[m]=fabs(imref->mri[i][j][k]-imres->mri[i][j][k]); m++;}

qsort(residu,nb_tot,sizeof(double),double_compare_function2);

sigma=1.4826*residu[(int)floor(nb_tot/2)];

seuil=1.96*sigma;
			
taille=x1[0]-x0[0]; 
N=taille*taille*taille;
difference=alloc_imatrix_3d(taille,taille,taille);


nb_tot=0;

for (topi=0;topi<D;topi++)
for (topj=0;topj<D;topj++)
for (topk=0;topk<D;topk++)
	if (masque_param[topi][topj][topk]!=0)
	{
	nb_tot++;
	ind=1.0*TOP_conv_ind(topi,topj,topk,nb_param)/3;
	N=0;
	for (i=x0[ind];i<x1[ind];i++)
	for (j=y0[ind];j<y1[ind];j++)
	for (k=z0[ind];k<z1[ind];k++)
		if ((imref->mri[i][j][k]!=0)&&(imres->mri[i][j][k]!=0))
		{
		difference[i-x0[ind]][j-y0[ind]][k-z0[ind]]=1.0*(imref->mri[i][j][k]-imres->mri[i][j][k]);
		N++;
		}
		
	moyenne=0;
	for (i=0;i<taille;i++)
	for (j=0;j<taille;j++)
	for (k=0;k<taille;k++)			
		moyenne=moyenne+difference[i][j][k];
	
	moyenne=1.0*moyenne/N;
	

	if (abs(moyenne)<1.0*seuil/sqrt(N))
			{masque_param[topi][topj][topk]=0;nb++;}

	}

percentage=100.0*nb/nb_tot;
printf("%f %% bloc en moins dans test_correlation\n",percentage);	

free(residu);
free_imatrix_3d(difference);
}


/*******************************************************************************
**      test_correlation2
**	
**      elimine les blocs ou les residus sont nuls
*******************************************************************************/

void test_correlation2(grphic3d *imref, grphic3d *imreca,int nb_param,int ***masque_param)
{
int D,wdth,hght,dpth;
int topi,topj,topk,i,j,k,ind,nb,nb_tot,N,l,m;
int *x0,*x1,*y0,*y1,*z0,*z1;
double percentage,moy1,var1,moy2,var2,seuil,sigma1,sigma2,cov;
	
wdth=imref->width;
hght=imref->height;
dpth=imref->depth;

D=TOP_D(nb_param);
x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

nb = nb_tot= 0;
printf("debut test_correlation\n");	
seuil=0.95;

for (topi=0;topi<D;topi++)
for (topj=0;topj<D;topj++)
for (topk=0;topk<D;topk++)
	if (masque_param[topi][topj][topk]!=0)
	{
	nb_tot++;
	ind=1.0*TOP_conv_ind(topi,topj,topk,nb_param)/3;
	
	// calcul de la moyenne
	l=0;
	moy1=0;
	for (i=x0[ind];i<x1[ind];i++)
	for (j=y0[ind];j<y1[ind];j++)
	for (k=z0[ind];k<z1[ind];k++)
		if (imref->mri[i][j][k]!=0)
		{
		moy1+=imref->mri[i][j][k];
		l++;
		}
	moy1=1.0*moy1/l;
	
	m=0;	
	moy2=0;
	for (i=x0[ind];i<x1[ind];i++)
	for (j=y0[ind];j<y1[ind];j++)
	for (k=z0[ind];k<z1[ind];k++)
		if (imreca->mri[i][j][k]!=0)
		{
		moy2+=imreca->mri[i][j][k];
		m++;
		}
	moy2=1.0*moy2/m;
		
	// calcul des variance
	var1=0;
	for (i=x0[ind];i<x1[ind];i++)
	for (j=y0[ind];j<y1[ind];j++)
	for (k=z0[ind];k<z1[ind];k++)
		if (imref->mri[i][j][k]!=0)
		{
		var1+=pow((imref->mri[i][j][k]-moy1),2.0);
		}
		
	var2=0;
	for (i=x0[ind];i<x1[ind];i++)
	for (j=y0[ind];j<y1[ind];j++)
	for (k=z0[ind];k<z1[ind];k++)
		if (imreca->mri[i][j][k]!=0)
		{
		var2+=pow((imreca->mri[i][j][k]-moy2),2.0);
		}
	
	sigma1=sqrt(1.0*var1/l);
	sigma2=sqrt(1.0*var2/m);
	

/*calcul de la covariance*/
cov=0;
N=0;
	for (i=x0[ind];i<x1[ind];i++)
	for (j=y0[ind];j<y1[ind];j++)
	for (k=z0[ind];k<z1[ind];k++)
		 	 if ((imref->mri[i][j][k])||(imreca->mri[i][j][k]))
			{
			cov+=(imref->mri[i][j][k]-moy1)*(imreca->mri[i][j][k]-moy2); 
			N++;
			}

cov=1.0*cov/(N*sigma1*sigma2);
cov=cov*cov;	
	
	if (cov>seuil)
			{masque_param[topi][topj][topk]=0;
		nb++;}

	}

percentage=100.0*nb/nb_tot;
printf("%f %% bloc en moins dans test_correlation\n",percentage);	

}

double dim_min_box(TOPTbox b)
{ double wx, wy, wz;
  wx = b.xM - b.xm;  wy = b.yM - b.ym; wz = b.zM - b.zm;
  if ((wx<=wy)&&(wx<=wz)) return wx;
  else if (wy<=wz) return wy;
  else return wz;
}



TOPTbox find_box(double x0, double y0, double z0, TOPTliste myListe, int *status)
{
TOPTbox myBox,*ptr;

*status=0;

ptr=myListe.tete;


while ((ptr!=NULL)&&(*status==0))
	{
	if ((x0>=ptr->xm)&&(x0<ptr->xM)&&(y0>=ptr->ym)&&(y0<ptr->yM)&&(z0>=ptr->zm)&&(z0<ptr->zM))
		{
		myBox=*ptr; 
		*status=1;
		}
	ptr=ptr->next;
	}

return myBox;
}


/******************************************************************************
** -- save_histo ---------------------------------------------------------------
** 
**   Enregistre les valeurs de l'histogramme dans un fichier texte
**
******************************************************************************/

void save_histo(void)
{ 
  grphic3d *im1;
  int pos3d,err,Nbpas;
	float xmax,xmin;
 	char tab[80],nomfichres[500],str[500];

  pos3d=GET_WND3D(TEXT0001);

	im1=ptr_img_3d(pos3d);
  
  sprintf(tab,"max=%f Choose the maximum :\n",im1->max_pixel*im1->rcoeff);
  xmax=GET_FLOAT(tab, 2,&err);
  if(err!=0) { printf("Attention ...\n"); }
  sprintf(tab,"min=%f Choose the minimum :\n",im1->min_pixel*im1->rcoeff);
  xmin=GET_FLOAT(tab, 2, &err);
  if(err!=0) {printf("Attention ...\n"); }
  Nbpas= GET_INT(TEXT0067, 2, &err);
  if(err!=0) { printf("Attention ...\n");}
	sprintf(str,"File ? in %s",_CurrentPath);
 	strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 
 save_histo_3d(pos3d,xmax,xmin,Nbpas,nomfichres);
}



/*******************************************
** --  histo_3d() ----------------------------
**
** Trace de l'histogramme de la ZI 3D
********************************************/
int	save_histo_3d(int im_1,float xmax,float xmin,int Nb,char *nomfichres)
{ 
  int 	roi3d_;
  int 	i,Nbpas;
  float	*x,*y;
  float moy,ecty,skew,kurt;
  float mediane,per5quantile,per95quantile;
  grphic3d *im1,*roi3d;
  FILE * fic;

  im1=ptr_img_3d(im_1);
  roi3d_=-im_1;
  roi3d=ptr_img_3d(roi3d_);

  y=CALLOC(Nb+2,float);
  x=CALLOC(Nb+2,float);

  for (i=0;i<Nb;i++) y[i]=0.0;
  if( imx_histogram_3d_p(im1,roi3d,xmin,xmax,Nb,x,y) == 0 ) return(0);

/*  Trace de l'histogramme */
#ifndef COMPILE_FOR_MEDIPY
	plot_histo(x,y,Nb,"Histogram of the ROI","ROI");
#endif

/*      Calcul statistique de l'histo   */
  imx_statistique_3d(im_1,roi3d_,xmin,xmax,&Nbpas,&moy,&ecty,&skew,&kurt
		  ,&mediane,&per5quantile,&per95quantile);

  PUTI(TEXT0097,Nbpas," _");
  PUTF(TEXT0098,moy," _");
  PUTF("+/-_",ecty," ");
  PUTF("Skewness >0 a droit <0 a gauche =_",skew," _");
  PUTF("Kurtosis >0 up <0 platy =_",kurt," ");
  PUTF("mediane= =_",mediane,"_");
  PUTF("quantile d'ordre 20 (5%)=_",per5quantile,"_");
  PUTF(",",per95quantile,"");

/* enregistrement dans un fichier */

/* Ouverture du fichier */
fic = fopen(nomfichres, "w+");

/* Ecriture des donnï¿½s dans le fichier */
for (i=0;i<Nb;i++)
	fprintf(fic, "%f \t %f \n", x[i], 1.0*y[i]/Nbpas);

/* Fermeture du fichier */
fclose(fic);


  FREE(x);
  FREE(y);
  return(1);
}

double my_funct(double x, double y, double z, TOPTliste liste1, TOPTliste liste2, int i, int j,int k )
{
double res;
TOPTbox myBaux;
int err;
double phi,x1,y1,z1;

myBaux=find_box(x,y,z,liste1,&err);
		
		if (err==0) { printf("Y a pas de boite qui marche !!!\n"); myBaux.ax=myBaux.ay=myBaux.az=0;}
		phi=1.0*sin(PI*1.0*(x-myBaux.xm)/(myBaux.xM-myBaux.xm))*sin(PI*1.0*(y-myBaux.ym)/(myBaux.yM-myBaux.ym))*sin(PI*1.0*(z-myBaux.zm)/(myBaux.zM-myBaux.zm));
		x1=x+myBaux.ax*phi;
		y1=y+myBaux.ay*phi;
		z1=z+myBaux.az*phi;
		//printf("%f %f %f \r",x1,y1,z1);
		myBaux=find_box(x1,y1,z1,liste2,&err);
		
		if (err==0) { printf("Y a pas de boite qui marche !!!\n"); myBaux.ax=myBaux.ay=myBaux.az=0;}
		phi=1.0*sin(PI*1.0*(x1-myBaux.xm)/(myBaux.xM-myBaux.xm))*sin(PI*1.0*(y1-myBaux.ym)/(myBaux.yM-myBaux.ym))*sin(PI*1.0*(z1-myBaux.zm)/(myBaux.zM-myBaux.zm));
		
		
		res=sqrt(pow(((x1+myBaux.ax*phi)-i),2.0)+pow(((y1+myBaux.ay*phi)-j),2.0)+pow(((z1+myBaux.az*phi)-k),2.0));
		//printf("erreur : %f \n",res);
	return (res);
}

double my_funct2(double x, double y, double z, TOPTliste liste1, int i, int j,int k )
{
double res;
TOPTbox myBaux;
int err;
double phi,x1,y1,z1;

myBaux=find_box(x,y,z,liste1,&err);
		
		if (err==0) { printf("Y a pas de boite qui marche !!!\n"); myBaux.ax=myBaux.ay=myBaux.az=0;}
		phi=1.0*sin(PI*1.0*(x-myBaux.xm)/(myBaux.xM-myBaux.xm))*sin(PI*1.0*(y-myBaux.ym)/(myBaux.yM-myBaux.ym))*sin(PI*1.0*(z-myBaux.zm)/(myBaux.zM-myBaux.zm));
		x1=x+myBaux.ax*phi;
		y1=y+myBaux.ay*phi;
		z1=z+myBaux.az*phi;
		//printf("%f %f %f \r",x1,y1,z1);
		
		
		res=sqrt(pow((x1-i),2.0)+pow((y1-j),2.0)+pow((z1-k),2.0));
		//printf("erreur : %f \n",res);
	return (res);
}

/*******************************************************************************
**   apply_inverse_transf_3d()                                  
**                                                                       
**    	appliquer a une image une transformation contenu dans un fichier
**    TODO ... anisotropie?        
*******************************************************************************/
void apply_inverse_transf_3d(void)
{

 field3d *champ;
 transf3d *transfo,*transfores;
 int	im_deb,im_res,e,inter_type;
 char 	nomfichier[256];
 grphic3d *imdeb,*imres, *imTemp=NULL;
 InterpolationFct interpol; 

  champ=NULL;
 
 //question sur les images
 im_deb=GET_PLACE3D(TEXT0030);
 im_res=GET_PLACE3D(TEXT0006);
 imdeb=ptr_img_3d(im_deb);imres=ptr_img_3d(im_res);


 // Creation de imtemp si necessaire
  if (imdeb==imres)
  {
   imTemp=cr_grphic3d(imdeb);
   if (!imTemp) { fprintf (stderr, "memory allocation failed in imx_apply_transf_3d_p\n"); return; }
   imx_copie_3d_p(imdeb, imTemp);
  }
  else imTemp=imdeb;



 //fichier contenant la transformation
 strcpy(nomfichier,GET_FILE("*.trf", &e));
 if ((test_fichier_transf_3d(nomfichier))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}

 //methode d'interpolation
 inter_type=imx_query_interpolation_type_3D(0);

 //chargement de la transformation
 transfo=load_transf_3d(nomfichier);
 if( transfo->typetrans != RIGID3D &&
	 transfo->typetrans != RIGIDZOOM3D &&
	 transfo->typetrans != AFFINE3D       )
 {
	 PUT_ERROR("Inversion of this transform not supported");
	 return;
 }

  transfores=cr_transf3d(transfo->width,transfo->height,transfo->depth,NULL);

	imres->dx=transfo->dx ; imres->dy=transfo->dy ; imres->dz=transfo->dz;
  imres->width=transfo->width ; imres->height=transfo->height ; imres->depth=transfo->depth ;


 	interpol=imx_choose_interpolation_fct(inter_type);

	inverse_transf_3d_p(transfo,transfores,0.0);

 //allocation du champ
  champ=	transf_to_field_3d(transfores,imTemp,imres);
	
 //calcul de l'image resultat
  interpol(imTemp,champ,imres); 
 
  // mise a jour de imres 
  imx_inimaxminpixel_3d_p(imres);
 	imx_reinitialise_visual_params_3d_p(imres);

 show_picture_3d(im_res);
 free_field3d(champ);
 free_transf3d(transfo);
 free_transf3d(transfores);
 if (imTemp!=imdeb) free_grphic3d(imTemp);

}

/*******************************************************************************
** 	inverse_transf_3d                                   
**                                                                   
**	Cette fonction permet de calculer l'inverse d'un champ avec la méthode   
**     itkInverseDeformationFieldImageFilter de Itk 
*******************************************************************************/
void	itkInverseDeformationFieldImageFilter_3d(void)
{
#ifdef  HAS_ITK
char  nomfich1[500],nomfichres[500],str[500];
transf3d *transfo1,*transfores;
int	e=0,wdth,hght,dpth,subsampling;
field3d *ch1,*chres;

 /*lecture du premier champ*/
 strcpy(nomfich1,GET_FILE("fichier trf a inverser",&e));
 if(e != 0)
     return ;
 put_file_extension(nomfich1,".trf",nomfich1);

  if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}


 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&e));
 put_file_extension(nomfichres,".trf",nomfichres);


    	/*lecture des champ*/
    	transfo1=load_transf_3d(nomfich1);
		
		wdth=transfo1->width;
		hght=transfo1->height;
		dpth=transfo1->depth;
		
	ch1=transf_to_field_3d(transfo1,NULL,NULL);
	
    		
    	subsampling=GET_INT("Subsampling factor", 16, &e);
	
	chres=cr_field3d(wdth,hght,dpth);
	
	itkInverseDeformationFieldImageFilter_3d_p(ch1,chres,subsampling);
	//itk_demons_registration_filter_3d_menu();
		
	transfores=field_to_transf_3d(chres,NULL,NULL);
	transfores->typetrans=CHAMP3D;
	transfores->dx=transfo1->dx;
	transfores->dy=transfo1->dy;
	transfores->dz=transfo1->dz;
	
	free_field3d(chres);
		
	/*enregistrement du champ resultat*/
    	save_transf_3d(transfores,nomfichres);
    	
	   /* Libere memoire allouee */
    free_transf3d(transfo1); free_transf3d(transfores);
#endif     
}	




/*******************************************************************************
** 	inverse_transf_3d                                   
**                                                                   
**	Cette fonction permet de calculer la transfo inverse       
*******************************************************************************/
void	inverse_transf_3d(void)
{
char  nomfich1[500],nomfichres[500],str[500];
transf3d *transfo1,*transfores;
int	e=0,wdth,hght,dpth,i,choix, wdth_res, hght_res, dpth_res,im_ref;
float dxres, dyres, dzres;
double *param,prec;
field3d /**ch1,*/*chres;
char *quest[4];
grphic3d *imref;
	
 for(i=0;i<4;i++)
    quest[i]=CALLOC(80,char);

 /*lecture du premier champ*/
 strcpy(nomfich1,GET_FILE("fichier trf a inverser",&e));
 if(e != 0)
     return ;
 put_file_extension(nomfich1,".trf",nomfich1);

  if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}


 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&e));
 put_file_extension(nomfichres,".trf",nomfichres);


 /*distance*/
  strcpy(quest[0],"Ne pas preciser");
  strcpy(quest[1],"Preciser manuellement");
  strcpy(quest[2],"Preciser a partir d'une image");
  strcpy(quest[3],"\0");

  choix=GETV_QCM("Caracteristiques de l'espace de depart",(const char **)quest);

  switch(choix)
  	{
	case 0:
	break;
	
	case 1:
	wdth_res=GET_INT("Width", 256, &e);
	hght_res=GET_INT("Height", 256, &e);
	dpth_res=GET_INT("Depth", 256, &e);

	dxres=GET_FLOAT("dx", 1, &e);
	dyres=GET_FLOAT("dy", 1, &e);
	dzres=GET_FLOAT("dz", 1, &e);

	break;
	
	case 2:
	
	im_ref= GET_PLACE3D("image de reference");
	imref=ptr_img_3d(im_ref);
	
	wdth_res=imref->width;
 	hght_res=imref->height;
 	dpth_res=imref->depth;

 	dxres=imref->dx;
 	dyres=imref->dy;
 	dzres=imref->dz;

	break;
	
	default:
	choix=0;
	break;
	}
	


    /*lecture des champ*/
    transfo1=load_transf_3d(nomfich1);
		
		wdth=transfo1->width;
		hght=transfo1->height;
		dpth=transfo1->depth;
		

	if (choix==0)
		{
		wdth_res=transfo1->width;
 		hght_res=transfo1->height;
 		dpth_res=transfo1->depth;

 		dxres=transfo1->dx;
 		dyres=transfo1->dy;
 		dzres=transfo1->dz;
		}


    transfores=cr_transf3d(wdth_res,hght_res,dpth_res,NULL);
		
if ((transfo1->typetrans==RIGID3D)||	(transfo1->typetrans==RIGIDZOOM3D)||(transfo1->typetrans==AFFINE3D))
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
		{	rigid_to_rigidz_3d(param);rigidz_to_affine_3d(param);}
	
	if (transfo1->typetrans==RIGIDZOOM3D)
		{	rigidz_to_affine_3d(param);}
	
	transfo1->param=param;
		    	 
	 	e=inv_affine_3d(transfo1,transfores);
		}
 else 
 	{
		
		
		if (transfo1->typetrans==BSPLINE3D)
			{
			chres=cr_field3d(wdth,hght,dpth);
		
			prec=GET_FLOAT("Precision (en voxel)", 0.1, &e);

			e=inv_bspline_3d(transfo1,chres,prec);
			}
		else
			{
			
			prec=GET_FLOAT("Precision (en voxel)", 0.1, &e);
	
		//ch1=transf_to_field_3d(transfo1,NULL,NULL);
		
		
		/*conversion en champ voxelique*/
		/*for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			ch1->raw[i][j][k].x=ch1->raw[i][j][k].x/transfo1->dx;
			ch1->raw[i][j][k].y=ch1->raw[i][j][k].y/transfo1->dy;
			ch1->raw[i][j][k].z=ch1->raw[i][j][k].z/transfo1->dz;
			}*/
			
		chres=cr_field3d(wdth,hght,dpth);
	
		e=inv_field_3d(transfo1,chres,prec);
		
		//free_field3d(ch1);
	
		/*conversion en champ milimetrique*/
		/*for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			chres->raw[i][j][k].x=chres->raw[i][j][k].x*transfo1->dx;
			chres->raw[i][j][k].y=chres->raw[i][j][k].y*transfo1->dy;
			chres->raw[i][j][k].z=chres->raw[i][j][k].z*transfo1->dz;
			}*/
			
		}
	
		
		
		transfores=field_to_transf_3d(chres,NULL,NULL);
	  transfores->typetrans=CHAMP3D;
		transfores->dx=dxres;
		transfores->dy=dyres;
		transfores->dz=dzres;
	
		free_field3d(chres);
	}   
		
		if (e==1)
			{
			
	 		/*enregistrement du champ resultat*/
    	save_transf_3d(transfores,nomfichres);
    	}
    /* Libere memoire allouee */
    free_transf3d(transfo1); free_transf3d(transfores);
     
}	



/*******************************************************************************
** 	inverse_transf_3d                                   
**                                                                   
**	Cette fonction permet de calculer la transfo inverse       
*******************************************************************************/
void	inverse_transf_3d_p(transf3d *transfo1,transf3d *transfores,double prec)
{
int	e=0,wdth,hght,dpth,i/*,j,k*/;
double *param;
field3d /**ch1,*/*chres;
		
		wdth=transfo1->width;
		hght=transfo1->height;
		dpth=transfo1->depth;
		
		
if ((transfo1->typetrans==RIGID3D)||	(transfo1->typetrans==RIGIDZOOM3D)||(transfo1->typetrans==AFFINE3D))
		{	
		/* copie des parametres */
    transfores->typetrans=AFFINE3D;
		transfores->dx=transfo1->dx;
		transfores->dy=transfo1->dy;
		transfores->dz=transfo1->dz;
		transfores->nb_param=15;
		transfores->param=CALLOC(transfores->nb_param,double);
		param=CALLOC(transfores->nb_param,double);
		
		for(i=0;i<transfo1->nb_param;i++)
			param[i]=transfo1->param[i];
		
	if (transfo1->typetrans==RIGID3D)
		{	rigid_to_rigidz_3d(param);rigidz_to_affine_3d(param);}
	
	if (transfo1->typetrans==RIGIDZOOM3D)
		{	rigidz_to_affine_3d(param);}
	
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
				
		//ch1=transf_to_field_3d(transfo1,NULL,NULL);
		
		
		/*conversion en champ voxelique*/
		/*for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			ch1->raw[i][j][k].x=ch1->raw[i][j][k].x/transfo1->dx;
			ch1->raw[i][j][k].y=ch1->raw[i][j][k].y/transfo1->dy;
			ch1->raw[i][j][k].z=ch1->raw[i][j][k].z/transfo1->dz;
			}*/
			
		chres=cr_field3d(wdth,hght,dpth);
	
		e=inv_field_3d(transfo1,chres,prec);
		
		//free_field3d(ch1);
	
		/*conversion en champ milimetrique*/
		/*for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			chres->raw[i][j][k].x=chres->raw[i][j][k].x*transfo1->dx;
			chres->raw[i][j][k].y=chres->raw[i][j][k].y*transfo1->dy;
			chres->raw[i][j][k].z=chres->raw[i][j][k].z*transfo1->dz;
			}*/
			
		}
	
		
		
		//transfores=field_to_transf_3d(chres,NULL,NULL);
		field_to_transf_3d_noalloc(chres,transfores,NULL,NULL);
	  transfores->typetrans=CHAMP3D;
		transfores->dx=transfo1->dx;
		transfores->dy=transfo1->dy;
		transfores->dz=transfo1->dz;
	
		free_field3d(chres);
	}   
		     
}	
/*******************************************************************************
** 	inv_affine_3d(transfo,transfores)                                   
**                                                                   
**	Cette fonction permet de calculer la transfo inverse       
*******************************************************************************/
int inv_affine_3d(transf3d *transfo,transf3d *transfores)
{
double a11,a12,a13,a21,a22,a23,a31,a32,a33; 
double tx,ty,tz;
int res=0;
mat3d A,invA;


		a11=transfo->param[0];a12=transfo->param[1];a13=transfo->param[2];
		a21=transfo->param[3];a22=transfo->param[4];a23=transfo->param[5];
		a31=transfo->param[6];a32=transfo->param[7];a33=transfo->param[8];
		
		tx=transfo->param[9];
		ty=transfo->param[10];
		tz=transfo->param[11];

		A.H[0][0]=a11;A.H[0][1]=a12;A.H[0][2]=a13;
		A.H[1][0]=a21;A.H[1][1]=a22;A.H[1][2]=a23;
		A.H[2][0]=a31;A.H[2][1]=a32;A.H[2][2]=a33;
		
		res=inversion_mat3d(&A,&invA);
		if (res!=-1) 
			{
				transfores->param[0]=invA.H[0][0];
				transfores->param[1]=invA.H[0][1];
				transfores->param[2]=invA.H[0][2];
				transfores->param[3]=invA.H[1][0];
				transfores->param[4]=invA.H[1][1];
				transfores->param[5]=invA.H[1][2];
				transfores->param[6]=invA.H[2][0];
				transfores->param[7]=invA.H[2][1];
				transfores->param[8]=invA.H[2][2];
		
				transfores->param[12]=transfo->param[12];
				transfores->param[13]=transfo->param[13];
				transfores->param[14]=transfo->param[14];
		
				transfores->param[9]=-(invA.H[0][0]*tx+invA.H[0][1]*ty+invA.H[0][2]*tz);
				transfores->param[10]=-(invA.H[1][0]*tx+invA.H[1][1]*ty+invA.H[1][2]*tz);
				transfores->param[11]=-(invA.H[2][0]*tx+invA.H[2][1]*ty+invA.H[2][2]*tz);
			
				res=1;
			}
			else
			{
			  PUT_ERROR("La transformation affine n'est pas inversible !");
			}
		

	
	
return(res);
}



/*******************************************************************************
** 	inv_bspline_3d(transfo,chres,prec)                                   
**                                                                   
**	Cette fonction permet de calculer la transfo inverse d'une transfo Bspline       
*******************************************************************************/
int inv_bspline_3d(transf3d *transfo, field3d *chres,double prec)
{

int wdth,hght,dpth,i,j,k,ix,iy,iz;
double *p;
int nb_param;
scal_func_t scal_func;
field3d *ch,*chmin,*chmax,*ch_tmp;
vector3d ***datares,***data;
dvector3d htmp;
float ***distance;
int tdebut,tfin/*,err*/;
int u,v,w,Ntot,Ntot2,N;
double xm,xM,ym,yM,zm,zM,dist,max,min,moyenne; 

PRECISION_NUMERIQUE_INV_BSPLINE=prec*0.001;
INV_BSPLINE_RESOL=transfo->resol;
INV_BSPLINE_COEFF2L=pow(2.0,INV_BSPLINE_RESOL);

datares=chres->raw;

wdth=transfo->width;
hght=transfo->height;
dpth=transfo->depth;

ch=cr_field3d(wdth,hght,dpth);
ch_tmp=cr_field3d(wdth+1,hght+1,dpth+1);
chmin=cr_field3d(wdth,hght,dpth);
chmax=cr_field3d(wdth,hght,dpth);
distance=alloc_matrix_3d(wdth,hght,dpth);


data=ch_tmp->raw;
			
	tdebut=time(NULL);


		
p = CALLOC(transfo->nb_param,double);  

if (p==NULL)
	{
	printf("Problème d'allocation dans inv_bspline_3d \n");
	exit(-1);
	}

scal_func=Bsplined1;


//nb_param = init_base_3d(wdth, hght, dpth, scal_func);

//		for (i = BASE3D.resol; i<transfo->resol; i++)
//		  nb_param = base_resol_up_3d(p,nb_param);
//		if (nb_param != transfo->nb_param)
//		  printf("ERREUR dans transf_to_field_3d\n");

nb_param=transfo->nb_param;


/* conversion des param de mm en voxel */
 for (i=0;i<nb_param/3;i++)
  		{
   		p[3*i]=transfo->param[3*i]/transfo->dx;p[3*i+1]=transfo->param[3*i+1]/transfo->dy;p[3*i+2]=transfo->param[3*i+2]/transfo->dz;
			}


/* initialisation de la transoformation inverse en recherchant le plus proche voisin */

//base_to_field_3d(transfo->nb_param, p, ch,NULL,NULL);

ParamBspline_to_field_3d(transfo->nb_param, p, ch,NULL,NULL);

 for (i=0;i<=wdth;i++)
   for (j=0;j<=hght;j++)
    for (k=0;k<=dpth;k++)
			{
			if ((i!=wdth)&&(j!=hght)&&(k!=dpth))
				data[i][j][k]=ch->raw[i][j][k];
			else
				{data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;
				}
			}
			
			
 for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
			{distance[i][j][k]=HUGE_VAL;
				chmin->raw[i][j][k].x=chmin->raw[i][j][k].y=chmin->raw[i][j][k].z=HUGE_VAL;
				chmax->raw[i][j][k].x=chmax->raw[i][j][k].y=chmax->raw[i][j][k].z=-HUGE_VAL;		
			}
 
 
  for (i=0;i<=wdth;i++)
   for (j=0;j<=hght;j++)
    for (k=0;k<=dpth;k++)
				{
				data[i][j][k].x=i+data[i][j][k].x;
				data[i][j][k].y=j+data[i][j][k].y;
				data[i][j][k].z=k+data[i][j][k].z;
				}
							
	for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
				{
				ch->raw[i][j][k].x=ch->raw[i][j][k].x+i;
				ch->raw[i][j][k].y=ch->raw[i][j][k].y+j;
				ch->raw[i][j][k].z=ch->raw[i][j][k].z+k;
				}
				
 for (i=1;i<=wdth;i++)
   for (j=1;j<=hght;j++)
    for (k=1;k<=dpth;k++)
			{
			
			
				
			xm=HUGE_VAL;xM=-HUGE_VAL;
			ym=HUGE_VAL;yM=-HUGE_VAL;
			zm=HUGE_VAL;zM=-HUGE_VAL;
			
			for (u=0;u<2;u++)
  	 	for (v=0;v<2;v++)
  	 	for (w=0;w<2;w++)
  	 		{
				xm=MINI(data[i-u][j-v][k-w].x,xm);
				ym=MINI(data[i-u][j-v][k-w].y,ym);
				zm=MINI(data[i-u][j-v][k-w].z,zm);
			
				xM=MAXI(data[i-u][j-v][k-w].x,xM);
				yM=MAXI(data[i-u][j-v][k-w].y,yM);
				zM=MAXI(data[i-u][j-v][k-w].z,zM);
				}
			
	
			xm=ceil(xm);
	
			ym=ceil(ym);
			
			zm=ceil(zm);
			
			xM=floor(xM)+1;
			yM=floor(yM)+1;
			zM=floor(zM)+1;
			
			xM=MINI(xM,wdth);
			yM=MINI(yM,hght);
			zM=MINI(zM,dpth);
			
			xm=MAXI(xm,0);
			ym=MAXI(ym,0);
			zm=MAXI(zm,0);

		
			
				for (ix=xm;ix<xM;ix++)
  	 		for (iy=ym;iy<yM;iy++)
  	 		for (iz=zm;iz<zM;iz++)
  				{
					
						for (u=0;u<2;u++)
  	 				for (v=0;v<2;v++)
  	 				for (w=0;w<2;w++)
  						{
							dist=sqrt((ix-data[i-u][j-v][k-w].x)*(ix-data[i-u][j-v][k-w].x)+(iy-data[i-u][j-v][k-w].y)*(iy-data[i-u][j-v][k-w].y)+(iz-data[i-u][j-v][k-w].z)*(iz-data[i-u][j-v][k-w].z));
							
							if (dist<distance[ix][iy][iz])
								{
								datares[ix][iy][iz].x=i-u-ix;
								datares[ix][iy][iz].y=j-v-iy;
								datares[ix][iy][iz].z=k-w-iz;
								distance[ix][iy][iz]=dist;
								}
							
							}
					
	
								
					chmin->raw[ix][iy][iz].x=MINI(chmin->raw[ix][iy][iz].x,i-1);
					chmin->raw[ix][iy][iz].y=MINI(chmin->raw[ix][iy][iz].y,j-1);
					chmin->raw[ix][iy][iz].z=MINI(chmin->raw[ix][iy][iz].z,k-1);
					
					chmax->raw[ix][iy][iz].x=MAXI(chmax->raw[ix][iy][iz].x,i);
					chmax->raw[ix][iy][iz].y=MAXI(chmax->raw[ix][iy][iz].y,j);
					chmax->raw[ix][iy][iz].z=MAXI(chmax->raw[ix][iy][iz].z,k);
	
	
					
	
					}
			
		
			}
		free_field3d(ch_tmp);
		

NB_BLOC_POINT_FIXE=0;
Ntot2=0;

/* raffinement du champ inverse */
for (i=0;i<wdth;i++)
	{
	//printf("Bloc %d sur wdth \r",i);
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
				if(sqrt(distance[i][j][k]>prec)&&(0.5*sqrt(3)>prec))
					{
					distance[i][j][k]=raffinement_bspline_inv_3d(p,nb_param,ch,chres,chmin,chmax,i,j,k,prec);
					Ntot2++;
					}

	}

/* verification a posteriori de la précision */
for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
		{
		
		xm=i+datares[i][j][k].x;
		ym=j+datares[i][j][k].y;
		zm=k+datares[i][j][k].z;
		
			
		eval_transfo_bspline1_3d(p, nb_param, chres->width, chres->height, chres->depth, xm, ym, zm, &htmp);
	
		distance[i][j][k]=sqrt( (htmp.x-i)*(htmp.x-i)+(htmp.y-j)*(htmp.y-j)+(htmp.z-k)*(htmp.z-k) );
		
			
		}

/* conversion du champ en mm */

for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
				{
				datares[i][j][k].x=datares[i][j][k].x*transfo->dx;
				datares[i][j][k].y=datares[i][j][k].y*transfo->dy;
				datares[i][j][k].z=datares[i][j][k].z*transfo->dz;
				}
		


  tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }

/* estimation de la précision de la transfo inverse */

max=0;
Ntot=0;N=0;
min=HUGE_VAL;
moyenne=0;
for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
				if(distance[i][j][k]>0.0)
				{
				if (isinf(distance[i][j][k]))
					printf("HUGE_VAL en %d %d %d\n",i,j,k);
				
				max=MAXI(distance[i][j][k],max);
				min=MINI(distance[i][j][k],min);
			
				Ntot++;
				
				if(distance[i][j][k]>prec)
						printf("%d %d  %d  distance = %f \n",i,j,k,distance[i][j][k]);
				else
					N++;
					
				moyenne+=distance[i][j][k];
				}

moyenne = moyenne/Ntot;

printf("Erreur max = %f \t Erreur min = %f \t Erreur moyenne %f \n",max,min,moyenne);
printf("Pourcentage de voxel avec une précision inferieur à %f :  %f pourcent \n",prec,100.0*N/Ntot);
printf("Nombre de voxel : %d \t Pourcentage de voxel pour lesquels le point fixe n'a pas fonctionne :  %f pourcent \n",NB_BLOC_POINT_FIXE,100.0*(double)NB_BLOC_POINT_FIXE/(double)Ntot2);


    

//		end_base_3d();
		free(p);
		free_matrix_3d(distance);
		free_field3d(ch);
		free_field3d(chmin);
		free_field3d(chmax);


		
	 	return(1);
}




/*******************************************************************************
** 	inv_field_3d(ch,chres)                                   
**                                                                   
**	Cette fonction permet de calculer la transfo inverse d'un champ dense       
*******************************************************************************/
int inv_field_3d(transf3d *transfo_field3d,field3d *chres, double prec)
{
transf3d* transfo_base3d;
field3d* chtmp;
int i,j,k;

// Conversion d'un champ dense en une transfo Bspline
transfo_base3d=ConvertTransfoField3dToBase3d(transfo_field3d);

//ProjectBase3dintoTopologyPreservingTransformation(transfo_base3d );


chtmp=cr_field3d(transfo_base3d->width,transfo_base3d->height,transfo_base3d->depth);

// Inversion de la transfo Bspline
inv_bspline_3d(transfo_base3d, chtmp, prec);


// Restriction du champ 
		for (i=0;i<transfo_field3d->width;i++)
		for (j=0;j<transfo_field3d->height;j++)
		for (k=0;k<transfo_field3d->depth;k++)
			{
			chres->raw[i][j][k].x=chtmp->raw[i][j][k].x;
			chres->raw[i][j][k].y=chtmp->raw[i][j][k].y;
			chres->raw[i][j][k].z=chtmp->raw[i][j][k].z;
			}
			
free_transf3d(transfo_base3d);
free_field3d(chtmp);

return(1);
}


/*******************************************************************************
** 	transfo* ConvertTransfoField3dToBase3d(ch);                                   
**                                                                   
**	Converti une transfo field3d en une transfo Bspline de degré 1  
** transfo est alloué dans la fonction
*******************************************************************************/
transf3d* ConvertTransfoField3dToBase3d(transf3d *transfo_field3d )
{
int wdth, hght, dpth, size_max, resol,nb_param,i,j,k,l,u,v,w;
double *param;
field3d *ch;
transf3d* transfo_base3d;

wdth=transfo_field3d->width;
hght=transfo_field3d->height;
dpth=transfo_field3d->depth;

size_max=MAXI(wdth,hght);
size_max=MAXI(size_max, dpth);

ch=transf_to_field_3d(transfo_field3d,NULL,NULL);

/*conversion en champ voxelique*/
		/*for (i=0;i<wdth;i++)
		for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			ch->raw[i][j][k].x=ch->raw[i][j][k].x/transfo_field3d->dx;
			ch->raw[i][j][k].y=ch->raw[i][j][k].y/transfo_field3d->dy;
			ch->raw[i][j][k].z=ch->raw[i][j][k].z/transfo_field3d->dz;
			}*/


resol=(int)ceil(log(size_max)/log(2.0));

size_max=pow(2.0,resol);
nb_param=pow(2.0,resol)-1;
nb_param=(int)3.0*pow(nb_param,3.0);

transfo_base3d=cr_transf3d(size_max,size_max,size_max,NULL);
transfo_base3d->nb_param=nb_param;
transfo_base3d->param=CALLOC(nb_param,double);
transfo_base3d->typetrans=BSPLINE3D;
transfo_base3d->resol=resol;
transfo_base3d->degre=1;
transfo_base3d->dx=transfo_field3d->dx;
transfo_base3d->dy=transfo_field3d->dy;
transfo_base3d->dz=transfo_field3d->dz;

param=transfo_base3d->param;

 for (i=0;i<nb_param;i++)
 		 param[i]=0.0;
		 
	for (i=1;i<wdth;i++)
	for (j=1;j<hght;j++)
	for (k=1;k<dpth;k++)
		{
		u=i-1;v=j-1;w=k-1;
		l=TOP_conv_ind(u,v,w,nb_param);
		param[l]=ch->raw[i][j][k].x;
		param[l+1]=ch->raw[i][j][k].y;
		param[l+2]=ch->raw[i][j][k].z;
		}			 

free_field3d(ch);
return(transfo_base3d);	 
}

/*******************************************************************************
** 	ProjectBase3dintoTopologyPreservingTransformation(ch);                                   
**                                                                   
**	Modifie les parametres d'une transfo Bspline pour qu'elle préserve la topologie
*******************************************************************************/
int ProjectBase3dintoTopologyPreservingTransformation(transf3d *transfo_base3d )
{
int wdth, hght, dpth, resol,nb_param,i,j,k,l,u,v,w,topD;
double *param, *p_norm, *p, *grad;
double Jmin=0.0, Jmax=10000000;
int *x00,*x11,*y00,*y11,*z00,*z11;
double xm,xM,ym,yM,zm,zM;
TSlpqr Slpqr[8];
scal_func_t scal_func=Bsplined1;
int nb_bloc, Ntot;
double cfmax,erreur_max,tmp;
int ***masque_param;
int tdebut,tfin;
 
tdebut=time(NULL);

if ((transfo_base3d->typetrans!=BSPLINE3D)||(transfo_base3d->degre!=1))
	return(0);
	
wdth=transfo_base3d->width;
hght=transfo_base3d->height;
dpth=transfo_base3d->depth;
resol=transfo_base3d->resol;
resol=transfo_base3d->resol;
nb_param=transfo_base3d->nb_param;
param=transfo_base3d->param;
topD=(int)pow(2.0,1.0*resol)-1;

p_norm=CALLOC(nb_param,double);
p=CALLOC(nb_param,double);
grad=CALLOC(nb_param,double);

masque_param=alloc_imatrix_3d(topD,topD,topD);
	
// Initialisation de Based3D
nb_param=init_base_3d(wdth,hght,dpth,scal_func);

for (l=1; l<resol; l++)	
 nb_param=base_resol_up_light_3d(p_norm,nb_param);

if (nb_param!=transfo_base3d->nb_param)
	printf("Probleme dans le nombre de parametre\n");


x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;

for (l=0; l<nb_param; l++)
	p[l]= param[l];

	 
for (l=0; l<nb_param/3; l++) 
	{		
	p_norm[3*l] =1.0*param[3*l]/wdth;
	p_norm[3*l+1]=1.0*param[3*l+1]/hght;
	p_norm[3*l+2]=1.0*param[3*l+2]/dpth;
	}

		 

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
	{
	Slpqr[i].Jm=Jmin;
	Slpqr[i].JM=Jmax;
	} 



	for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
	for (k=0;k<topD;k++)
		masque_param[i][j][k]=0;


Ntot=topD*topD*topD;

/* On met a zero les parametres qui ne preserve pas la topologie */
nb_bloc =1;
while (nb_bloc>0)
{
nb_bloc =0;

	for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
	for (k=0;k<topD;k++)
		if (masque_param[i][j][k]==0)
		{
		
		 	
		//---------------------------------------------------------------
  		//----- initialisation Slpqr            -------------------------
  		//---------------------------------------------------------------
		l=TOP_conv_ind(i,j,k,nb_param);
  		
		printf("bloc %d sur %d \r",l, Ntot);
		xm = 1.0*x00[l/3]/wdth; xM = 1.0*x11[l/3]/wdth; ym = 1.0*y00[l/3]/hght; 
  		yM = 1.0*y11[l/3]/hght; zm = 1.0*z00[l/3]/dpth; zM =1.0*z11[l/3]/dpth;

  		Slpqr[0].b.xm = xm; Slpqr[0].b.xM = (xm+xM)/2;
  		Slpqr[1].b.xm = xm; Slpqr[1].b.xM = (xm+xM)/2;
  		Slpqr[2].b.xm = xm; Slpqr[2].b.xM = (xm+xM)/2;
  		Slpqr[3].b.xm = xm; Slpqr[3].b.xM = (xm+xM)/2;
  		Slpqr[4].b.xM = xM; Slpqr[4].b.xm = (xm+xM)/2;
  		Slpqr[5].b.xM = xM; Slpqr[5].b.xm = (xm+xM)/2;
  		Slpqr[6].b.xM = xM; Slpqr[6].b.xm = (xm+xM)/2;
  		Slpqr[7].b.xM = xM; Slpqr[7].b.xm = (xm+xM)/2;

 		Slpqr[0].b.ym = ym; Slpqr[0].b.yM = (ym+yM)/2;
  		Slpqr[1].b.ym = ym; Slpqr[1].b.yM = (ym+yM)/2;
  		Slpqr[4].b.ym = ym; Slpqr[4].b.yM = (ym+yM)/2;
  		Slpqr[5].b.ym = ym; Slpqr[5].b.yM = (ym+yM)/2;
  		Slpqr[2].b.yM = yM; Slpqr[2].b.ym = (ym+yM)/2;
  		Slpqr[3].b.yM = yM; Slpqr[3].b.ym = (ym+yM)/2;
  		Slpqr[6].b.yM = yM; Slpqr[6].b.ym = (ym+yM)/2;
  		Slpqr[7].b.yM = yM; Slpqr[7].b.ym = (ym+yM)/2;

  		Slpqr[0].b.zm = zm; Slpqr[0].b.zM = (zm+zM)/2;
  		Slpqr[2].b.zm = zm; Slpqr[2].b.zM = (zm+zM)/2;
  		Slpqr[4].b.zm = zm; Slpqr[4].b.zM = (zm+zM)/2;
  		Slpqr[6].b.zm = zm; Slpqr[6].b.zM = (zm+zM)/2;
  		Slpqr[1].b.zM = zM; Slpqr[1].b.zm = (zm+zM)/2;
  		Slpqr[3].b.zM = zM; Slpqr[3].b.zm = (zm+zM)/2;
  		Slpqr[5].b.zM = zM; Slpqr[5].b.zm = (zm+zM)/2;
  		Slpqr[7].b.zM = zM; Slpqr[7].b.zm = (zm+zM)/2;

	
		
		/*cfmax = TOP_linemin_maxpas_3d(nb_param,p_norm,grad,i,j,k,Jmin,Jmax,Slpqr);

		
		if (cfmax>=1)
			{
			p_norm[l]=grad[l];
			p_norm[l+1]=grad[l+1];
			p_norm[l+2]=grad[l+2];
			p[l]  =p_norm[l]*wdth  ;
			p[l+1]=p_norm[l+1]*hght;
			p[l+2]=p_norm[l+2]*dpth;
			
			masque_param[l]=masque_param[l+1]=masque_param[l+2]=1;
			nb_bloc++;
			}
		else
			{
			if (cfmax>0.01)
				{
				p_norm[l]=grad[l]*cfmax;
				p_norm[l+1]=grad[l+1]*cfmax;
				p_norm[l+2]=grad[l+2]*cfmax;   
				p[l]  =p_norm[l]*wdth  ;
				p[l+1]=p_norm[l+1]*hght;
				p[l+2]=p_norm[l+2]*dpth;
				
				nb_bloc++;
				}
			}
		
		*/
		
		
		/*if (TOP_verification_locale(nb_param,param,p_norm,Jmin,Jmax,i,j,k,Slpqr)==0)
			{
			nb_bloc++; 
			
			px=param[l];
			py=param[l+1];
			pz=param[l+2];
			u=1;
			while ((TOP_verification_locale(nb_param,param,p_norm,Jmin,Jmax,i,j,k,Slpqr)==0)&&(u<=nb_pas))
				{
				
				alpha=(1.0-(double)u/(double)nb_pas);
				param[l]=px*alpha;
				param[l+1]=py*alpha;
				param[l+2]=pz*alpha;
				
				p_norm[l]  =param[l]/wdth  ;
				p_norm[l+1]=param[l+1]/hght;
				p_norm[l+2]=param[l+2]/dpth;
				
				u++;
				}
				
			if (TOP_verification_locale(nb_param,param,p_norm,Jmin,Jmax,i,j,k,Slpqr)==0)	
				printf("Le probleme n'a pas ete resolu sur le bloc %d %d %d \n",i,j,k);
			else
				printf("probleme a ete resolu sur le bloc %d %d %d \n",i,j,k);
				
			
			}*/
	
		if (TOP_verification_locale(nb_param,param,p_norm,Jmin,Jmax,i,j,k,Slpqr)==0)
			{
			nb_bloc++; 
			param[l]=param[l+1]=param[l+2]=	p_norm[l]=p_norm[l+1]=p_norm[l+2]=0.0;
			
			for (u=MAXI(0,i-1);u<MINI(i+2,topD);u++)
			for (v=MAXI(0,j-1);v<MINI(j+2,topD);v++)
			for (w=MAXI(0,k-1);w<MINI(k+2,topD);w++)
				masque_param[u][v][w]=0;
			
			
			}
		else
			{
			masque_param[i][j][k]=1;
			}
	
		}			 

printf("nombre de bloc ne conservant pas la topologie : %d \n",nb_bloc);

}


/* on met a 0 les elements de masque param ou il faut optimiser */

	for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
	for (k=0;k<topD;k++)
		{
		l=TOP_conv_ind(i,j,k,nb_param);
  		
		if ((param[l]!=p[l])||(param[l+1]!=p[l+1])||(param[l+2]!=p[l+2]))
			masque_param[i][j][k]=0;
		
		}
		
		
/* on estime les parametres les plus proches possible de la solution */

nb_bloc =1;
while (nb_bloc>0)
{
nb_bloc =0;

	for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
	for (k=0;k<topD;k++)
		if (masque_param[i][j][k]==0)
		{
		
		 	
		//---------------------------------------------------------------
  		//----- initialisation Slpqr            -------------------------
  		//---------------------------------------------------------------
		l=TOP_conv_ind(i,j,k,nb_param);
  		
		printf("bloc %d sur %d \r",l, Ntot);
		xm = 1.0*x00[l/3]/wdth; xM = 1.0*x11[l/3]/wdth; ym = 1.0*y00[l/3]/hght; 
  		yM = 1.0*y11[l/3]/hght; zm = 1.0*z00[l/3]/dpth; zM =1.0*z11[l/3]/dpth;

  		Slpqr[0].b.xm = xm; Slpqr[0].b.xM = (xm+xM)/2;
  		Slpqr[1].b.xm = xm; Slpqr[1].b.xM = (xm+xM)/2;
  		Slpqr[2].b.xm = xm; Slpqr[2].b.xM = (xm+xM)/2;
  		Slpqr[3].b.xm = xm; Slpqr[3].b.xM = (xm+xM)/2;
  		Slpqr[4].b.xM = xM; Slpqr[4].b.xm = (xm+xM)/2;
  		Slpqr[5].b.xM = xM; Slpqr[5].b.xm = (xm+xM)/2;
  		Slpqr[6].b.xM = xM; Slpqr[6].b.xm = (xm+xM)/2;
  		Slpqr[7].b.xM = xM; Slpqr[7].b.xm = (xm+xM)/2;

 		Slpqr[0].b.ym = ym; Slpqr[0].b.yM = (ym+yM)/2;
  		Slpqr[1].b.ym = ym; Slpqr[1].b.yM = (ym+yM)/2;
  		Slpqr[4].b.ym = ym; Slpqr[4].b.yM = (ym+yM)/2;
  		Slpqr[5].b.ym = ym; Slpqr[5].b.yM = (ym+yM)/2;
  		Slpqr[2].b.yM = yM; Slpqr[2].b.ym = (ym+yM)/2;
  		Slpqr[3].b.yM = yM; Slpqr[3].b.ym = (ym+yM)/2;
  		Slpqr[6].b.yM = yM; Slpqr[6].b.ym = (ym+yM)/2;
  		Slpqr[7].b.yM = yM; Slpqr[7].b.ym = (ym+yM)/2;

  		Slpqr[0].b.zm = zm; Slpqr[0].b.zM = (zm+zM)/2;
  		Slpqr[2].b.zm = zm; Slpqr[2].b.zM = (zm+zM)/2;
  		Slpqr[4].b.zm = zm; Slpqr[4].b.zM = (zm+zM)/2;
  		Slpqr[6].b.zm = zm; Slpqr[6].b.zM = (zm+zM)/2;
  		Slpqr[1].b.zM = zM; Slpqr[1].b.zm = (zm+zM)/2;
  		Slpqr[3].b.zM = zM; Slpqr[3].b.zm = (zm+zM)/2;
  		Slpqr[5].b.zM = zM; Slpqr[5].b.zm = (zm+zM)/2;
  		Slpqr[7].b.zM = zM; Slpqr[7].b.zm = (zm+zM)/2;

	
		grad[l]=(param[l]-p[l])/wdth;
		grad[l+1]=(param[l+1]-p[l+1])/hght;
		grad[l+2]=(param[l+2]-p[l+2])/dpth;
		
		cfmax = TOP_linemin_maxpas_3d(nb_param,p_norm,grad,i,j,k,Jmin,Jmax,Slpqr);

		
		if (cfmax>=1)
			{
			p_norm[l]=p_norm[l]-grad[l];
			p_norm[l+1]=p_norm[l+1]-grad[l+1];
			p_norm[l+2]=p_norm[l+2]-grad[l+2];
			
			param[l]  =p_norm[l]*wdth  ;
			param[l+1]=p_norm[l+1]*hght;
			param[l+2]=p_norm[l+2]*dpth;
			
			masque_param[i][j][k]=1;
			nb_bloc++;
			}
		else
			{
			if (cfmax>0.01)
				{
				p_norm[l]=p_norm[l]-0.99*grad[l]*cfmax;
				p_norm[l+1]=p_norm[l+1]-0.99*grad[l+1]*cfmax;
				p_norm[l+2]=p_norm[l+2]-0.99*grad[l+2]*cfmax;   
				
				param[l]  =p_norm[l]*wdth  ;
				param[l+1]=p_norm[l+1]*hght;
				param[l+2]=p_norm[l+2]*dpth;
				
				nb_bloc++;
				
				if (sqrt((param[l]-p[l])*(param[l]-p[l])+(param[l+1]-p[l+1])*(param[l+1]-p[l+1])+(param[l+2]-p[l+2])*(param[l+2]-p[l+2]))<0.001)
					masque_param[i][j][k]=1;
				
				}
			}
		
		
	
		}			 

printf("nombre de bloc ayant ete traite avec succes  : %d \n",nb_bloc);

}

erreur_max=0.0;

	for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
	for (k=0;k<topD;k++)
		{
		tmp=sqrt((param[l]-p[l])*(param[l]-p[l])+(param[l+1]-p[l+1])*(param[l+1]-p[l+1])+(param[l+2]-p[l+2])*(param[l+2]-p[l+2]));
		
		if (tmp>erreur_max)
			erreur_max=tmp;
		}



  tfin=time(NULL);
    {
     int h,m,s,ttotal;
     ttotal=tfin-tdebut;
     h=ttotal/3600;
     m=(ttotal-h*3600)/60;
     s=(ttotal-h*3600-m*60);
     aff_log("TEMPS PARTIEL: %dh %dm %ds \n",h,m,s);
    }
    
printf("Erreur max lors de la projection : %f \n",erreur_max);

end_base_3d();

free(p_norm);free(p);free(grad);free_imatrix_3d(masque_param);
	
return(1);	 
}

/*******************************************************************************
** 	deform_anim_3d(void)                                   
**                                                                   
**	Cette fonction permet de generer une sequence d'image illustrant l'evolution de la deformation       
*******************************************************************************/
void deform_anim_3d(void) 
{
int im_1,nb,l,i,j,k;
grphic3d *im1,*imres;
char  nomfich1[500],nomfichres[500];
transf3d *transfo1;
int	e=0,wdth,hght,dpth,inter_type;
field3d* champ,*champ_tmp;
InterpolationFct inter_fct;
grphic3d *maskdeb,*maskres;
 
im_1=GET_PLACE3D(TEXT0224);

 /*lecture du premier champ*/
 strcpy(nomfich1,GET_FILE("fichier trf",&e));
 if(e != 0)
     return ;
 put_file_extension(nomfich1,".trf",nomfich1);

  if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}

inter_type=imx_query_interpolation_type_3D(0);

 	nb= GET_INT("nombre d'images", 10, &e);

 
strcpy(nomfichres,GET_FILE(TEXT0420,&e));
	put_file_extension(nomfichres,".ipb",nomfichres);
	
	if(nomfichres == NULL || strlen(nomfichres) <= 0 || e) 
	  PUT_WARN( TEXT0421);


im1=ptr_img_3d(im_1);
imres=cr_grphic3d(im1);imx_copie_param_3d_p(im1,imres);

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

champ=cr_field3d(wdth,hght,dpth);
champ_tmp=cr_field3d(wdth,hght,dpth);

inter_fct = imx_choose_interpolation_fct(inter_type);	
/*lecture des champ*/

transfo1=load_transf_3d(nomfich1);
champ=transf_to_field_3d(transfo1,NULL,im1);

save_mri_ipb_3d_p(nomfichres,im1);

	 
  
maskdeb=cr_grphic3d(im1);imx_copie_3d_p(im1,maskdeb);
maskres=cr_grphic3d(imres);imx_copie_3d_p(imres,maskres);  
  
maskdeb->max_pixel=1;
maskdeb->min_pixel=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	if (im1->mri[i][j][k]!=0) maskdeb->mri[i][j][k]=1;
	}

for (l=1;l<=nb;l++)
{
	for(i=0;i<wdth;i++)
	for(j=0;j<hght;j++)
	for(k=0;k<dpth;k++)
		{
		champ_tmp->raw[i][j][k].x=1.0*l*champ->raw[i][j][k].x/nb;
		champ_tmp->raw[i][j][k].y=1.0*l*champ->raw[i][j][k].y/nb;
		champ_tmp->raw[i][j][k].z=1.0*l*champ->raw[i][j][k].z/nb;
		}
inter_fct(im1,champ_tmp,imres); 
inter_labeled_3d(maskdeb,champ_tmp,maskres);

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
		if (maskres->mri[i][j][k]==0)
			imres->mri[i][j][k]=0;
			
save_mri_ipb_3d_p(nomfichres,imres);
}

free_transf3d(transfo1);free_grphic3d(imres);free_field3d(champ);free_field3d(champ_tmp);
free_grphic3d(maskdeb);free_grphic3d(maskres);
}


/******************************************************************************
 **    filtre_champ_diffusion_anisotrope_3d() :
 **
 */
/*! filtrage de diffusion anisotrope tel que preconise dans
**  "Magnetic Resonance Image Tissue Claissification Using a Partial Volume Model"
**  Shattuck et al.
**  avec kd=5.0 et t0=1/8
**
************************************************************************/
void filtre_champ_diffusion_anisotrope_3d()
{
  int im_src, nb_iter=0;
  int err=0;
	char  nomfich1[500],nomfichres[500],str[500];
	double Kd;
	
 	/*lecture du champ*/
 	strcpy(nomfich1,GET_FILE("fichier trf a Filtrer",&err));
 	if(err != 0)
     return ;
 	put_file_extension(nomfich1,".trf",nomfich1);

  if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}


  //choix des images
  im_src=GET_PLACE3D("Image Source");
 
  //choix du nombre d'iteration
  nb_iter=GET_INT("Nombre d'iterations", 0, &err);

  //choix du nombre d'iteration
  Kd=GET_FLOAT("Kd", 10, &err);


		/*nom du fichier resultat*/
 		sprintf(str,"File ? in %s",_CurrentPath);
 		strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 		put_file_extension(nomfichres,".trf",nomfichres);

  imx_filtre_champ_diffusion_anisotrope_3d(im_src, nomfich1, nomfichres, nb_iter,Kd);

}




/******************************************************************************
 **    imx_filtre_champ_diffusion_anisotrope_3d() :
 **
 */
/*! filtrage de diffusion anisotrope tel que preconise dans
**  "Magnetic Resonance Image Tissue Claissification Using a Partial Volume Model"
**  Shattuck et al.
**  avec kd=5.0 et t0=1/8
**
**  \param im_src :image source
**  \param im_dst : image resultat (E/S)
**  \param taille_filtre : nb_iter : nombre d'iterations
**
****************************************************************************/
int imx_filtre_champ_diffusion_anisotrope_3d(int im_src, char* nomfich1, char* nomfichres, int nb_iter,double Kd)
{
  grphic3d *imsrc=NULL;
  int res=0,wdth,hght,dpth;
	transf3d *transfo,*transfores;
	field3d *ch,*chres;

  imsrc=ptr_img_3d(im_src);
 
  /*lecture des champ*/
    transfo=load_transf_3d(nomfich1);
		
		wdth=transfo->width;
		hght=transfo->height;
		dpth=transfo->depth;
		

    transfores=cr_transf3d(wdth,hght,dpth,NULL);
		
		ch=transf_to_field_3d(transfo,NULL,NULL);
		chres=cr_field3d(wdth,hght,dpth);

  	res=imx_filtre_champ_diffusion_anisotrope_3d_p(imsrc,ch,chres,nb_iter,Kd);

		transfores=field_to_transf_3d(chres,NULL,NULL);
	  transfores->typetrans=CHAMP3D;
		transfores->dx=transfo->dx;
		transfores->dy=transfo->dy;
		transfores->dz=transfo->dz;
	
		free_field3d(chres);
		free_field3d(ch);
    save_transf_3d(transfores,nomfichres);
    free_transf3d(transfo); free_transf3d(transfores);

  return res;
}

/******************************************************************************
 **    imx_filtre_champ_diffusion_anisotrope_3d_p() :
 **
 */
/*! filtrage de diffusion anisotrope tel que preconise dans
**  "Magnetic Resonance Image Tissue Claissification Using a Partial Volume Model"
**  Shattuck et al.
**  avec kd=5.0 et t0=1/8
**
**  \param imsrc :image source
**  \param imdst : image resultat (E/S)
**  \param nb_iter : nombre d'iterations
**
****************************************************************************/
int imx_filtre_champ_diffusion_anisotrope_3d_p(grphic3d *imsrc,field3d* ch,field3d* chres,int nb_iter,double kd)
{
 unsigned int i,j,k;
 TDimension width, height, depth;
 double *tmp_res_n0x=NULL, *tmp_res_n1x=NULL,*tmp_res_n0y=NULL, *tmp_res_n1y=NULL,*tmp_res_n0z=NULL, *tmp_res_n1z=NULL,*imagen0=NULL,*imagen1=NULL,*norm_grad=NULL;
 int iter, res=0;
 double /*kd=5.0, */t0=1.0/8.0, val, term_f=0.0,/* max, iVal,*/ rcoeff;

 if ((!imsrc))
  { fprintf(stderr,"les deux images doivent etre non nulles dans imx_filtre_diffusion_anisotrope_3d_p\n"); return 1; }

 if ((imsrc->width!=ch->width)||(imsrc->height!=ch->height)||(imsrc->depth!=ch->depth))
  { fprintf(stderr,"Le champ et l'image doivent avoir la meme taille\n"); return 1; }

 width=imsrc->width; height=imsrc->height; depth=imsrc->depth;

printf("Filtrage par diffusion anisotrope : Nbiter : %d  Kd %f \n",nb_iter,kd);
 //allocations memoire
 tmp_res_n0x=CALLOC(width*height*depth, double);
 tmp_res_n1x=CALLOC(width*height*depth, double);
 tmp_res_n0y=CALLOC(width*height*depth, double);
 tmp_res_n1y=CALLOC(width*height*depth, double);
 tmp_res_n0z=CALLOC(width*height*depth, double);
 tmp_res_n1z=CALLOC(width*height*depth, double);
 
 imagen0=CALLOC(width*height*depth, double);
 imagen1=CALLOC(width*height*depth, double);
 
 norm_grad=CALLOC(width*height*depth, double);

 //on remplit avec les valeurs du champ
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    {
		tmp_res_n0x[i*height*depth+j*depth+k]=ch->raw[i][j][k].x;
		tmp_res_n0y[i*height*depth+j*depth+k]=ch->raw[i][j][k].y;
		tmp_res_n0z[i*height*depth+j*depth+k]=ch->raw[i][j][k].z;
		}

 //on remplit l'image double
 rcoeff=imsrc->rcoeff;
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    imagen0[i*height*depth+j*depth+k]=(double)rcoeff*imsrc->mri[i][j][k];
	
 //filtrage de l'image
 for (iter=0;iter<nb_iter;iter++)
 {
  //calcul de la norme du gradient de l'image
  res=imx_norme_sqr_quad_gradient_3d_p(imagen0, norm_grad, width, height, depth);

  //calcul de l'image des coefficients de diffusion
  for (i=0;i<width;i++)
   for (j=0;j<height;j++)
    for (k=0;k<depth;k++)
    {
     val=exp(-norm_grad[i*height*depth+j*depth+k]/(kd*kd));
     norm_grad[i*height*depth+j*depth+k]=val;
    }

  //calcul de l'image filtree
  for (i=0;i<width;i++)
   for (j=0;j<height;j++)
    for (k=0;k<depth;k++)
    {
     	///// suivant x
			
			//calcul du terme de diffusion
     	term_f=0.0;
     	val=tmp_res_n0x[i*height*depth+j*depth+k];
     	if ((i>0)&&(i<(width-1)))
     	{
     	 term_f+=norm_grad[(i-1)*height*depth+j*depth+k]*(tmp_res_n0x[(i+1)*height*depth+j*depth+k]-val);
     	 term_f+=norm_grad[(i+1)*height*depth+j*depth+k]*(tmp_res_n0x[(i-1)*height*depth+j*depth+k]-val);
     	}
     	if ((j>0)&&(j<(height-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+(j-1)*depth+k]*(tmp_res_n0x[i*height*depth+(j+1)*depth+k]-val);
     	 term_f+=norm_grad[i*height*depth+(j+1)*depth+k]*(tmp_res_n0x[i*height*depth+(j-1)*depth+k]-val);
     	}
     	if ((k>0)&&(k<(depth-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+j*depth+k-1]*(tmp_res_n0x[i*height*depth+j*depth+k+1]-val);
     	 term_f+=norm_grad[i*height*depth+j*depth+k+1]*(tmp_res_n0x[i*height*depth+j*depth+k-1]-val);
     	}
     	term_f*=t0;
     	//update de l'image res
     	tmp_res_n1x[i*height*depth+j*depth+k]=val+term_f;
    	
			
			///// suivant y
			
			//calcul du terme de diffusion
     	term_f=0.0;
     	val=tmp_res_n0y[i*height*depth+j*depth+k];
     	if ((i>0)&&(i<(width-1)))
     	{
     	 term_f+=norm_grad[(i-1)*height*depth+j*depth+k]*(tmp_res_n0y[(i+1)*height*depth+j*depth+k]-val);
     	 term_f+=norm_grad[(i+1)*height*depth+j*depth+k]*(tmp_res_n0y[(i-1)*height*depth+j*depth+k]-val);
     	}
     	if ((j>0)&&(j<(height-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+(j-1)*depth+k]*(tmp_res_n0y[i*height*depth+(j+1)*depth+k]-val);
     	 term_f+=norm_grad[i*height*depth+(j+1)*depth+k]*(tmp_res_n0y[i*height*depth+(j-1)*depth+k]-val);
     	}
     	if ((k>0)&&(k<(depth-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+j*depth+k-1]*(tmp_res_n0y[i*height*depth+j*depth+k+1]-val);
     	 term_f+=norm_grad[i*height*depth+j*depth+k+1]*(tmp_res_n0y[i*height*depth+j*depth+k-1]-val);
     	}
     	term_f*=t0;
     	//update de l'image res
     	tmp_res_n1y[i*height*depth+j*depth+k]=val+term_f;
    	
			
			///// suivant z
			
			//calcul du terme de diffusion
     	term_f=0.0;
     	val=tmp_res_n0z[i*height*depth+j*depth+k];
     	if ((i>0)&&(i<(width-1)))
     	{
     	 term_f+=norm_grad[(i-1)*height*depth+j*depth+k]*(tmp_res_n0z[(i+1)*height*depth+j*depth+k]-val);
     	 term_f+=norm_grad[(i+1)*height*depth+j*depth+k]*(tmp_res_n0z[(i-1)*height*depth+j*depth+k]-val);
     	}
     	if ((j>0)&&(j<(height-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+(j-1)*depth+k]*(tmp_res_n0z[i*height*depth+(j+1)*depth+k]-val);
     	 term_f+=norm_grad[i*height*depth+(j+1)*depth+k]*(tmp_res_n0z[i*height*depth+(j-1)*depth+k]-val);
     	}
     	if ((k>0)&&(k<(depth-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+j*depth+k-1]*(tmp_res_n0z[i*height*depth+j*depth+k+1]-val);
     	 term_f+=norm_grad[i*height*depth+j*depth+k+1]*(tmp_res_n0z[i*height*depth+j*depth+k-1]-val);
     	}
     	term_f*=t0;
     	//update de l'image res
     	tmp_res_n1z[i*height*depth+j*depth+k]=val+term_f;
    	
			
			///// sur l'image de reference
			
			//calcul du terme de diffusion
     	term_f=0.0;
     	val=imagen0[i*height*depth+j*depth+k];
     	if ((i>0)&&(i<(width-1)))
     	{
     	 term_f+=norm_grad[(i-1)*height*depth+j*depth+k]*(imagen0[(i+1)*height*depth+j*depth+k]-val);
     	 term_f+=norm_grad[(i+1)*height*depth+j*depth+k]*(imagen0[(i-1)*height*depth+j*depth+k]-val);
     	}
     	if ((j>0)&&(j<(height-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+(j-1)*depth+k]*(imagen0[i*height*depth+(j+1)*depth+k]-val);
     	 term_f+=norm_grad[i*height*depth+(j+1)*depth+k]*(imagen0[i*height*depth+(j-1)*depth+k]-val);
     	}
     	if ((k>0)&&(k<(depth-1)))
     	{
     	 term_f+=norm_grad[i*height*depth+j*depth+k-1]*(imagen0[i*height*depth+j*depth+k+1]-val);
     	 term_f+=norm_grad[i*height*depth+j*depth+k+1]*(imagen0[i*height*depth+j*depth+k-1]-val);
     	}
     	term_f*=t0;
     	//update de l'image res
     	imagen1[i*height*depth+j*depth+k]=val+term_f;
    	
		}

  //on reboucle si on n'a pas atteint le nombre d'iterations voulues
  if (iter<(nb_iter-1))
  {
   for (i=0;i<width;i++)
    for (j=0;j<height;j++)
     for (k=0;k<depth;k++)
      {
       tmp_res_n0x[i*height*depth+j*depth+k]=tmp_res_n1x[i*height*depth+j*depth+k];
       tmp_res_n0y[i*height*depth+j*depth+k]=tmp_res_n1y[i*height*depth+j*depth+k];
       tmp_res_n0z[i*height*depth+j*depth+k]=tmp_res_n1z[i*height*depth+j*depth+k];
       imagen0[i*height*depth+j*depth+k]=imagen1[i*height*depth+j*depth+k];
			}
  }
 }

 
 for (i=0;i<width;i++)
  for (j=0;j<height;j++)
   for (k=0;k<depth;k++)
    {
		chres->raw[i][j][k].x=tmp_res_n1x[i*height*depth+j*depth+k];
		chres->raw[i][j][k].y=tmp_res_n1y[i*height*depth+j*depth+k];
		chres->raw[i][j][k].z=tmp_res_n1z[i*height*depth+j*depth+k];
		}


 FREE(tmp_res_n0x); FREE(tmp_res_n1x); 
 FREE(tmp_res_n0y); FREE(tmp_res_n1y); 
 FREE(tmp_res_n0z); FREE(tmp_res_n1z); 
 
 free(imagen0);
 free(imagen1);
 
 FREE(norm_grad);

 return res;
}

/******************************************************************************
*********
*********  ----------------- filtre_champ_gaussien_3d
*********
********************************************************************************/ 

void filtre_champ_gaussien_3d(void)
{
  int wnd=0;
  int e=0;
  char *quest[5];
  int i;
  int val;
  double sigma;
	char  nomfich1[500],nomfichres[500],str[500];
	
 	/*lecture du champ*/
 	strcpy(nomfich1,GET_FILE("fichier trf a Filtrer",&e));
 	if(e != 0)
     return ;
 	put_file_extension(nomfich1,".trf",nomfich1);

  if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}



	sigma=(double)GET_FLOAT("Sigma ?",1,&e);
  
  for(i=0;i<5;i++)
     quest[i]=CALLOC(80,char);

  strcpy(quest[0],"3 x 3 x 3");
  strcpy(quest[1],"5 x 5 x 5");
  strcpy(quest[2],"7 x 7 x 7");
  strcpy(quest[3],"9 x 9 x 9");
  strcpy(quest[4],"\0"); 

  val=GETV_QCM("Mask ?",(const char **)quest);


  for(i=0;i<5;i++)
     free(quest[i]);

  if(val==0)
    wnd=3;
  else if(val==1)
    wnd=5;
  else if(val==2)
    wnd=7;
  else if(val==3)
    wnd=9;

	/*nom du fichier resultat*/
 	sprintf(str,"File ? in %s",_CurrentPath);
 	strcpy(nomfichres,SAVE_FILE(str,NULL,&e));
 	put_file_extension(nomfichres,".trf",nomfichres);


  filtre_champ_gaussien_3d_p(nomfich1,nomfichres,wnd,sigma);


}


/******************************************************************************
*********
*********  ----------------- filtre_champ_gaussien_3d_p
*********
********************************************************************************/ 

int filtre_champ_gaussien_3d_p(char* nomfich1, char* nomfichres, int wnd,double sigma)
{
  int res=0,wdth,hght,dpth;
	transf3d *transfo,*transfores;
	field3d *ch,*chres;

 
  /*lecture des champ*/
    transfo=load_transf_3d(nomfich1);
		
		wdth=transfo->width;
		hght=transfo->height;
		dpth=transfo->depth;
		

    transfores=cr_transf3d(wdth,hght,dpth,NULL);
		
		ch=transf_to_field_3d(transfo,NULL,NULL);
		chres=cr_field3d(wdth,hght,dpth);

  	imx_filtre_champ_gaussien_3d_p(ch,chres,wnd,sigma);

		transfores=field_to_transf_3d(chres,NULL,NULL);
	  transfores->typetrans=CHAMP3D;
		transfores->dx=transfo->dx;
		transfores->dy=transfo->dy;
		transfores->dz=transfo->dz;
	
		free_field3d(chres);
		free_field3d(ch);
    save_transf_3d(transfores,nomfichres);
    free_transf3d(transfo); free_transf3d(transfores);

  return res;
}

/******************************************************************************
*********
*********  ----------------- imx_filtre_champ_gaussien_3d_p
*********
********************************************************************************/
void imx_filtre_champ_gaussien_3d_p(field3d *imdeb, field3d *imres, int wnd, double sigma)
{
 field3d *inter;
 int i,j,k,p;
 double *filter,sum=0.0,coeff;
 double norm=0.0;
 int width,height,depth;

wnd = MAXI((int) floor(sigma+0.5)*3,3);

 width=imdeb->width;
 height=imdeb->height;
 depth=imdeb->depth;

 filter=alloc_dvector(wnd);
 inter=cr_field3d(width,height,depth);
 
 coeff=pow((double) 2.0*PI,(double)0.5);
 for(i=0;i<wnd;i++)
/* JP 3/99   filter[i]=(1/(2.0*PI*sigma))*exp(-1.0*((double)(i-wnd/2))*((double)(i-wnd/2))/(2.0*sigma*sigma));*/
/* 2.50662828 correspond a sqrt(2PI)   */
    filter[i]=(1/(coeff*sigma))*exp(-1.0*((double)(i-wnd/2))*((double)(i-wnd/2))/(2.0*sigma*sigma));
 for(i=0;i<wnd;i++)
    norm+=filter[i];

 for(i=0;i<wnd;i++)
    filter[i]/=norm;

  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	    {
			inter->raw[i][j][k].x=0;
			inter->raw[i][j][k].y=0;
			inter->raw[i][j][k].z=0;
			}
			
			
/*  for(k=wnd/2;k<depth-wnd/2;k++)
  for(j=wnd/2;j<height-wnd/2;j++)
	for(i=wnd/2;i<width-wnd/2;i++)
	   {
      sum=0.0;
     	for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(imdeb->raw[i-p][j][k].x))*filter[p+wnd/2];
			inter->raw[i][j][k].x=sum;
	   
		 sum=0.0;
     	for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(imdeb->raw[i-p][j][k].y))*filter[p+wnd/2];
			inter->raw[i][j][k].y=sum;
	   
		 sum=0.0;
     	for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(imdeb->raw[i-p][j][k].z))*filter[p+wnd/2];
			inter->raw[i][j][k].z=sum;
	   
		 
		 }

  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	    {
			imres->raw[i][j][k].x=0;
			imres->raw[i][j][k].y=0;
			imres->raw[i][j][k].z=0;
			}

  for(k=wnd/2;k<depth-wnd/2;k++)
     for(i=wnd/2;i<width-wnd/2;i++)
        for(j=wnd/2;j<height-wnd/2;j++)
	   {
      sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(inter->raw[i][j-p][k].x))*filter[p+wnd/2];
			imres->raw[i][j][k].x=sum;
	   
		 sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(inter->raw[i][j-p][k].y))*filter[p+wnd/2];
			imres->raw[i][j][k].y=sum;
	   
		 sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(inter->raw[i][j-p][k].z))*filter[p+wnd/2];
			imres->raw[i][j][k].z=sum;
	   
		 }


  for(i=wnd/2;i<width-wnd/2;i++)
     for(j=wnd/2;j<height-wnd/2;j++)
        for(k=wnd/2;k<depth-wnd/2;k++)
	   {
      sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(imres->raw[i][j][k-p].x))*filter[p+wnd/2];
			inter->raw[i][j][k].x=sum;

			sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(imres->raw[i][j][k-p].y))*filter[p+wnd/2];
			inter->raw[i][j][k].y=sum;

 			sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       sum+=((double)(imres->raw[i][j][k-p].z))*filter[p+wnd/2];
			inter->raw[i][j][k].z=sum;


	   }
*/

for(k=0;k<depth;k++)
  for(j=0;j<height;j++)
	for(i=0;i<width;i++)
	   {
      sum=0.0;
     	for(p=-wnd/2;p<=wnd/2;p++)
	       if ((i-p>=0)&&(i-p<=width-1))
				 sum+=((double)(imdeb->raw[i-p][j][k].x))*filter[p+wnd/2];
			inter->raw[i][j][k].x=sum;
	   
		 sum=0.0;
     	for(p=-wnd/2;p<=wnd/2;p++)
	       if ((i-p>=0)&&(i-p<=width-1))
				 sum+=((double)(imdeb->raw[i-p][j][k].y))*filter[p+wnd/2];
			inter->raw[i][j][k].y=sum;
	   
		 sum=0.0;
     	for(p=-wnd/2;p<=wnd/2;p++)
	       if ((i-p>=0)&&(i-p<=width-1))
				 sum+=((double)(imdeb->raw[i-p][j][k].z))*filter[p+wnd/2];
			inter->raw[i][j][k].z=sum;
	   
		 
		 }

  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	    {
			imres->raw[i][j][k].x=0;
			imres->raw[i][j][k].y=0;
			imres->raw[i][j][k].z=0;
			}

  for(k=0;k<depth;k++)
     for(i=0;i<width;i++)
        for(j=0;j<height;j++)
	   {
      sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       if ((j-p>=0)&&(j-p<=height-1))
				 sum+=((double)(inter->raw[i][j-p][k].x))*filter[p+wnd/2];
			imres->raw[i][j][k].x=sum;
	   
		 sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	       if ((j-p>=0)&&(j-p<=height-1))
				  sum+=((double)(inter->raw[i][j-p][k].y))*filter[p+wnd/2];
			imres->raw[i][j][k].y=sum;
	   
		 sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	        if ((j-p>=0)&&(j-p<=height-1))
				 sum+=((double)(inter->raw[i][j-p][k].z))*filter[p+wnd/2];
			imres->raw[i][j][k].z=sum;
	   
		 }


  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
        for(k=0;k<depth;k++)
	   {
      sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	        if ((k-p>=0)&&(k-p<=depth-1))
				 sum+=((double)(imres->raw[i][j][k-p].x))*filter[p+wnd/2];
			inter->raw[i][j][k].x=sum;

			sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	        if ((k-p>=0)&&(k-p<=depth-1))
			sum+=((double)(imres->raw[i][j][k-p].y))*filter[p+wnd/2];
			inter->raw[i][j][k].y=sum;

 			sum=0.0;
      for(p=-wnd/2;p<=wnd/2;p++)
	         if ((k-p>=0)&&(k-p<=depth-1))
			 sum+=((double)(imres->raw[i][j][k-p].z))*filter[p+wnd/2];
			inter->raw[i][j][k].z=sum;


	   }

  for(i=0;i<width;i++)
     for(j=0;j<height;j++)
	for(k=0;k<depth;k++)
	    {
			imres->raw[i][j][k].x=inter->raw[i][j][k].x;
			imres->raw[i][j][k].y=inter->raw[i][j][k].y;
			imres->raw[i][j][k].z=inter->raw[i][j][k].z;
			}


 free_dvector(filter,wnd);
 free_field3d(inter);

}


/******************************************************************************
*********
*********  ----------------- imx_filtre_champ_gaussien_3d_p
*********
********************************************************************************/
void imx_filtre_param_gaussien_3d_p(double *param_deb, double *param_res, int nb_param, double sigma)
{
int size,i,j,k,l,wnd=5;
field3d *ch_deb, *ch_res;

size=TOP_D(nb_param);

ch_deb=cr_field3d(size,size,size);
ch_res=cr_field3d(size,size,size);

l=0;
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	ch_deb->raw[i][j][k].x=param_deb[3*l];
	ch_deb->raw[i][j][k].y=param_deb[3*l+1];
	ch_deb->raw[i][j][k].z=param_deb[3*l+2];
	l++;
	}

imx_filtre_champ_gaussien_3d_p(ch_deb, ch_res, wnd,  sigma);

l=0;
for (i=0;i<size;i++)
for (j=0;j<size;j++)
for (k=0;k<size;k++)
	{
	param_res[3*l]=ch_res->raw[i][j][k].x;
	param_res[3*l+1]=ch_res->raw[i][j][k].y;
	param_res[3*l+2]=ch_res->raw[i][j][k].z;
	l++;
	}


free_field3d(ch_deb);
free_field3d(ch_res);
}


/******************************************************************************
*********
*********  ----------------- eval_deplacement_bspline1_3d
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

void eval_deplacement_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, double x, double y, double z, dvector3d* u)
{
int resol,coeff2l,topi,topj,topk,l,topMax;
double x0,y0,z0,xm,ym,zm,norm;
double ax000,ax001,ax010,ax100,ax110,ax101,ax011,ax111;
double ay000,ay001,ay010,ay100,ay110,ay101,ay011,ay111;
double az000,az001,az010,az100,az110,az101,az011,az111;
double dx0,dy0,dz0,dxm,dym,dzm;


if (INV_BSPLINE_RESOL!=-1)
	resol=INV_BSPLINE_RESOL;
else
	resol=floor(log (pow(nb_param/3.0,1.0/3.0)+1)/log(2.0)+0.5);

if (INV_BSPLINE_COEFF2L!=-1)
	coeff2l=INV_BSPLINE_COEFF2L;
else
	coeff2l=pow(2.0,resol);
	
	
topMax=coeff2l-1;


topi=floor(x/(double)wdth*(double)coeff2l);

if (x==wdth)
	topi=topi-1;

topj=floor(y/(double)hght*(double)coeff2l);

if (y==hght)
	topj=topj-1;

topk=floor(z/(double)dpth*(double)coeff2l);

if (z==dpth)
	topk=topk-1;

x0=(double)topi*(double)wdth/(double)coeff2l;
y0=(double)topj*(double)hght/(double)coeff2l;
z0=(double)topk*(double)dpth/(double)coeff2l;
xm=((double)topi+1.0)*(double)wdth/(double)coeff2l;
ym=((double)topj+1.0)*(double)hght/(double)coeff2l;
zm=((double)topk+1.0)*(double)dpth/(double)coeff2l;

norm=(xm-x0)*(ym-y0)*(zm-z0);


/* recopie des paramètres relatif au bloc considéré */
	
	if ((topi<topMax)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj,topk,topMax)/3;
		ax111=param[3*l];ay111=param[3*l+1];az111=param[3*l+2];
		}
	else
		{
		ax111=0.0;ay111=0.0;az111=0.0;
		}
		
	if ((topi>0)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi-1,topj,topk,topMax)/3;
		ax011=param[3*l];ay011=param[3*l+1];az011=param[3*l+2];
		}
	else
		{
		ax011=0.0;ay011=0.0;az011=0.0;
		}

	if ((topj>0)&&(topi<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj-1,topk,topMax)/3;
		ax101=param[3*l];ay101=param[3*l+1];az101=param[3*l+2];
		}
	else
		{
		ax101=0.0;ay101=0.0;az101=0.0;
		}


	if ((topk>0)&&(topi<topMax)&&(topj<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj,topk-1,topMax)/3;
		ax110=param[3*l];ay110=param[3*l+1];az110=param[3*l+2];
		}
	else
		{
		ax110=0.0;ay110=0.0;az110=0.0;
		}

	if ((topi>0)&&(topj>0)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi-1,topj-1,topk,topMax)/3;
		ax001=param[3*l];ay001=param[3*l+1];az001=param[3*l+2];
		}
	else
		{
		ax001=0.0;ay001=0.0;az001=0.0;
		}


	if ((topi>0)&&(topk>0)&&(topj<topMax))
		{
		l=TOP_conv_ind_topD(topi-1,topj,topk-1,topMax)/3;
		ax010=param[3*l];ay010=param[3*l+1];az010=param[3*l+2];
		}
	else
		{
		ax010=0.0;ay010=0.0;az010=0.0;
		}


	if ((topj>0)&&(topk>0)&&(topi<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj-1,topk-1,topMax)/3;
		ax100=param[3*l];ay100=param[3*l+1];az100=param[3*l+2];
		}
	else
		{
		ax100=0.0;ay100=0.0;az100=0.0;
		}


	if ((topi>0)&&(topj>0)&&(topk>0))
		{
		l=TOP_conv_ind_topD(topi-1,topj-1,topk-1,topMax)/3;
		ax000=param[3*l];ay000=param[3*l+1];az000=param[3*l+2];
		}
	else
		{
		ax000=0.0;ay000=0.0;az000=0.0;
		}

dxm=xm-x;
dym=ym-y;
dzm=zm-z;
dx0=x-x0;
dy0=y-y0;
dz0=z-z0;


u->x=
 ax000*dxm*dym*dzm+
 ax100*dx0*dym*dzm+
 ax010*dxm*dy0*dzm+
 ax001*dxm*dym*dz0+
 ax110*dx0*dy0*dzm+
 ax011*dxm*dy0*dz0+
 ax101*dx0*dym*dz0+
 ax111*dx0*dy0*dz0;

u->x=u->x/norm;

 
u->y=
 ay000*dxm*dym*dzm+
 ay100*dx0*dym*dzm+
 ay010*dxm*dy0*dzm+
 ay001*dxm*dym*dz0+
 ay110*dx0*dy0*dzm+
 ay011*dxm*dy0*dz0+
 ay101*dx0*dym*dz0+
 ay111*dx0*dy0*dz0;


u->y=u->y/norm;


u->z=
 az000*dxm*dym*dzm+
 az100*dx0*dym*dzm+
 az010*dxm*dy0*dzm+
 az001*dxm*dym*dz0+
 az110*dx0*dy0*dzm+
 az011*dxm*dy0*dz0+
 az101*dx0*dym*dz0+
 az111*dx0*dy0*dz0;

u->z=u->z/norm;


/* Forme factorisee pour gagner peut-être du temps CPU */
/*u->x=
 dzm*(
 		dym*(ax000*dxm+ax100*dx0)+
 		dy0*(ax010*dxm+ax110*dx0))+
 dz0*(
 		dym*(ax001*dxm+ax101*dx0)+
		dy0*(ax011*dxm+ax111*dx0));


u->x=u->x/norm;

 
u->y=
dzm*(
 		dym*(ay000*dxm+ay100*dx0)+
 		dy0*(ay010*dxm+ay110*dx0))+
 dz0*(
 		dym*(ay001*dxm+ay101*dx0)+
		dy0*(ay011*dxm+ay111*dx0));



u->y=u->y/norm;


u->z=
dzm*(
 		dym*(az000*dxm+az100*dx0)+
 		dy0*(az010*dxm+az110*dx0))+
 dz0*(
 		dym*(az001*dxm+az101*dx0)+
		dy0*(az011*dxm+az111*dx0));

u->z=u->z/norm;
*/

}

/******************************************************************************
*********
*********  ----------------- eval_transfo_bspline1_3d
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

void eval_transfo_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, double x, double y, double z, dvector3d* u)
{
eval_deplacement_bspline1_3d(param, nb_param, wdth, hght, dpth, x, y,  z,  u);
u->x=u->x+x;
u->y=u->y+y;
u->z=u->z+z;
}

/******************************************************************************
*********
*********  ----------------- eval_jacobien_transfo_bspline1_3d
*********
*********		Evaluation en un point à coordonnées réelles de la matrice jacobienne 
*********    d'un champ de déformation de Bspline de degré 1
********************************************************************************/

void eval_jacobien_transfo_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, double x, double y, double z, double** J)
{
int resol,coeff2l,topi,topj,topk,l,topMax;
double x0,y0,z0,xm,ym,zm,norm;
double ax000,ax001,ax010,ax100,ax110,ax101,ax011,ax111;
double ay000,ay001,ay010,ay100,ay110,ay101,ay011,ay111;
double az000,az001,az010,az100,az110,az101,az011,az111;
double dx0,dy0,dz0,dxm,dym,dzm;

if (INV_BSPLINE_RESOL!=-1)
	resol=INV_BSPLINE_RESOL;
else
	resol=floor(log (pow(nb_param/3.0,1.0/3.0)+1)/log(2.0)+0.5);

if (INV_BSPLINE_COEFF2L!=-1)
	coeff2l=INV_BSPLINE_COEFF2L;
else
	coeff2l=pow(2.0,resol);
	
topMax=coeff2l-1;


topi=floor(x/wdth*coeff2l);
topj=floor(y/hght*coeff2l);
topk=floor(z/dpth*coeff2l);

x0=topi*wdth/coeff2l;
y0=topj*hght/coeff2l;
z0=topk*dpth/coeff2l;
xm=(topi+1)*wdth/coeff2l;
ym=(topj+1)*hght/coeff2l;
zm=(topk+1)*dpth/coeff2l;

norm=(xm-x0)*(ym-y0)*(zm-z0);

/* recopie des paramètres relatif au bloc considéré */
	
	if ((topi<topMax)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
		ax111=param[3*l];ay111=param[3*l+1];az111=param[3*l+2];
		}
	else
		{
		ax111=0.0;ay111=0.0;az111=0.0;
		}
		
	if ((topi>0)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi-1,topj,topk,nb_param)/3;
		ax011=param[3*l];ay011=param[3*l+1];az011=param[3*l+2];
		}
	else
		{
		ax011=0.0;ay011=0.0;az011=0.0;
		}

	if ((topj>0)&&(topi<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi,topj-1,topk,nb_param)/3;
		ax101=param[3*l];ay101=param[3*l+1];az101=param[3*l+2];
		}
	else
		{
		ax101=0.0;ay101=0.0;az101=0.0;
		}


	if ((topk>0)&&(topi<topMax)&&(topj<topMax))
		{
		l=TOP_conv_ind(topi,topj,topk-1,nb_param)/3;
		ax110=param[3*l];ay110=param[3*l+1];az110=param[3*l+2];
		}
	else
		{
		ax110=0.0;ay110=0.0;az110=0.0;
		}

	if ((topi>0)&&(topj>0)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi-1,topj-1,topk,nb_param)/3;
		ax001=param[3*l];ay001=param[3*l+1];az001=param[3*l+2];
		}
	else
		{
		ax001=0.0;ay001=0.0;az001=0.0;
		}


	if ((topi>0)&&(topk>0)&&(topj<topMax))
		{
		l=TOP_conv_ind(topi-1,topj,topk-1,nb_param)/3;
		ax010=param[3*l];ay010=param[3*l+1];az010=param[3*l+2];
		}
	else
		{
		ax010=0.0;ay010=0.0;az010=0.0;
		}


	if ((topj>0)&&(topk>0)&&(topi<topMax))
		{
		l=TOP_conv_ind(topi,topj-1,topk-1,nb_param)/3;
		ax100=param[3*l];ay100=param[3*l+1];az100=param[3*l+2];
		}
	else
		{
		ax100=0.0;ay100=0.0;az100=0.0;
		}


	if ((topi>0)&&(topj>0)&&(topk>0))
		{
		l=TOP_conv_ind(topi-1,topj-1,topk-1,nb_param)/3;
		ax000=param[3*l];ay000=param[3*l+1];az000=param[3*l+2];
		}
	else
		{
		ax000=0.0;ay000=0.0;az000=0.0;
		}

dxm=xm-x;
dym=ym-y;
dzm=zm-z;
dx0=x-x0;
dy0=y-y0;
dz0=z-z0;

/* dhx/dx */
J[0][0]=1.0+(
				-	ax000*dym*dzm
				+	ax100*dym*dzm
 				-	ax010*dy0*dzm
 				-	ax001*dym*dz0
 				+	ax110*dy0*dzm
 				-	ax011*dy0*dz0
 				+	ax101*dym*dz0
 				+	ax111*dy0*dz0)/norm;



/* dhx/dy */
J[0][1]=(
				-	ax000*dxm*dzm
				-	ax100*dx0*dzm
 				+	ax010*dxm*dzm
 				-	ax001*dxm*dz0
 				+	ax110*dx0*dzm
 				+	ax011*dxm*dz0
 				-	ax101*dx0*dz0
 				+	ax111*dx0*dz0)/norm;


/* dhx/dz */
J[0][2]=(
				-	ax000*dxm*dym
				-	ax100*dx0*dym
 				-	ax010*dxm*dy0
 				+	ax001*dxm*dym
				-	ax110*dx0*dy0
 				+	ax011*dxm*dy0
 				+	ax101*dx0*dym
 				+	ax111*dx0*dy0)/norm;




/* dhy/dx */
J[1][0]=(
				-	ay000*dym*dzm
				+	ay100*dym*dzm
 				-	ay010*dy0*dzm
 				-	ay001*dym*dz0
 				+	ay110*dy0*dzm
 				-	ay011*dy0*dz0
 				+	ay101*dym*dz0
 				+	ay111*dy0*dz0)/norm;


/* dhy/dy */
J[1][1]=1.0+(
				-	ay000*dxm*dzm
				-	ay100*dx0*dzm
 				+	ay010*dxm*dzm
 				-	ay001*dxm*dz0
 				+	ay110*dx0*dzm
 				+	ay011*dxm*dz0
 				-	ay101*dx0*dz0
 				+	ay111*dx0*dz0)/norm;


/* dhy/dz */
J[1][2]=(
				-	ay000*dxm*dym
				-	ay100*dx0*dym
 				-	ay010*dxm*dy0
 				+	ay001*dxm*dym
				-	ay110*dx0*dy0
 				+	ay011*dxm*dy0
 				+	ay101*dx0*dym
 				+	ay111*dx0*dy0)/norm;



/* dhz/dx */
J[2][0]=(
				-	az000*dym*dzm
				+	az100*dym*dzm
 				-	az010*dy0*dzm
 				-	az001*dym*dz0
 				+	az110*dy0*dzm
 				-	az011*dy0*dz0
 				+	az101*dym*dz0
 				+	az111*dy0*dz0)/norm;


/* dhz/dy */
J[2][1]=(
				-	az000*dxm*dzm
				-	az100*dx0*dzm
 				+	az010*dxm*dzm
 				-	az001*dxm*dz0
 				+	az110*dx0*dzm
 				+	az011*dxm*dz0
 				-	az101*dx0*dz0
 				+	az111*dx0*dz0)/norm;


/* dhz/dz */
J[2][2]=1.0+(
				-	az000*dxm*dym
				-	az100*dx0*dym
 				-	az010*dxm*dy0
 				+	az001*dxm*dym
				-	az110*dx0*dy0
 				+	az011*dxm*dy0
 				+	az101*dx0*dym
 				+	az111*dx0*dy0)/norm;


}



/******************************************************************************
*********
*********  ----------------- eval_jacobien_transfo_bspline1_BoiteInvBspline_3d
*********
*********		Evaluation en un point à coordonnées réelles de la matrice jacobienne 
*********    d'un champ de déformation de Bspline de degré 1
********************************************************************************/

void eval_jacobien_transfo_bspline1_BoiteInvBspline_3d(BoiteInvBspline *ParamInvBspline, double x, double y, double z, double** J)
{
double dx0,dy0,dz0,dxm,dym,dzm;


dxm=ParamInvBspline->xm-x;
dym=ParamInvBspline->ym-y;
dzm=ParamInvBspline->zm-z;
dx0=x-ParamInvBspline->x0;
dy0=y-ParamInvBspline->y0;
dz0=z-ParamInvBspline->z0;

/* dhx/dx */
J[0][0]=1.0+(
				-	ParamInvBspline->ax000*dym*dzm
				+	ParamInvBspline->ax100*dym*dzm
 				-	ParamInvBspline->ax010*dy0*dzm
 				-	ParamInvBspline->ax001*dym*dz0
 				+	ParamInvBspline->ax110*dy0*dzm
 				-	ParamInvBspline->ax011*dy0*dz0
 				+	ParamInvBspline->ax101*dym*dz0
 				+	ParamInvBspline->ax111*dy0*dz0)/ParamInvBspline->norm;



/* dhx/dy */
J[0][1]=(
				-	ParamInvBspline->ax000*dxm*dzm
				-	ParamInvBspline->ax100*dx0*dzm
 				+	ParamInvBspline->ax010*dxm*dzm
 				-	ParamInvBspline->ax001*dxm*dz0
 				+	ParamInvBspline->ax110*dx0*dzm
 				+	ParamInvBspline->ax011*dxm*dz0
 				-	ParamInvBspline->ax101*dx0*dz0
 				+	ParamInvBspline->ax111*dx0*dz0)/ParamInvBspline->norm;


/* dhx/dz */
J[0][2]=(
				-	ParamInvBspline->ax000*dxm*dym
				-	ParamInvBspline->ax100*dx0*dym
 				-	ParamInvBspline->ax010*dxm*dy0
 				+	ParamInvBspline->ax001*dxm*dym
				-	ParamInvBspline->ax110*dx0*dy0
 				+	ParamInvBspline->ax011*dxm*dy0
 				+	ParamInvBspline->ax101*dx0*dym
 				+	ParamInvBspline->ax111*dx0*dy0)/ParamInvBspline->norm;




/* dhy/dx */
J[1][0]=(
				-	ParamInvBspline->ay000*dym*dzm
				+	ParamInvBspline->ay100*dym*dzm
 				-	ParamInvBspline->ay010*dy0*dzm
 				-	ParamInvBspline->ay001*dym*dz0
 				+	ParamInvBspline->ay110*dy0*dzm
 				-	ParamInvBspline->ay011*dy0*dz0
 				+	ParamInvBspline->ay101*dym*dz0
 				+	ParamInvBspline->ay111*dy0*dz0)/ParamInvBspline->norm;


/* dhy/dy */
J[1][1]=1.0+(
				-	ParamInvBspline->ay000*dxm*dzm
				-	ParamInvBspline->ay100*dx0*dzm
 				+	ParamInvBspline->ay010*dxm*dzm
 				-	ParamInvBspline->ay001*dxm*dz0
 				+	ParamInvBspline->ay110*dx0*dzm
 				+	ParamInvBspline->ay011*dxm*dz0
 				-	ParamInvBspline->ay101*dx0*dz0
 				+	ParamInvBspline->ay111*dx0*dz0)/ParamInvBspline->norm;


/* dhy/dz */
J[1][2]=(
				-	ParamInvBspline->ay000*dxm*dym
				-	ParamInvBspline->ay100*dx0*dym
 				-	ParamInvBspline->ay010*dxm*dy0
 				+	ParamInvBspline->ay001*dxm*dym
				-	ParamInvBspline->ay110*dx0*dy0
 				+	ParamInvBspline->ay011*dxm*dy0
 				+	ParamInvBspline->ay101*dx0*dym
 				+	ParamInvBspline->ay111*dx0*dy0)/ParamInvBspline->norm;



/* dhz/dx */
J[2][0]=(
				-	ParamInvBspline->az000*dym*dzm
				+	ParamInvBspline->az100*dym*dzm
 				-	ParamInvBspline->az010*dy0*dzm
 				-	ParamInvBspline->az001*dym*dz0
 				+	ParamInvBspline->az110*dy0*dzm
 				-	ParamInvBspline->az011*dy0*dz0
 				+	ParamInvBspline->az101*dym*dz0
 				+	ParamInvBspline->az111*dy0*dz0)/ParamInvBspline->norm;


/* dhz/dy */
J[2][1]=(
				-	ParamInvBspline->az000*dxm*dzm
				-	ParamInvBspline->az100*dx0*dzm
 				+	ParamInvBspline->az010*dxm*dzm
 				-	ParamInvBspline->az001*dxm*dz0
 				+	ParamInvBspline->az110*dx0*dzm
 				+	ParamInvBspline->az011*dxm*dz0
 				-	ParamInvBspline->az101*dx0*dz0
 				+	ParamInvBspline->az111*dx0*dz0)/ParamInvBspline->norm;


/* dhz/dz */
J[2][2]=1.0+(
				-	ParamInvBspline->az000*dxm*dym
				-	ParamInvBspline->az100*dx0*dym
 				-	ParamInvBspline->az010*dxm*dy0
 				+	ParamInvBspline->az001*dxm*dym
				-	ParamInvBspline->az110*dx0*dy0
 				+	ParamInvBspline->az011*dxm*dy0
 				+	ParamInvBspline->az101*dx0*dym
 				+	ParamInvBspline->az111*dx0*dy0)/ParamInvBspline->norm;


}

/******************************************************************************
*********
*********  ----------------- raffinement_bspline_inv_3d
*********
*********		Raffinement de l'inversion du champ Bspline 
*********
*********				
********************************************************************************/
double raffinement_bspline_inv_3d(double *param,int nb_param,field3d *chdepart,field3d *chres,field3d *chmin,field3d *chmax,int i,int j,int k,double prec)
{
TOPTliste myQ,myL;
TOPTbox myB,auxB;
dvector3d hmid;
double distance,distance2,prec_tmp/*,prec_num=0.0001*/,reduc;
int wdth,hght,dpth,umin,umax,vmin,vmax,wmin,wmax,nb_boite=0;
int resol,coeff2l,imin,imax,jmin,jmax,kmin,kmax,ii,jj,kk;
int u,v,w,u1,v1,w1,nb_iter=0;
double  xmid,ymid,zmid;

 prec_tmp=prec;
 myQ.nelem = 0;
 myL.nelem = 0;
 
if (INV_BSPLINE_RESOL!=-1)
	resol=INV_BSPLINE_RESOL;
else
	resol=floor(log (pow(nb_param/3.0,1.0/3.0)+1)/log(2.0)+0.5);

if (INV_BSPLINE_COEFF2L!=-1)
	coeff2l=INV_BSPLINE_COEFF2L;
else
	coeff2l=pow(2.0,resol);
	

//printf("Bloc %d %d %d \n",i,j,k);


 		
 wdth=chdepart->width;
 hght=chdepart->height;
 dpth=chdepart->depth;
 
 umin=MAXI(chmin->raw[i][j][k].x,0);
 vmin=MAXI(chmin->raw[i][j][k].y,0);
 wmin=MAXI(chmin->raw[i][j][k].z,0);
 
 umax=MINI(chmax->raw[i][j][k].x,wdth);
 vmax=MINI(chmax->raw[i][j][k].y,hght);
 wmax=MINI(chmax->raw[i][j][k].z,dpth);
 


/*auxB.xm=umin;auxB.xM=umax;
auxB.ym=vmin;auxB.yM=vmax;
auxB.zm=wmin;auxB.zM=wmax;

hmid.x=2.0*i-chdepart->raw[i][j][k].x;
hmid.y=2.0*j-chdepart->raw[i][j][k].y;
hmid.z=2.0*k-chdepart->raw[i][j][k].z;


if (inv_bspline_find_solution_with_fixed_point_algo_global(param, nb_param, wdth, hght, dpth, prec, (double)i, (double)j, (double) k, &auxB, &hmid, &distance))
		 	{
			chres->raw[i][j][k].x=hmid.x-i;
			chres->raw[i][j][k].y=hmid.y-j;
			chres->raw[i][j][k].z=hmid.z-k;
			return(distance);
			}
			
			
NB_BLOC_POINT_FIXE++;	
*/





 //---- On découpe la boite pour que chaque sous-boite soi incluse sur un Slpqr

 if (umin<wdth)
 	imin=floor(1.0*umin/wdth*coeff2l);
 else
	imin=floor(1.0*umin/wdth*coeff2l)-1;

 
 if (umax<wdth)
	imax=floor(1.0*umax/wdth*coeff2l)+1;
 else
  imax=floor(1.0*umax/wdth*coeff2l);
 
 
 if (vmin<hght)
 	jmin=floor(1.0*vmin/hght*coeff2l);
 else
 jmin=floor(1.0*vmin/hght*coeff2l)-1;
 
 if (vmax<hght)
 	jmax=floor(1.0*vmax/hght*coeff2l)+1;
 else
 	jmax=floor(1.0*vmax/hght*coeff2l);
 
 
 if (wmin<dpth)
 	kmin=floor(1.0*wmin/dpth*coeff2l);
 else
 	kmin=floor(1.0*wmin/dpth*coeff2l)-1;
 
 if (wmax<dpth)
 	kmax=floor(1.0*wmax/dpth*coeff2l)+1;
 else
 	kmax=floor(1.0*wmax/dpth*coeff2l);
 
 //if ((i==5)&&(j==124)&&(k==128)) 
//	printf("tutu\n"); 
 
 for (ii=imin;ii<imax;ii++)
 for (jj=jmin;jj<jmax;jj++)
 for (kk=kmin;kk<kmax;kk++)
 		{
		u=(double)ii*(double)wdth/(double)coeff2l;
		u1=(double)(ii+1.0)*(double)wdth/(double)coeff2l;
		
		v=(double)jj*(double)hght/(double)coeff2l;
		v1=(double)(jj+1.0)*(double)hght/(double)coeff2l;
		
		w=(double)kk*(double)dpth/(double)coeff2l;
		w1=(double)(kk+1.0)*(double)dpth/(double)coeff2l;
		
		u=MAXI(u,umin);
		v=MAXI(v,vmin);
		w=MAXI(w,wmin);
		
		u1=MINI(u1,umax);
		v1=MINI(v1,vmax);
		w1=MINI(w1,wmax);
		
		
		auxB.xm=u;
		auxB.xM=u1;
		auxB.ym=v;
		auxB.yM=v1;
		auxB.zm=w;
		auxB.zM=w1;
		
		auxB.hx[0].x=chdepart->raw[u][v][w].x;
		auxB.hx[0].y=chdepart->raw[u][v][w].y;
		auxB.hx[0].z=chdepart->raw[u][v][w].z;
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u, v, w, &auxB.hx[0]);
		
		if (u1<wdth) 
		{
		auxB.hx[1].x=chdepart->raw[u1][v][w].x;
		auxB.hx[1].y=chdepart->raw[u1][v][w].y;
		auxB.hx[1].z=chdepart->raw[u1][v][w].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u1, v, w, &auxB.hx[1]);
		
		if (v1<hght) 
		{
		auxB.hx[2].x=chdepart->raw[u][v1][w].x;
		auxB.hx[2].y=chdepart->raw[u][v1][w].y;
		auxB.hx[2].z=chdepart->raw[u][v1][w].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u, v1, w, &auxB.hx[2]);
		
		
		if ((u1<wdth)&&(v1<hght)) 
		{
		auxB.hx[3].x=chdepart->raw[u1][v1][w].x;
		auxB.hx[3].y=chdepart->raw[u1][v1][w].y;
		auxB.hx[3].z=chdepart->raw[u1][v1][w].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u1, v1, w, &auxB.hx[3]);

		
		if ((w1<dpth)) 
		{
		auxB.hx[4].x=chdepart->raw[u][v][w1].x;
		auxB.hx[4].y=chdepart->raw[u][v][w1].y;
		auxB.hx[4].z=chdepart->raw[u][v][w1].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u, v, w1, &auxB.hx[4]);

		if ((u1<wdth)&&(w1<dpth)) 
		{
		auxB.hx[5].x=chdepart->raw[u1][v][w1].x;
		auxB.hx[5].y=chdepart->raw[u1][v][w1].y;
		auxB.hx[5].z=chdepart->raw[u1][v][w1].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u1, v, w1, &auxB.hx[5]);

		if ((v1<hght)&&(w1<dpth)) 
		{
		auxB.hx[6].x=chdepart->raw[u][v1][w1].x;
		auxB.hx[6].y=chdepart->raw[u][v1][w1].y;
		auxB.hx[6].z=chdepart->raw[u][v1][w1].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u, v1, w1, &auxB.hx[6]);

		if ((u1<wdth)&&(v1<hght)&&(w1<dpth)) 
		{
		auxB.hx[7].x=chdepart->raw[u1][v1][w1].x;
		auxB.hx[7].y=chdepart->raw[u1][v1][w1].y;
		auxB.hx[7].z=chdepart->raw[u1][v1][w1].z;
		}
		else
		eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, u1, v1, w1, &auxB.hx[7]);
		
		
		if ((u<u1)&&(v<v1)&&(w<w1))   // On fait ce test car il y a des boites qui peuvent etre de largeur nulle
			{
			
			/*if (inv_bspline_find_solution_with_fixed_point_algo(param, nb_param, wdth, hght, dpth, prec, (double)i, (double)j, (double) k, &auxB, &hmid, &distance))
		 	{

			
			chres->raw[i][j][k].x=hmid.x-i;
			chres->raw[i][j][k].y=hmid.y-j;
			chres->raw[i][j][k].z=hmid.z-k;
			// on vide la liste chainée
			while (myQ.nelem>0)
				depile_end(&myQ);
				
			//printf("Le point fixe a fonctionne pour le bloc %d %d %d\n",i,j,k);	
			return(distance);
			}
		*/
			empile_front(&myQ,auxB);
			
			}
		
		}
			
		//printf("Le point fixe n'a fonctionne pour le bloc %d %d %d\n",i,j,k);	
		
NB_BLOC_POINT_FIXE++;	

do
{
	
while (myQ.nelem>0)
	{
	
	 	myB = depile_end(&myQ); 
		
		
		if ((myB.xm>myB.xM)||(myB.ym>myB.yM)||(myB.zm>myB.zM))
			printf("Il y a un gros soucis !!!!!!!!!! Inversion des bornes !!!! \n");
		/*else
		{
		// Aucune solution n'a ete trouve. Le champ ne preservait pas la topologie !!! On renvoie une valeur nulle pour ce point .... 
		chres->raw[i][j][k].x=0;
		chres->raw[i][j][k].y=0;
		chres->raw[i][j][k].z=0;
 		
		eval_transfo_bspline1_3d(param, nb_param, chres->width, chres->height, chres->depth, i, j, k, &hmid);
		distance=sqrt((hmid.x-i)*(hmid.x-i)+(hmid.y-j)*(hmid.y-j)+(hmid.z-k)*(hmid.z-k));
		
		printf("le champ de deformation n'est pas inversible au point (%d,%d,%d) \n",i,j,k);
		
		return(distance);
		}*/
	 nb_boite++;
		
	
	if (acknowledge_invbspline(param, nb_param, chres->width, chres->height, chres->depth,i,j,k,myB))
	{
				
	 distance2=0.5*sqrt((myB.xM-myB.xm)*(myB.xM-myB.xm)+(myB.yM-myB.ym)*(myB.yM-myB.ym)+(myB.zM-myB.zm)*(myB.zM-myB.zm));
	 xmid= (myB.xm+ myB.xM)/2.0;ymid= (myB.ym+ myB.yM)/2.0;zmid= (myB.zm+ myB.zM)/2.0;
	 eval_transfo_bspline1_3d(param, nb_param, chres->width, chres->height, chres->depth, xmid, ymid, zmid, &hmid);
	 distance=sqrt((hmid.x-i)*(hmid.x-i)+(hmid.y-j)*(hmid.y-j)+(hmid.z-k)*(hmid.z-k));
	 
		if ((distance>prec_tmp)||(distance2>prec_tmp))
			{
			//-- On reduit la taille de la boite	
				reduc=inv_bspline_reduction_boite_3d(param,  nb_param, chres->width, chres->height, chres->depth,  prec_tmp,  i,  j,  k, &myB);
			
				if (reduc==-1) // la solution est trouvee lors du grignotage
				{
				empile_front(&myL,myB);
			
				}
				else
				{
					if (reduc<=1)
					{
					xmid= (myB.xm+ myB.xM)/2.0;ymid= (myB.ym+ myB.yM)/2.0;zmid= (myB.zm+ myB.zM)/2.0;
					eval_transfo_bspline1_3d(param, nb_param, chres->width, chres->height, chres->depth, xmid, ymid, zmid, &hmid);
					inv_bspline_casser_boite_3d(param,  nb_param, chres->width, chres->height, chres->depth, &myB, &myQ,xmid,ymid,zmid,hmid);
					}
					
				}		
			
			
			}
			else
			{
			//-- La boite est suffisamment précise, on peut quitter la fonction
				
			empile_front(&myL,myB);
				/*chres->raw[i][j][k].x=xmid-i;
				chres->raw[i][j][k].y=ymid-j;
				chres->raw[i][j][k].z=zmid-k;
				
				// on vide la liste chainée
				while (myQ.nelem>0)
					depile_end(&myQ);
									
				return(distance);*/
			}
		
			
	}
			
	}

myB.xm=HUGE_VAL;myB.ym=HUGE_VAL;myB.zm=HUGE_VAL;
myB.xM=-HUGE_VAL;myB.yM=-HUGE_VAL;myB.zM=-HUGE_VAL;


if ((myL.nelem==0)||(nb_iter>10))
		{
		// Aucune solution n'a ete trouve. Le champ ne preservait pas la topologie !!! On renvoie une valeur nulle pour ce point .... 
		chres->raw[i][j][k].x=0;
		chres->raw[i][j][k].y=0;
		chres->raw[i][j][k].z=0;
 		
		eval_transfo_bspline1_3d(param, nb_param, chres->width, chres->height, chres->depth, i, j, k, &hmid);
		distance=sqrt((hmid.x-i)*(hmid.x-i)+(hmid.y-j)*(hmid.y-j)+(hmid.z-k)*(hmid.z-k));
		
		printf("le champ de deformation n'est pas inversible au point (%d,%d,%d). \n",i,j,k);
		printf("la precision requise a ce point est de %f \n",prec_tmp);
		
		while (myL.nelem>0)
		 	depile_end(&myL);
			
		return(distance);
		}
		
while (myL.nelem>0)
	{
	 	auxB = depile_end(&myL);
		myB.xm=MINI(myB.xm,auxB.xm);
		myB.ym=MINI(myB.ym,auxB.ym);
		myB.zm=MINI(myB.zm,auxB.zm);
		
		myB.xM=MAXI(myB.xM,auxB.xM);
		myB.yM=MAXI(myB.yM,auxB.yM);
		myB.zM=MAXI(myB.zM,auxB.zM); 
		
		empile_front(&myQ,auxB);
	}

 distance2=0.5*sqrt((myB.xM-myB.xm)*(myB.xM-myB.xm)+(myB.yM-myB.ym)*(myB.yM-myB.ym)+(myB.zM-myB.zm)*(myB.zM-myB.zm));
 xmid= (myB.xm+ myB.xM)/2.0;ymid= (myB.ym+ myB.yM)/2.0;zmid= (myB.zm+ myB.zM)/2.0;
 eval_transfo_bspline1_3d(param, nb_param, chres->width, chres->height, chres->depth, xmid, ymid, zmid, &hmid);
 distance=sqrt((hmid.x-i)*(hmid.x-i)+(hmid.y-j)*(hmid.y-j)+(hmid.z-k)*(hmid.z-k));

 if ((distance2<prec)&&(distance<prec))
 	{
	chres->raw[i][j][k].x=xmid-i;
	chres->raw[i][j][k].y=ymid-j;
	chres->raw[i][j][k].z=zmid-k;
	 
	 // on vide la liste chainée
	while (myQ.nelem>0)
	  depile_end(&myQ);
				
	 return(distance);
	}


prec_tmp=0.5*prec_tmp;
nb_iter++;


} while(1);

}

/******************************************************************************
*********
*********  ----------------- acknowledge
*********
*********		Verifie si c'est une boite admissible
*********
*********				
********************************************************************************/
int acknowledge_invbspline(double *param, int nb_param, int width, int height, int depth, double topi, double topj, double topk, TOPTbox myB)
{
int l;
double xmin,xmax,ymin,ymax,zmin,zmax;

//update_bloc_transfo_bspline1_3d(param, nb_param, width, height, depth, &myB);


xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;
	
	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,myB.hx[l].x);ymin=MINI(ymin,myB.hx[l].y);zmin=MINI(zmin,myB.hx[l].z);
	xmax=MAXI(xmax,myB.hx[l].x);ymax=MAXI(ymax,myB.hx[l].y);zmax=MAXI(zmax,myB.hx[l].z);
	} 

if ((xmin<=topi)&&(xmax>=topi)&&(ymin<=topj)&&(ymax>=topj)&&(zmin<=topk)&&(zmax>=topk))
	return(1);
else
	return(0);
	
}


/******************************************************************************
*********
*********  ----------------- inv_bspline_casser_boite_3d
*********
*********		Casse une boite en 8 sous-boites
*********
*********				
********************************************************************************/
void inv_bspline_casser_boite_3d(double *param, int nb_param, int width, int height, int depth, TOPTbox *myB, TOPTliste *myQ, double xmid, double ymid, double zmid, dvector3d hmid)
{
 dvector3d hcube[27];
 TOPTbox auxB;
update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  myB);

	hcube[0]=myB->hx[0];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->ym, myB->zm, &hcube[1]);
	hcube[2]=myB->hx[1];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, ymid, myB->zm, &hcube[3]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, ymid, myB->zm, &hcube[4]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, ymid, myB->zm, &hcube[5]);
	hcube[6]=myB->hx[2];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->yM, myB->zm, &hcube[7]);
	hcube[8]=myB->hx[3];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, myB->ym, zmid, &hcube[9]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->ym, zmid, &hcube[10]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, myB->ym, zmid, &hcube[11]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm,ymid, zmid, &hcube[12]);
	hcube[13]=hmid;
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, ymid, zmid, &hcube[14]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, myB->yM, zmid, &hcube[15]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->yM, zmid, &hcube[16]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, myB->yM, zmid, &hcube[17]);
	hcube[18]=myB->hx[4];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->ym, myB->zM, &hcube[19]);
	hcube[20]=myB->hx[5];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, ymid, myB->zM, &hcube[21]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, ymid, myB->zM, &hcube[22]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, ymid, myB->zM, &hcube[23]);
	hcube[24]=myB->hx[6];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->yM, myB->zM, &hcube[25]);
	hcube[26]=myB->hx[7];	
		
		
	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[0];
	 auxB.hx[1]=hcube[1];
	 auxB.hx[2]=hcube[3];
	 auxB.hx[3]=hcube[4];
	 auxB.hx[4]=hcube[9];
	 auxB.hx[5]=hcube[10];
	 auxB.hx[6]=hcube[12];
	 auxB.hx[7]=hcube[13];

	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[9];
	 auxB.hx[1]=hcube[10];
	 auxB.hx[2]=hcube[12];
	 auxB.hx[3]=hcube[13];
	 auxB.hx[4]=hcube[18];
	 auxB.hx[5]=hcube[19];
	 auxB.hx[6]=hcube[21];
	 auxB.hx[7]=hcube[22];
	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[3];
	 auxB.hx[1]=hcube[4];
	 auxB.hx[2]=hcube[6];
	 auxB.hx[3]=hcube[7];
	 auxB.hx[4]=hcube[12];
	 auxB.hx[5]=hcube[13];
	 auxB.hx[6]=hcube[15];
	 auxB.hx[7]=hcube[16];
	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[12];
	 auxB.hx[1]=hcube[13];
	 auxB.hx[2]=hcube[15];
	 auxB.hx[3]=hcube[16];
	 auxB.hx[4]=hcube[21];
	 auxB.hx[5]=hcube[22];
	 auxB.hx[6]=hcube[24];
	 auxB.hx[7]=hcube[25];
	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[1];
	 auxB.hx[1]=hcube[2];
	 auxB.hx[2]=hcube[4];
	 auxB.hx[3]=hcube[5];
	 auxB.hx[4]=hcube[10];
	 auxB.hx[5]=hcube[11];
	 auxB.hx[6]=hcube[13];
	 auxB.hx[7]=hcube[14];
	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[10];
	 auxB.hx[1]=hcube[11];
	 auxB.hx[2]=hcube[13];
	 auxB.hx[3]=hcube[14];
	 auxB.hx[4]=hcube[19];
	 auxB.hx[5]=hcube[20];
	 auxB.hx[6]=hcube[22];
	 auxB.hx[7]=hcube[23];
	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[4];
	 auxB.hx[1]=hcube[5];
	 auxB.hx[2]=hcube[7];
	 auxB.hx[3]=hcube[8];
	 auxB.hx[4]=hcube[13];
	 auxB.hx[5]=hcube[14];
	 auxB.hx[6]=hcube[16];
	 auxB.hx[7]=hcube[17];
	 //update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[13];
	 auxB.hx[1]=hcube[14];
	 auxB.hx[2]=hcube[16];
	 auxB.hx[3]=hcube[17];
	 auxB.hx[4]=hcube[22];
	 auxB.hx[5]=hcube[23];
	 auxB.hx[6]=hcube[25];
	 auxB.hx[7]=hcube[26];
	//update_bloc_transfo_bspline1_3d(param,  nb_param,  width, height, depth,  &auxB);
	 empile_front(myQ,auxB);
}

//-------------------------------------------------------------------------------------------------------------------------------------------
void inv_bspline_casser_boite_3d2(double *param, int nb_param, int width, int height, int depth, TOPTbox *myB, TOPTliste *myQ, double xmid, double ymid, double zmid, dvector3d hmid)
{
 dvector3d hcube[27];
 TOPTbox auxB;

	hcube[0]=myB->hx[0];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->ym, myB->zm, &hcube[1]);
	hcube[2]=myB->hx[1];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, ymid, myB->zm, &hcube[3]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, ymid, myB->zm, &hcube[4]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, ymid, myB->zm, &hcube[5]);
	hcube[6]=myB->hx[2];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->yM, myB->zm, &hcube[7]);
	hcube[8]=myB->hx[3];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, myB->ym, zmid, &hcube[9]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->ym, zmid, &hcube[10]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, myB->ym, zmid, &hcube[11]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm,ymid, zmid, &hcube[12]);
	hcube[13]=hmid;
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, ymid, zmid, &hcube[14]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, myB->yM, zmid, &hcube[15]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->yM, zmid, &hcube[16]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, myB->yM, zmid, &hcube[17]);
	hcube[18]=myB->hx[4];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->ym, myB->zM, &hcube[19]);
	hcube[20]=myB->hx[5];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xm, ymid, myB->zM, &hcube[21]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, ymid, myB->zM, &hcube[22]);
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, myB->xM, ymid, myB->zM, &hcube[23]);
	hcube[24]=myB->hx[6];
	eval_transfo_bspline1_3d(param, nb_param, width, height, depth, xmid, myB->yM, myB->zM, &hcube[25]);
	hcube[26]=myB->hx[7];	
		
		
	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[0];
	 auxB.hx[1]=hcube[1];
	 auxB.hx[2]=hcube[3];
	 auxB.hx[3]=hcube[4];
	 auxB.hx[4]=hcube[9];
	 auxB.hx[5]=hcube[10];
	 auxB.hx[6]=hcube[12];
	 auxB.hx[7]=hcube[13];
	 empile_front(myQ,auxB);

	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[9];
	 auxB.hx[1]=hcube[10];
	 auxB.hx[2]=hcube[12];
	 auxB.hx[3]=hcube[13];
	 auxB.hx[4]=hcube[18];
	 auxB.hx[5]=hcube[19];
	 auxB.hx[6]=hcube[21];
	 auxB.hx[7]=hcube[22];
	 empile_front(myQ,auxB);

	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[3];
	 auxB.hx[1]=hcube[4];
	 auxB.hx[2]=hcube[6];
	 auxB.hx[3]=hcube[7];
	 auxB.hx[4]=hcube[12];
	 auxB.hx[5]=hcube[13];
	 auxB.hx[6]=hcube[15];
	 auxB.hx[7]=hcube[16];
	 empile_front(myQ,auxB);

	 auxB.xm=myB->xm;
	 auxB.xM=xmid;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[12];
	 auxB.hx[1]=hcube[13];
	 auxB.hx[2]=hcube[15];
	 auxB.hx[3]=hcube[16];
	 auxB.hx[4]=hcube[21];
	 auxB.hx[5]=hcube[22];
	 auxB.hx[6]=hcube[24];
	 auxB.hx[7]=hcube[25];
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[1];
	 auxB.hx[1]=hcube[2];
	 auxB.hx[2]=hcube[4];
	 auxB.hx[3]=hcube[5];
	 auxB.hx[4]=hcube[10];
	 auxB.hx[5]=hcube[11];
	 auxB.hx[6]=hcube[13];
	 auxB.hx[7]=hcube[14];
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=myB->ym;
	 auxB.yM=ymid;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[10];
	 auxB.hx[1]=hcube[11];
	 auxB.hx[2]=hcube[13];
	 auxB.hx[3]=hcube[14];
	 auxB.hx[4]=hcube[19];
	 auxB.hx[5]=hcube[20];
	 auxB.hx[6]=hcube[22];
	 auxB.hx[7]=hcube[23];
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=myB->zm;
	 auxB.zM=zmid;
	 auxB.hx[0]=hcube[4];
	 auxB.hx[1]=hcube[5];
	 auxB.hx[2]=hcube[7];
	 auxB.hx[3]=hcube[8];
	 auxB.hx[4]=hcube[13];
	 auxB.hx[5]=hcube[14];
	 auxB.hx[6]=hcube[16];
	 auxB.hx[7]=hcube[17];
	 empile_front(myQ,auxB);

	 auxB.xm=xmid;
	 auxB.xM=myB->xM;
	 auxB.ym=ymid;
	 auxB.yM=myB->yM;
	 auxB.zm=zmid;
	 auxB.zM=myB->zM;
	 auxB.hx[0]=hcube[13];
	 auxB.hx[1]=hcube[14];
	 auxB.hx[2]=hcube[16];
	 auxB.hx[3]=hcube[17];
	 auxB.hx[4]=hcube[22];
	 auxB.hx[5]=hcube[23];
	 auxB.hx[6]=hcube[25];
	 auxB.hx[7]=hcube[26];
	 empile_front(myQ,auxB);
}



/******************************************************************************
*********
*********  ----------------- inv_bsppline_reduction_boite_3d
*********
*********		Grignotte une boite pour mieux encadrer la solution
*********
*********				
********************************************************************************/



 double inv_bspline_reduction_boite_3d(double *param, int nb_param, int wdth, int hght, int dpth, double precis, double i, double j, double k, TOPTbox* myB)
 {
int resol,coeff2l,topi,topj,topk,l,topMax,nb_iter;
double x0,y0,z0,xm,ym,zm,norm,pourcentage,result,tmp;
double ax000,ax001,ax010,ax100,ax110,ax101,ax011,ax111;
double ay000,ay001,ay010,ay100,ay110,ay101,ay011,ay111;
double az000,az001,az010,az100,az110,az101,az011,az111;
BoiteInvBspline ParamInvBspline;
double xmid,ymid,zmid;
double distance,distance2,pourcentage_reduction=0.9;
dvector3d hmid;

//---- Initialisation des paramètres caractéristiques du bloc myB considéré  -----

if (INV_BSPLINE_RESOL!=-1)
	resol=INV_BSPLINE_RESOL;
else
	resol=floor(log (pow(nb_param/3.0,1.0/3.0)+1)/log(2.0)+0.5);

if (INV_BSPLINE_COEFF2L!=-1)
	coeff2l=INV_BSPLINE_COEFF2L;
else
	coeff2l=pow(2.0,resol);
	
topMax=coeff2l-1;


topi=floor(myB->xm/wdth*coeff2l);

if (myB->xm==wdth)
	topi=topi-1;

topj=floor(myB->ym/hght*coeff2l);

if (myB->ym==hght)
	topj=topj-1;

topk=floor(myB->zm/dpth*coeff2l);

if (myB->zm==dpth)
	topk=topk-1;

x0=(double)topi*wdth/coeff2l;
y0=(double)topj*hght/coeff2l;
z0=(double)topk*dpth/coeff2l;
xm=(double)(topi+1.0)*wdth/coeff2l;
ym=(double)(topj+1.0)*hght/coeff2l;
zm=(double)(topk+1.0)*dpth/coeff2l;

norm=(xm-x0)*(ym-y0)*(zm-z0);


// recopie des paramètres relatif au bloc considéré 
	
	if ((topi<topMax)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
		ax111=param[3*l];ay111=param[3*l+1];az111=param[3*l+2];
		}
	else
		{
		ax111=0.0;ay111=0.0;az111=0.0;
		}
		
	if ((topi>0)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi-1,topj,topk,nb_param)/3;
		ax011=param[3*l];ay011=param[3*l+1];az011=param[3*l+2];
		}
	else
		{
		ax011=0.0;ay011=0.0;az011=0.0;
		}

	if ((topj>0)&&(topi<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi,topj-1,topk,nb_param)/3;
		ax101=param[3*l];ay101=param[3*l+1];az101=param[3*l+2];
		}
	else
		{
		ax101=0.0;ay101=0.0;az101=0.0;
		}


	if ((topk>0)&&(topi<topMax)&&(topj<topMax))
		{
		l=TOP_conv_ind(topi,topj,topk-1,nb_param)/3;
		ax110=param[3*l];ay110=param[3*l+1];az110=param[3*l+2];
		}
	else
		{
		ax110=0.0;ay110=0.0;az110=0.0;
		}

	if ((topi>0)&&(topj>0)&&(topk<topMax))
		{
		l=TOP_conv_ind(topi-1,topj-1,topk,nb_param)/3;
		ax001=param[3*l];ay001=param[3*l+1];az001=param[3*l+2];
		}
	else
		{
		ax001=0.0;ay001=0.0;az001=0.0;
		}


	if ((topi>0)&&(topk>0)&&(topj<topMax))
		{
		l=TOP_conv_ind(topi-1,topj,topk-1,nb_param)/3;
		ax010=param[3*l];ay010=param[3*l+1];az010=param[3*l+2];
		}
	else
		{
		ax010=0.0;ay010=0.0;az010=0.0;
		}


	if ((topj>0)&&(topk>0)&&(topi<topMax))
		{
		l=TOP_conv_ind(topi,topj-1,topk-1,nb_param)/3;
		ax100=param[3*l];ay100=param[3*l+1];az100=param[3*l+2];
		}
	else
		{
		ax100=0.0;ay100=0.0;az100=0.0;
		}


	if ((topi>0)&&(topj>0)&&(topk>0))
		{
		l=TOP_conv_ind(topi-1,topj-1,topk-1,nb_param)/3;
		ax000=param[3*l];ay000=param[3*l+1];az000=param[3*l+2];
		}
	else
		{
		ax000=0.0;ay000=0.0;az000=0.0;
		}



//----- Mise à jour de ParamInvBspline
ParamInvBspline.resol=resol;
ParamInvBspline.topi=topi;
ParamInvBspline.topj=topj;
ParamInvBspline.topk=topk;
ParamInvBspline.x0=x0;
ParamInvBspline.y0=y0;
ParamInvBspline.z0=z0;
ParamInvBspline.xm=xm;
ParamInvBspline.ym=ym;
ParamInvBspline.zm=zm;

ParamInvBspline.ax000=ax000;
ParamInvBspline.ax001=ax001;
ParamInvBspline.ax010=ax010;
ParamInvBspline.ax100=ax100;
ParamInvBspline.ax110=ax110;
ParamInvBspline.ax101=ax101;
ParamInvBspline.ax011=ax011;
ParamInvBspline.ax111=ax111;
									
ParamInvBspline.ay000=ay000;
ParamInvBspline.ay001=ay001;
ParamInvBspline.ay010=ay010;
ParamInvBspline.ay100=ay100;
ParamInvBspline.ay110=ay110;
ParamInvBspline.ay101=ay101;
ParamInvBspline.ay011=ay011;
ParamInvBspline.ay111=ay111;
	
ParamInvBspline.az000=az000;
ParamInvBspline.az001=az001;
ParamInvBspline.az010=az010;
ParamInvBspline.az100=az100;
ParamInvBspline.az110=az110;
ParamInvBspline.az101=az101;
ParamInvBspline.az011=az011;
ParamInvBspline.az111=az111;

ParamInvBspline.norm=(xm-x0)*(ym-y0)*(zm-z0);


if ((i==132)&&(j==183)&&(k==123))
	nb_iter=0;

nb_iter=0;
result=1;
do
{
pourcentage=grignote_boite_suivant_x(&ParamInvBspline,myB,i,j,k);
tmp=grignote_boite_suivant_y(&ParamInvBspline,myB,i,j,k);
pourcentage=tmp*pourcentage;
tmp=grignote_boite_suivant_z(&ParamInvBspline,myB,i,j,k);
pourcentage=tmp*pourcentage;
nb_iter++;
result=result*pourcentage;

if (pourcentage<pourcentage_reduction)
	{
	distance2=0.5*sqrt((myB->xM-myB->xm)*(myB->xM-myB->xm)+(myB->yM-myB->ym)*(myB->yM-myB->ym)+(myB->zM-myB->zm)*(myB->zM-myB->zm));

	xmid= (myB->xm+ myB->xM)/2.0;ymid= (myB->ym+ myB->yM)/2.0;zmid= (myB->zm+ myB->zM)/2.0;
	eval_transfo_bspline1_BoiteInvBspline_3d(&ParamInvBspline,  xmid, ymid, zmid, &hmid);
	distance=sqrt((hmid.x-i)*(hmid.x-i)+(hmid.y-j)*(hmid.y-j)+(hmid.z-k)*(hmid.z-k));	
	}
else
	distance=100*precis;	
	
}
while((pourcentage<pourcentage_reduction)&&((distance>precis)||(distance2>precis))&&(pourcentage>0)&&(nb_iter<20));

//if (pourcentage<1)
//printf("Le grignotage a fontionné  en %d iteration sur la boite %d %d %d  et la reduction est de %f \n",nb_iter,(int)i,(int)j,(int)k,result);   


if ((distance<precis)&&(distance2<precis))
	result=-1;

if (pourcentage==0)
	result=2;

	
 return(result);
 }

/******************************************************************************
*********
*********  ----------------- comp_double
*********
*********		Comparaison de double pour la fonction de tri qsort
*********				
********************************************************************************/
int comp_double(const void *num1, const void *num2)
{
  if ((double *)num1 <  (double *)num2) return -1;
  if ((double *)num1 == (double *)num2) return  0;
  if ((double *)num1 >  (double *)num2) return  1;
  
  return(1); //pour eviter un warning a la compilation
}


/******************************************************************************
*********
*********  ----------------- grignote_boite_suivant_x
*********
*********		Grignotte une boite suivant x
*********		valeur de retour : pourcentage de réduction de la boite
*********				
********************************************************************************/
double grignote_boite_suivant_x(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k)
{
double size_depart,size_arrive,pourcentage;
int is_xm,is_xM, test;

	
size_depart=myB->xM-myB->xm;	 

is_xm = test_if_x_suivant_hx_is_updatable(ParamInvBspline, myB, i, myB->xm);
is_xM = test_if_x_suivant_hx_is_updatable(ParamInvBspline, myB, i, myB->xM);
test=is_xm*is_xM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_xm!=0)||(is_xM!=0))
 if (update_x_suivant_hx(ParamInvBspline, myB, i,  is_xm,  is_xM)==0)
 	return(0);

is_xm = test_if_x_suivant_hy_is_updatable(ParamInvBspline, myB, j, myB->xm);
is_xM = test_if_x_suivant_hy_is_updatable(ParamInvBspline, myB, j, myB->xM);
test=is_xm*is_xM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_xm!=0)||(is_xM!=0))
	if(update_x_suivant_hy(ParamInvBspline, myB, j,  is_xm,  is_xM)==0)
	return(0);

is_xm = test_if_x_suivant_hz_is_updatable(ParamInvBspline, myB, k, myB->xm);
is_xM = test_if_x_suivant_hz_is_updatable(ParamInvBspline, myB, k, myB->xM);
test=is_xm*is_xM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_xm!=0)||(is_xM!=0))
	if(update_x_suivant_hz(ParamInvBspline, myB, k,  is_xm,  is_xM)==0)
	return(0);

	
size_arrive=myB->xM-myB->xm;

if (myB->xm>myB->xM)
	printf("Il y a un problème pour x!!!\n");

pourcentage=1.0-(size_depart-size_arrive)/size_depart;


return(pourcentage);
}



//-------------------------------------------------------------------------------------------------------------------
double grignote_boite_suivant_x2(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k)
{
double size_depart,size_arrive,pourcentage,x1,x2,x3;
double tabx[12];
int l,stop;

size_depart=myB->xM-myB->xm;

//---- on considère les valeurs 4 coins des faces d'abscisse x

grignote_eval_proposition_x(ParamInvBspline,myB,myB->ym, myB->zm,i,j,k, &x1, &x2, &x3);
tabx[0]=x1;tabx[1]=x2;tabx[2]=x3;

grignote_eval_proposition_x(ParamInvBspline,myB,myB->ym, myB->zM,i,j,k,  &x1, &x2, &x3);
tabx[3]=x1;tabx[4]=x2;tabx[5]=x3;

grignote_eval_proposition_x(ParamInvBspline,myB,myB->yM, myB->zm,i,j,k,  &x1, &x2, &x3);
tabx[6]=x1;tabx[7]=x2;tabx[8]=x3;

grignote_eval_proposition_x(ParamInvBspline,myB,myB->yM, myB->zM,i,j,k, &x1, &x2, &x3);
tabx[9]=x1;tabx[10]=x2;tabx[11]=x3;

qsort(tabx,12,sizeof(double),comp_double);
//tri_rapide_double(tabx, -1, 12);

l=0;stop=0;

while ((stop==0)&&(l<12))
	{
	if  ((tabx[l]>-1)&&(tabx[l]<myB->xM)&&(fabs(tabx[l]-myB->xm)>2.0*PRECISION_NUMERIQUE_INV_BSPLINE))
		{
		if (grignote_boite_update_xm(ParamInvBspline, myB, tabx[l] ,  i, j, k)==0)
			stop=1;
		
		}
	l++;	
	}
	

l=11;stop=0;

while ((stop==0)&&(l>=0))
	{
	if ((tabx[l]>-1)&&(tabx[l]>myB->xm)&&(fabs(tabx[l]-myB->xM)>2.0*PRECISION_NUMERIQUE_INV_BSPLINE))
		{
		if (grignote_boite_update_xM(ParamInvBspline, myB, tabx[l] ,  i, j, k)==0)
			stop=1;
		
		}
	l--;	
	}
	
/*if (!isinf(xm)&&(xm>myB->xm))
	{
	grignote_boite_update_xm(ParamInvBspline, myB, xm ,  i, j, k);
	}

	

if (!isinf(xM)&&(xM<myB->xM))
	{
	grignote_boite_update_xM(ParamInvBspline, myB, xM ,  i, j, k);
	}*/
	
size_arrive=myB->xM-myB->xm;

if (myB->xm>myB->xM)
	printf("Il y a un problème pour x!!!\n");

pourcentage=1.0-(size_depart-size_arrive)/size_depart;
return(pourcentage);
}


/******************************************************************************
*********
*********  ----------------- grignote_boite_update_xm
*********
*********		Met a jour xm si la valeur x propose est admissible
*********				
********************************************************************************/
int grignote_boite_update_xm(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k)
{
dvector3d h[8];
double xmin,xmax,ymin,ymax,zmin,zmax;
int l,retour=0;

//Pour ce prémunir des erreurs de precision 
x=x-PRECISION_NUMERIQUE_INV_BSPLINE;

/* on verifie si xm est effectivement une borne inf */
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->ym, myB->zm, &h[0]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->yM, myB->zm, &h[1]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->ym, myB->zM, &h[2]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->yM, myB->zM, &h[3]);
	 
	 
h[4]=myB->hx[0];
h[5]=myB->hx[2];
h[6]=myB->hx[4];
h[7]=myB->hx[6];
		 
xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;

	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,h[l].x);ymin=MINI(ymin,h[l].y);zmin=MINI(zmin,h[l].z);
	xmax=MAXI(xmax,h[l].x);ymax=MAXI(ymax,h[l].y);zmax=MAXI(zmax,h[l].z);
	} 

//xmin=xmin+PRECISION_NUMERIQUE_INV_BSPLINE;ymin=ymin+PRECISION_NUMERIQUE_INV_BSPLINE;zmin=zmin+PRECISION_NUMERIQUE_INV_BSPLINE;
//xmax=xmax-PRECISION_NUMERIQUE_INV_BSPLINE;ymax=ymax-PRECISION_NUMERIQUE_INV_BSPLINE;zmax=zmax-PRECISION_NUMERIQUE_INV_BSPLINE;
	
	if (!((xmin<=i)&&(xmax>=i)&&(ymin<=j)&&(ymax>=j)&&(zmin<=k)&&(zmax>=k)))
	{myB->xm=x;
	 myB->hx[0] =h[0];
 	 myB->hx[2] =h[1];
 	 myB->hx[4] =h[2];
 	 myB->hx[6] =h[3];
	 retour=1;
	}

return(retour);
				
}


/******************************************************************************
*********
*********  ----------------- grignote_boite_update_xM
*********
*********		Met a jour xM si la valeur x propose est admissible
*********				
********************************************************************************/
int grignote_boite_update_xM(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x , double i,double j,double k)
{
dvector3d h[8];
double xmin,xmax,ymin,ymax,zmin,zmax;
int l,retour=0;

//Pour ce prémunir des erreurs de precision 
x=x+PRECISION_NUMERIQUE_INV_BSPLINE;

/* on verifie si xM est effectivement une borne inf */
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->ym, myB->zm, &h[0]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->yM, myB->zm, &h[1]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->ym, myB->zM, &h[2]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  x, myB->yM, myB->zM, &h[3]);
	 
h[4]=myB->hx[1];
h[5]=myB->hx[3];
h[6]=myB->hx[5];
h[7]=myB->hx[7];
	 
xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;

	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,h[l].x);ymin=MINI(ymin,h[l].y);zmin=MINI(zmin,h[l].z);
	xmax=MAXI(xmax,h[l].x);ymax=MAXI(ymax,h[l].y);zmax=MAXI(zmax,h[l].z);
	} 
		
//xmin=xmin+PRECISION_NUMERIQUE_INV_BSPLINE;ymin=ymin+PRECISION_NUMERIQUE_INV_BSPLINE;zmin=zmin+PRECISION_NUMERIQUE_INV_BSPLINE;
//xmax=xmax-PRECISION_NUMERIQUE_INV_BSPLINE;ymax=ymax-PRECISION_NUMERIQUE_INV_BSPLINE;zmax=zmax-PRECISION_NUMERIQUE_INV_BSPLINE;

	if (!((xmin<=i)&&(xmax>=i)&&(ymin<=j)&&(ymax>=j)&&(zmin<=k)&&(zmax>=k)))
	{myB->xM=x;
	 myB->hx[1]=h[0];
	 myB->hx[3]=h[1];
	 myB->hx[5]=h[2];
	 myB->hx[7]=h[3];
	 retour=1;
	}


return(retour);
				
}


/******************************************************************************
*********
*********  ----------------- grignote_boite_suivant_y
*********
*********		Grignotte une boite suivant y
*********		valeur de retour : pourcentage de réduction de la boite
*********				
********************************************************************************/
double grignote_boite_suivant_y(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k)
{
double size_depart,size_arrive,pourcentage;
int is_ym,is_yM, test;

size_depart=myB->yM-myB->ym;	 


is_ym = test_if_y_suivant_hx_is_updatable(ParamInvBspline, myB, i, myB->ym);
is_yM = test_if_y_suivant_hx_is_updatable(ParamInvBspline, myB, i, myB->yM);

test=is_ym*is_yM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_ym!=0)||(is_yM!=0))
	if(update_y_suivant_hx(ParamInvBspline, myB, i,  is_ym,  is_yM)==0)
	return(0);

is_ym = test_if_y_suivant_hy_is_updatable(ParamInvBspline, myB, j, myB->ym);
is_yM = test_if_y_suivant_hy_is_updatable(ParamInvBspline, myB, j, myB->yM);
test=is_ym*is_yM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_ym!=0)||(is_yM!=0))
	if(update_y_suivant_hy(ParamInvBspline, myB, j,  is_ym,  is_yM)==0)
	return(0);

is_ym = test_if_y_suivant_hz_is_updatable(ParamInvBspline, myB, k, myB->ym);
is_yM = test_if_y_suivant_hz_is_updatable(ParamInvBspline, myB, k, myB->yM);
test=is_ym*is_yM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_ym!=0)||(is_yM!=0))
	if(update_y_suivant_hz(ParamInvBspline, myB, k,  is_ym,  is_yM)==0)
	return(0);
	
	
size_arrive=myB->yM-myB->ym;

if (myB->ym>myB->yM)
	printf("Il y a un problème pour y!!!\n");

pourcentage=1.0-(size_depart-size_arrive)/size_depart;
return(pourcentage);
}

//---------------------------------------------------------------------------------

double grignote_boite_suivant_y2(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k)
{
double size_depart,size_arrive,pourcentage,x1,x2,x3;
double tabx[12];
int l,stop;

size_depart=myB->yM-myB->ym;

//---- on considère les valeurs 4 coins des faces d'abscisse y

grignote_eval_proposition_y(ParamInvBspline,myB,myB->xm, myB->zm,i,j,k, &x1, &x2, &x3);
tabx[0]=x1;tabx[1]=x2;tabx[2]=x3;

grignote_eval_proposition_y(ParamInvBspline,myB,myB->xm, myB->zM,i,j,k, &x1, &x2, &x3);
tabx[3]=x1;tabx[4]=x2;tabx[5]=x3;

grignote_eval_proposition_y(ParamInvBspline,myB,myB->xM, myB->zm,i,j,k, &x1, &x2, &x3);
tabx[6]=x1;tabx[7]=x2;tabx[8]=x3;

grignote_eval_proposition_y(ParamInvBspline,myB,myB->xM, myB->zM,i,j,k, &x1, &x2, &x3);
tabx[9]=x1;tabx[10]=x2;tabx[11]=x3;


//tri_rapide_double(tabx, -1, 12);
qsort(tabx,12,sizeof(double),comp_double);

l=0;stop=0;

while ((stop==0)&&(l<12))
	{
	if ((tabx[l]>-1)&&(tabx[l]<myB->yM)&&(fabs(tabx[l]-myB->ym)>2.0*PRECISION_NUMERIQUE_INV_BSPLINE))
		{
		if (grignote_boite_update_ym(ParamInvBspline, myB, tabx[l] ,  i, j, k)==0)
			stop=1;
		
		}
	l++;	
	}
	

l=11;stop=0;

while ((stop==0)&&(l>=0))
	{
	if ((tabx[l]>-1)&&(tabx[l]>myB->ym)&&(fabs(tabx[l]-myB->yM)>2.0*PRECISION_NUMERIQUE_INV_BSPLINE))
		{
		if (grignote_boite_update_yM(ParamInvBspline, myB, tabx[l] ,  i, j, k)==0)
			stop=1;
		
		}
	l--;	
	}
		
/*if (!isinf(ym)&&(ym>myB->ym))
	{
	grignote_boite_update_ym(ParamInvBspline, myB, ym ,  i, j, k);	
	}
	
	
	
	
	if (!isinf(yM)&&(yM<myB->yM))
	{
	grignote_boite_update_yM(ParamInvBspline, myB, yM ,  i, j, k);	
	}*/


if (myB->ym>myB->yM)
	printf("Il y a un problème pour y!!!\n");
	
		
size_arrive=myB->yM-myB->ym;

pourcentage=1.0-(size_depart-size_arrive)/size_depart;
return(pourcentage);
}

/******************************************************************************
*********
*********  ----------------- grignote_boite_update_ym
*********
*********		Met a jour ym si la valeur y propose est admissible
*********				
********************************************************************************/
int grignote_boite_update_ym(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y , double i,double j,double k)
{
dvector3d h[8];
double xmin,xmax,ymin,ymax,zmin,zmax;
int l,retour=0;

//Pour ce prémunir des erreurs de precision 
y=y-PRECISION_NUMERIQUE_INV_BSPLINE;

/* on verifie si xm est effectivement une borne inf */
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xm, y, myB->zm, &h[0]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xM, y, myB->zm, &h[1]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xm, y, myB->zM, &h[2]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xM, y, myB->zM, &h[3]);
	 
h[4]=myB->hx[0];
h[5]=myB->hx[1];
h[6]=myB->hx[4];
h[7]=myB->hx[5];
	 
xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;

	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,h[l].x);ymin=MINI(ymin,h[l].y);zmin=MINI(zmin,h[l].z);
	xmax=MAXI(xmax,h[l].x);ymax=MAXI(ymax,h[l].y);zmax=MAXI(zmax,h[l].z);
	} 

//xmin=xmin+PRECISION_NUMERIQUE_INV_BSPLINE;ymin=ymin+PRECISION_NUMERIQUE_INV_BSPLINE;zmin=zmin+PRECISION_NUMERIQUE_INV_BSPLINE;
//xmax=xmax-PRECISION_NUMERIQUE_INV_BSPLINE;ymax=ymax-PRECISION_NUMERIQUE_INV_BSPLINE;zmax=zmax-PRECISION_NUMERIQUE_INV_BSPLINE;
				
	if (!((xmin<=i)&&(xmax>=i)&&(ymin<=j)&&(ymax>=j)&&(zmin<=k)&&(zmax>=k)))
	{myB->ym=y;
	myB->hx[0]=h[0];
	myB->hx[1]=h[1];
	myB->hx[4]=h[2];
	myB->hx[5]=h[3];
	retour=1;
	}

return(retour);
				
}


/******************************************************************************
*********
*********  ----------------- grignote_boite_update_yM
*********
*********		Met a jour yM si la valeur y propose est admissible
*********				
********************************************************************************/
int grignote_boite_update_yM(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y , double i,double j,double k)
{
dvector3d h[8];
double xmin,xmax,ymin,ymax,zmin,zmax;
int l,retour=0;

//Pour ce prémunir des erreurs de precision 
y=y+PRECISION_NUMERIQUE_INV_BSPLINE;

/* on verifie si xM est effectivement une borne inf */
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xm,  y, myB->zm, &h[0]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xM,  y, myB->zm, &h[1]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xm,  y, myB->zM, &h[2]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xM,  y, myB->zM, &h[3]);
	 
h[4]=myB->hx[2];
h[5]=myB->hx[3];
h[6]=myB->hx[6];
h[7]=myB->hx[7];
	 
xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;


	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,h[l].x);ymin=MINI(ymin,h[l].y);zmin=MINI(zmin,h[l].z);
	xmax=MAXI(xmax,h[l].x);ymax=MAXI(ymax,h[l].y);zmax=MAXI(zmax,h[l].z);
	} 

//xmin=xmin+PRECISION_NUMERIQUE_INV_BSPLINE;ymin=ymin+PRECISION_NUMERIQUE_INV_BSPLINE;zmin=zmin+PRECISION_NUMERIQUE_INV_BSPLINE;
//xmax=xmax-PRECISION_NUMERIQUE_INV_BSPLINE;ymax=ymax-PRECISION_NUMERIQUE_INV_BSPLINE;zmax=zmax-PRECISION_NUMERIQUE_INV_BSPLINE;

	if (!((xmin<=i)&&(xmax>=i)&&(ymin<=j)&&(ymax>=j)&&(zmin<=k)&&(zmax>=k)))
	{myB->yM=y;
	myB->hx[2]=h[0];
	myB->hx[3]=h[1];
	myB->hx[6]=h[2];
	myB->hx[7]=h[3];
	retour=1;
	}
return(retour);
				
}

/******************************************************************************
*********
*********  ----------------- grignote_boite_suivant_z
*********
*********		Grignotte une boite suivant z
*********		valeur de retour : pourcentage de réduction de la boite
*********				
********************************************************************************/
double grignote_boite_suivant_z(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k)
{
double size_depart,size_arrive,pourcentage;
int is_zm,is_zM, test;

size_depart=myB->zM-myB->zm;	 

is_zm = test_if_z_suivant_hx_is_updatable(ParamInvBspline, myB, i, myB->zm);
is_zM = test_if_z_suivant_hx_is_updatable(ParamInvBspline, myB, i, myB->zM);
test=is_zm*is_zM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);


if ((is_zm!=0)||(is_zM!=0))
	if(update_z_suivant_hx(ParamInvBspline, myB, i,  is_zm,  is_zM)==0)
	return(0);

is_zm = test_if_z_suivant_hy_is_updatable(ParamInvBspline, myB, j, myB->zm);
is_zM = test_if_z_suivant_hy_is_updatable(ParamInvBspline, myB, j, myB->zM);
test=is_zm*is_zM;

if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_zm!=0)||(is_zM!=0))
	if(update_z_suivant_hy(ParamInvBspline, myB, j,  is_zm,  is_zM)==0)
	return(0);

is_zm = test_if_z_suivant_hz_is_updatable(ParamInvBspline, myB, k, myB->zm);
is_zM = test_if_z_suivant_hz_is_updatable(ParamInvBspline, myB, k, myB->zM);
test=is_zm*is_zM;


if (test==1) // la solution n'est pas sur cette boite
	return(0);

if ((is_zm!=0)||(is_zM!=0))
	if(update_z_suivant_hz(ParamInvBspline, myB, k,  is_zm,  is_zM)==0)
	return(0);

	
size_arrive=myB->zM-myB->zm;

if (myB->zm>myB->zM)
	printf("Il y a un problème pour z!!!\n");

pourcentage=1.0-(size_depart-size_arrive)/size_depart;
return(pourcentage);
}

//--------------------------------------------------------------------------------
double grignote_boite_suivant_z2(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i,double j,double k)
{
double size_depart,size_arrive,pourcentage,x1,x2,x3;
double tabx[12];
int l,stop;

size_depart=myB->zM-myB->zm;

//---- on considère les valeurs 4 coins des faces d'abscisse z

grignote_eval_proposition_z(ParamInvBspline,myB,myB->xm, myB->ym,i,j,k, &x1, &x2, &x3);
tabx[0]=x1;tabx[1]=x2;tabx[2]=x3;

grignote_eval_proposition_z(ParamInvBspline,myB,myB->xm, myB->yM,i,j,k, &x1, &x2, &x3);
tabx[3]=x1;tabx[4]=x2;tabx[5]=x3;

grignote_eval_proposition_z(ParamInvBspline,myB,myB->xM, myB->ym,i,j,k, &x1, &x2, &x3);
tabx[6]=x1;tabx[7]=x2;tabx[8]=x3;

grignote_eval_proposition_z(ParamInvBspline,myB,myB->xM, myB->yM,i,j,k, &x1, &x2, &x3);
tabx[9]=x1;tabx[10]=x2;tabx[11]=x3;


qsort(tabx,12,sizeof(double),comp_double);
//tri_rapide_double(tabx, -1, 12);

l=0;stop=0;

while ((stop==0)&&(l<12))
	{
	if ((tabx[l]>-1)&&(tabx[l]<myB->zM)&&(fabs(tabx[l]-myB->zm)>2.0*PRECISION_NUMERIQUE_INV_BSPLINE))
		{
		if (grignote_boite_update_zm(ParamInvBspline, myB, tabx[l] ,  i, j, k)==0)
			stop=1;
		
		}
	l++;	
	}
	

l=11;stop=0;

while ((stop==0)&&(l>=0))
	{
	if ((tabx[l]>-1)&&(tabx[l]>myB->zm)&&(fabs(tabx[l]-myB->zM)>2.0*PRECISION_NUMERIQUE_INV_BSPLINE))
		{
		if (grignote_boite_update_zM(ParamInvBspline, myB, tabx[l] ,  i, j, k)==0)
			stop=1;
		
		}
	l--;	
	}
	

if (myB->zm>myB->zM)
	printf("Il y a un problème pour z!!!\n");

/*if (!isinf(zm)&&(zm>myB->zm))
	{
	grignote_boite_update_zm(ParamInvBspline, myB, zm ,  i, j, k);
	}
	
	

if (!isinf(zM)&&(zM<myB->zM))
	{
	grignote_boite_update_zM(ParamInvBspline, myB, zM ,  i, j, k);
	}*/


size_arrive=myB->zM-myB->zm;

pourcentage=1.0-(size_depart-size_arrive)/size_depart;
return(pourcentage);
}

/******************************************************************************
*********
*********  ----------------- grignote_boite_update_zm
*********
*********		Met a jour zm si la valeur z propose est admissible
*********				
********************************************************************************/
int grignote_boite_update_zm(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double z , double i,double j,double k)
{
dvector3d h[8];
double xmin,xmax,ymin,ymax,zmin,zmax;
int l,retour=0;

//Pour ce prémunir des erreurs de precision 
z=z-PRECISION_NUMERIQUE_INV_BSPLINE;

/* on verifie si xm est effectivement une borne inf */
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline, myB->xm, myB->ym, z, &h[0]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline, myB->xm, myB->yM, z, &h[1]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline, myB->xM, myB->ym, z, &h[2]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline, myB->xM, myB->yM, z, &h[3]);
	 
h[4]=myB->hx[0];
h[5]=myB->hx[2];
h[6]=myB->hx[1];
h[7]=myB->hx[3];

xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;


	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,h[l].x);ymin=MINI(ymin,h[l].y);zmin=MINI(zmin,h[l].z);
	xmax=MAXI(xmax,h[l].x);ymax=MAXI(ymax,h[l].y);zmax=MAXI(zmax,h[l].z);
	} 
	
//xmin=xmin+PRECISION_NUMERIQUE_INV_BSPLINE;ymin=ymin+PRECISION_NUMERIQUE_INV_BSPLINE;zmin=zmin+PRECISION_NUMERIQUE_INV_BSPLINE;
//xmax=xmax-PRECISION_NUMERIQUE_INV_BSPLINE;ymax=ymax-PRECISION_NUMERIQUE_INV_BSPLINE;zmax=zmax-PRECISION_NUMERIQUE_INV_BSPLINE;



	if (!((xmin<=i)&&(xmax>=i)&&(ymin<=j)&&(ymax>=j)&&(zmin<=k)&&(zmax>=k)))
	{myB->zm=z;
	myB->hx[0]=h[0];
	myB->hx[2]=h[1];
	myB->hx[1]=h[2];
	myB->hx[3]=h[3];
	retour=1;
	}

return(retour);

				
}


/******************************************************************************
*********
*********  ----------------- grignote_boite_update_zM
*********
*********		Met a jour zM si la valeur z propose est admissible
*********				
********************************************************************************/
int grignote_boite_update_zM(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double z , double i,double j,double k)
{
dvector3d h[8];
double xmin,xmax,ymin,ymax,zmin,zmax;
int l,retour=0;

//Pour ce prémunir des erreurs de precision 
z=z+PRECISION_NUMERIQUE_INV_BSPLINE;

/* on verifie si xM est effectivement une borne inf */
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xm, myB->ym, z, &h[0]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xm, myB->yM, z, &h[1]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xM, myB->ym, z, &h[2]);
eval_transfo_bspline1_BoiteInvBspline_3d(ParamInvBspline,  myB->xM, myB->yM, z, &h[3]);
	 
h[4]=myB->hx[4];
h[5]=myB->hx[6];
h[6]=myB->hx[5];
h[7]=myB->hx[7];
	 
xmin=HUGE_VAL;ymin=HUGE_VAL;zmin=HUGE_VAL;
xmax=-HUGE_VAL;ymax=-HUGE_VAL;zmax=-HUGE_VAL;

	for (l=0;l<8;l++)
	{
	xmin=MINI(xmin,h[l].x);ymin=MINI(ymin,h[l].y);zmin=MINI(zmin,h[l].z);
	xmax=MAXI(xmax,h[l].x);ymax=MAXI(ymax,h[l].y);zmax=MAXI(zmax,h[l].z);
	} 
	
//xmin=xmin+PRECISION_NUMERIQUE_INV_BSPLINE;ymin=ymin+PRECISION_NUMERIQUE_INV_BSPLINE;zmin=zmin+PRECISION_NUMERIQUE_INV_BSPLINE;
//xmax=xmax-PRECISION_NUMERIQUE_INV_BSPLINE;ymax=ymax-PRECISION_NUMERIQUE_INV_BSPLINE;zmax=zmax-PRECISION_NUMERIQUE_INV_BSPLINE;

	if (!((xmin<=i)&&(xmax>=i)&&(ymin<=j)&&(ymax>=j)&&(zmin<=k)&&(zmax>=k)))
	{myB->zM=z;
	myB->hx[4]=h[0];
	myB->hx[6]=h[1];
	myB->hx[5]=h[2];
	myB->hx[7]=h[3];
	retour=1;
	}
return(retour);
				
}


/******************************************************************************
*********
*********  ----------------- grignote_eval_proposition_x
*********
*********		Propose un intervale xmin et xmax correspondant à un couple (y,z)
*********		
********************************************************************************/

void grignote_eval_proposition_x(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z,double i,double j,double k, double *x1, double *x2, double *x3)
{
double new_borne, borne_min, borne_max,tmp;

borne_min=HUGE_VAL;
borne_max=-HUGE_VAL;

//----- equation x+ux=i

tmp=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)
								+(ParamInvBspline->ax100-ParamInvBspline->ax000)*(ParamInvBspline->ym- y)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ax101-ParamInvBspline->ax001)*(ParamInvBspline->ym- y)*( z-ParamInvBspline->z0)
								+(ParamInvBspline->ax110-ParamInvBspline->ax010)*( y-ParamInvBspline->y0)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ax111-ParamInvBspline->ax011)*( y-ParamInvBspline->y0)*( z-ParamInvBspline->z0));
								
								
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*i
						-(ParamInvBspline->ax000*ParamInvBspline->xm-ParamInvBspline->ax100*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax001*ParamInvBspline->xm-ParamInvBspline->ax101*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ax010*ParamInvBspline->xm-ParamInvBspline->ax110*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax011*ParamInvBspline->xm-ParamInvBspline->ax111*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(z-ParamInvBspline->z0))/tmp;
	
					
	if ((new_borne>myB->xm)&&(new_borne<myB->xM))
		{
		//borne_max=MAXI(borne_max,new_borne);
		//borne_min=MINI(borne_min,new_borne);
		*x1=new_borne;
		}
	else
		*x1=-1;
	}
else
	*x1=-1;



//----- equation y+uy=j

tmp= (
								(ParamInvBspline->ay100-ParamInvBspline->ay000)*(ParamInvBspline->ym- y)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ay101-ParamInvBspline->ay001)*(ParamInvBspline->ym- y)*( z-ParamInvBspline->z0)
								+(ParamInvBspline->ay110-ParamInvBspline->ay010)*( y-ParamInvBspline->y0)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ay111-ParamInvBspline->ay011)*( y-ParamInvBspline->y0)*( z-ParamInvBspline->z0));
								
if (tmp!=0)
	{
		new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(j-y)
						-(ParamInvBspline->ay000*ParamInvBspline->xm-ParamInvBspline->ay100*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay001*ParamInvBspline->xm-ParamInvBspline->ay101*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ay010*ParamInvBspline->xm-ParamInvBspline->ay110*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay011*ParamInvBspline->xm-ParamInvBspline->ay111*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(z-ParamInvBspline->z0))/tmp;
	
				
		if ((new_borne>myB->xm)&&(new_borne<myB->xM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*x2=new_borne;
			}
		else
			*x2=-1;
	
	}
else
	*x2=-1;
	
	//----- equation z+uz=k

tmp= (
								(ParamInvBspline->az100-ParamInvBspline->az000)*(ParamInvBspline->ym- y)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->az101-ParamInvBspline->az001)*(ParamInvBspline->ym- y)*( z-ParamInvBspline->z0)
								+(ParamInvBspline->az110-ParamInvBspline->az010)*( y-ParamInvBspline->y0)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->az111-ParamInvBspline->az011)*( y-ParamInvBspline->y0)*( z-ParamInvBspline->z0));

if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(k-z)
						-(ParamInvBspline->az000*ParamInvBspline->xm-ParamInvBspline->az100*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az001*ParamInvBspline->xm-ParamInvBspline->az101*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->az010*ParamInvBspline->xm-ParamInvBspline->az110*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az011*ParamInvBspline->xm-ParamInvBspline->az111*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(z-ParamInvBspline->z0))/tmp;
					
		if ((new_borne>myB->xm)&&(new_borne<myB->xM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*x3=new_borne;
			}
			else
			*x3=-1;
	}
	else
	*x3=-1;
	
// *xmin=borne_min;
// *xmax=borne_max;

}


/******************************************************************************
*********
*********  ----------------- grignote_eval_proposition_y
*********
*********		Propose un intervale ymin et ymax correspondant à un couple (x,z)
*********		
********************************************************************************/

void grignote_eval_proposition_y(BoiteInvBspline *ParamInvBspline, TOPTbox* myB, double x, double z,double i,double j,double k, double *y1, double *y2, double *y3)
{
double new_borne, borne_min, borne_max,tmp;

borne_min=HUGE_VAL;
borne_max=-HUGE_VAL;

//----- equation x+ux=i

tmp=		(
	 	(ParamInvBspline->ax010-ParamInvBspline->ax000)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ax011-ParamInvBspline->ax001)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
	 +(ParamInvBspline->ax110-ParamInvBspline->ax100)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ax111-ParamInvBspline->ax101)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0));
	 


if (tmp!=0)
	{

	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(i-x)
						-(ParamInvBspline->ax000*ParamInvBspline->ym-ParamInvBspline->ax010*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax001*ParamInvBspline->ym-ParamInvBspline->ax011*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ax100*ParamInvBspline->ym-ParamInvBspline->ax110*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax101*ParamInvBspline->ym-ParamInvBspline->ax111*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0))/tmp;

					
		if ((new_borne>myB->ym)&&(new_borne<myB->yM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*y1=new_borne;
			}
		else
			*y1=-1;
	}
	else
	*y1=-1;



//----- equation y+uy=j
tmp=		((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)
	 +(ParamInvBspline->ay010-ParamInvBspline->ay000)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ay011-ParamInvBspline->ay001)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
	 +(ParamInvBspline->ay110-ParamInvBspline->ay100)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ay111-ParamInvBspline->ay101)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0));
	 
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*j
						-(ParamInvBspline->ay000*ParamInvBspline->ym-ParamInvBspline->ay010*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay001*ParamInvBspline->ym-ParamInvBspline->ay011*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ay100*ParamInvBspline->ym-ParamInvBspline->ay110*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay101*ParamInvBspline->ym-ParamInvBspline->ay111*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0))/tmp;

					
	if ((new_borne>myB->ym)&&(new_borne<myB->yM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*y2=new_borne;
			}
			else 
			*y2=-1;
	}
	else
	*y2=-1;


//----- equation z+uz=k
tmp=		(
	 	(ParamInvBspline->az010-ParamInvBspline->az000)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->az011-ParamInvBspline->az001)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
	 +(ParamInvBspline->az110-ParamInvBspline->az100)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->az111-ParamInvBspline->az101)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0));
	 
	 
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(k-z)
						-(ParamInvBspline->az000*ParamInvBspline->ym-ParamInvBspline->az010*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az001*ParamInvBspline->ym-ParamInvBspline->az011*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->az100*ParamInvBspline->ym-ParamInvBspline->az110*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az101*ParamInvBspline->ym-ParamInvBspline->az111*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0))/tmp;

					


	if ((new_borne>myB->ym)&&(new_borne<myB->yM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*y3=new_borne;
			}
			else
			*y3=-1;
	}
	else 
	*y3=-1;

// *ymin=borne_min;
// *ymax=borne_max;

}

/******************************************************************************
*********
*********  ----------------- grignote_eval_proposition_z
*********
*********		Propose un intervale ymin et ymax correspondant à un couple (x,y)
*********		
********************************************************************************/

void grignote_eval_proposition_z(BoiteInvBspline *ParamInvBspline, TOPTbox* myB, double x, double y,double i,double j,double k, double *z1, double *z2, double *z3)
{
double new_borne, borne_min, borne_max,tmp;

borne_min=HUGE_VAL;
borne_max=-HUGE_VAL;



//----- equation x+ux=i
tmp=(
		(ParamInvBspline->ax001-ParamInvBspline->ax000)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ax011-ParamInvBspline->ax010)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
	 +(ParamInvBspline->ax101-ParamInvBspline->ax100)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ax111-ParamInvBspline->ax110)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0));
	 


if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(i-x)
						-(ParamInvBspline->ax000*ParamInvBspline->zm-ParamInvBspline->ax001*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ax010*ParamInvBspline->zm-ParamInvBspline->ax011*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
						-(ParamInvBspline->ax100*ParamInvBspline->zm-ParamInvBspline->ax101*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ax110*ParamInvBspline->zm-ParamInvBspline->ax111*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0))/tmp;
	
					
	if ((new_borne>myB->zm)&&(new_borne<myB->zM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*z1=new_borne;
			}
			else
			*z1=-1;
	}
	else
	*z1=-1;



//----- equation y+uy=j
tmp=(
		(ParamInvBspline->ay001-ParamInvBspline->ay000)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ay011-ParamInvBspline->ay010)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
	 +(ParamInvBspline->ay101-ParamInvBspline->ay100)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ay111-ParamInvBspline->ay110)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0));
								
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(j-y)
						-(ParamInvBspline->ay000*ParamInvBspline->zm-ParamInvBspline->ay001*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ay010*ParamInvBspline->zm-ParamInvBspline->ay011*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
						-(ParamInvBspline->ay100*ParamInvBspline->zm-ParamInvBspline->ay101*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ay110*ParamInvBspline->zm-ParamInvBspline->ay111*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0))/tmp;
	
					
	if ((new_borne>myB->zm)&&(new_borne<myB->zM))
			{
			//borne_max=MAXI(borne_max,new_borne);
			//borne_min=MINI(borne_min,new_borne);
			*z2=new_borne;
			}
			else
			*z2=-1;
	}
	else
	*z2=-1;
	

//----- equation z+uz=k
tmp=(
		(ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)
	 +(ParamInvBspline->az001-ParamInvBspline->az000)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->az011-ParamInvBspline->az010)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
	 +(ParamInvBspline->az101-ParamInvBspline->az100)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->az111-ParamInvBspline->az110)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0));
	 
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*k
						-(ParamInvBspline->az000*ParamInvBspline->zm-ParamInvBspline->az001*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->az010*ParamInvBspline->zm-ParamInvBspline->az011*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
						-(ParamInvBspline->az100*ParamInvBspline->zm-ParamInvBspline->az101*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->az110*ParamInvBspline->zm-ParamInvBspline->az111*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0))/tmp;
												

			if ((new_borne>myB->zm)&&(new_borne<myB->zM))
				{
				//borne_max=MAXI(borne_max,new_borne);
				//borne_min=MINI(borne_min,new_borne);
				*z3=new_borne;
				}
				else
				*z3=-1;
	}
	else
	*z3=-1;
	
// *zmin=borne_min;
// *zmax=borne_max;

}


/******************************************************************************
*********
*********  ----------------- update_bloc_transfo_bspline1_3d
*********
*********		Evaluation du champ de déformation aux 8 sommets du bloc
*********   
********************************************************************************/

void update_bloc_transfo_bspline1_3d(double *param, int nb_param, int wdth, int hght, int dpth, TOPTbox* myB)
{
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xm, myB->ym, myB->zm, &myB->hx[0]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xM, myB->ym, myB->zm, &myB->hx[1]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xm, myB->yM, myB->zm, &myB->hx[2]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xM, myB->yM, myB->zm, &myB->hx[3]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xm, myB->ym, myB->zM, &myB->hx[4]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xM, myB->ym, myB->zM, &myB->hx[5]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xm, myB->yM, myB->zM, &myB->hx[6]);
eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, myB->xM, myB->yM, myB->zM, &myB->hx[7]);
}


/******************************************************************************
*********
*********  ----------------- eval_deplacement_bspline1_BoiteInvBspline_3d
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

void eval_deplacement_bspline1_BoiteInvBspline_3d(BoiteInvBspline *ParamInvBspline, double x, double y, double z, dvector3d* u)
{
double norm, dxm,dym,dzm,dx0,dy0,dz0;

norm=(ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0);

dxm=ParamInvBspline->xm-x;
dym=ParamInvBspline->ym-y;
dzm=ParamInvBspline->zm-z;
dx0=x-ParamInvBspline->x0;
dy0=y-ParamInvBspline->y0;
dz0=z-ParamInvBspline->z0;


u->x=
 ParamInvBspline->ax000*dxm*dym*dzm+
 ParamInvBspline->ax100*dx0*dym*dzm+
 ParamInvBspline->ax010*dxm*dy0*dzm+
 ParamInvBspline->ax001*dxm*dym*dz0+
 ParamInvBspline->ax110*dx0*dy0*dzm+
 ParamInvBspline->ax011*dxm*dy0*dz0+
 ParamInvBspline->ax101*dx0*dym*dz0+
 ParamInvBspline->ax111*dx0*dy0*dz0;

u->x=u->x/norm;

 
u->y=
 ParamInvBspline->ay000*dxm*dym*dzm+
 ParamInvBspline->ay100*dx0*dym*dzm+
 ParamInvBspline->ay010*dxm*dy0*dzm+
 ParamInvBspline->ay001*dxm*dym*dz0+
 ParamInvBspline->ay110*dx0*dy0*dzm+
 ParamInvBspline->ay011*dxm*dy0*dz0+
 ParamInvBspline->ay101*dx0*dym*dz0+
 ParamInvBspline->ay111*dx0*dy0*dz0;


u->y=u->y/norm;


u->z=
 ParamInvBspline->az000*dxm*dym*dzm+
 ParamInvBspline->az100*dx0*dym*dzm+
 ParamInvBspline->az010*dxm*dy0*dzm+
 ParamInvBspline->az001*dxm*dym*dz0+
 ParamInvBspline->az110*dx0*dy0*dzm+
 ParamInvBspline->az011*dxm*dy0*dz0+
 ParamInvBspline->az101*dx0*dym*dz0+
 ParamInvBspline->az111*dx0*dy0*dz0;

u->z=u->z/norm;



}

/******************************************************************************
*********
*********  ----------------- eval_transfo_bspline1_BoiteInvBspline_3d
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

void eval_transfo_bspline1_BoiteInvBspline_3d(BoiteInvBspline *ParamInvBspline, double x, double y, double z, dvector3d* u)
{
eval_deplacement_bspline1_BoiteInvBspline_3d(ParamInvBspline, x, y,  z,  u);
u->x=u->x+x;
u->y=u->y+y;
u->z=u->z+z;
}




/******************************************************************************
*********
*********  ----------------- eval_deplacement_bspline1_BoiteInvBspline_3d_x
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

double eval_deplacement_bspline1_BoiteInvBspline_3d_x(BoiteInvBspline *ParamInvBspline, double x, double y, double z)
{
double  dxm,dym,dzm,dx0,dy0,dz0, res;

//norm=(ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0);


dxm=ParamInvBspline->xm-x;
dym=ParamInvBspline->ym-y;
dzm=ParamInvBspline->zm-z;
dx0=x-ParamInvBspline->x0;
dy0=y-ParamInvBspline->y0;
dz0=z-ParamInvBspline->z0;


res=
 ParamInvBspline->ax000*dxm*dym*dzm+
 ParamInvBspline->ax100*dx0*dym*dzm+
 ParamInvBspline->ax010*dxm*dy0*dzm+
 ParamInvBspline->ax001*dxm*dym*dz0+
 ParamInvBspline->ax110*dx0*dy0*dzm+
 ParamInvBspline->ax011*dxm*dy0*dz0+
 ParamInvBspline->ax101*dx0*dym*dz0+
 ParamInvBspline->ax111*dx0*dy0*dz0;

res=res/ParamInvBspline->norm;

return(res);
}
 

/******************************************************************************
*********
*********  ----------------- eval_transfo_bspline1_BoiteInvBspline_3d_x
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

double eval_transfo_bspline1_BoiteInvBspline_3d_x(BoiteInvBspline *ParamInvBspline, double x, double y, double z)
{
return(eval_deplacement_bspline1_BoiteInvBspline_3d_x(ParamInvBspline, x, y,  z)+x);
}



/******************************************************************************
*********
*********  ----------------- eval_deplacement_bspline1_BoiteInvBspline_3d_y
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

double eval_deplacement_bspline1_BoiteInvBspline_3d_y(BoiteInvBspline *ParamInvBspline, double x, double y, double z)
{
double  dxm,dym,dzm,dx0,dy0,dz0, res;

//norm=(ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0);

dxm=ParamInvBspline->xm-x;
dym=ParamInvBspline->ym-y;
dzm=ParamInvBspline->zm-z;
dx0=x-ParamInvBspline->x0;
dy0=y-ParamInvBspline->y0;
dz0=z-ParamInvBspline->z0;


res=
 ParamInvBspline->ay000*dxm*dym*dzm+
 ParamInvBspline->ay100*dx0*dym*dzm+
 ParamInvBspline->ay010*dxm*dy0*dzm+
 ParamInvBspline->ay001*dxm*dym*dz0+
 ParamInvBspline->ay110*dx0*dy0*dzm+
 ParamInvBspline->ay011*dxm*dy0*dz0+
 ParamInvBspline->ay101*dx0*dym*dz0+
 ParamInvBspline->ay111*dx0*dy0*dz0;


res=res/ParamInvBspline->norm;

return(res);
}
 

/******************************************************************************
*********
*********  ----------------- eval_transfo_bspline1_BoiteInvBspline_3d_y
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

double eval_transfo_bspline1_BoiteInvBspline_3d_y(BoiteInvBspline *ParamInvBspline, double x, double y, double z)
{
return(eval_deplacement_bspline1_BoiteInvBspline_3d_y(ParamInvBspline, x, y,  z)+y);
}

/******************************************************************************
*********
*********  ----------------- eval_deplacement_bspline1_BoiteInvBspline_3d_z
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

double eval_deplacement_bspline1_BoiteInvBspline_3d_z(BoiteInvBspline *ParamInvBspline, double x, double y, double z)
{
double  dxm,dym,dzm,dx0,dy0,dz0, res;

//norm=(ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0);

dxm=ParamInvBspline->xm-x;
dym=ParamInvBspline->ym-y;
dzm=ParamInvBspline->zm-z;
dx0=x-ParamInvBspline->x0;
dy0=y-ParamInvBspline->y0;
dz0=z-ParamInvBspline->z0;


res=
  ParamInvBspline->az000*dxm*dym*dzm+
 ParamInvBspline->az100*dx0*dym*dzm+
 ParamInvBspline->az010*dxm*dy0*dzm+
 ParamInvBspline->az001*dxm*dym*dz0+
 ParamInvBspline->az110*dx0*dy0*dzm+
 ParamInvBspline->az011*dxm*dy0*dz0+
 ParamInvBspline->az101*dx0*dym*dz0+
 ParamInvBspline->az111*dx0*dy0*dz0;



res=res/ParamInvBspline->norm;

return(res);
}
 

/******************************************************************************
*********
*********  ----------------- eval_transfo_bspline1_BoiteInvBspline_3d_z
*********
*********		Evaluation en un point à coordonnées réelles d'un champ de déformation
*********    Bspline de degré 1
********************************************************************************/

double eval_transfo_bspline1_BoiteInvBspline_3d_z(BoiteInvBspline *ParamInvBspline, double x, double y, double z)
{
return(eval_deplacement_bspline1_BoiteInvBspline_3d_z(ParamInvBspline, x, y,  z)+z);
}



/*******************************************************************************
**     base_to_field_3d(nb_param,param,champ)                          
*/
/*!                                                                   
**     transformation liee a la base de fonctions ondelettes
**     \param nb_param : nbre d'element dans param
**     \param param : tableau de parametres
**     \param champ : (E/S)
**     \param imref,imreca : pas utilise pour le moment (a voir pour les transfos anisotropiques)
**     \retval : 1
*******************************************************************************/
int ParamBspline_to_field_3d(int nb_param, double *param, field3d *champ, grphic3d *imref, grphic3d *imreca)
{   
  int i,j,k,l;
  int wdth,hght,dpth;
  double dxreca,dyreca,dzreca;
  double dxref,dyref,dzref;
 	//double *param_voxel;
	dvector3d u;
	
	
		
  vector3d ***data;
 
  wdth=champ->width;hght=champ->height;dpth=champ->depth;
  data=champ->raw;
  /*mise a zero du champ*/
  /*for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     {data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;}*/
    
 
  if (imreca)
  {
    dxreca = imreca->dx;
    dyreca = imreca->dy;
    dzreca = imreca->dz;
  }
  else
  {
    dxreca = 1.0;
    dyreca = 1.0;
    dzreca = 1.0;
  }
	if (imref)
	{
    if (wdth != (int)imref->width || hght != (int)imref->height || dpth != (int)imref->depth)
      return 0;
		dxref = imref->dx;	  
		dyref = imref->dy;	  
		dzref = imref->dz;	  
	}
	else
	{
		   dxref = 1.0;
		   dyref = 1.0;
		   dzref = 1.0;
	}

if ((dxref!=dxreca)||(dyref!=dyreca)||(dzref!=dzreca))
		{   
		PUT_ERROR("ERREUR dans base_to_field_3d : le champ est l'image n'ont pas les meme dx,dy,dz\n");
 	   printf("ERREUR dans base_to_field_3d : le champ est l'image n'ont pas les meme dx,dy,dz\n");
		 return(-1);
  	}


/*param_voxel=(double *)malloc(nb_param*sizeof(double));
 
 if (param_voxel==NULL)
 	{
	printf("Erreur d'allocation dans ParamBspline_to_field_3d \n");
	exit(-1);
	}*/
 
 
 for (l=0;l<nb_param/3;l++)
		{
		param[3*l]=param[3*l]/dxreca;
		param[3*l+1]=param[3*l+1]/dyreca;
		param[3*l+2]=param[3*l+2]/dzreca;
		}

 for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		{
		
		eval_deplacement_bspline1_3d(param, nb_param, wdth, hght, dpth, i, j,  k,  &u);
 
		data[i][j][k].x=u.x;
		data[i][j][k].y=u.y;
		data[i][j][k].z=u.z;
		}

for (l=0;l<nb_param/3;l++)
		{
		param[3*l]=param[3*l]*dxreca;
		param[3*l+1]=param[3*l+1]*dyreca;
		param[3*l+2]=param[3*l+2]*dzreca;
		}
 
// free(param_voxel);
 
 
 
  return(1);
}

/******************************************************************************
*********
*********  ----------------- inv_bspline_find_solution_with_fixed_point_algo
*********
*********		Essaie de trouver l'inverse sur une boite avec la méthode
*********		du point fixe
*********    Retourne 1 si une solution avec une précision prec est trouve, 0 sinon
********************************************************************************/


 int inv_bspline_find_solution_with_fixed_point_algo(double *param, int nb_param, int wdth, int hght, int dpth, double precis, double i, double j, double k, TOPTbox* myB, dvector3d* res, double *distance)
 {
 dvector3d xcourant,  xnew;
 int nb_iter = 0,u,v;
 double **J,F[3],*prod,**invJ;
int resol,coeff2l,topi,topj,topk,l,topMax;
double x0,y0,z0,xm,ym,zm,norm;
double ax000,ax001,ax010,ax100,ax110,ax101,ax011,ax111;
double ay000,ay001,ay010,ay100,ay110,ay101,ay011,ay111;
double az000,az001,az010,az100,az110,az101,az011,az111;
BoiteInvBspline ParamInvBspline;
mat3d A,invA;

  J=alloc_dmatrix(3,3);
 invJ=alloc_dmatrix(3,3);



//---- Initialisation des paramètres caractéristiques du bloc myB considéré  -----

if (INV_BSPLINE_RESOL!=-1)
	resol=INV_BSPLINE_RESOL;
else
	resol=floor(log (pow(nb_param/3.0,1.0/3.0)+1)/log(2.0)+0.5);

if (INV_BSPLINE_COEFF2L!=-1)
	coeff2l=INV_BSPLINE_COEFF2L;
else
	coeff2l=pow(2.0,resol);
	
topMax=coeff2l-1;


topi=floor(myB->xm/wdth*coeff2l);

if (myB->xm==wdth)
	topi=topi-1;

topj=floor(myB->ym/hght*coeff2l);

if (myB->ym==hght)
	topj=topj-1;

topk=floor(myB->zm/dpth*coeff2l);

if (myB->zm==dpth)
	topk=topk-1;

x0=(double)topi*wdth/coeff2l;
y0=(double)topj*hght/coeff2l;
z0=(double)topk*dpth/coeff2l;
xm=(double)(topi+1.0)*wdth/coeff2l;
ym=(double)(topj+1.0)*hght/coeff2l;
zm=(double)(topk+1.0)*dpth/coeff2l;

norm=(xm-x0)*(ym-y0)*(zm-z0);


/* recopie des paramètres relatif au bloc considéré */
	
	if ((topi<topMax)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj,topk,topMax)/3;
		ax111=param[3*l];ay111=param[3*l+1];az111=param[3*l+2];
		}
	else
		{
		ax111=0.0;ay111=0.0;az111=0.0;
		}
		
	if ((topi>0)&&(topj<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi-1,topj,topk,topMax)/3;
		ax011=param[3*l];ay011=param[3*l+1];az011=param[3*l+2];
		}
	else
		{
		ax011=0.0;ay011=0.0;az011=0.0;
		}

	if ((topj>0)&&(topi<topMax)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj-1,topk,topMax)/3;
		ax101=param[3*l];ay101=param[3*l+1];az101=param[3*l+2];
		}
	else
		{
		ax101=0.0;ay101=0.0;az101=0.0;
		}


	if ((topk>0)&&(topi<topMax)&&(topj<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj,topk-1,topMax)/3;
		ax110=param[3*l];ay110=param[3*l+1];az110=param[3*l+2];
		}
	else
		{
		ax110=0.0;ay110=0.0;az110=0.0;
		}

	if ((topi>0)&&(topj>0)&&(topk<topMax))
		{
		l=TOP_conv_ind_topD(topi-1,topj-1,topk,topMax)/3;
		ax001=param[3*l];ay001=param[3*l+1];az001=param[3*l+2];
		}
	else
		{
		ax001=0.0;ay001=0.0;az001=0.0;
		}


	if ((topi>0)&&(topk>0)&&(topj<topMax))
		{
		l=TOP_conv_ind_topD(topi-1,topj,topk-1,topMax)/3;
		ax010=param[3*l];ay010=param[3*l+1];az010=param[3*l+2];
		}
	else
		{
		ax010=0.0;ay010=0.0;az010=0.0;
		}


	if ((topj>0)&&(topk>0)&&(topi<topMax))
		{
		l=TOP_conv_ind_topD(topi,topj-1,topk-1,topMax)/3;
		ax100=param[3*l];ay100=param[3*l+1];az100=param[3*l+2];
		}
	else
		{
		ax100=0.0;ay100=0.0;az100=0.0;
		}


	if ((topi>0)&&(topj>0)&&(topk>0))
		{
		l=TOP_conv_ind_topD(topi-1,topj-1,topk-1,topMax)/3;
		ax000=param[3*l];ay000=param[3*l+1];az000=param[3*l+2];
		}
	else
		{
		ax000=0.0;ay000=0.0;az000=0.0;
		}



//----- Mise à jour de ParamInvBspline
ParamInvBspline.resol=resol;
ParamInvBspline.topi=topi;
ParamInvBspline.topj=topj;
ParamInvBspline.topk=topk;
ParamInvBspline.x0=x0;
ParamInvBspline.y0=y0;
ParamInvBspline.z0=z0;
ParamInvBspline.xm=xm;
ParamInvBspline.ym=ym;
ParamInvBspline.zm=zm;

ParamInvBspline.ax000=ax000;
ParamInvBspline.ax001=ax001;
ParamInvBspline.ax010=ax010;
ParamInvBspline.ax100=ax100;
ParamInvBspline.ax110=ax110;
ParamInvBspline.ax101=ax101;
ParamInvBspline.ax011=ax011;
ParamInvBspline.ax111=ax111;
									
ParamInvBspline.ay000=ay000;
ParamInvBspline.ay001=ay001;
ParamInvBspline.ay010=ay010;
ParamInvBspline.ay100=ay100;
ParamInvBspline.ay110=ay110;
ParamInvBspline.ay101=ay101;
ParamInvBspline.ay011=ay011;
ParamInvBspline.ay111=ay111;
	
ParamInvBspline.az000=az000;
ParamInvBspline.az001=az001;
ParamInvBspline.az010=az010;
ParamInvBspline.az100=az100;
ParamInvBspline.az110=az110;
ParamInvBspline.az101=az101;
ParamInvBspline.az011=az011;
ParamInvBspline.az111=az111;

ParamInvBspline.norm=(xm-x0)*(ym-y0)*(zm-z0);



 // On initialise la methode avec le milieu de la boite
 xcourant.x=(myB->xM+myB->xm)/2.0;
 xcourant.y=(myB->yM+myB->ym)/2.0;
 xcourant.z=(myB->zM+myB->zm)/2.0;
 
 
 do
 {

 // on teste si le point est dans la boite, sinon on quitte le programme

/*if ((xcourant.x<myB->xm)||(xcourant.x>myB->xM))
  	{
	if (deblocx<=0)
		{
		xcourant.x=MINI(xcourant.x,myB->xM-PRECISION_NUMERIQUE_INV_BSPLINE);
		xcourant.x=MAXI(xcourant.x,myB->xm+PRECISION_NUMERIQUE_INV_BSPLINE);
		
		deblocx=2;
		}
	else
		{free_dmatrix(J, 3, 3);
		return(0);}
		
	}

deblocx--;

if ((xcourant.y<myB->ym)||(xcourant.y>myB->yM))
  	{
	if (deblocy<=0)
		{
		xcourant.y=MINI(xcourant.y,myB->yM-PRECISION_NUMERIQUE_INV_BSPLINE);
		xcourant.y=MAXI(xcourant.y,myB->ym+PRECISION_NUMERIQUE_INV_BSPLINE);
		
		deblocy=2;
		}
	else
		{free_dmatrix(J, 3, 3);
		return(0);}		
	}

deblocy--;

if ((xcourant.z<myB->zm)||(xcourant.z>myB->zM))
  	{
	if (deblocz<=0)
		{
		xcourant.z=MINI(xcourant.z,myB->zM-PRECISION_NUMERIQUE_INV_BSPLINE);
		xcourant.z=MAXI(xcourant.z,myB->zm+PRECISION_NUMERIQUE_INV_BSPLINE);
		
		deblocz=2;
		}
	else
		{free_dmatrix(J, 3, 3);
		return(0);}
	}	

deblocz--;
*/

if ((xcourant.x<myB->xm)||(xcourant.x>myB->xM)||(xcourant.y<myB->ym)||(xcourant.y>myB->yM)||(xcourant.z<myB->zm)||(xcourant.z>myB->zM))
  	{
	free_dmatrix(J, 3, 3);
	free_dmatrix(invJ, 3, 3);
	return(0);
	}

// eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, xcourant.x, xcourant.y, xcourant.z, &xnew);
eval_transfo_bspline1_BoiteInvBspline_3d(&ParamInvBspline,   xcourant.x, xcourant.y, xcourant.z, &xnew);

// on teste si le nouveau point est une solution a la precision prec
*distance=sqrt((xnew.x-i)*(xnew.x-i)+(xnew.y-j)*(xnew.y-j)+(xnew.z-k)*(xnew.z-k));

	if (*distance<precis)
		{
		//if (nb_iter>2)
		//	printf("ca converge...\n");
		res->x=xcourant.x;
		res->y=xcourant.y;
		res->z=xcourant.z;
		free_dmatrix(J, 3, 3);
		free_dmatrix(invJ, 3, 3);
		return(1);
		
		}
F[0]=xnew.x-i;
F[1]=xnew.y-j;
F[2]=xnew.z-k;
	
//eval_jacobien_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, xcourant.x, xcourant.y, xcourant.z, J);
eval_jacobien_transfo_bspline1_BoiteInvBspline_3d(&ParamInvBspline, xcourant.x, xcourant.y, xcourant.z, J);



for (u=0;u<3;u++)
for (v=0;v<3;v++)
	A.H[u][v]=J[u][v];

inversion_mat3d(&A, &invA);

for (u=0;u<3;u++)
for (v=0;v<3;v++)
	invJ[u][v]=invA.H[u][v];


//invJ=matrix_inversion(J, 3);


prod = matrix_multiply_vector(invJ, 3,3, F);


// on calcul le prochain point		
// xcourant.x= xcourant.x-xnew.x+i;
// xcourant.y= xcourant.y-xnew.y+j;
// xcourant.z= xcourant.z-xnew.z+k;

xcourant.x= xcourant.x-prod[0];
xcourant.y= xcourant.y-prod[1];
xcourant.z= xcourant.z-prod[2];


free_dvector(prod,3); 
//free_dmatrix(invJ, 3, 3);
 nb_iter++;
 }
 while (nb_iter<100);
 
 printf("Nombre d'iteration depassé dans inv_bspline_find_solution_with_fixed_point_algo pour %f %f %f \n",i,j,k);
 free_dmatrix(J, 3, 3);
 free_dmatrix(invJ, 3, 3);

return(0); 
 }


/******************************************************************************
*********
*********  ----------------- inv_bspline_find_solution_with_fixed_point_algo_global
*********
*********		Essaie de trouver l'inverse sur une boite avec la méthode
*********		du point fixe
*********    Retourne 1 si une solution avec une précision prec est trouve, 0 sinon
********************************************************************************/


 int inv_bspline_find_solution_with_fixed_point_algo_global(double *param, int nb_param, int wdth, int hght, int dpth, double precis, double i, double j, double k, TOPTbox* myB, dvector3d* res, double *distance)
 {
 dvector3d xcourant,  xnew;
 int nb_iter = 0;
 double **J,F[3],*prod,**invJ;


  J=alloc_dmatrix(3,3);


 // On initialise la methode avec le milieu de la boite
 xcourant.x=(myB->xM+myB->xm)/2.0;
 xcourant.y=(myB->yM+myB->ym)/2.0;
 xcourant.z=(myB->zM+myB->zm)/2.0;
 
myB->xm=0;myB->xM=wdth;
myB->ym=0;myB->yM=hght;
myB->zm=0;myB->zM=dpth;

 
 xcourant.x=res->x;
 xcourant.y=res->y;
 xcourant.z=res->z;

 do
 {

 // on teste si le point est dans la boite, sinon on quitte le programme

/*if ((xcourant.x<myB->xm)||(xcourant.x>myB->xM))
  	{
	if (deblocx<=0)
		{
		xcourant.x=MINI(xcourant.x,myB->xM-PRECISION_NUMERIQUE_INV_BSPLINE);
		xcourant.x=MAXI(xcourant.x,myB->xm+PRECISION_NUMERIQUE_INV_BSPLINE);
		
		deblocx=2;
		}
	else
		{free_dmatrix(J, 3, 3);
		return(0);}
		
	}

deblocx--;

if ((xcourant.y<myB->ym)||(xcourant.y>myB->yM))
  	{
	if (deblocy<=0)
		{
		xcourant.y=MINI(xcourant.y,myB->yM-PRECISION_NUMERIQUE_INV_BSPLINE);
		xcourant.y=MAXI(xcourant.y,myB->ym+PRECISION_NUMERIQUE_INV_BSPLINE);
		
		deblocy=2;
		}
	else
		{free_dmatrix(J, 3, 3);
		return(0);}		
	}

deblocy--;

if ((xcourant.z<myB->zm)||(xcourant.z>myB->zM))
  	{
	if (deblocz<=0)
		{
		xcourant.z=MINI(xcourant.z,myB->zM-PRECISION_NUMERIQUE_INV_BSPLINE);
		xcourant.z=MAXI(xcourant.z,myB->zm+PRECISION_NUMERIQUE_INV_BSPLINE);
		
		deblocz=2;
		}
	else
		{free_dmatrix(J, 3, 3);
		return(0);}
	}	

deblocz--;
*/

if ((xcourant.x<myB->xm)||(xcourant.x>myB->xM)||(xcourant.y<myB->ym)||(xcourant.y>myB->yM)||(xcourant.z<myB->zm)||(xcourant.z>myB->zM))
  	{
	free_dmatrix(J, 3, 3);
	return(0);
	}

 eval_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, xcourant.x, xcourant.y, xcourant.z, &xnew);

// on teste si le nouveau point est une solution a la precision prec
*distance=sqrt((xnew.x-i)*(xnew.x-i)+(xnew.y-j)*(xnew.y-j)+(xnew.z-k)*(xnew.z-k));

	if (*distance<precis)
		{
		//if (nb_iter>2)
		//	printf("ca converge...\n");
		res->x=xcourant.x;
		res->y=xcourant.y;
		res->z=xcourant.z;
		free_dmatrix(J, 3, 3);
		return(1);
		
		}
F[0]=xnew.x-i;
F[1]=xnew.y-j;
F[2]=xnew.z-k;
	
eval_jacobien_transfo_bspline1_3d(param, nb_param, wdth, hght, dpth, xcourant.x, xcourant.y, xcourant.z, J);

invJ=matrix_inversion(J, 3);
prod = matrix_multiply_vector(invJ, 3,3, F);


// on calcul le prochain point		
// xcourant.x= xcourant.x-xnew.x+i;
// xcourant.y= xcourant.y-xnew.y+j;
// xcourant.z= xcourant.z-xnew.z+k;

xcourant.x= xcourant.x-prod[0];
xcourant.y= xcourant.y-prod[1];
xcourant.z= xcourant.z-prod[2];


free_dvector(prod,3); 
free_dmatrix(invJ, 3, 3);
 nb_iter++;
 }
 while (nb_iter<100);
 
 printf("Nombre d'iteration depassé dans inv_bspline_find_solution_with_fixed_point_algo pour %f %f %f \n",i,j,k);
 free_dmatrix(J, 3, 3);
return(0); 
 }

/******************************************************************************
*********
*********  ----------------- test_if_x_suivant_hx_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_x_suivant_hx_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double x)
{
double infhx, suphx,tmp;

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,  x,  myB->ym,  myB->zm);
infhx=tmp;suphx=tmp;

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,  x,  myB->yM,  myB->zm);
infhx=MINI(infhx,tmp);
suphx=MAXI(suphx,tmp);

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,  x,  myB->ym,  myB->zM);
infhx=MINI(infhx,tmp);
suphx=MAXI(suphx,tmp);

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,  x,  myB->yM,  myB->zM);
infhx=MINI(infhx,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphx=MAXI(suphx,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhx>i)
	return(1);

if (suphx<i)
	return(-1);
	
return(0);	
}

/******************************************************************************
*********
*********  ----------------- test_if_x_suivant_hy_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_x_suivant_hy_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double j, double x)
{
double infhy, suphy,tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  x,  myB->ym,  myB->zm);
infhy=tmp;suphy=tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  x,  myB->yM,  myB->zm);
infhy=MINI(infhy,tmp);
suphy=MAXI(suphy,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  x,  myB->ym,  myB->zM);
infhy=MINI(infhy,tmp);
suphy=MAXI(suphy,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  x,  myB->yM,  myB->zM);
infhy=MINI(infhy,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphy=MAXI(suphy,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhy>j)
	return(1);

if (suphy<j)
	return(-1);
	
return(0);

}

/******************************************************************************
*********
*********  ----------------- test_if_x_suivant_hz_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_x_suivant_hz_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double k, double x)
{
double infhz, suphz,tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  x,  myB->ym,  myB->zm);
infhz=tmp;suphz=tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  x,  myB->yM,  myB->zm);
infhz=MINI(infhz,tmp);
suphz=MAXI(suphz,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  x,  myB->ym,  myB->zM);
infhz=MINI(infhz,tmp);
suphz=MAXI(suphz,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  x,  myB->yM,  myB->zM);
infhz=MINI(infhz,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphz=MAXI(suphz,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhz>k)
	return(1);

if (suphz<k)
	return(-1);
	
return(0);
}


/******************************************************************************
*********
*********  ----------------- test_if_y_suivant_hx_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_y_suivant_hx_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double y)
{
double infhx, suphx,tmp;

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xm,y, myB->zm);
infhx=tmp;suphx=tmp;

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xM,y, myB->zm);
infhx=MINI(infhx,tmp);
suphx=MAXI(suphx,tmp);

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xm,y, myB->zM);
infhx=MINI(infhx,tmp);
suphx=MAXI(suphx,tmp);

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xM, y, myB->zM);
infhx=MINI(infhx,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphx=MAXI(suphx,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhx>i)
	return(1);

if (suphx<i)
	return(-1);
	
return(0);
}

/******************************************************************************
*********
*********  ----------------- test_if_y_suivant_hy_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_y_suivant_hy_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double j, double y)
{
double infhy, suphy,tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline, myB->xm,y, myB->zm );
infhy=tmp;suphy=tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline, myB->xM,y, myB->zm );
infhy=MINI(infhy,tmp);
suphy=MAXI(suphy,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline, myB->xm,y, myB->zM );
infhy=MINI(infhy,tmp);
suphy=MAXI(suphy,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline, myB->xM, y, myB->zM);
infhy=MINI(infhy,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphy=MAXI(suphy,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhy>j)
	return(1);

if (suphy<j)
	return(-1);
	
return(0);
}

/******************************************************************************
*********
*********  ----------------- test_if_y_suivant_hz_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_y_suivant_hz_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double k, double y)
{
double infhz, suphz,tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xm,y, myB->zm);
infhz=tmp;suphz=tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xM,y, myB->zm);
infhz=MINI(infhz,tmp);
suphz=MAXI(suphz,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xm,y, myB->zM);
infhz=MINI(infhz,tmp);
suphz=MAXI(suphz,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xM, y, myB->zM);
infhz=MINI(infhz,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphz=MAXI(suphz,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhz>k)
	return(1);

if (suphz<k)
	return(-1);
	
return(0);	
}


/******************************************************************************
*********
*********  ----------------- test_if_z_suivant_hx_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_z_suivant_hx_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, double z)
{
double infhx, suphx,tmp;

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xm,myB->ym,z);
infhx=tmp;suphx=tmp;

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xM,myB->ym,z);
infhx=MINI(infhx,tmp);
suphx=MAXI(suphx,tmp);

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xm,myB->yM,z);
infhx=MINI(infhx,tmp);
suphx=MAXI(suphx,tmp);

tmp = eval_transfo_bspline1_BoiteInvBspline_3d_x(ParamInvBspline,   myB->xM,myB->yM,z);
infhx=MINI(infhx,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphx=MAXI(suphx,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhx>i)
	return(1);

if (suphx<i)
	return(-1);
	
return(0);
}

/******************************************************************************
*********
*********  ----------------- test_if_z_suivant_hy_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_z_suivant_hy_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double j, double z)
{
double infhy, suphy,tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  myB->xm,myB->ym,z);
infhy=tmp;suphy=tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  myB->xM,myB->ym,z);
infhy=MINI(infhy,tmp);
suphy=MAXI(suphy,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  myB->xm,myB->yM,z);
infhy=MINI(infhy,tmp);
suphy=MAXI(suphy,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_y(ParamInvBspline,  myB->xM,myB->yM,z);
infhy=MINI(infhy,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphy=MAXI(suphy,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhy>j)
	return(1);

if (suphy<j)
	return(-1);
	
return(0);
}

/******************************************************************************
*********
*********  ----------------- test_if_z_suivant_hz_is_updatable
*********
*********		Test si x peut-etre mis a jour en considerant hx
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int test_if_z_suivant_hz_is_updatable(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double k, double z)
{
double infhz, suphz,tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xm,myB->ym,z);
infhz=tmp;suphz=tmp;

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xM,myB->ym,z);
infhz=MINI(infhz,tmp);
suphz=MAXI(suphz,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xm,myB->yM,z);
infhz=MINI(infhz,tmp);
suphz=MAXI(suphz,tmp);

tmp=eval_transfo_bspline1_BoiteInvBspline_3d_z(ParamInvBspline,  myB->xM,myB->yM,z);
infhz=MINI(infhz,tmp)-PRECISION_NUMERIQUE_INV_BSPLINE;
suphz=MAXI(suphz,tmp)+PRECISION_NUMERIQUE_INV_BSPLINE;

if (infhz>k)
	return(1);

if (suphz<k)
	return(-1);
	
return(0);
}




/******************************************************************************
*********
*********  ----------------- solve_x_suivant_hx
*********
*********		Trouve la solution x de l'equation hx(x,y,z)=i 
*********				
********************************************************************************/
double solve_x_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z, double i)
{
double tmp, new_borne;

tmp=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)
								+(ParamInvBspline->ax100-ParamInvBspline->ax000)*(ParamInvBspline->ym- y)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ax101-ParamInvBspline->ax001)*(ParamInvBspline->ym- y)*( z-ParamInvBspline->z0)
								+(ParamInvBspline->ax110-ParamInvBspline->ax010)*( y-ParamInvBspline->y0)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ax111-ParamInvBspline->ax011)*( y-ParamInvBspline->y0)*( z-ParamInvBspline->z0));
								
								
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*i
						-(ParamInvBspline->ax000*ParamInvBspline->xm-ParamInvBspline->ax100*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax001*ParamInvBspline->xm-ParamInvBspline->ax101*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ax010*ParamInvBspline->xm-ParamInvBspline->ax110*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax011*ParamInvBspline->xm-ParamInvBspline->ax111*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(z-ParamInvBspline->z0))/tmp;
	
					
	return(new_borne);

	}
else
	return(-1);

}


/******************************************************************************
*********
*********  ----------------- solve_x_suivant_hy
*********
*********		Trouve la solution x de l'equation hy(x,y,z)=j 
*********				
********************************************************************************/
double solve_x_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z, double j)
{
double tmp, new_borne;


tmp= (
								(ParamInvBspline->ay100-ParamInvBspline->ay000)*(ParamInvBspline->ym- y)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ay101-ParamInvBspline->ay001)*(ParamInvBspline->ym- y)*( z-ParamInvBspline->z0)
								+(ParamInvBspline->ay110-ParamInvBspline->ay010)*( y-ParamInvBspline->y0)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->ay111-ParamInvBspline->ay011)*( y-ParamInvBspline->y0)*( z-ParamInvBspline->z0));
								
if (tmp!=0)
	{
		new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(j-y)
						-(ParamInvBspline->ay000*ParamInvBspline->xm-ParamInvBspline->ay100*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay001*ParamInvBspline->xm-ParamInvBspline->ay101*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ay010*ParamInvBspline->xm-ParamInvBspline->ay110*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay011*ParamInvBspline->xm-ParamInvBspline->ay111*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(z-ParamInvBspline->z0))/tmp;
	
					
		return(new_borne);

	}
else
	return(-1);

}

/******************************************************************************
*********
*********  ----------------- solve_x_suivant_hz
*********
*********		Trouve la solution x de l'equation hz(x,y,z)=k 
*********				
********************************************************************************/
double solve_x_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double y, double z, double k)
{
double tmp, new_borne;

tmp= (
								(ParamInvBspline->az100-ParamInvBspline->az000)*(ParamInvBspline->ym- y)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->az101-ParamInvBspline->az001)*(ParamInvBspline->ym- y)*( z-ParamInvBspline->z0)
								+(ParamInvBspline->az110-ParamInvBspline->az010)*( y-ParamInvBspline->y0)*(ParamInvBspline->zm- z)
								+(ParamInvBspline->az111-ParamInvBspline->az011)*( y-ParamInvBspline->y0)*( z-ParamInvBspline->z0));

if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(k-z)
						-(ParamInvBspline->az000*ParamInvBspline->xm-ParamInvBspline->az100*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az001*ParamInvBspline->xm-ParamInvBspline->az101*ParamInvBspline->x0)*(ParamInvBspline->ym-y)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->az010*ParamInvBspline->xm-ParamInvBspline->az110*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az011*ParamInvBspline->xm-ParamInvBspline->az111*ParamInvBspline->x0)*(y-ParamInvBspline->y0)*(z-ParamInvBspline->z0))/tmp;
					
	
					
	return(new_borne);

	}
else
	return(-1);

}


/******************************************************************************
*********
*********  ----------------- solve_y_suivant_hx
*********
*********		Trouve la solution y de l'equation hx(x,y,z)=i 
*********				
********************************************************************************/
double solve_y_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z, double i)
{
double tmp, new_borne;


tmp=		(
	 	(ParamInvBspline->ax010-ParamInvBspline->ax000)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ax011-ParamInvBspline->ax001)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
	 +(ParamInvBspline->ax110-ParamInvBspline->ax100)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ax111-ParamInvBspline->ax101)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0));
	 


if (tmp!=0)
	{

	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(i-x)
						-(ParamInvBspline->ax000*ParamInvBspline->ym-ParamInvBspline->ax010*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax001*ParamInvBspline->ym-ParamInvBspline->ax011*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ax100*ParamInvBspline->ym-ParamInvBspline->ax110*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ax101*ParamInvBspline->ym-ParamInvBspline->ax111*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0))/tmp;

					
	return(new_borne);

	}
else
	return(-1);

}




/******************************************************************************
*********
*********  ----------------- solve_y_suivant_hy
*********
*********		Trouve la solution y de l'equation hy(x,y,z)=j 
*********				
********************************************************************************/
double solve_y_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z, double j)
{
double tmp, new_borne;


tmp=		((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)
	 +(ParamInvBspline->ay010-ParamInvBspline->ay000)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ay011-ParamInvBspline->ay001)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
	 +(ParamInvBspline->ay110-ParamInvBspline->ay100)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->ay111-ParamInvBspline->ay101)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0));
	 
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*j
						-(ParamInvBspline->ay000*ParamInvBspline->ym-ParamInvBspline->ay010*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay001*ParamInvBspline->ym-ParamInvBspline->ay011*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->ay100*ParamInvBspline->ym-ParamInvBspline->ay110*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->ay101*ParamInvBspline->ym-ParamInvBspline->ay111*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0))/tmp;

						
					
		return(new_borne);

	}
else
	return(-1);

}

/******************************************************************************
*********
*********  ----------------- solve_y_suivant_hz
*********
*********		Trouve la solution x de l'equation hz(x,y,z)=k 
*********				
********************************************************************************/
double solve_y_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double z, double k)
{
double tmp, new_borne;

tmp=		(
	 	(ParamInvBspline->az010-ParamInvBspline->az000)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->az011-ParamInvBspline->az001)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
	 +(ParamInvBspline->az110-ParamInvBspline->az100)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
	 +(ParamInvBspline->az111-ParamInvBspline->az101)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0));
	 
	 
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(k-z)
						-(ParamInvBspline->az000*ParamInvBspline->ym-ParamInvBspline->az010*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az001*ParamInvBspline->ym-ParamInvBspline->az011*ParamInvBspline->y0)*(ParamInvBspline->xm-x)*(z-ParamInvBspline->z0)
						-(ParamInvBspline->az100*ParamInvBspline->ym-ParamInvBspline->az110*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(ParamInvBspline->zm-z)
						-(ParamInvBspline->az101*ParamInvBspline->ym-ParamInvBspline->az111*ParamInvBspline->y0)*(x-ParamInvBspline->x0)*(z-ParamInvBspline->z0))/tmp;

	
					
	return(new_borne);

	}
else
	return(-1);

}


/******************************************************************************
*********
*********  ----------------- solve_z_suivant_hx
*********
*********		Trouve la solution z de l'equation hx(x,y,z)=i 
*********				
********************************************************************************/
double solve_z_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y, double i)
{
double tmp, new_borne;

tmp=(
		(ParamInvBspline->ax001-ParamInvBspline->ax000)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ax011-ParamInvBspline->ax010)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
	 +(ParamInvBspline->ax101-ParamInvBspline->ax100)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ax111-ParamInvBspline->ax110)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0));
	 


if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(i-x)
						-(ParamInvBspline->ax000*ParamInvBspline->zm-ParamInvBspline->ax001*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ax010*ParamInvBspline->zm-ParamInvBspline->ax011*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
						-(ParamInvBspline->ax100*ParamInvBspline->zm-ParamInvBspline->ax101*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ax110*ParamInvBspline->zm-ParamInvBspline->ax111*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0))/tmp;
	
					
	return(new_borne);

	}
else
	return(-1);

}


/******************************************************************************
*********
*********  ----------------- solve_z_suivant_hy
*********
*********		Trouve la solution z de l'equation hy(x,y,z)=j 
*********				
********************************************************************************/
double solve_z_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y, double j)
{
double tmp, new_borne;


tmp=(
		(ParamInvBspline->ay001-ParamInvBspline->ay000)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ay011-ParamInvBspline->ay010)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
	 +(ParamInvBspline->ay101-ParamInvBspline->ay100)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->ay111-ParamInvBspline->ay110)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0));
								
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*(j-y)
						-(ParamInvBspline->ay000*ParamInvBspline->zm-ParamInvBspline->ay001*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ay010*ParamInvBspline->zm-ParamInvBspline->ay011*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
						-(ParamInvBspline->ay100*ParamInvBspline->zm-ParamInvBspline->ay101*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->ay110*ParamInvBspline->zm-ParamInvBspline->ay111*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0))/tmp;
	
					
		return(new_borne);

	}
else
	return(-1);

}

/******************************************************************************
*********
*********  ----------------- solve_z_suivant_hz
*********
*********		Trouve la solution z de l'equation hz(x,y,z)=k 
*********				
********************************************************************************/
double solve_z_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double x, double y, double k)
{
double tmp, new_borne;

tmp=(

		(ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)
	 +(ParamInvBspline->az001-ParamInvBspline->az000)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->az011-ParamInvBspline->az010)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
	 +(ParamInvBspline->az101-ParamInvBspline->az100)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
	 +(ParamInvBspline->az111-ParamInvBspline->az110)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0));
	 
if (tmp!=0)
	{
	new_borne=((ParamInvBspline->xm-ParamInvBspline->x0)*(ParamInvBspline->ym-ParamInvBspline->y0)*(ParamInvBspline->zm-ParamInvBspline->z0)*k
						-(ParamInvBspline->az000*ParamInvBspline->zm-ParamInvBspline->az001*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->az010*ParamInvBspline->zm-ParamInvBspline->az011*ParamInvBspline->z0)*(ParamInvBspline->xm-x)*(y-ParamInvBspline->y0)
						-(ParamInvBspline->az100*ParamInvBspline->zm-ParamInvBspline->az101*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(ParamInvBspline->ym-y)
						-(ParamInvBspline->az110*ParamInvBspline->zm-ParamInvBspline->az111*ParamInvBspline->z0)*(x-ParamInvBspline->x0)*(y-ParamInvBspline->y0))/tmp;
												
					
	
					
	return(new_borne);

	}
else
	return(-1);

}



/******************************************************************************
*********
*********  ----------------- update_x_suivant_hx
*********
*********		Mais a jour xm et/ou xM en regardant les solutions de hx=i
*********		valeur de retour : 1 si oui, 0 sinon
*********				
**************************************************************************void update_x_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM)
******/
int update_x_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_xm, int is_xM)
{
double xmin=myB->xM, xmax=myB->xm, x;
double xm=myB->xm-PRECISION_NUMERIQUE_INV_BSPLINE,xM=myB->xM+PRECISION_NUMERIQUE_INV_BSPLINE;

int is_solution=0;

x = solve_x_suivant_hx(ParamInvBspline, myB, myB->ym, myB->zm,  i);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}


x = solve_x_suivant_hx(ParamInvBspline, myB, myB->yM, myB->zm,  i);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}


x = solve_x_suivant_hx(ParamInvBspline, myB, myB->ym, myB->zM,  i);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}
	

x = solve_x_suivant_hx(ParamInvBspline, myB, myB->yM, myB->zM,  i);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}

if (is_solution>0)
{
if (is_xm)
	myB->xm=MAXI(xmin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->xm);
	
if (is_xM)
	myB->xM=MINI(xmax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->xM);
}

return(is_solution);
	
}


/******************************************************************************
*********
*********  ----------------- update_x_suivant_hy
*********
*********		Mais a jour xm et/ou xM en regardant les solutions de hy=j
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_x_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double j, int is_xm, int is_xM)
{
double xmin=myB->xM, xmax=myB->xm, x;
double xm=myB->xm-PRECISION_NUMERIQUE_INV_BSPLINE,xM=myB->xM+PRECISION_NUMERIQUE_INV_BSPLINE;

int is_solution=0;

x = solve_x_suivant_hy(ParamInvBspline, myB, myB->ym, myB->zm,  j);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}


x = solve_x_suivant_hy(ParamInvBspline, myB, myB->yM, myB->zm,  j);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}


x = solve_x_suivant_hy(ParamInvBspline, myB, myB->ym, myB->zM,  j);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}
	

x = solve_x_suivant_hy(ParamInvBspline, myB, myB->yM, myB->zM,  j);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}

if (is_solution>0)
{
if (is_xm)
	myB->xm=MAXI(xmin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->xm);
	
if (is_xM)
	myB->xM=MINI(xmax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->xM);
}

return(is_solution);	
}

/******************************************************************************
*********
*********  ----------------- update_x_suivant_hz
*********
*********		Mais a jour xm et/ou xM en regardant les solutions de hz=k
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_x_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double k, int is_xm, int is_xM)
{
double xmin=myB->xM, xmax=myB->xm, x;
double xm=myB->xm-PRECISION_NUMERIQUE_INV_BSPLINE,xM=myB->xM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

x = solve_x_suivant_hz(ParamInvBspline, myB, myB->ym, myB->zm,  k);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}


x = solve_x_suivant_hz(ParamInvBspline, myB, myB->yM, myB->zm,  k);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}


x = solve_x_suivant_hz(ParamInvBspline, myB, myB->ym, myB->zM,  k);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}
	

x = solve_x_suivant_hz(ParamInvBspline, myB, myB->yM, myB->zM,  k);

if ((x>xm)&&(x<xM))
	{
	xmin=MINI(xmin, x);
	xmax=MAXI(xmax,x);
	is_solution++;
	}

if (is_solution>0)
{
if (is_xm)
	myB->xm=MAXI(xmin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->xm);
	
if (is_xM)
	myB->xM=MINI(xmax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->xM);

}	
return(is_solution);	
}


/******************************************************************************
*********
*********  ----------------- update_y_suivant_hx
*********
*********		Mais a jour ym et/ou yM en regardant les solutions de hx=i
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_y_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_ym, int is_yM)
{
double ymin=myB->yM, ymax=myB->ym, y;
double ym=myB->ym-PRECISION_NUMERIQUE_INV_BSPLINE,yM=myB->yM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

y = solve_y_suivant_hx(ParamInvBspline, myB, myB->xm, myB->zm,  i);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


y = solve_y_suivant_hx(ParamInvBspline, myB, myB->xM, myB->zm,  i);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


y = solve_y_suivant_hx(ParamInvBspline, myB, myB->xm, myB->zM,  i);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}
	

y = solve_y_suivant_hx(ParamInvBspline, myB, myB->xM, myB->zM,  i);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}

if (is_solution>0)
{
if (is_ym)
	myB->ym=MAXI(ymin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->ym);
	
if (is_yM)
	myB->yM=MINI(ymax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->yM);

}
return(is_solution);	
}


/******************************************************************************
*********
*********  ----------------- update_y_suivant_hy
*********
*********		Mais a jour ym et/ou yM en regardant les solutions de hy=j
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_y_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double j, int is_ym, int is_yM)
{
double ymin=myB->yM, ymax=myB->ym, y;
double ym=myB->ym-PRECISION_NUMERIQUE_INV_BSPLINE,yM=myB->yM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

y = solve_y_suivant_hy(ParamInvBspline, myB, myB->xm, myB->zm,  j);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


y = solve_y_suivant_hy(ParamInvBspline, myB, myB->xM, myB->zm,  j);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


y = solve_y_suivant_hy(ParamInvBspline, myB, myB->xm, myB->zM,  j);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}
	

y = solve_y_suivant_hy(ParamInvBspline, myB, myB->xM, myB->zM,  j);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


if (is_solution>0)
{
if (is_ym)
	myB->ym=MAXI(ymin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->ym);
	
if (is_yM)
	myB->yM=MINI(ymax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->yM);
}
return(is_solution);	
}

/******************************************************************************
*********
*********  ----------------- update_y_suivant_hz
*********
*********		Mais a jour ym et/ou yM en regardant les solutions de hz=k
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_y_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double k, int is_ym, int is_yM)
{
double ymin=myB->yM, ymax=myB->ym, y;
double ym=myB->ym-PRECISION_NUMERIQUE_INV_BSPLINE,yM=myB->yM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

y = solve_y_suivant_hz(ParamInvBspline, myB, myB->xm, myB->zm,  k);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


y = solve_y_suivant_hz(ParamInvBspline, myB, myB->xM, myB->zm,  k);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}


y = solve_y_suivant_hz(ParamInvBspline, myB, myB->xm, myB->zM,  k);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}
	

y = solve_y_suivant_hz(ParamInvBspline, myB, myB->xM, myB->zM,  k);

if ((y>ym)&&(y<yM))
	{
	ymin=MINI(ymin, y);
	ymax=MAXI(ymax,y);
	is_solution++;
	}

if (is_solution>0)
{
if (is_ym)
	myB->ym=MAXI(ymin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->ym);
	
if (is_yM)
	myB->yM=MINI(ymax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->yM);
}	
return(is_solution);	
}


/******************************************************************************
*********
*********  ----------------- update_z_suivant_hx
*********
*********		Mais a jour zm et/ou zM en regardant les solutions de hx=i
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_z_suivant_hx(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double i, int is_zm, int is_zM)
{
double zmin=myB->zM, zmax=myB->zm, z;
double zm=myB->zm-PRECISION_NUMERIQUE_INV_BSPLINE,zM=myB->zM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

z = solve_z_suivant_hx(ParamInvBspline, myB, myB->xm, myB->ym,  i);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}


z = solve_z_suivant_hx(ParamInvBspline, myB, myB->xM, myB->ym,  i);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}


z = solve_z_suivant_hx(ParamInvBspline, myB, myB->xm, myB->yM,  i);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}
	

z = solve_z_suivant_hx(ParamInvBspline, myB, myB->xM, myB->yM,  i);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}

if (is_solution>0)
{
if (is_zm)
	myB->zm=MAXI(zmin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->zm);
	
if (is_zM)
	myB->zM=MINI(zmax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->zM);
}	
return(is_solution);	
}


/******************************************************************************
*********
*********  ----------------- update_z_suivant_hy
*********
*********		Mais a jour zm et/ou zM en regardant les solutions de hy=j
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_z_suivant_hy(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double j, int is_zm, int is_zM)
{
double zmin=myB->zM, zmax=myB->zm, z;
double zm=myB->zm-PRECISION_NUMERIQUE_INV_BSPLINE,zM=myB->zM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

z = solve_z_suivant_hy(ParamInvBspline, myB, myB->xm, myB->ym,  j);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}


z = solve_z_suivant_hy(ParamInvBspline, myB, myB->xM, myB->ym,  j);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}


z = solve_z_suivant_hy(ParamInvBspline, myB, myB->xm, myB->yM,  j);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}
	

z = solve_z_suivant_hy(ParamInvBspline, myB, myB->xM, myB->yM,  j);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}

if (is_solution>0)
{
if (is_zm)
	myB->zm=MAXI(zmin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->zm);
	
if (is_zM)
	myB->zM=MINI(zmax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->zM);
}	
return(is_solution);	
}

/******************************************************************************
*********
*********  ----------------- update_z_suivant_hz
*********
*********		Mais a jour xm et/ou xM en regardant les solutions de hz=k
*********		valeur de retour : 1 si oui, 0 sinon
*********				
********************************************************************************/
int update_z_suivant_hz(BoiteInvBspline *ParamInvBspline, TOPTbox* myB,double k, int is_zm, int is_zM)
{
double zmin=myB->zM, zmax=myB->zm, z;
double zm=myB->zm-PRECISION_NUMERIQUE_INV_BSPLINE,zM=myB->zM+PRECISION_NUMERIQUE_INV_BSPLINE;
int is_solution=0;

z = solve_z_suivant_hz(ParamInvBspline, myB, myB->xm, myB->ym,  k);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}


z = solve_z_suivant_hz(ParamInvBspline, myB, myB->xM, myB->ym,  k);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}


z = solve_z_suivant_hz(ParamInvBspline, myB, myB->xm, myB->yM,  k);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}
	

z = solve_z_suivant_hz(ParamInvBspline, myB, myB->xM, myB->yM,  k);

if ((z>zm)&&(z<zM))
	{
	zmin=MINI(zmin, z);
	zmax=MAXI(zmax,z);
	is_solution++;
	}

if (is_solution>0)
{
if (is_zm)
	myB->zm=MAXI(zmin-PRECISION_NUMERIQUE_INV_BSPLINE,myB->zm);
	
if (is_zM)
	myB->zM=MINI(zmax+PRECISION_NUMERIQUE_INV_BSPLINE,myB->zM);
}	

return(is_solution);
	
}


/**************************************************************************
**      base_resol_up_light_3d(param,nb_param)
*/	
/*!      Passage de la base de fonction a la resolution superieure (avec moins d'allocation memoire)
**	Modifie les parametres en consequence :
**		nouvelle taille en resultat de la fonction
**		nouvelles valeurs dans le tableau param
**	\retval le nombre de parametres
**  \remark modifie la variable globale BASE
**************************************************************************/

int	base_resol_up_light_3d(double *param, int nb_param)
{
 int	i,j,k,l;
 int	pasx,pasy,pasz,wdth,hght,dpth;
 int	resol,nbint,tfilt,nb_func,nbfx,nbfy,nbfz,supx,supy,supz,sup;
 int	*x0,*x1,*y0,*y1,*z0,*z1;
 int	(*scal_func)();
 
  wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;
  
 /*parametre de la base dans les variables locales*/
  resol=BASE3D.resol;resol++; /*passage a la resolution superieure*/
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;
  scal_func=BASE3D.scal_func;
  tfilt=BASE3D.tfilt;
  
     nbfx=tfilt+2*nbfx-2;nbfy=tfilt+2*nbfy-2;nbfz=tfilt+2*nbfz-2;
 
   /*nombre d'intervalle de decoupage*/
    nbint=(int)(pow(2.0,(double)resol));
   /*pas suivant x y et z*/ 
    pasx=wdth/nbint;pasy=hght/nbint;pasz=dpth/nbint;
 
 
 if (scal_func==Bsplined0) sup = 1;
 if (scal_func==Bsplined1) sup = 2;
 if (scal_func==Bsplined2) sup = 3;
 
   /*mise a jour de l'echantillonage de la fonction d'echelle*/
    supx=sup;
    supy=sup;
    supz=sup;
   /*verification que ca correspond avec ce qui est deja calcule*/
    if (nbfx!=nbint+1-supx || nbfy!=nbint+1-supy || nbfz!=nbint+1-supz) printf("ERREUR 1 dans base_resol_up_3d\n");
    nb_func=nbfx*nbfy*nbfz;
  /*reallocation des element*/
    x0=REALLOC(x0,nb_func,int);
    y0=REALLOC(y0,nb_func,int);
    z0=REALLOC(z0,nb_func,int);  
    x1=REALLOC(x1,nb_func,int);
    y1=REALLOC(y1,nb_func,int);
    z1=REALLOC(z1,nb_func,int); 
  
    if (x0==NULL || x1==NULL || y0==NULL || y1==NULL || z0==NULL || z1==NULL )
     {printf("ERREUR d'allocation dans base_resol_up_3d\n");exit(1);}
  /*remplissage de la base*/
   l=0;
    for (i=0;i<nbfx;i++)
     for (j=0;j<nbfy;j++) 
      for (k=0;k<nbfz;k++)
       {
        x0[l]=i*pasx;x1[l]=i*pasx+supx*pasx;
	y0[l]=j*pasy;y1[l]=j*pasy+supy*pasy;
	z0[l]=k*pasz;z1[l]=k*pasz+supz*pasz;
	l++;
       }
   BASE3D.x0=x0;BASE3D.x1=x1;BASE3D.y0=y0;BASE3D.y1=y1;BASE3D.z0=z0;BASE3D.z1=z1;
   BASE3D.resol=resol;BASE3D.nb_func=nb_func;
   BASE3D.nbfx=nbfx;BASE3D.nbfy=nbfy;BASE3D.nbfz=nbfz;
   
  
 return(3*nb_func);
}

/**************************************************************************
**      fusion_transf_3d()
**  Fusionne 2 champs de deformation defini sur une certaine ROI
**************************************************************************/
void fusion_transf_3d()
{
 int im_1,im_2,err;
 char nomfichier1[PATH_LEN],nomfichier2[PATH_LEN],nomfichres[PATH_LEN],str[PATH_LEN];
 
  /*fichier contenant la transformation*/
  {
	sprintf(nomfichier1,"%s",GET_FILE("Premier fichier trf",&err));
	
    if (err)
      return ;
	  
	put_file_extension(nomfichier1,".trf",nomfichier1);

    if (test_fichier_transf_3d(nomfichier1) != 1)
    {
      PUT_ERROR("This file contains no 3D transformation");
      return ;
    }
  }
 

 im_1=GET_PLACE3D("ROI associee au champ 1");
 
  /*fichier contenant la transformation*/
  {
	sprintf(nomfichier2,"%s",GET_FILE("Second fichier trf",&err));
	
    if (err)
      return ;
	  
	put_file_extension(nomfichier2,".trf",nomfichier2);

    if (test_fichier_transf_3d(nomfichier2) != 1)
    {
      PUT_ERROR("This file contains no 3D transformation");
      return ;
    }
  }
  
 im_2= GET_PLACE3D("ROI associee au champ 2");


 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 put_file_extension(nomfichres,".trf",nomfichres);


imx_fusion_transf_3d(im_1, nomfichier1, im_2, nomfichier2, nomfichres);

}

/**************************************************************************
**      imx_fusion_transf_3d()
**  Fusionne 2 champs de deformation defini sur une certaine ROI
**************************************************************************/
void imx_fusion_transf_3d(int im_1, char* nomfichier1, int im_2, char* nomfichier2, char* nomfichres)
{
grphic3d *im1,*im2;
field3d *ch1, *ch2, *chres;
transf3d *transfo1,*transfo2,*transfores;


im1=ptr_img_3d(im_1);
im2=ptr_img_3d(im_2);

transfo1=load_transf_3d(nomfichier1);
transfo2=load_transf_3d(nomfichier2);

if ((transfo1->width!=transfo2->width)||(transfo1->height!=transfo2->height)||(transfo1->depth!=transfo2->depth))
	{
	printf("Les transformations n'ont pas la meme taille !!! \n");
	return;
	}
	
if ((transfo1->width!=im1->width)||(transfo1->height!=im1->height)||(transfo1->depth!=im1->depth))
	{
	printf("Les images et les champs n'ont pas la meme taille !!! \n");
	return;
	}


if ((transfo2->width!=im2->width)||(transfo2->height!=im2->height)||(transfo2->depth!=im2->depth))
	{
	printf("Les images et les champs n'ont pas la meme taille !!! \n");
	return;
	}

ch1=transf_to_field_3d(transfo1,NULL,NULL);
ch2=transf_to_field_3d(transfo2,NULL,NULL);


chres=cr_field3d(transfo1->width,transfo1->height,transfo1->depth);
    

imx_fusion_transf_3d_p(im1, ch1, im2, ch2, chres);

transfores=field_to_transf_3d(chres,NULL,NULL);
transfores->dx=transfo1->dx;
transfores->dy=transfo1->dy;
transfores->dz=transfo1->dz;



save_transf_3d(transfores,nomfichres);
  

free_field3d(ch1);free_field3d(ch2);free_field3d(chres);
free_transf3d(transfo1);free_transf3d(transfo2);free_transf3d(transfores);
   
}

/**************************************************************************
**      imx_fusion_transf_3d_p()
**  Fusionne 2 champs de deformation defini sur une certaine ROI
**************************************************************************/
void imx_fusion_transf_3d_p(grphic3d *im1, field3d *ch1, grphic3d *im2, field3d *ch2, field3d *chres)
{
grphic3d *distance1,*distance2;
int i,j,k,wdth,hght,dpth;
double d1,d2,dtot;


wdth=im1->width;
hght=im1->height;
dpth=im1->depth;


distance1=cr_grphic3d(im1);imx_copie_param_3d_p(im1,distance1);
distance2=cr_grphic3d(im2);imx_copie_param_3d_p(im2,distance2);


imx_chamfer_distance_3d_p(im1,distance1);
imx_chamfer_distance_3d_p(im2,distance2);

for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
for (k=0; k<dpth; k++)
	if ((distance1->mri[i][j][k]!=0)||(distance2->mri[i][j][k]!=0))
	{
	d1=distance1->mri[i][j][k]*distance1->rcoeff;
	d2=distance2->mri[i][j][k]*distance2->rcoeff;
	dtot=d1+d2;
	d1=d1/dtot;
	d2=d2/dtot;
	
	chres->raw[i][j][k].x=d2*ch1->raw[i][j][k].x+d1*ch2->raw[i][j][k].x;
	chres->raw[i][j][k].y=d2*ch1->raw[i][j][k].y+d1*ch2->raw[i][j][k].y;
	chres->raw[i][j][k].z=d2*ch1->raw[i][j][k].z+d1*ch2->raw[i][j][k].z;
	
	}
	else
	{
	/* le champ1 est priviligie par rapport au champ2 dans les zones d'ambiguite */
	chres->raw[i][j][k].x=ch1->raw[i][j][k].x;
	chres->raw[i][j][k].y=ch1->raw[i][j][k].y;
	chres->raw[i][j][k].z=ch1->raw[i][j][k].z;
	}
	
	
	

}


/**************************************************************************
**      convert_IPB_to_trf()
**  converti 3 fichier ipb en un trf
**************************************************************************/
void convert_IPB_to_trf()
{
 int im_x,im_y,im_z,err;
 char nomfichres[PATH_LEN],str[PATH_LEN];
 
 im_x=GET_PLACE3D("champs suivant x (en voxel)");
 im_y=GET_PLACE3D("champs suivant y (en voxel)");
 im_z=GET_PLACE3D("champs suivant z (en voxel)");
 

 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 put_file_extension(nomfichres,".trf",nomfichres);


imx_convert_IPB_to_trf(im_x, im_y, im_z, nomfichres);

}

/**************************************************************************
**      imx_convert_IPB_to_trf()
**   converti 3 fichier ipb en un trf
**************************************************************************/
void imx_convert_IPB_to_trf(int im_x, int im_y, int im_z, char* nomfichres) 
{
grphic3d *imx,*imy,*imz;
field3d *chres;
transf3d *transfores;


imx=ptr_img_3d(im_x);
imy=ptr_img_3d(im_y);
imz=ptr_img_3d(im_z);


if ((imx->width!=imy->width)||(imx->height!=imy->height)||(imx->depth!=imy->depth)||(imz->width!=imy->width)||(imz->height!=imy->height)||(imz->depth!=imy->depth)||
(imz->width!=imx->width)||(imz->height!=imx->height)||(imz->depth!=imx->depth))
	{
	printf("Les images n'ont pas la meme taille !!! \n");
	return;
	}
	


chres=cr_field3d(imx->width,imx->height,imx->depth);
    

imx_convert_IPB_to_trf_p(imx, imy, imz, chres);



transfores=field_to_transf_3d(chres,NULL,NULL);
transfores->dx=imx->dx;
transfores->dy=imx->dy;
transfores->dz=imx->dz;



save_transf_3d(transfores,nomfichres);
  

free_field3d(chres);
free_transf3d(transfores);
   
}

/**************************************************************************
**      imx_convert_IPB_to_trf_p()
**   converti 3 fichier ipb en un trf
**************************************************************************/
void imx_convert_IPB_to_trf_p(grphic3d *imx, grphic3d *imy, grphic3d *imz, field3d* chres)
{
int wdth,hght,dpth,i,j,k;

wdth=imx->width;
hght=imx->height;
dpth=imx->depth;

for(i=0;i<wdth;i++)
for(j=0;j<hght;j++)
for(k=0;k<dpth;k++)
	{
	chres->raw[i][j][k].x=imx->mri[i][j][k]*imx->rcoeff*imx->dx;
	chres->raw[i][j][k].y=imy->mri[i][j][k]*imy->rcoeff*imx->dy;
	chres->raw[i][j][k].z=imz->mri[i][j][k]*imz->rcoeff*imx->dz;
	}

}

/******************************************************************************
** -- ComputeMoment3D ---------------------------------------------------------------
** 
**   Calcul les moments 3D d'une image binaire ou entre 0et 1
**
******************************************************************************/

void ComputeMoment3D(void)
{ 
  int im_deb,im_res,err;
  int p,q,r;
  float rayon;
  
im_deb=GET_PLACE3D("Image entre 0 et 1");
im_res=GET_PLACE3D("Image resultat");
p= GET_INT("p ?", 0, &err);
q= GET_INT("q ?", 0, &err);
r= GET_INT("r ?", 0, &err);
rayon= GET_FLOAT("Rayon", 2, &err);

imx_ComputeMoment3D( im_deb,  im_res,  p,  q, r, rayon);

show_picture_3d(im_res);
  
 
}

/******************************************************************************
** -- imx_ComputeMoment3D ---------------------------------------------------------------
** 
**   Calcul les moments 3D d'une image binaire ou entre 0et 1
**
******************************************************************************/

void imx_ComputeMoment3D(int im_deb, int im_res, int p, int q, int r, float rayon)
{
  grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_ComputeMoment3D_p(imdeb,  imres,  p,  q, r, rayon);


}


/******************************************************************************
** -- imx_ComputeMoment3D_p ---------------------------------------------------------------
** 
**   Calcul les moments 3D d'une image binaire ou entre 0 et 1
**
******************************************************************************/

void imx_ComputeMoment3D_p(grphic3d * imdeb, grphic3d * imres, int p, int q, int r, float rayon)
{
 int i,j,k,l,m,n,wdth,hght,dpth;
 double ***res;
 double ***masque;
 grphic3d *tmp;
 double pow1, pow2, max1, min1, max2, min2, max3, min3;
 double norm;

wdth=imdeb->width;
hght=imdeb->height;
dpth=imdeb->depth;
 

 tmp=cr_grphic3d(imdeb);
 imx_dilat_3d_p(imdeb,tmp,4,2);//4 pour voisinage sphérique sauf que je ne sais pas le rayon utilisé, 2 = nombre d'itération ca sert à quoi???
 res = alloc_dmatrix_3d(wdth, hght, dpth);

 //préparation du masque de convolution
 masque = alloc_dmatrix_3d(2*rayon+1, 2*rayon+1, 2*rayon+1);
 norm=0;
   for (i=0;i<2*rayon+1;i++)
    {
    pow1 = pow(i-rayon,p);
    for (j=0;j<2*rayon+1;j++)
      {
      pow2 = pow(j-rayon,q);
      for (k=0;k<2*rayon+1;k++)
	{
	      	if(pow(i-rayon,2)+pow(j-rayon,2)+pow(k-rayon,2) > pow(rayon,2))
			 masque[i][j][k] = 0;
		else {   masque[i][j][k] = pow1 * pow2 * pow(k-rayon,r);
			 norm = norm+abs(masque[i][j][k]);
		     }
	}
      } 
    }
 //revoir s'il y a une fonction qui le calcul directement
   for (i=0;i<2*rayon+1;i++)
	for (j=0;j<2*rayon+1;j++)
		for (k=0;k<2*rayon+1;k++)
			masque[i][j][k] = masque[i][j][k]/norm;

 //calcul du GMI
 imx_copie_param_3d_p(imdeb,imres);
 
 for (i=0;i<wdth;i++)
 { 
    max1 = MAXI(0,i-rayon);
    min1 = MINI(wdth,i+rayon+1);
    for (j=0;j<hght;j++)
    {
	max2 = MAXI(0,j-rayon);
        min2 = MINI(hght,j+rayon+1);
      	for (k=0;k<dpth;k++)
      	{
		max3 = MAXI(0,k-rayon);
       		min3 = MINI(dpth,k+rayon+1);
		res[i][j][k] = 0;
        	if (imdeb->mri[i][j][k] != 0) //(tmp->mri[i][j][k] != 0 )
         	{
		for (l=max1;l<min1;l++)// faire l'arrondi de rayon et voir probleme de rayon = 1
		  for (m=max2;m<min2;m++)
      		    for (n=max3;n<min3;n++)
                	res[i][j][k] = res[i][j][k]+ masque[(int)(rayon-(i-l))][(int)(rayon-(j-m))][(int)(rayon-(k-n))]* 
				(imdeb->mri[l][m][n]*imdeb->rcoeff);
		}  

       	}
    }
  }
 
imx_convert_dmatrix3d_grphic3d(res,imres);
imx_inimaxminpixel_3d_p(imres);
free_dmatrix_3d(res);
free_grphic3d(tmp);

}

/******************************************************************************
** -- ComputeRotationInvariantMoment3D ---------------------------------------------------------------
** 
**   Calcul les moments Invariants par rotation 3D d'une image binaire ou entre 0 et 1
**
******************************************************************************/

void ComputeRotationInvariantMoment3D(void)
{ 
  int im_deb,im_res,err;
  int num;
  float rayon;
  
  
im_deb=GET_PLACE3D("Image entre 0 et 1");
im_res=GET_PLACE3D("Image resultat");
num= GET_INT("num ?", 0, &err);
rayon= GET_FLOAT("Rayon", 2, &err);

imx_ComputeRotationInvariantMoment3D( im_deb,  im_res,  num, rayon);

show_picture_3d(im_res);
 
}
/******************************************************************************
** -- imx_ComputeRotationInvariantMoment3D ---------------------------------------------------------------
** 
**   Calcul les moments Invariants par rotation 3D d'une image binaire ou entre 0 et 1
**
******************************************************************************/

void imx_ComputeRotationInvariantMoment3D(int im_deb, int im_res, int num, float rayon)
{
  grphic3d *imdeb,*imres;

  imdeb=ptr_img_3d(im_deb);
  imres=ptr_img_3d(im_res);

  imx_ComputeRotationInvariantMoment3D_p(imdeb, imres, num, rayon);


}

/******************************************************************************
** -- imx_ComputeRotationInvariantMoment3D_p ----------------------------------
** 
**   Calcul les moments invariants par rotation 3D d'une image binaire ou entre 0 et 1
**
******************************************************************************/

void imx_ComputeRotationInvariantMoment3D_p(grphic3d * imdeb, grphic3d * imres, int num, float rayon)
{

 grphic3d *tmp1, *tmp2, *tmp3 ;
////je dois ecrire les equations pour QUE LA PROGRAMMATION SOIT CLAIRE
//travailler sur une copi de imres puis utiliser la fonction imx_copi_3d_p dans imres, ceci évite le problème quand src=dst dans l'affichage

tmp1=cr_grphic3d(imdeb);
tmp2=cr_grphic3d(imdeb);
tmp3=cr_grphic3d(imdeb);





  switch(num){
  case 1: 
	imx_ComputeMoment3D_p(imdeb, imres, 0, 0, 0, rayon);break;
  case 2: 
	imx_ComputeMoment3D_p(imdeb, tmp1, 2, 0, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 0, 2, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp3, 0, 0, 2, rayon);
	imx_add_3d_p(tmp1,tmp2,imres);
	imx_add_3d_p(tmp3,imres,imres);break;
  case 3: 
	imx_ComputeMoment3D_p(imdeb, tmp1, 2, 0, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 0, 2, 0, rayon);
	imx_mul_3d_p(tmp1,tmp2,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 2, 0, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 0, 0, 2, rayon);
	imx_mul_3d_p(tmp1,tmp2,tmp3);

	imx_add_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 0, 2, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 0, 0, 2, rayon);
	imx_mul_3d_p(tmp1,tmp2,tmp3);

	imx_add_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 1, 0, 1, rayon);
	imx_mul_3d_p(tmp1,tmp1,tmp3);

	imx_subabs_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 1, 1, 0, rayon);
	imx_mul_3d_p(tmp1,tmp1,tmp3);

	imx_subabs_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 0, 1, 1, rayon);
	imx_mul_3d_p(tmp1,tmp1,tmp3);

	imx_subabs_3d_p(imres,tmp3,imres);break;
  case 4: 
	imx_ComputeMoment3D_p(imdeb, tmp1, 2, 0, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 0, 2, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp3, 0, 0, 2, rayon);
	imx_mul_3d_p(tmp1,tmp2,imres);
	imx_mul_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 0, 0, 2, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 1, 1, 0, rayon);
	imx_mul_3d_p(tmp1,tmp2,tmp3);
	imx_mul_3d_p(tmp3,tmp2,tmp3);

	imx_subabs_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 0, 2, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 1, 0, 1, rayon);
	imx_mul_3d_p(tmp1,tmp2,tmp3);
	imx_mul_3d_p(tmp3,tmp2,tmp3);

	imx_subabs_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 2, 0, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 0, 1, 1, rayon);
	imx_mul_3d_p(tmp1,tmp2,tmp3);
	imx_mul_3d_p(tmp3,tmp2,tmp3);

	imx_subabs_3d_p(imres,tmp3,imres);

	imx_ComputeMoment3D_p(imdeb, tmp1, 1, 1, 0, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp2, 1, 0, 1, rayon);
	imx_ComputeMoment3D_p(imdeb, tmp3, 0, 1, 1, rayon);
	imx_mul_3d_p(tmp2,tmp3,tmp3);
	imx_mul_3d_p(tmp3,tmp1,tmp3);
	imx_mul_coe_3d_p(tmp3,2,tmp3);

	imx_add_3d_p(imres,tmp3,imres);	break;
 default: printf("ce n'est pas un choix valable\n");break;
  }




free_grphic3d(tmp1);
free_grphic3d(tmp2);
free_grphic3d(tmp3);



	
}

 
void stat_fich_histo(void)
{ 

char nomfichier[500];
int err;
FILE * fic;
float ng, pb;
float prob5,prob95,prob50;
float ng5,ng95,ng50;
bool trv1,trv2,trv3;
 
sprintf(nomfichier,"%s",GET_FILE("*.*",&err));
if (err)
   return ;
   


prob5  = 0.0;  
prob95 = 0.0;
prob50 = 0.0;
trv1 = FALSE;
trv2 = FALSE;
trv3 = FALSE;

fic = fopen(nomfichier, "r");

while( fscanf(fic,"%f %f",&ng,&pb) !=EOF )
{
	//printf("niveau %f : proba : %f \n",ng,pb);
	prob5 = prob5 + pb;
	prob95 = prob95 + pb;
	prob50 = prob50 + pb;
	if (prob5>=0.05 && trv1 == FALSE)   
	{
		trv1 = TRUE;
		ng5 = ng;
	}

	if (prob95>=0.95 && trv2 == FALSE)   
	{
		trv2 = TRUE;
		ng95 = ng;
	}

	if (prob50>=0.5 && trv3== FALSE)   
	{
		trv3 = TRUE;
		ng50 = ng;
	}
	
}
 
fclose(fic);


  PUTF("mediane= =_",ng50,"_");
  PUTF("quantile d'ordre 20 (5%)=_",ng5,"_");
  PUTF(",",ng95,"");



 
}

/*******************************************************************************
**        Compute_energie_groupwise_variance                                                
**                                                                        
*******************************************************************************/

void Compute_energie_groupwise_variance()
{
#ifndef COMPILE_FOR_MEDIPY

char * file,s[250];
char * nb_img_ini;
int e, ans, l, wdth, hght, dpth,max=0,step=0, i,j,k;
int nb_img=0, im_mask, nb=0;
double var=0.0, moy=0.0, Etot=0.0,tmp;
grphic3d** serie_groupwise;
grphic3d* mask;


file=strdup(GET_FILE("Serie d'images a recaler (*.ipb)", &e)); 
	if (e)
	{
		free(file);
		return;
	}    

nb_img_ini = getheader_interfile(file,"number of images=", ENTIER,0);	
nb_img=atoi(nb_img_ini);


 sprintf(s,"Voulez vous consider un mask \n");
 ans=GET_EXIST(s,&e);
 
 /*question sur les images a recaler*/
 if (ans==1)
 {
 im_mask=GET_PLACE3D("Mask");
 mask=ptr_img_3d(im_mask);
 }
 else 
 {
 mask=NULL;
 }
 
serie_groupwise=CALLOC(nb_img,ptr_grphic3d);

wdth = getmri_width(file, 1, MODE_3D);
hght = getmri_height(file, 1, MODE_3D);
dpth = getmri_depth(file, 1, MODE_3D);

 
for (l = 0; l<nb_img;l++) {
    serie_groupwise[l] = cr_grphic3d_modif(wdth, hght, dpth, 1., 1., 1);		
    load_mri_3d(file,  serie_groupwise[l], l+1, max, step);
  }


// Calcul de la variance moyenne sur le mask
if (mask)
{
for (i=0; i< wdth; i++)
for (j=0; j< hght; j++)
for (k=0; k< dpth; k++)
	if (mask->mri[i][j][k]>0)
		{
		moy=0.0;
		for (l = 0; l<nb_img;l++) 
			moy+=serie_groupwise[l]->mri[i][j][k]*serie_groupwise[l]->rcoeff;

		moy=moy/nb_img;
		
		var=0.0;
		for (l = 0; l<nb_img;l++) 
			{
			tmp=(serie_groupwise[l]->mri[i][j][k]*serie_groupwise[l]->rcoeff-moy);
			var+=tmp*tmp;
			}
		
		var=var/nb_img;
		
		Etot+=var;
		nb++;
		}
}
else
{
for (i=0; i< wdth; i++)
for (j=0; j< hght; j++)
for (k=0; k< dpth; k++)
		{
		moy=0.0;
		for (l = 0; l<nb_img;l++) 
			moy+=serie_groupwise[l]->mri[i][j][k]*serie_groupwise[l]->rcoeff;

		moy=moy/nb_img;
		
		var=0.0;
		for (l = 0; l<nb_img;l++) 
			{
			tmp=(serie_groupwise[l]->mri[i][j][k]*serie_groupwise[l]->rcoeff-moy);
			var+=tmp*tmp;
			}
		
		var=var/nb_img;
		
		Etot+=var;
		nb++;
		}
}

Etot=Etot/nb;

printf("variance moyenne : %f \n", Etot);


free(file);

// Liberation memoire des elements de la structure ParamRecalageBspline
for (l = 0; l<nb_img;l++) 
	{
	free_grphic3d(serie_groupwise[l]);
	}

free(serie_groupwise);
#endif

}

/*******************************************************************************
**        NLMWeightOnePoint3D                                                
**                                                                        
*******************************************************************************/

void NLMWeightOnePoint3D(grphic3d* Image, double *** weights, float NLMsmooth, int NLMhwnx, int NLMhwny, int NLMhwnz, int NLMhwvsx, int NLMhwvsy, int NLMhwvsz,  int x, int y, int z)
{
double sum = 0;
int xx,yy,zz;
int sx,sy,sz;
int px,py,pz;
double constante= Image->rcoeff * Image->rcoeff /NLMsmooth;
int minsx=MAXI(-x,-NLMhwvsx);
int minsy=MAXI(-y,-NLMhwvsy);
int minsz=MAXI(-z,-NLMhwvsz);

int maxsx=MINI(Image->width-x-1,NLMhwvsx);
int maxsy=MINI(Image->height-y-1,NLMhwvsy);
int maxsz=MINI(Image->depth-z-1,NLMhwvsz);

  //initialisation de weights
  for(xx=0; xx!=2*NLMhwvsx+1; xx++)
    for(yy=0; yy!=2*NLMhwvsy+1; yy++)
      for(zz=0; zz!=2*NLMhwvsz+1; zz++)
        weights[xx][yy][zz] = 0;



  //Voxel selection in search volume
  for( sx=minsx;sx<=maxsx;sx++){		    
    xx = x+sx;
      for( sy=minsy;sy<=maxsy;sy++){
        yy = y+sy;
          for( sz=minsz;sz<=maxsz;sz++){
          double diff=0;	    
  	      double weight = 0;
	      double dist = 0;
			    
		  int minpx=MAXI(-x,-NLMhwnx);
          int minpy=MAXI(-y,-NLMhwny);
          int minpz=MAXI(-z,-NLMhwnz);
          int maxpx=MINI(Image->width-x-1,NLMhwnx);
          int maxpy=MINI(Image->height-y-1,NLMhwny);
          int maxpz=MINI(Image->depth-z-1,NLMhwnz);
          int xpx,ypy,zpz,xxpx,yypy,zzpz;		    
          
		  zz = z+sz;
		      
 	      minpx=MAXI(minpx,-xx);
		  
		  minpy=MAXI(minpy,-yy);
		  
		  minpz=MAXI(minpz,-zz);
		  
		  maxpx=MINI(maxpx,Image->width-xx-1);
		  
		  maxpy=MINI(maxpy,Image->height-yy-1);
		
		  maxpz=MINI(maxpz,Image->depth-zz-1);
		
				
	          //distance computation between patches
              for( px=minpx;px<=maxpx;px++){
	        xpx = x + px;
	        xxpx= xx+ px;
	          for( py=minpy;py<=maxpy;py++){
		    ypy = y + py;
		    yypy= yy+ py;
  	                for( pz=minpz;pz<=maxpz;pz++){
			  zpz = z + pz;
			  zzpz= zz+ pz;
  		            diff = Image->mri[xpx][ypy][zpz] - Image->mri[xxpx][yypy][zzpz];
  			    dist += diff*diff;
		          }
			}
	              }

  	      weight = exp(-dist *constante);
              weights[sx+NLMhwvsx][sy+NLMhwvsy][sz+NLMhwvsx] = weight;
              sum += weight;                           


			  
			
	    }
	  }
	}

  if(sum>0)
    for(xx=0; xx!=2*NLMhwvsx+1; xx++)
      for(yy=0; yy!=2*NLMhwvsy+1; yy++)
        for(zz=0; zz!=2*NLMhwvsz+1; zz++)
          weights[xx][yy][zz] /= sum;

}

/*******************************************************************************
**        NLMSmoothComputation3D                                                
**                                                                        
*******************************************************************************/

float NLMSmoothComputation3D(grphic3d* Image, int NLMhwnx, int NLMhwny, int NLMhwnz, float NLMbeta, int padding)
{
  uint count = 0;
  double sigma2 = 0;
  int x,y,z;
  float NLMsmooth;
  for(x=1;x<Image->width-1;x++)
    for(y=1;y<Image->height-1;y++)
      for(z=1;z<Image->depth-1;z++)
        if(Image->mri[x][y][z] * Image->rcoeff >padding)
        {
	float ei = sqrt(6/7.0)*Image->rcoeff*(Image->mri[x][y][z] -(Image->mri[x+1][y][z]+Image->mri[x-1][y][z]+Image->mri[x][y+1][z]+Image->mri[x][y-1][z]
				+Image->mri[x][y][z+1]+Image->mri[x][y][z-1])/6.0);
	sigma2 += ei*ei;
	count ++;
        }

  sigma2 = sigma2 / count;
  NLMsmooth = 2 * NLMbeta * sigma2 * (2*NLMhwnx+1) * (2*NLMhwny+1) * (2*NLMhwnz+1);
  return NLMsmooth;

}
