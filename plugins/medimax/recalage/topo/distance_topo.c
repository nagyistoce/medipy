/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>

#include "noyau/imagix.h"

#include "recalage/mtch_3d.h"
#include "recalage/chps_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/validation_topo_3d.h"
#include "recalage/topo/distance_topo.h"
#include "recalage/topo/normalisation_topo.h"
#include "recalage/topo/divers_topo.h"
#include "math/imx_matrix.h"



/*******************************************************************************
**  Energie_globale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_t dist)
**
**	Calcule l'energie globale selon le critere dist entre deux image
*******************************************************************************/

double Energie_globale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_locale_t dist,reg_func_locale_t regularisation,int ***masque_param)
{
grphic3d *mask=NULL;

if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d))
	mask=imref;

if (dist==Energie_IM_locale_3d)
{
double E;

E=Calcul_IM_global_3d(imref, imreca,  nb_param, param);

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*regularisation_globale_3d(nb_param, param,param_norm, masque_param,mask,regularisation);
}

return(E);
}
else
{
double E,tmp;
int topi,topj,topk,topD,resol,l,i;
int *x00,*x11,*y00,*y11,*z00,*z11;
int width,height,depth;
double xm,xM,ym,yM,zm,zM;
TSlpqr Slpqr[8];
   
width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
resol=BASE3D.resol; 
x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
topD=(int)pow(2.0,1.0*resol)-1;

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
	{
	Slpqr[i].Jm=-HUGE_VAL;
	Slpqr[i].JM=HUGE_VAL;
	} 
  
E=0.0;
	for (topi=0; topi<topD; topi=topi+2)
  	for (topj=0; topj<topD; topj=topj+2)
  		for (topk=0; topk<topD; topk=topk+2)
  			if (masque_param[topi][topj][topk]!=0)
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

				
				tmp=dist(imref,imreca,nb_param,param,param_norm,topi,topj,topk,Slpqr,regularisation);
									
				E=E+tmp*tmp;
				
					
			  }
E=sqrt(E);
//E=E/(width*height*depth);
return(E); 
}

}

/*******************************************************************************
**  Energie_quad_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_quad_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
	
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
        
			auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
			E=E+auxE*auxE;
				
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
       
			 auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
				E=E+auxE*auxE;  
				   
       }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL;
       }
      }
   }
   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


//E=sqrt(E)*imreca->rcoeff;
E=sqrt(E);

free_field3d(champ);

return(E);

}

/*******************************************************************************
**  Energie_quad_sym_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie symetrique locale quadratique relative a une boite
*******************************************************************************/
  
double Energie_quad_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double Ereg;

  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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
    dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
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
        
		auxE=auxE-imref->mri[i][j][k];
		if (auxE!=0)
			{xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
			fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
			E=E+(1.0+fabs(J))*auxE*auxE;
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
        auxE=auxE-imref->mri[i][j][k];
	
		if (auxE!=0)
			{xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
			fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
			E=E+(1.0+fabs(J))*auxE*auxE;
			}     
       }
       else
       {
        /*auxE=0-imref->mri[i][j][k];
				if (auxE!=0)
					{xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
					fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
					E=E+(1.0+fabs(J))*auxE*auxE;
					}*/
					E=HUGE_VAL;
       }
      }
   }
}   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);
E=E/2.0;
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=sqrt(E);
free_field3d(champ);
return(E);
} 
/*******************************************************************************
**  Energie_Lp_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_Lp_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
double p=ORDRE_NORME_LP;
field3d *champ;
vector3d ***data;
double Ereg;

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
        
		auxE=fabs(auxE-imref->mri[i+bx0][j+by0][k+bz0]);
		auxE=1.0*pow(auxE,p);
		E=E+auxE;
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
        auxE=fabs(auxE-imref->mri[i+bx0][j+by0][k+bz0]);
				auxE=1.0*pow(auxE,p);
				E=E+auxE;     
       }
       else
       {
	      /*auxE=fabs(imref->mri[i+bx0][j+by0][k+bz0]);
				auxE=1.0*pow(auxE,p);
				E=E+auxE;*/
				E=HUGE_VAL; 
       }
      }
   }

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}
   
E=sqrt(E);
free_field3d(champ);
return(E);
}

/*******************************************************************************
**  Energie_Lp_sym_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_Lp_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double p=ORDRE_NORME_LP;
double Ereg;
 
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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
    dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
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
        
		auxE=fabs(auxE-imref->mri[i][j][k]);
		auxE=1.0*pow(auxE,p);
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
		E=E+(1.0+fabs(J))*auxE;
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
	    auxE=fabs(auxE-imref->mri[i][j][k]);
		auxE=1.0*pow(auxE,p);
		
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
		E=E+(1.0+fabs(J))*auxE;     
       }
       else
       {
        /*auxE=fabs(0-imref->mri[i][j][k]);
				auxE=1.0*pow(auxE,p);
				xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
				fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
				E=E+(1.0+fabs(J))*auxE; */
				E=HUGE_VAL;
       }
      }
   }
}   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);
E=E/2.0;
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);
free_field3d(champ);
return(E);
} 


/*******************************************************************************
**  Energie_L1L2_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_L1L2_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
double Ereg;

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
        
		auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		auxE=auxE*auxE+TOPO_EPSILON;
		E=E+sqrt(auxE);
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
        auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
				auxE=auxE*auxE+TOPO_EPSILON;
				E=E+sqrt(auxE);     
       }
       else
       {
	      /*auxE=imref->mri[i+bx0][j+by0][k+bz0];
				auxE=auxE*auxE+TOPO_EPSILON;
				E=E+sqrt(auxE); */
				E=HUGE_VAL;
       }
      }
   }

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}
   
free_field3d(champ);
return(sqrt(E));
}

/*******************************************************************************
**  Energie_L1L2_sym_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_L1L2_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double Ereg;
 
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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
    dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
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

		auxE=auxE-imref->mri[i][j][k];
		auxE=auxE*auxE+TOPO_EPSILON;
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
		E=E+(1.0+fabs(J))*sqrt(auxE);
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

    auxE=auxE-imref->mri[i][j][k];
		auxE=auxE*auxE+TOPO_EPSILON;
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
		E=E+(1.0+fabs(J))*sqrt(auxE);
       }
       else
       {
    		/*auxE=imref->mri[i][j][k];
				auxE=auxE*auxE+TOPO_EPSILON;
				xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
				fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
				E=E+(1.0+fabs(J))*sqrt(auxE);*/
				E=HUGE_VAL;
       }
      }
   }
}   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);
E=E/2.0;
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);
free_field3d(champ);
return(E);
} 

/*******************************************************************************
**  Energie_L1norm_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_L1norm_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
double Ereg;

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

    if ((imref->mri[i+bx0][j+by0][k+bz0]+auxE)!=0)   
			E=E+1.0*fabs(auxE-imref->mri[i+bx0][j+by0][k+bz0])/(auxE+imref->mri[i+bx0][j+by0][k+bz0]);
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
    
		if ((imref->mri[i+bx0][j+by0][k+bz0]+auxE)!=0)
			E=E+1.0*fabs(auxE-imref->mri[i+bx0][j+by0][k+bz0])/(auxE+imref->mri[i+bx0][j+by0][k+bz0]);
       }
       else
       {
	   /*if ((imref->mri[i+bx0][j+by0][k+bz0]!=0))   
				E=E+1.0;*/
				E=HUGE_VAL;
       }
      }
   }
   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);
E=E/2.0;
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

free_field3d(champ);
return(E);
}

/*******************************************************************************
**  Energie_quad_topo_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique avec penalisation du jaconien relative a une boite
*******************************************************************************/

double Energie_quad_topo_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin,stop=1;
double Jmin,Jmax;
double Ereg;

Jmin=Slpqr[0].Jm;
Jmax=Slpqr[0].JM;

  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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
	for (i=ideb;((i<ifin)&&(stop==1));i++)
	for (j=jdeb;((j<jfin)&&(stop==1));j++)
	for (k=kdeb;((k<kfin)&&(stop==1));k++)
	   {
    dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
    i2=(int)dx;j2=(int)dy;k2=(int)dz;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
		
		if ((J<=Jmin)||(J>=Jmax))
			stop=0;
			
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
        
		auxE=auxE-imref->mri[i][j][k];
		auxE=auxE*auxE;
		E=E+auxE;
	  
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

    auxE=(auxE-imref->mri[i][j][k]);
		auxE=auxE*auxE;
		
			E=E+auxE;
		     
       }
       else
       {
        /*auxE=(0-imref->mri[i][j][k]);
				auxE=auxE*auxE;
				E=E+auxE;
				*/
				E=HUGE_VAL;
		   }
      }
   }
} 
if (stop==1)  
{
if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=sqrt(E);
}
else 
E=HUGE_VAL;



free_field3d(champ);
return(E);
} 
/*******************************************************************************
**  Energie_quad_sym_topo_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique avec penalisation du jaconien relative a une boite
*******************************************************************************/

double Energie_quad_sym_topo_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin,stop=1;
double Jmin,Jmax;
double Ereg;
 
Jmin=Slpqr[0].Jm;
Jmax=Slpqr[0].JM;

  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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
	for (i=ideb;((i<ifin)&&(stop==1));i++)
	for (j=jdeb;((j<jfin)&&(stop==1));j++)
	for (k=kdeb;((k<kfin)&&(stop==1));k++)
	   {
    dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
    i2=(int)dx;j2=(int)dy;k2=(int)dz;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
		xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
		fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
		
		if ((J<=Jmin)||(J>=Jmax))
			stop=0;
			
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
        
		auxE=auxE-imref->mri[i][j][k];
		auxE=auxE*auxE;
		E=E+(1.0+fabs(J))*auxE;
	  
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

	    auxE=(auxE-imref->mri[i][j][k]);
		auxE=auxE*auxE;
		
			E=E+(1.0+fabs(J))*auxE;
		     
       }
       else
       {
        /*auxE=(0-imref->mri[i][j][k]);
				auxE=auxE*auxE;
		
				E=E+(1.0+fabs(J))*auxE;*/
				E=HUGE_VAL;
		   }
      }
   }
} 
if (stop==1)  
{
E=E/2.0;

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=sqrt(E);
}

else 
E=HUGE_VAL;




free_field3d(champ);
return(E);
} 

/*******************************************************************************
**  Energie_geman_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale robuste relative a une boite
*******************************************************************************/

double Energie_geman_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE,u2;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
double Ereg;

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
        
		auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		u2=auxE*auxE;
		E=E+u2/(TOPO_ALPHA_ROBUST+u2);
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

	        auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		u2=auxE*auxE;
		E=E+u2/(TOPO_ALPHA_ROBUST+u2);
	     }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
		
				u2=auxE*auxE;
				E=E+u2/(TOPO_ALPHA_ROBUST+u2);*/
				E=HUGE_VAL;
		   }
      }
   }

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}
   
E=sqrt(E);
free_field3d(champ);

if (isnan(E))
	printf("Attention y a un bug ! TOPO_ALPHA_ROBUST = %f \n",TOPO_ALPHA_ROBUST);

return(E);
}


/*******************************************************************************
**  Energie_geman_sym_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale robuste relative a une boite
*******************************************************************************/

double Energie_geman_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{


int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0;
double dx,dy,dz,auxE,u2;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double Ereg;
 
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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
    dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
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
	 
		auxE=auxE-imref->mri[i][j][k];
		u2=auxE*auxE;
		if (auxE!=0)
			{xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
			fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
			E=E+(1.0+fabs(J))*u2/(TOPO_ALPHA_ROBUST+u2);
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
	
	auxE=auxE-imref->mri[i][j][k];
				u2=auxE*auxE;
			
		if (auxE!=0)
			{xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
			fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
			E=E+(1.0+fabs(J))*u2/(TOPO_ALPHA_ROBUST+u2);
			}     
       }
       else
       {
        /*auxE=0-imref->mri[i][j][k];
				if (auxE!=0)
					{xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
					fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
					E=E+(1.0+fabs(J))*auxE*auxE;
					}*/
					E=HUGE_VAL;
       }
      }
   }
}   

E=E/2.0;
if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}
E=sqrt(E);
free_field3d(champ);

if (isnan(E))
	printf("Attention y a un bug ! TOPO_ALPHA_ROBUST = %f \n",TOPO_ALPHA_ROBUST);

return(E);
}


/*******************************************************************************
**  Energie_globale_sym_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_t dist)
**
**	Calcule l'energie globale symetrique (en utilisant la transfo inverse) selon le critere L2 entre deux image
*******************************************************************************/

double Energie_globale_sym_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_locale_t dist,int ***masque_param,reg_func_locale_t regularisation)
{
/*grphic3d *imtreca,*imtref;

field3d *champ,*champ_inv;
int wdth,hght,dpth,i,j,k;
double E;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

imtreca=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imtreca);
imtref=cr_grphic3d(imref);imx_copie_param_3d_p(imref,imtref);

champ=cr_field3d(wdth,hght,dpth);
champ_inv=cr_field3d(wdth,hght,dpth);
	 
base_to_field_3d(nb_param,param,champ,NULL,NULL);
inter_qsinc3_3d(imreca,champ,imtref);

//inv_field_3d(champ,champ_inv);
inter_qsinc3_3d(imref,champ_inv,imtreca);

E=0;

	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)	 
		{
	
E=E+0.5*(1.0*(imreca->mri[i][j][k]-imtreca->mri[i][j][k])*(imreca->mri[i][j][k]-imtreca->mri[i][j][k])
				+1.0*(imref->mri[i][j][k]-imtref->mri[i][j][k])*(imref->mri[i][j][k]-imtref->mri[i][j][k]));		
		}



free_field3d(champ);
free_field3d(champ_inv);
free_grphic3d(imtreca);
free_grphic3d(imtref);

E=sqrt(E);
return(E); */
return(-1);
}

/*******************************************************************************
**  			update_TOPO_ALPHA_ROBUST_global
**
**	Calcul global du parametre d'échelle dans la fonction de Geman Mcclure
*******************************************************************************/

void update_TOPO_ALPHA_ROBUST_global(grphic3d *imref, grphic3d *imreca, int nb_param, double *param)
{
grphic3d *imres;

field3d *champ;
int wdth,hght,dpth,i,j,k,Ntot,l;
double *residu;

wdth=imreca->width;
hght=imreca->height;
dpth=imreca->depth;

imres=cr_grphic3d(imref);imx_copie_param_3d_p(imref,imres);

champ=cr_field3d(wdth,hght,dpth);
	 
base_to_field_3d(nb_param,param,champ,NULL,NULL);
inter_qsinc3_3d(imreca,champ,imres);

Ntot=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if ((imreca->mri[i][j][k]!=0)&&(imres->mri[i][j][k]!=0))
		  Ntot++;
			
residu=(double *) malloc (Ntot*sizeof(double));

l=0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if ((imreca->mri[i][j][k]!=0)&&(imres->mri[i][j][k]!=0))
			{residu[l]=fabs(imref->mri[i][j][k]-imres->mri[i][j][k]); l++;}

qsort(residu,Ntot,sizeof(double),double_compare_function2);

TOPO_ALPHA_ROBUST=1.4826*residu[(int)floor(Ntot/2)];
TOPO_ALPHA_ROBUST=TOPO_ALPHA_ROBUST*TOPO_ALPHA_ROBUST;

printf("Parametre TOPO_ALPHA_ROBUST = %f \n",TOPO_ALPHA_ROBUST);

free(residu);
free_field3d(champ);
free_grphic3d(imres);

}


/*******************************************************************************
**  			update_TOPO_ALPHA_ROBUST_local
**
**	Calcul local du parametre d'échelle dans la fonction de Geman Mcclure
*******************************************************************************/

void update_TOPO_ALPHA_ROBUST_local(grphic3d *imref, grphic3d *imreca, int nb_param, double *param,int topi, int topj, int topk)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz,ind;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
double *residu,*residu_tmp;
field3d *champ;
vector3d ***data;


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

residu_tmp=(double *) malloc (topDx*topDy*topDz*sizeof(double));

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

ind=0;
// On balaye chacun des pixels de la boite
for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
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
        
		auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		
		if ((imref->mri[i+bx0][j+by0][k+bz0]!=0)&&(auxE!=0))
		residu_tmp[ind]=fabs(auxE);
		else
		residu_tmp[ind]=0;
		
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

       auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		
		if ((imref->mri[i+bx0][j+by0][k+bz0]!=0)&&(auxE!=0))
		residu_tmp[ind]=fabs(auxE);
		else
		residu_tmp[ind]=0;
		
	     }
       else
       {
	      auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
				
				if ((imref->mri[i+bx0][j+by0][k+bz0]!=0))
				residu_tmp[ind]=fabs(auxE);
				else
				residu_tmp[ind]=0;
		
		   }
      }
   ind ++;
	 }
	 
	 
l=0;
for (i=0; i<ind;i++)
	if (residu_tmp[i]>0)
	 	l++;
		
residu=(double *) malloc (l*sizeof(double));

l=0;
for (i=0; i<ind;i++)
	if (residu_tmp[i]>0)
	 	{residu[l]=residu_tmp[i];l++;}



qsort(residu,l,sizeof(double),double_compare_function2);

TOPO_ALPHA_ROBUST=1.4826*residu[(int)floor(l/2)];
TOPO_ALPHA_ROBUST=TOPO_ALPHA_ROBUST*TOPO_ALPHA_ROBUST;

if (l==0)
	TOPO_ALPHA_ROBUST=1;
 	
free(residu);free(residu_tmp);
free_field3d(champ);

}


/*******************************************************************************
**  			regularisation_energie_membrane_local
**
**	Calcul local du terme de regularisation representant l'energie de membrane
*******************************************************************************/

double regularisation_energie_membrane_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
int uxdeb,uydeb,uzdeb;
double Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,reg;
field3d *champ;
vector3d ***data;
int ideb,jdeb,kdeb,condition;


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

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;


reg=0;
	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	   {
    
		if (mask!=NULL)
			condition=mask->mri[ux+bx0][uy+by0][uz+bz0];
		else
			condition=1;
			
		
		if (condition>0)
		{
			if (ux<ideb)
			uxdeb=0;
			else
			uxdeb=ideb;
			
			if (uy<jdeb)
			uydeb=0;
			else
			uydeb=jdeb;
		
			if (uz<kdeb)
			uzdeb=0;
			else
			uzdeb=kdeb;
				
		 
		 	Jxx=		1.0*(data[uxdeb+1][uy][uz].x-data[uxdeb][uy][uz].x);
		 	//Jxx=	1.0+	1.0*(data[uxdeb+1][uy][uz].x-data[uxdeb][uy][uz].x);
		 	Jxy=		1.0*(data[ux][uydeb+1][uz].x-data[ux][uydeb][uz].x);
		 	Jxz=		1.0*(data[ux][uy][uzdeb+1].x-data[ux][uy][uzdeb].x);
		 	Jyx=		1.0*(data[uxdeb+1][uy][uz].y-data[uxdeb][uy][uz].y);
		 	Jyy=		1.0*(data[ux][uydeb+1][uz].y-data[ux][uydeb][uz].y);
		 	//Jyy=	1.0+	1.0*(data[ux][uydeb+1][uz].y-data[ux][uydeb][uz].y);
		 	Jyz=		1.0*(data[ux][uy][uzdeb+1].y-data[ux][uy][uzdeb].y);
		 	Jzx=		1.0*(data[uxdeb+1][uy][uz].z-data[uxdeb][uy][uz].z);
		 	Jzy=		1.0*(data[ux][uydeb+1][uz].z-data[ux][uydeb][uz].z);
		 	Jzz=		1.0*(data[ux][uy][uzdeb+1].z-data[ux][uy][uzdeb].z);
		 	//Jzz=	1.0+	1.0*(data[ux][uy][uzdeb+1].z-data[ux][uy][uzdeb].z);
		 
		 	reg+=Jxx*Jxx+Jxy*Jxy+Jxz*Jxz+Jyx*Jyx+Jyy*Jyy+Jyz*Jyz+Jzx*Jzx+Jzy*Jzy+Jzz*Jzz;
		 
		 	//J=Jxx*(Jyy*Jzz-Jzy*Jyz)-Jyx*(Jxy*Jzz-Jxz*Jzy)+Jzx*(Jxy*Jyz-Jyy*Jxz);
			} 
		 }

//reg=reg-3*topDx*topDy*topDz;


free_field3d(champ);

reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique
return(reg);
}

/*******************************************************************************
**  			regularisation_energie_membrane_Lp_local
**
**	Calcul local du terme de regularisation representant l'energie de membrane (avec norme Lp)
*******************************************************************************/

double regularisation_energie_membrane_Lp_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
int uxdeb,uydeb,uzdeb;
double Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz,reg,tmp;
field3d *champ;
vector3d ***data;
int ideb,jdeb,kdeb,condition;


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

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;


reg=0;
	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	   {
    
		
		if (mask!=NULL)
			condition=mask->mri[ux+bx0][uy+by0][uz+bz0];
		else
			condition=1;
			
		
		if (condition>0)
		{
			if (ux<ideb)
				uxdeb=0;
			else
				uxdeb=ideb;
			
			if (uy<jdeb)
				uydeb=0;
			else
				uydeb=jdeb;
		
			if (uz<kdeb)
				uzdeb=0;
			else
				uzdeb=kdeb;
				
		 
			 Jxx=		1.0*(data[uxdeb+1][uy][uz].x-data[uxdeb][uy][uz].x);
			 Jxy=		1.0*(data[ux][uydeb+1][uz].x-data[ux][uydeb][uz].x);
			 Jxz=		1.0*(data[ux][uy][uzdeb+1].x-data[ux][uy][uzdeb].x);
			 Jyx=		1.0*(data[uxdeb+1][uy][uz].y-data[uxdeb][uy][uz].y);
			 Jyy=		1.0*(data[ux][uydeb+1][uz].y-data[ux][uydeb][uz].y);
			 Jyz=		1.0*(data[ux][uy][uzdeb+1].y-data[ux][uy][uzdeb].y);
			 Jzx=		1.0*(data[uxdeb+1][uy][uz].z-data[uxdeb][uy][uz].z);
			 Jzy=		1.0*(data[ux][uydeb+1][uz].z-data[ux][uydeb][uz].z);
			 Jzz=		1.0*(data[ux][uy][uzdeb+1].z-data[ux][uy][uzdeb].z);
		 
			 tmp=pow(Jxx,ORDRE_REGULARISATION_LP);
			 if (isnan(tmp))
				 	tmp=0;
				reg+=tmp;	 
		 
		 
			 tmp=pow(Jxy,ORDRE_REGULARISATION_LP);
			 if (isnan(tmp))
				 	tmp=0;
				reg+=tmp;	 
		 
		 	tmp=pow(Jxz,ORDRE_REGULARISATION_LP);
		 	if (isnan(tmp))
				 	tmp=0;
				reg+=tmp;	 
		 
		 	tmp=pow(Jyx,ORDRE_REGULARISATION_LP);
		 	if (isnan(tmp))
				 	tmp=0;
				reg+=tmp;	 
		 
		 	tmp=pow(Jyy,ORDRE_REGULARISATION_LP);
		 	if (isnan(tmp))
			 	tmp=0;
			reg+=tmp;	 
		 
		 tmp=pow(Jyz,ORDRE_REGULARISATION_LP);
		 if (isnan(tmp))
			 	tmp=0;
			reg+=tmp;	 
		 
		 tmp=pow(Jzx,ORDRE_REGULARISATION_LP);
		 if (isnan(tmp))
			 	tmp=0;
			reg+=tmp;	 
		 
		 tmp=pow(Jzy,ORDRE_REGULARISATION_LP);
		 if (isnan(tmp))
			 	tmp=0;
			reg+=tmp;	 
		 
		 tmp=pow(Jzz,ORDRE_REGULARISATION_LP);
		 if (isnan(tmp))
			 	tmp=0;
			reg+=tmp;	 
		 
		} 
		 }



free_field3d(champ);
reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,ORDRE_REGULARISATION_LP)); // Pour assurer la compatibilite avec le recalage symetrique

return(reg);
}


/*******************************************************************************
**  			regularisation_energie_membrane_jacobien_local
**
**	Calcul local du terme de regularisation representant l'energie de membrane
*******************************************************************************/

double regularisation_energie_membrane_jacobien_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,*dfx,*dfy,*dfz, dfdx, dfdy,dfdz,d2fdxdy, d2fdxdz, d2fdydz, px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
double Jacobien,Jx,Jy,Jz,reg;
field3d *champ,*champdx, *champdy, *champdz,*champdxy,*champdxz,*champdyz;
mat3d mat;

vector3d ***data,***datadx,***datady,***datadz,***datadxy,***datadxz,***datadyz;
int ideb,jdeb,kdeb,condition;


resol=BASE3D.resol;
width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
topD=(int)pow(2.0,1.0*resol)-1;


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;

champ=cr_field3d(topDx,topDy,topDz);
champdx=cr_field3d(topDx,topDy,topDz);
champdy=cr_field3d(topDx,topDy,topDz);
champdz=cr_field3d(topDx,topDy,topDz);


champdxy=cr_field3d(topDx,topDy,topDz);
champdxz=cr_field3d(topDx,topDy,topDz);
champdyz=cr_field3d(topDx,topDy,topDz);


data=champ->raw;
datadx=champdx->raw;
datady=champdy->raw;
datadz=champdz->raw;

datadxy=champdxy->raw;
datadxz=champdxz->raw;
datadyz=champdyz->raw;


for (i=0;i<topDx;i++)
	for(j=0;j<topDy;j++)
	 	for (k=0;k<topDz;k++)
			{
			data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;
			datadx[i][j][k].x=datadx[i][j][k].y=datadx[i][j][k].z=0.0;
			datady[i][j][k].x=datady[i][j][k].y=datady[i][j][k].z=0.0;
			datadz[i][j][k].x=datadz[i][j][k].y=datadz[i][j][k].z=0.0;
			datadxy[i][j][k].x=datadxy[i][j][k].y=datadxy[i][j][k].z=0.0;
			datadyz[i][j][k].x=datadyz[i][j][k].y=datadyz[i][j][k].z=0.0;
			datadxz[i][j][k].x=datadxz[i][j][k].y=datadxz[i][j][k].z=0.0;
			}
			
			
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
      			dfdx=dfx[auxi-x00]*fy[auxj-y00]*fz[auxk-z00];
			dfdy=fx[auxi-x00]*dfy[auxj-y00]*fz[auxk-z00];
			dfdz=fx[auxi-x00]*fy[auxj-y00]*dfz[auxk-z00];
      			
			d2fdxdy=dfx[auxi-x00]*dfy[auxj-y00]*fz[auxk-z00];
			d2fdydz=fx[auxi-x00]*dfy[auxj-y00]*dfz[auxk-z00];
			d2fdxdz=dfx[auxi-x00]*fy[auxj-y00]*dfz[auxk-z00];
			
			
			ux=auxi-bx0;uy=auxj-by0;uz=auxk-bz0;
												
      			data[ux][uy][uz].x=data[ux][uy][uz].x+px*f;
      			data[ux][uy][uz].y=data[ux][uy][uz].y+py*f;
      			data[ux][uy][uz].z=data[ux][uy][uz].z+pz*f;
			
			datadx[ux][uy][uz].x=datadx[ux][uy][uz].x+px*dfdx;
      			datadx[ux][uy][uz].y=datadx[ux][uy][uz].y+py*dfdx;
      			datadx[ux][uy][uz].z=datadx[ux][uy][uz].z+pz*dfdx;
			
			datady[ux][uy][uz].x=datady[ux][uy][uz].x+px*dfdy;
      			datady[ux][uy][uz].y=datady[ux][uy][uz].y+py*dfdy;
      			datady[ux][uy][uz].z=datady[ux][uy][uz].z+pz*dfdy;
			
			
			datadz[ux][uy][uz].x=datadz[ux][uy][uz].x+px*dfdz;
      			datadz[ux][uy][uz].y=datadz[ux][uy][uz].y+py*dfdz;
      			datadz[ux][uy][uz].z=datadz[ux][uy][uz].z+pz*dfdz;
			
			
			datadxy[ux][uy][uz].x=datadxy[ux][uy][uz].x+px*d2fdxdy;
      			datadxy[ux][uy][uz].y=datadxy[ux][uy][uz].y+py*d2fdxdy;
      			datadxy[ux][uy][uz].z=datadxy[ux][uy][uz].z+pz*d2fdxdy;
			
			datadyz[ux][uy][uz].x=datadyz[ux][uy][uz].x+px*d2fdydz;
      			datadyz[ux][uy][uz].y=datadyz[ux][uy][uz].y+py*d2fdydz;
      			datadyz[ux][uy][uz].z=datadyz[ux][uy][uz].z+pz*d2fdydz;
			
			
			datadxz[ux][uy][uz].x=datadxz[ux][uy][uz].x+px*d2fdxdz;
      			datadxz[ux][uy][uz].y=datadxz[ux][uy][uz].y+py*d2fdxdz;
      			datadxz[ux][uy][uz].z=datadxz[ux][uy][uz].z+pz*d2fdxdz;
			
			
			
			}
  
		}

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;


reg=0;
	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	   {
    
		if (mask!=NULL)
			condition=mask->mri[ux+bx0][uy+by0][uz+bz0];
		else
			condition=1;
			
		
		if (condition>0)
		{
			
		// calcul du jacobien
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		
		Jacobien=determinant_mat3d(&mat);
		
		
		
		
		// calcul de la derivee du jacobien suivant x
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datadxy[ux][uy][uz].x;
		mat.H[1][1]= 		datadxy[ux][uy][uz].y;
		mat.H[2][1]= 		datadxy[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		
		Jx=determinant_mat3d(&mat);
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadxz[ux][uy][uz].x;
		mat.H[1][2]= 		datadxz[ux][uy][uz].y;
		mat.H[2][2]= 		datadxz[ux][uy][uz].z;
		
		Jx=Jx+determinant_mat3d(&mat);
		
	
		// calcul de la derivee du jacobien suivant y
		
		mat.H[0][0]= 		datadxy[ux][uy][uz].x;
		mat.H[1][0]= 		datadxy[ux][uy][uz].y;
		mat.H[2][0]= 		datadxy[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		Jy=determinant_mat3d(&mat);
	
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadyz[ux][uy][uz].x;
		mat.H[1][2]= 		datadyz[ux][uy][uz].y;
		mat.H[2][2]= 		datadyz[ux][uy][uz].z;
		
		Jy=Jy+determinant_mat3d(&mat);
		
		
		
		// calcul de la derivee du jacobien suivant z
		
		mat.H[0][0]= 	 	datadxz[ux][uy][uz].x;
		mat.H[1][0]= 		datadxz[ux][uy][uz].y;
		mat.H[2][0]= 		datadxz[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		Jz=determinant_mat3d(&mat);
	
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datadyz[ux][uy][uz].x;
		mat.H[1][1]= 		datadyz[ux][uy][uz].y;
		mat.H[2][1]= 		datadyz[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;	 
		 
		Jz=Jz+determinant_mat3d(&mat);
		
		if (Jacobien>0)
		reg=reg+(Jx*Jx+Jy*Jy+Jz*Jz)/(Jacobien*Jacobien);
		
		
		if (isnan(reg))
			printf("y a un bug reg\n");
		}

	}

free_field3d(champ);
free_field3d(champdx);
free_field3d(champdy);
free_field3d(champdz);
free_field3d(champdxy);
free_field3d(champdxz);
free_field3d(champdyz);



reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique
return(reg);
}


/*******************************************************************************
**  			regularisation_log_jacobien_local
**
**	Calcul local du terme de regularisation representant l'energie de membrane
*******************************************************************************/

double regularisation_log_jacobien_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,wdth,hght,dpth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int i,j,k;
int ideb,ifin,jdeb,jfin,kdeb,kfin;	
int  D=TOP_D(nb_param);
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2,mymono[21];
double auxdbl21a[21], ordonnee[21];                      // variables auxiliaire
double J=1,xt,yt,zt;
double reg=0;

resol=BASE3D.resol;
wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

topD=(int)pow(2.0,1.0*resol)-1;

  //-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
   for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;}
  

  for (aux0=0; aux0<8; aux0++)
  	Slpqr[aux0].ind = 7 - aux0;

  TOP_init_slpqr(Slpqr  ,topi-1,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+1,topi-1,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+2,topi-1,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+3,topi-1,topj  ,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+4,topi  ,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+5,topi  ,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+6,topi  ,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+7,topi  ,topj  ,topk  ,D,nb_param,param_norm);


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;

for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    
    	
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
	xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
	fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
	J=log(fabs(J));
	reg+= pow(J,2.0);
	
	 }
 }


reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique
return(reg);
}


/*******************************************************************************
**  			regularisation_log_jacobien_centre_local
**
*******************************************************************************/

double regularisation_log_jacobien_centre_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,wdth,hght,dpth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int i,j,k, auxi, auxj, auxk;
int ideb,ifin,jdeb,jfin,kdeb,kfin;	
int  D=TOP_D(nb_param);
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2,mymono[21];
double auxdbl21a[21], ordonnee[21];                      // variables auxiliaire
double J=1,xt,yt,zt,Jmean,Jtmp;
double reg=0;

resol=BASE3D.resol;
wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

topD=(int)pow(2.0,1.0*resol)-1;

  //-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
   for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;}
  

  for (aux0=0; aux0<8; aux0++)
  	Slpqr[aux0].ind = 7 - aux0;

  TOP_init_slpqr(Slpqr  ,topi-1,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+1,topi-1,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+2,topi-1,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+3,topi-1,topj  ,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+4,topi  ,topj-1,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+5,topi  ,topj-1,topk  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+6,topi  ,topj  ,topk-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+7,topi  ,topj  ,topk  ,D,nb_param,param_norm);


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
topDx=bx1-bx0;
topDy=by1-by0;
topDz=bz1-bz0;

for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    
    	
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
	xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
	fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
	
	J=log(fabs(J));
	
	Jmean=0;
	
	for (auxi=i-1;auxi<i+2;auxi++)
	for (auxj=j-1;auxj<j+2;auxj++)
	for (auxk=k-1;auxk<k+2;auxk++)
		{
		xt=(double)(auxi)/wdth;yt=(double)(auxj)/hght;zt=(double)(auxk)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &Jtmp, mymono);
		Jmean+=log(fabs(Jtmp));
		}
	
	Jmean=Jmean/27.0;
	J=J-Jmean;
	reg+= pow(J,2.0);
	
	 }
 }


reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique
return(reg);
}

/*******************************************************************************
**  			regularisation_dist_identite_local
**
*******************************************************************************/

double regularisation_dist_identite_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
double reg;
field3d *champ;
vector3d ***data;
int ideb,jdeb,kdeb;


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

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;


reg=0;
	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	   {
    
		
			reg+= data[ux][uy][uz].x*data[ux][uy][uz].x+data[ux][uy][uz].y*data[ux][uy][uz].y+data[ux][uy][uz].z*data[ux][uy][uz].z;
			 
		 }


free_field3d(champ);

reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique
return(reg);
}



/*******************************************************************************
**  Energie_globale_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param, double *param_norm, dist_func_t dist)
**
**	Calcule l'energie globale selon le critere dist entre deux image
*******************************************************************************/

double regularisation_globale_3d(int nb_param, double *param, double *param_norm,int ***masque_param,grphic3d *mask, reg_func_locale_t regularisation)
{

if ((regularisation==regularisation_energie_membrane_local)||(regularisation==regularisation_energie_membrane_Lp_local)||(regularisation==regularisation_energie_membrane_jacobien_local))
{
double E;
int topi,topj,topk,topD,resol,l;
  
resol=BASE3D.resol; 
topD=(int)pow(2.0,1.0*resol)-1;

E=0.0;
	for (topi=0; topi<topD; topi=topi+2)
  	for (topj=0; topj<topD; topj=topj+2)
  		for (topk=0; topk<topD; topk=topk+2)
  			if (masque_param[topi][topj][topk]!=0)
				{
				l = TOP_conv_ind(topi,topj,topk,nb_param);
				E+=regularisation(nb_param, param,param_norm, topi,  topj,  topk, mask, NULL);
				}


E=sqrt(E);

return(E); 
}
else
{
double E;
int topi,topj,topk,topD,resol,l,i;
int *x00,*x11,*y00,*y11,*z00,*z11;
int width,height,depth;
double xm,xM,ym,yM,zm,zM;
TSlpqr Slpqr[8];
   
width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
resol=BASE3D.resol; 
x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
topD=(int)pow(2.0,1.0*resol)-1;

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
	{
	Slpqr[i].Jm=-HUGE_VAL;
	Slpqr[i].JM=HUGE_VAL;
	} 
  
E=0.0;
	for (topi=0; topi<topD; topi=topi+2)
  	for (topj=0; topj<topD; topj=topj+2)
  		for (topk=0; topk<topD; topk=topk+2)
  			if (masque_param[topi][topj][topk]!=0)
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

				E+=regularisation(nb_param, param,param_norm, topi,  topj,  topk, mask,Slpqr);
			  }
E=sqrt(E);
return(E); 
}

}


/*******************************************************************************
**  Energie_quad_locale_pondere_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_quad_locale_pondere_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3,Ntot;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh,Ptot,poids;
field3d *champ;
vector3d ***data;

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
		

		    
		poids=imref->mask->mri[i+bx0][j+by0][k+bz0];
		
		Ptot+=poids;
		   
		auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		Ntot++;		
		E=E+poids*auxE*auxE;
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

        auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
		
		
			
		poids=imref->mask->mri[i+bx0][j+by0][k+bz0];
		
		Ptot+=poids;
		Ntot++;
		
			E=E+poids*auxE*auxE;
       }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL;
       }
      }
   }
   
if (Ptot>0)
E=1.0*Ntot*E/Ptot;
else 
E=0;

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);

free_field3d(champ);

return(E);
}


/*******************************************************************************
**  Energie_ML_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule de la vraisemblance relative a une boite entre une image labelisée et une image réelle
*******************************************************************************/

double Energie_ML_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3,Ntot;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh,Ptot;
field3d *champ;
vector3d ***data;
int label;



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
		 label=imref->mask->mri[i+bx0][j+by0][k+bz0]-1;
			if ((label>-1)&&(MELANGE_GAUSSIENNE_2D.param[label].l>0))		
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
		    
		auxE=auxE-MELANGE_GAUSSIENNE_2D.param[label].my;
		auxE=auxE*auxE*MELANGE_GAUSSIENNE_2D.param[label].sx-2*auxE*(imref->mri[i+bx0][j+by0][k+bz0]-MELANGE_GAUSSIENNE_2D.param[label].mx)*MELANGE_GAUSSIENNE_2D.param[label].sxy;
    auxE=auxE/(MELANGE_GAUSSIENNE_2D.param[label].sx*MELANGE_GAUSSIENNE_2D.param[label].sy-MELANGE_GAUSSIENNE_2D.param[label].sxy*MELANGE_GAUSSIENNE_2D.param[label].sxy);
		  
			if ((isnan(auxE))||(auxE>1000000000))
					printf("Y a un bug \n");
			
			E+=auxE;
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
	
			auxE=auxE-MELANGE_GAUSSIENNE_2D.param[label].my;
	
			auxE=auxE*auxE*MELANGE_GAUSSIENNE_2D.param[label].sx-2*auxE*(imref->mri[i+bx0][j+by0][k+bz0]-MELANGE_GAUSSIENNE_2D.param[label].mx)*MELANGE_GAUSSIENNE_2D.param[label].sxy;
    	auxE=auxE/(MELANGE_GAUSSIENNE_2D.param[label].sx*MELANGE_GAUSSIENNE_2D.param[label].sy-MELANGE_GAUSSIENNE_2D.param[label].sxy*MELANGE_GAUSSIENNE_2D.param[label].sxy);
	 
			if ((isnan(auxE))||(auxE>1000000000))
						printf("Y a un bug \n");
		
		E+=auxE;
		}
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL;
				printf("Attention E vaut HUGE_VAL\n");
       }
      }
   }
   }


	
if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


//E=sqrt(E);

if (isnan(E))
	printf("Y a un bug \n");

free_field3d(champ);
return(E);
}


/*******************************************************************************
**  Calcul_IM_global_3d
**
*******************************************************************************/

double Calcul_IM_global_3d(grphic3d *imref, grphic3d *imreca, int nb_param, double *param)
{
grphic3d *imres;
int wdth, hght, dpth;
field3d *champ;

wdth=imref->width;
hght=imref->height;
dpth=imref->depth;


/* calcul de l'image déformée */
imres=cr_grphic3d(imreca);imx_copie_param_3d_p(imreca,imres);
champ=cr_field3d(wdth,hght,dpth);
base_to_field_3d(nb_param,param,champ,NULL,NULL);
inter_linear_3d(imreca,champ,imres); 

update_HistoJoint_global(imref,imres);


free_field3d(champ);
free_grphic3d(imres);

return (HISTOJOINT.IM);
}



/*******************************************************************************
**  Energie_IM_locale_3d
**
*******************************************************************************/

double Energie_IM_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
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
double va,vb,vc,vd,ve,vf,vg,vh,Ptot,IM,Ereg;
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


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
	HISTOJOINT.aux_hij[i][j]=HISTOJOINT.hij[i][j];



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

IM=0;

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
		if ((i0<HISTOJOINT.nbI)&&(j0<HISTOJOINT.nbJ))HISTOJOINT.aux_hij[i0][j0]+=(i0-x+1.0)*(j0-y+1.0);
		if (i0<nbI) HISTOJOINT.aux_hij[i0+1][j0]+=(x-i0)*(j0-y+1.0);
		if (j0<nbJ) HISTOJOINT.aux_hij[i0][j0+1]+=(i0-x+1.0)*(y-j0);
		if ((i0<nbI)&&(j0<nbJ)) HISTOJOINT.aux_hij[i0+1][j0+1]+=(x-i0)*(y-j0);
		}
		
			}
      else
      {
      /* if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
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
		if ((i0<HISTOJOINT.nbI)&&(j0<HISTOJOINT.nbJ)) HISTOJOINT.aux_hij[i0][j0]+=(i0-x+1.0)*(j0-y+1.0);
		if (i0<nbI) HISTOJOINT.aux_hij[i0+1][j0]+=(x-i0)*(j0-y+1.0);
		if (j0<nbJ) HISTOJOINT.aux_hij[i0][j0+1]+=(i0-x+1.0)*(y-j0);
		if ((i0<nbI)&&(j0<nbJ)) HISTOJOINT.aux_hij[i0+1][j0+1]+=(x-i0)*(y-j0);

		}	
			}*/
			IM=HUGE_VAL;
	
      }
   }


free_field3d(champ);

//calcul de l'IM
for (i=0;i<HISTOJOINT.nbI;i++)
	HISTOJOINT.aux_hi[i]=0.0;




for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
	HISTOJOINT.aux_hi[i]+=HISTOJOINT.aux_hij[i][j];


for (j=0;j<HISTOJOINT.nbJ;j++)
	HISTOJOINT.aux_hj[j]=0.0;

for (j=0;j<HISTOJOINT.nbJ;j++)
for (i=0;i<HISTOJOINT.nbI;i++)
	HISTOJOINT.aux_hj[j]+=HISTOJOINT.aux_hij[i][j];


Ntot=0;
for (i=0;i<HISTOJOINT.nbI;i++)
	Ntot+=HISTOJOINT.aux_hi[i];

HISTOJOINT.Ntot=Ntot;


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
if (HISTOJOINT.aux_hij[i][j]>0.00000001) /* pour éviter des bug numériques */
	IM-=HISTOJOINT.aux_hij[i][j]*log(HISTOJOINT.Ntot*HISTOJOINT.aux_hij[i][j]/HISTOJOINT.aux_hi[i]/HISTOJOINT.aux_hj[j]);
	
	
	
IM=IM/HISTOJOINT.Ntot;

	
if (IM!=0)	
	IM=-1/IM;

if (isnan(IM))
	IM=HUGE_VAL;

HISTOJOINT.aux_IM=IM;

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

IM+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

return(IM);
}


/*******************************************************************************
**  Energie_ICP_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie entre deux nuages de points 
** imref est une image binaire
** imreca est une carte de distance
*******************************************************************************/

double Energie_ICP_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
grphic3d *mask=NULL;
mask=imref;

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
        
				
			E=E+auxE*auxE;
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
       
			E=E+auxE*auxE;
				   
       }
       else
       {
	      	E=HUGE_VAL;
       }
      }
		}

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param,param_norm, topi,  topj,  topk,mask,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


free_field3d(champ);

E=sqrt(E);

return(E);

}

/*******************************************************************************
**  Energie_ICP_sym_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie entre deux nuages de points 
** imref est une image binaire
** imreca est une carte de distance
*******************************************************************************/

double Energie_ICP_sym_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
grphic3d *mask=NULL;
mask=imref;

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
        
				
			E=E+auxE*auxE;
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
       
			E=E+auxE*auxE;
				   
       }
       else
       {
	      	E=HUGE_VAL;
       }
      }
		}


// On balaye chacun des pixels de la boite de l'image complémentaire
for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	   if (imref->mask->mri[i+bx0][j+by0][k+bz0]>0)
		 {
    dx=(double)i+data[i][j][k].x+bx0;
    dy=(double)j+data[i][j][k].y+by0;
    dz=(double)k+data[i][j][k].z+bz0;
    i2=(int)dx;j2=(int)dy;k2=(int)dz;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=(double)imreca->mask->mri[i2][j2][k2];
        vb=(double)imreca->mask->mri[i2][j2][k3];
        vc=(double)imreca->mask->mri[i2][j3][k2];
        vd=(double)imreca->mask->mri[i2][j3][k3];
        ve=(double)imreca->mask->mri[i3][j2][k2];
        vf=(double)imreca->mask->mri[i3][j2][k3];
        vg=(double)imreca->mask->mri[i3][j3][k2];
        vh=(double)imreca->mask->mri[i3][j3][k3];

        auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        
				
			E=E+auxE*auxE;
      }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
	
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
				if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imreca->mask->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imreca->mask->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imreca->mask->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imreca->mask->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imreca->mask->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imreca->mask->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imreca->mask->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imreca->mask->mri[i3][j3][k3];

        auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
       
			E=E+auxE*auxE;
				   
       }
       else
       {
	      	E=HUGE_VAL;
       }
      }
		}
if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param,param_norm, topi,  topj,  topk,mask,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


free_field3d(champ);

E=sqrt(E);

return(E);

}

/*******************************************************************************
**  Energie_atrophie_jacobien_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie relatif a la simulation d'atrophie avec comparaison des cartes de jacobien (norme L2)
*******************************************************************************/
  
double Energie_atrophie_jacobien_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz;
int i,j,k;
double E=0;
double auxE;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double Ereg;
 
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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


	

topD=(int)pow(2.0,1.0*resol)-1;


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
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
  
ideb=Slpqr[aux0].b.xm*width;
ifin=Slpqr[aux0].b.xM*width;
jdeb=Slpqr[aux0].b.ym*height;
jfin=Slpqr[aux0].b.yM*height;
kdeb=Slpqr[aux0].b.zm*depth;
kfin=Slpqr[aux0].b.zM*depth;

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	   if (imref->mri[i][j][k]>0)
	{
	
        
			xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
			fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
			auxE=J-(double)imref->mri[i][j][k]*imref->rcoeff;
			E=E+auxE*auxE;
			
	    }
}   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);
E=E/2.0;
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=sqrt(E);
return(E);
} 

/*******************************************************************************
**  Energie_atrophie_log_jacobien_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie relatif a la simulation d'atrophie avec comparaison des cartes de jacobien (norme L2)
*******************************************************************************/
  
double Energie_atrophie_log_jacobien_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz;
int i,j,k;
double E=0;
double auxE;
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2;
double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaire
int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
int nb_func,nbfx,nbfy,nbfz;
double J=1,mymono[21],xt,yt,zt;
int ideb,ifin,jdeb,jfin,kdeb,kfin;
double Ereg;
 
 
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;resol=BASE3D.resol;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;

	
  
       
  aux0 = TOP_conv_ind(topi,topj,topk,nb_param);
  

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


	

topD=(int)pow(2.0,1.0*resol)-1;


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
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
  
ideb=Slpqr[aux0].b.xm*width;
ifin=Slpqr[aux0].b.xM*width;
jdeb=Slpqr[aux0].b.ym*height;
jfin=Slpqr[aux0].b.yM*height;
kdeb=Slpqr[aux0].b.zm*depth;
kfin=Slpqr[aux0].b.zM*depth;

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	   if (imref->mri[i][j][k]>0)
		 {
	
        
			xt=(double)(i)/width;yt=(double)(j)/height;zt=(double)(k)/depth;
			fast_eval_fun(auxdbl21b, xt, yt, zt, &J, mymono);
			auxE=log(J/(double)imref->mri[i][j][k]/imref->rcoeff);
			auxE=auxE*auxE;
			
			if ((!isnan(auxE))&&(!isinf(auxE)))
				E=E+auxE;
			
	    }
}   

E=E*100000;

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);
E=E/2.0;
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=sqrt(E);
return(E);
} 


/*******************************************************************************
**  Energie_quad_locale_symetrique_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique symetrique relative a une boite
*******************************************************************************/

double Energie_quad_locale_symetrique_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dxreca,dyreca,dzreca,dxref,dyref,dzref,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
double valreca, valref;
field3d *champ;
vector3d ***data;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;



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
	   {
    dxreca=(double)i+data[i][j][k].x+bx0;
    dyreca=(double)j+data[i][j][k].y+by0;
    dzreca=(double)k+data[i][j][k].z+bz0;
 
 
    i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
		
		/* calcul de la valeur interpolée dans imreca */
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
	
        va=(double)imreca->mri[i2][j2][k2];
        vb=(double)imreca->mri[i2][j2][k3];
        vc=(double)imreca->mri[i2][j3][k2];
        vd=(double)imreca->mri[i2][j3][k3];
        ve=(double)imreca->mri[i3][j2][k2];
        vf=(double)imreca->mri[i3][j2][k3];
        vg=(double)imreca->mri[i3][j3][k2];
        vh=(double)imreca->mri[i3][j3][k3];

        valreca=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
		
				
      }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
	
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
				if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imreca->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imreca->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imreca->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imreca->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imreca->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imreca->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imreca->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imreca->mri[i3][j3][k3];

        valreca=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
		   
				
		
       }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL; valreca=0;
				
       }
      }
			
	dxref=(double)i-lambda_ref*data[i][j][k].x+bx0;
    	dyref=(double)j-lambda_ref*data[i][j][k].y+by0;
    	dzref=(double)k-lambda_ref*data[i][j][k].z+bz0;
 
 
    i2=(int)dxref;j2=(int)dyref;k2=(int)dzref;
    i3=i2+1;j3=j2+1;k3=k2+1;
    

			
		/* calcul de la valeur interpolée dans imref */
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
	
      
		    va=(double)imref->mri[i2][j2][k2];
        vb=(double)imref->mri[i2][j2][k3];
        vc=(double)imref->mri[i2][j3][k2];
        vd=(double)imref->mri[i2][j3][k3];
        ve=(double)imref->mri[i3][j2][k2];
        vf=(double)imref->mri[i3][j2][k3];
        vg=(double)imref->mri[i3][j3][k2];
        vh=(double)imref->mri[i3][j3][k3];

        valref=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
		
					
      }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
	
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
				
				
				if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imref->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imref->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imref->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imref->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imref->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imref->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imref->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imref->mri[i3][j3][k3];

        valref=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;

	 
       }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL;valref=0;
				
       }
      }
		
		auxE=(valreca*imreca->rcoeff-valref*imref->rcoeff);
		E=E+auxE*auxE;
	
		
				
   }
  
//E=0.5*E;	 

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk, NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);

free_field3d(champ);

return(E);

}

/*******************************************************************************
**  Energie_quad_locale_symetrique_coupe_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique symetrique relative a une boite
*******************************************************************************/

double Energie_quad_locale_symetrique_coupe_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dxreca,dyreca,dzreca,dxref,dyref,dzref,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
double valreca, valref;
field3d *champ;
vector3d ***data;
int Ntot=0;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;

grphic3d *mask=NULL;
mask=imref; // Il faudrait calculer le mask intelligement

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
	   {
    dxreca=(double)i+data[i][j][k].x+bx0;
    dyreca=(double)j+data[i][j][k].y+by0;
    dzreca=(double)k+data[i][j][k].z+bz0;
 
 
    i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
		
		/* calcul de la valeur interpolée dans imreca */
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
	
        va=(double)imreca->mri[i2][j2][k2];
        vb=(double)imreca->mri[i2][j2][k3];
        vc=(double)imreca->mri[i2][j3][k2];
        vd=(double)imreca->mri[i2][j3][k3];
        ve=(double)imreca->mri[i3][j2][k2];
        vf=(double)imreca->mri[i3][j2][k3];
        vg=(double)imreca->mri[i3][j3][k2];
        vh=(double)imreca->mri[i3][j3][k3];

        valreca=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
		
				
      }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
	
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
				if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imreca->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imreca->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imreca->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imreca->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imreca->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imreca->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imreca->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imreca->mri[i3][j3][k3];

        valreca=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
		   
				
		
       }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL; valreca=0;
				
       }
      }
			
	dxref=(double)i-lambda_ref*data[i][j][k].x+bx0;
    	dyref=(double)j-lambda_ref*data[i][j][k].y+by0;
    	dzref=(double)k-lambda_ref*data[i][j][k].z+bz0;
 
 
    i2=(int)dxref;j2=(int)dyref;k2=(int)dzref;
    i3=i2+1;j3=j2+1;k3=k2+1;
    

			
		/* calcul de la valeur interpolée dans imref */
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
	
      
		    va=(double)imref->mri[i2][j2][k2];
        vb=(double)imref->mri[i2][j2][k3];
        vc=(double)imref->mri[i2][j3][k2];
        vd=(double)imref->mri[i2][j3][k3];
        ve=(double)imref->mri[i3][j2][k2];
        vf=(double)imref->mri[i3][j2][k3];
        vg=(double)imref->mri[i3][j3][k2];
        vh=(double)imref->mri[i3][j3][k3];

        valref=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
		
					
      }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
	
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
				
				
				if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imref->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)imref->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)imref->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)imref->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)imref->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)imref->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)imref->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)imref->mri[i3][j3][k3];

        valref=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
		+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
		+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;

	 
       }
       else
       {
	      /*auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
					E=E+auxE*auxE; */
				E=HUGE_VAL;valref=0;
				
       }
      }
		
		if ((valreca>0)&&(valref>0))
			{auxE=(valreca-valref);
			E=E+auxE*auxE;
			Ntot++;
			}
   }
  
if (Ntot>0)
//E=0.5*E/(double)Ntot;	 
E=E/(double)Ntot;	 
else
E=0;

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,mask,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);

free_field3d(champ);

return(E);

}



/*******************************************************************************
**  Energie_quad_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite
*******************************************************************************/

double Energie_quad_sous_champ_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz;
int ux,uy,uz;
int i2,j2,k2;
field3d *champ;
vector3d ***data;
	
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
		{
		i2=i+bx0;
		j2=j+by0;
		k2=k+bz0;
		
		if (imref->mri[i2][j2][k2]>0)
			{
			dx=data[i][j][k].x-_ParamRecalageBspline.chref->raw[i2][j2][k2].x;
			dy=data[i][j][k].y-_ParamRecalageBspline.chref->raw[i2][j2][k2].y;
			dz=data[i][j][k].z-_ParamRecalageBspline.chref->raw[i2][j2][k2].z;
		
			E+=dx*dx+dy*dy+dz*dz;
			}
   		}
   

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);

free_field3d(champ);

return(E);

}


/*******************************************************************************
**  Energie_DTI_quad_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite pour une série DTI
*******************************************************************************/

double Energie_DTI_quad_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{


#ifdef COMPILE_FOR_MEDIPY
printf("Cette fonctionnalite n'est pas encore exportee sous medipy!\n");
return(-1);
#else
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dx,dy,dz,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
field3d *champ;
vector3d ***data;
	

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


for (l=0;l<_ParamRecalageBspline.item_ref->nb_directions;l++)
{
	// On balaye chacun des pixels de la boite
	for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
	for (k=0;k<topDz;k++)
	{
    		dx=(double)i+data[i][j][k].x+bx0;
    		dy=(double)j+data[i][j][k].y+by0;
    		dz=(double)k+data[i][j][k].z+bz0;
    		i2=(int)dx;j2=(int)dy;k2=(int)dz;
    		i3=i2+1;j3=j2+1;k3=k2+1;
    
    		if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      		{
        		xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        		va=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j2][k2];
        		vb=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j2][k3];
        		vc=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j3][k2];
        		vd=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j3][k3];
        		ve=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j2][k2];
        		vf=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j2][k3];
        		vg=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j3][k2];
        		vh=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j3][k3];

       		 	auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
			+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
			+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        
			auxE=auxE-_ParamRecalageBspline.item_ref->img_3d[l]->mri[i+bx0][j+by0][k+bz0];
			
				
			E=E+auxE*auxE;	
      		}
      		else
      		{
       			if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       			{
 		       	xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
			if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j2][k2];
        		if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j2][k3];
        		if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j3][k2];
        		if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i2][j3][k3];
        		if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j2][k2];
        		if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j2][k3];
        		if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j3][k2];
        		if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.item_reca->img_3d[l]->mri[i3][j3][k3];
	
        		auxE=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
			+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
			+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
       
			 auxE=auxE-_ParamRecalageBspline.item_ref->img_3d[l]->mri[i+bx0][j+by0][k+bz0];
			E=E+auxE*auxE;  
       			}
       			else
      			{
	   		E=HUGE_VAL;
       			}
      		}
   	}	
   
}

if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);

free_field3d(champ);

return(E);
#endif //fin COMPILE_FOR_MEDIPY
}




/*******************************************************************************
**  Energie_groupwise_variance_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie locale quadratique relative a une boite pour une série DTI
*******************************************************************************/

double Energie_groupwise_variance_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dxreca,dyreca,dzreca,dxref,dyref,dzref,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
double valreca, valref, valref2,tmp;
field3d *champ;
vector3d ***data;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
int numref=_ParamRecalageBspline.nb_reca_groupwise-1;



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
	   	{
    	/* calcul de la valeur interpolée dans reca */
   		dxreca=(double)i+data[i][j][k].x+bx0;
    	dyreca=(double)j+data[i][j][k].y+by0;
    	dzreca=(double)k+data[i][j][k].z+bz0;
 
     	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
		
		if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       		{
 		   	xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
			
			if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k2];
        	if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k3];
        	if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k2];
        	if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k3];
        	if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k2];
        	if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k3];
        	if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k2];
        	if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k3];
	
        	valreca=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
			+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
			+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
      		}
       	else
      		{
	   		free_field3d(champ);
			return(HUGE_VAL);
			}
      			
		/* calcul de la valeur interpolée dans ref et ref2 */
    	dxref=(double)i-lambda_ref*data[i][j][k].x+bx0;
    	dyref=(double)j-lambda_ref*data[i][j][k].y+by0;
    	dzref=(double)k-lambda_ref*data[i][j][k].z+bz0;
 
     	i2=(int)dxref;j2=(int)dyref;k2=(int)dzref;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    			
		if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       			{
 		       	xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
	      		
			/*if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j2][k2];
        		if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j2][k3];
        		if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j3][k2];
        		if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j3][k3];
        		if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j2][k2];
        		if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j2][k3];
        		if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j3][k2];
        		if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j3][k3];*/
			
			if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.dbl_ref_groupwise[i2][j2][k2];
        		if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.dbl_ref_groupwise[i2][j2][k3];
        		if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.dbl_ref_groupwise[i2][j3][k2];
        		if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.dbl_ref_groupwise[i2][j3][k3];
        		if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.dbl_ref_groupwise[i3][j2][k2];
        		if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.dbl_ref_groupwise[i3][j2][k3];
        		if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.dbl_ref_groupwise[i3][j3][k2];
        		if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.dbl_ref_groupwise[i3][j3][k3];
        		
			valref=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
				+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
				+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
      			
				
			/*if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.ref2_groupwise->mri[i2][j2][k2];
        		if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.ref2_groupwise->mri[i2][j2][k3];
        		if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.ref2_groupwise->mri[i2][j3][k2];
        		if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.ref2_groupwise->mri[i2][j3][k3];
        		if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.ref2_groupwise->mri[i3][j2][k2];
        		if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.ref2_groupwise->mri[i3][j2][k3];
        		if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.ref2_groupwise->mri[i3][j3][k2];
        		if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.ref2_groupwise->mri[i3][j3][k3];*/
	
			/*if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i2][j2][k2];
        		if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i2][j2][k3];
        		if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i2][j3][k2];
        		if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i2][j3][k3];
        		if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i3][j2][k2];
        		if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i3][j2][k3];
        		if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i3][j3][k2];
        		if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.dbl_ref2_groupwise[i3][j3][k3];
	
				valref2=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
				+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
				+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;*/
				
				
				valref2=0.0;
			
			for (l=0; l<_ParamRecalageBspline.nb_tot_groupwise; l++)
			if (l!=numref)
			{
			if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i2][j2][k2];
        		if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i2][j2][k3];
        		if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i2][j3][k2];
        		if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i2][j3][k3];
        		if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i3][j2][k2];
        		if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i3][j2][k3];
        		if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i3][j3][k2];
        		if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.serie_groupwise[l]->mri[i3][j3][k3];
	
				tmp=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
				+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
				+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
				
				tmp=tmp*_ParamRecalageBspline.serie_groupwise[l]->rcoeff;
				
			valref2+=tmp*tmp;
			}
			
		
				}
       		else
      			{
	   			free_field3d(champ);
				return(HUGE_VAL);
				}	
		
		
		valreca=valreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
	
		auxE=(valref+valreca)/_ParamRecalageBspline.nb_tot_groupwise;
	
		tmp=(valref2+valreca*valreca)/_ParamRecalageBspline.nb_tot_groupwise-auxE*auxE;
		
		
		
		//printf("[ %d , %d , %d ] : %f \n",i,j,k,tmp);
		
		//if (tmp>100)
			//printf("Ya un bug\n");
		
		if (tmp<0)
			tmp=0;
		
			
		E=E+tmp;	
   		}
  

//E=2.0*E; // Pour donner les memes resultats que pour le recalage symetrique avec 2 images


if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk, NULL,Slpqr);
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=2.0*sqrt(E);

free_field3d(champ);

return(E);

}

/*******************************************************************************
**  Energie_groupwise_variance_locale_nonsym_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
*******************************************************************************/

double Energie_groupwise_variance_locale_nonsym_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
double E=0,Ereg;
double dxreca,dyreca,dzreca,auxE;
int ux,uy,uz;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double va,vb,vc,vd,ve,vf,vg,vh;
double valreca, valref, valref2,tmp;
field3d *champ;
vector3d ***data;



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
	   	{
    	/* calcul de la valeur interpolée dans reca */
   		dxreca=(double)i+data[i][j][k].x+bx0;
    	dyreca=(double)j+data[i][j][k].y+by0;
    	dzreca=(double)k+data[i][j][k].z+bz0;
 
     	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
		
		if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       		{
 		   	xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
			
			if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k2];
        	if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k3];
        	if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k2];
        	if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k3];
        	if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k2];
        	if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k3];
        	if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k2];
        	if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k3];
	
        	valreca=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
			+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
			+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
      		}
       	else
      		{
	   		free_field3d(champ);
			return(HUGE_VAL);
			}
      			
				
		valref=_ParamRecalageBspline.dbl_ref_groupwise[i+bx0][j+by0][k+bz0];
		valref2=_ParamRecalageBspline.dbl_ref2_groupwise[i+bx0][j+by0][k+bz0];
		
      	//valref=valref*_ParamRecalageBspline.ref_groupwise->rcoeff;
		//valref2=valref2*_ParamRecalageBspline.ref2_groupwise->rcoeff;
		valreca=valreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
	
	
		auxE=(valref+valreca)/_ParamRecalageBspline.nb_tot_groupwise;
	
		tmp=(valref2+valreca*valreca)/_ParamRecalageBspline.nb_tot_groupwise-auxE*auxE;
		
		//printf("[ %d , %d , %d ] : %f \n",i,j,k,tmp);
		
		//if (tmp>100)
			//printf("Ya un bug\n");
		
		if (tmp<0)
			tmp=0;
		
		E=E+tmp;	
   		}
  

//E=2.0*E; // Pour donner les memes resultats que pour le recalage symetrique avec 2 images


if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk, NULL,Slpqr);
E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}

E=2.0*sqrt(E);

free_field3d(champ);

return(E);

}



/*******************************************************************************
**  Energie_patch_locale_3d(imref,imreca,imres,param,champ_ini,topi,topj,topk)
**
**	Calcule l'energie patch quadratique relative a une boite
*******************************************************************************/

double Energie_patch_locale_3d(grphic3d *imref, grphic3d *imreca,int nb_param, double *param, double *param_norm,int topi, int topj, int topk, TSlpqr *Slpqr,reg_func_locale_t regularisation)
{

int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz;
int i,j,k;
double E=0,Ereg;
double auxE;

dvector3d u;
int xref,yref,zref;
double xreca,yreca,zreca;
int wx,wy,wz;
int patchx,patchy,patchz;
int hwn = 2;
int nbp = hwn*2 + 1;
int hwp = 1;


 int wx_ref,wy_ref,wz_ref;
 int wx_reca,wy_reca,wz_reca;
 int xpx,ypy,zpz,wxpx,wypy,wzpz;
 int x2px,y2py,z2pz,wx2px,wy2py,wz2pz;
 int indice;

 //precalcul du champ de deplacement
 field3d *champ;
 vector3d ***data;
	
	

int xmin,ymin,zmin,xmax,ymax,zmax;

int interpolation_choice = 0; // 0 NN, 1 TL
int costfunction = 1; // 0 : L2, 1, patch
//int patchdef = 0; //0 non deforme, 1 deforme // non utilise pour le moment

int chpx,chpy,chpz;
	
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

nbp = nbp*nbp*nbp;


 //taille du champ :
 xmin = bx0 - hwn - hwp;
 if(xmin < 0)
   xmin = 0;
 ymin = by0 - hwn - hwp;
 if(ymin < 0)
  ymin = 0;
 zmin = bz0 - hwn - hwp;
 if(zmin < 0)
   zmin = 0;
 xmax = bx1 + hwn + hwp;
 if(xmax >= width)
   xmax = width-1;
 ymax = by1 + hwn + hwp;
 if(ymax >= height)
   ymax = height-1;
 zmax = bz1 + hwn + hwp;
 if(zmax >= depth)
   zmax = depth-1;

 //int chpx = xmax-xmin+1;
 //int chpy = ymax-ymin+1;
 //int chpz = zmax-zmin+1;

 //debug. calcul sur toute l'image
 chpx = width;
 chpy = height;
 chpz = depth;

 champ=cr_field3d(chpx,chpy,chpz);
 data=champ->raw;
 for (i=0;i<chpx;i++)
   for(j=0;j<chpy;j++)
     for (k=0;k<chpz;k++)
       data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;

 for(i=0;i<chpx;i++)
   for(j=0;j<chpy;j++)
     for(k=0;k<chpz;k++)
       {
	 //eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, i+xmin, j+ymin, k+zmin, &u);
         eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, i, j, k, &u);      
	 data[i][j][k].x=u.x;
	 data[i][j][k].y=u.y;
	 data[i][j][k].z=u.z;
       }


// On balaye chacun des pixels de la boite

for (i=0;i<topDx;i++)
  for (j=0;j<topDy;j++)
    for (k=0;k<topDz;k++)
    {
     //critere patch
     xref = i+bx0;
     yref = j+by0;
     zref = k+bz0;
     eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, xref, yref, zref, &u);   
     xreca = xref + u.x;
     yreca = yref + u.y;
     zreca = zref + u.z;
     

     //L2
     if(costfunction==0)
       {
       if(interpolation_choice==0)
         {
         xreca = (int)(xreca);
         yreca = (int)(yreca);
         zreca = (int)(zreca);
	 }
       if( (xreca>=0) & (xreca<width) & (yreca>=0) & (yreca<height) & (zreca>=0) & (zreca<depth) )
         {
           auxE = image_trilinear_value(imreca, xreca, yreca, zreca)  - imref->mri[xref][yref][zref];
           //auxE = imreca->mri[(int)(xreca)][(int)(yreca)][(int)(zreca)] - imref->mri[xref][yref][zref];
           E = E + auxE*auxE;
         }
       }   
     
     //souci dans le calcul du champ
     /*
     if( ( fabs(data[xref-xmin][yref-ymin][zref-zmin].x - u.x) > 0.01 ) || (fabs(data[xref-xmin][yref-ymin][zref-zmin].y - u.y) > 0.01) || (fabs(data[xref-xmin][yref-ymin][zref-zmin].z - u.z) > 0.01) )
       { 
	 printf("%d %d %d, %f %f %f \n",xref,yref,zref,u.x,u.y,u.z);
	 printf("%f %f %f\n",data[xref-xmin][yref-ymin][zref-zmin].x, data[xref-xmin][yref-ymin][zref-zmin].y, data[xref-xmin][yref-ymin][zref-zmin].z);
       }
     */

     /*
     if( ( fabs(data[xref][yref][zref].x - u.x) > 0.01 ) || (fabs(data[xref][yref][zref].y - u.y) > 0.01) || (fabs(data[xref][yref][zref].z - u.z) > 0.01) )
       { 
	 printf("%d %d %d, %f %f %f \n",xref,yref,zref,u.x,u.y,u.z);
	 printf("%f %f %f\n",data[xref][yref][zref].x, data[xref][yref][zref].y, data[xref][yref][zref].z);
       }
     */
     if(costfunction==1)
     {
     double wref[nbp];
     double wreca[nbp];
     double sum_ref = 0;
     double sum_reca = 0;

     int count = 0;

	 for(indice=0;indice<nbp;indice++)
       {
	 wref[indice] = 0;
	 wreca[indice] = 0;
       }
     
     //parcours du voisinage
     for(wx=-hwn;wx<hwn+1;wx++)
       for(wy=-hwn;wy<hwn+1;wy++)
         for(wz=-hwn;wz<hwn+1;wz++)
	   {
	     //position du voisin dans ima_ref
	     wx_ref = xref + wx;
	     wy_ref = yref + wy;
	     wz_ref = zref + wz;

	     if( (wx_ref-xmin>=0) & (wx_ref<width) & (wy_ref-ymin>=0) & (wy_ref<height) & (wz_ref-zmin>=0) & (wz_ref<depth) )
	       {
		 double poids_ref = 0;
		 double dist_ref = 0;
		 double tmpd_ref = 0;
		 double poids_reca = 0;
		 double dist_reca = 0;
		 double tmpd_reca = 0;
		
		 	       
		 //eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, wx_ref, wy_ref, wz_ref, &u);
		 //position du voisin dans ima_reca
		 //wx_reca = wx_ref + u.x;
		 //wy_reca = wy_ref + u.y;
		 //wz_reca = wz_ref + u.z;
		 
		 //wx_reca = wx_ref + data[wx_ref-xmin][wy_ref-ymin][wz_ref-zmin].x;
		 //wy_reca = wy_ref + data[wx_ref-xmin][wy_ref-ymin][wz_ref-zmin].y;
		 //wz_reca = wz_ref + data[wx_ref-xmin][wy_ref-ymin][wz_ref-zmin].z;
		 wx_reca = wx_ref + data[wx_ref][wy_ref][wz_ref].x;
		 wy_reca = wy_ref + data[wx_ref][wy_ref][wz_ref].y;
		 wz_reca = wz_ref + data[wx_ref][wy_ref][wz_ref].z;
		 

		 
		 
		 //parcours des patches
		 for(patchx=-hwp;patchx<hwp+1;patchx++)
		   for(patchy=-hwp;patchy<hwp+1;patchy++)
		     for(patchz=-hwp;patchz<hwp+1;patchz++)
		       {
			 xpx = xref + patchx;
			 ypy = yref + patchy;
			 zpz = zref + patchz;
			 wxpx = wx_ref + patchx;
			 wypy = wy_ref + patchy;
			 wzpz = wz_ref + patchz;
			 			 

			 if( (xpx-xmin>=0) & (xpx<width) & (ypy-ymin>=0) & (ypy<height) & (zpz-zmin>=0) & (zpz<depth) )
			   if( (wxpx-xmin>=0) & (wxpx<width) & (wypy-ymin>=0) & (wypy<height) & (wzpz-zmin>=0) & (wzpz<depth) )
			     {
			       //sans deformation du patch :
			       
				// x2px = xreca + patchx;
				// y2py = yreca + patchy;
				// z2pz = zreca + patchz;
				// wx2px = wx_reca + patchx;
				// wy2py = wy_reca + patchy;
				// wz2pz = wz_reca + patchz;
			       
			       
			       //avec deformation du patch :
			       
			       //eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, xpx, ypy, zpz, &u);
			       //x2px = xpx + u.x;
			       //y2py = ypy + u.y;
			       //z2pz = zpz + u.z;
			       //eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, wxpx, wypy, wzpz, &u);
			       //wx2px = wxpx + u.x;
			       //wy2py = wypy + u.y;
			       //wz2pz = wzpz + u.z;
			       
                               
			       //x2px = xpx + data[xpx-xmin][ypy-ymin][zpz-zmin].x;
			       //y2py = ypy + data[xpx-xmin][ypy-ymin][zpz-zmin].y;
			       //z2pz = zpz + data[xpx-xmin][ypy-ymin][zpz-zmin].z;
			       //wx2px = wxpx + data[wxpx-xmin][wypy-ymin][wzpz-zmin].x;
			       //wy2py = wypy + data[wxpx-xmin][wypy-ymin][wzpz-zmin].y;
			       //wz2pz = wzpz + data[wxpx-xmin][wypy-ymin][wzpz-zmin].z;
                               
                               x2px = xpx + data[xpx][ypy][zpz].x;
			       y2py = ypy + data[xpx][ypy][zpz].y;
			       z2pz = zpz + data[xpx][ypy][zpz].z;
			       wx2px = wxpx + data[wxpx][wypy][wzpz].x;
			       wy2py = wypy + data[wxpx][wypy][wzpz].y;
			       wz2pz = wzpz + data[wxpx][wypy][wzpz].z;

			 
         			 //NN
                               if(interpolation_choice==0)
                               {
                                 x2px = (int)x2px;
                                 y2py = (int)y2py;
                                 z2pz = (int)z2pz;
                                 wx2px = (int)wx2px;
                                 wy2py = (int)wy2py;
                                 wz2pz = (int)wz2pz;
			       
			       if( (x2px>=0) & (x2px<width) & (y2py>=0) & (y2py<height) & (z2pz>=0) & (z2pz<depth) )
				 if( (wx2px>=0) & (wx2px<width) & (wy2py>=0) & (wy2py<height) & (wz2pz>=0) & (wz2pz<depth) )
				   {
				     tmpd_ref = imref->mri[xpx][ypy][zpz] - imref->mri[wxpx][wypy][wzpz] ;
				     dist_ref += tmpd_ref*tmpd_ref;			    
				     tmpd_reca = imreca->mri[x2px][y2py][z2pz] - imreca->mri[wx2px][wy2py][wz2pz] ;   
				     dist_reca += tmpd_reca*tmpd_reca;		
				   }
			       }

		        	 // trilineaire
                         if(interpolation_choice==1)
			 {
			 if( (x2px>=0) & (x2px<width) & (y2py>=0) & (y2py<height) & (z2pz>=0) & (z2pz<depth) )
				 if( (wx2px>=0) & (wx2px<width) & (wy2py>=0) & (wy2py<height) & (wz2pz>=0) & (wz2pz<depth) )
			    {
			    
			    tmpd_ref = imref->mri[xpx][ypy][zpz] - imref->mri[wxpx][wypy][wzpz] ;
			    dist_ref += tmpd_ref*tmpd_ref;
			    //necessite interpolation
			    tmpd_reca = image_trilinear_value(imreca, x2px, y2py, z2pz) - image_trilinear_value(imreca, wx2px,wy2py, wz2pz);		      
			    dist_reca += tmpd_reca*tmpd_reca;		
			    }
			 }

                              }



		
		       } //fin de parcours des patchs

		 //TRUCS pour accelerer
		 //remplacer trilineaire par nearest neighbour (les indices sont deja des entiers !)
		 //prendre le patch autour du point central (en nearest)
		 //faire l'interpolation en dehors de la boucle 
		 //calcul du champ de déplacement en dehors de la boucle
		 
		 poids_ref = exp(-dist_ref/PATCH_REF);
		 sum_ref += poids_ref;
		 poids_reca = exp(-dist_reca/PATCH_RECA);
		 sum_reca += poids_reca;
		 
		 wref[count] = poids_ref;
		 wreca[count] = poids_reca;
		 count = count+1;
	       }
	   }

     if( (sum_ref > 0) && (sum_reca > 0) )
       {
	 for(indice=0;indice<nbp;indice++)
	   {
	     wref[indice] = wref[indice] / sum_ref;
	     wreca[indice] = wreca[indice] / sum_reca;

	     auxE = wref[indice] - wreca[indice];
	     E = E + auxE*auxE;
	   }
       }
     }
     //critere quadratique         
     /*
     xref = i+bx0;
     yref = j+by0;
     zref = k+bz0;     
     eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, xref, yref, zref, &u);     
     xreca = xref + u.x;
     yreca = yref + u.y;
     zreca = zref + u.z;
     */

     /*
     i2=(int)xreca;
     j2=(int)yreca;
     k2=(int)zreca;
     i3=i2+1;
     j3=j2+1;
     k3=k2+1;
    
    if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
      {
        xrel=xreca-(double)i2;
	yrel=yreca-(double)j2;
	zrel=zreca-(double)k2;
	
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

	// Caster ou pas ? :)
	//auxE=(int)(auxE+0.5);
        
	auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
	E=E+auxE*auxE;
				
      }
      else
      {
       if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       {
	
        xrel=xreca-(double)i2;yrel=yreca-(double)j2;zrel=zreca-(double)k2;
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
	//auxE=(int)(auxE+0.5);
       
	auxE=auxE-imref->mri[i+bx0][j+by0][k+bz0];
	E=E+auxE*auxE;  
				   
       }
       else
       {
	      //auxE=0-imref->mri[i+bx0][j+by0][k+bz0];
	      //E=E+auxE*auxE; 
				E=HUGE_VAL;
       }
      }
     */
      
    }
   


if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
{
Ereg=regularisation(nb_param, param, param_norm, topi,  topj,  topk,NULL,Slpqr);
//printf("E/sigma2 = %f   Ereg = %f   E/sigma2/Ereg= %f \n",E/TOPO_ALPHA_ROBUST,Ereg,E/TOPO_ALPHA_ROBUST/Ereg);

E+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*Ereg;
}


E=sqrt(E);

free_field3d(champ);

return(E);

}
double image_trilinear_value(grphic3d *im, double x, double y, double z)
{

  double res = 0;
  int i2=(int)x;
  int j2=(int)y;
  int k2=(int)z;
  int i3=i2+1;
  int j3=j2+1;
  int k3=k2+1;
  int width=im->width;
  int height=im->height;
  int depth=im->depth;
  double xrel, yrel, zrel;
  double va,vb,vc,vd,ve,vf,vg,vh;

  if (i2>=0 && i3<width && j2>=0 && j3<height && k2>=0 && k3<depth)
    {
      xrel=x-(double)i2;
      yrel=y-(double)j2;
      zrel=z-(double)k2;
      
      va=(double)im->mri[i2][j2][k2];
      vb=(double)im->mri[i2][j2][k3];
      vc=(double)im->mri[i2][j3][k2];
      vd=(double)im->mri[i2][j3][k3];
      ve=(double)im->mri[i3][j2][k2];
      vf=(double)im->mri[i3][j2][k3];
      vg=(double)im->mri[i3][j3][k2];
      vh=(double)im->mri[i3][j3][k3];
      
      res=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
	+(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
	+(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;				
    }
  else
    {
      if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
	{
	  
	  xrel=x-(double)i2;
	  yrel=y-(double)j2;
	  zrel=z-(double)k2;
	  if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)im->mri[i2][j2][k2];
	  if (i2==-1 || j2==-1 || k3==depth) vb=0.0; else vb=(double)im->mri[i2][j2][k3];
	  if (i2==-1 || j3==height || k2==-1) vc=0.0; else vc=(double)im->mri[i2][j3][k2];
	  if (i2==-1 || j3==height || k3==depth) vd=0.0; else vd=(double)im->mri[i2][j3][k3];
	  if (i3==width || j2==-1 || k2==-1) ve=0.0; else ve=(double)im->mri[i3][j2][k2];
	  if (i3==width || j2==-1 || k3==depth) vf=0.0; else vf=(double)im->mri[i3][j2][k3];
	  if (i3==width || j3==height || k2==-1) vg=0.0; else vg=(double)im->mri[i3][j3][k2];
	  if (i3==width || j3==height || k3==depth) vh=0.0; else vh=(double)im->mri[i3][j3][k3];
	  
	  res=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
	    +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
	    +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
	}
      else
	{
	  res = 0;
	  //E=HUGE_VAL;
	}
    }
  return res;
  
}

/*******************************************************************************
**  			regularisation_patch_imagebased_local
**
**	Calcul local du terme de regularisation avec des patch
*******************************************************************************/

double regularisation_patch_imagebased_local(int nb_param, double *param,double *param_norm, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
double reg;
field3d *champ;
vector3d ***data;
dvector3d u;
int ideb,jdeb,kdeb, ifin,jfin,kfin;
double ***poids, dx,dy,dz;
int nhx,nhy,nhz,spx,spy,spz;

 	
nhx=2*TOPO_NLMhwvsx+1;
nhy=2*TOPO_NLMhwvsy+1;
nhz=2*TOPO_NLMhwvsz+1;

spx=2*TOPO_NLMhwnx+1;
spy=2*TOPO_NLMhwny+1;
spz=2*TOPO_NLMhwnz+1;

poids = alloc_dmatrix_3d(nhx, nhy, nhz);

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



reg=0;
	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	{

ideb=MAXI(0,ux-TOPO_NLMhwnx+bx0);
jdeb=MAXI(0,uy-TOPO_NLMhwny+by0);
kdeb=MAXI(0,uz-TOPO_NLMhwnz+bz0);

ifin=MINI(width,ux+TOPO_NLMhwnx+bx0);
jfin=MINI(height,uy+TOPO_NLMhwny+by0);
kfin=MINI(depth,uz+TOPO_NLMhwnz+bz0);

auxi=ux+bx0;
auxj=uy+by0;
auxk=uz+bz0;

    	
 NLMWeightOnePoint3D(_ParamRecalageBspline.ref_groupwise, poids, TOPO_PATCH_COEFF,  
					TOPO_NLMhwnx,  TOPO_NLMhwny,  TOPO_NLMhwnz,  TOPO_NLMhwvsx,  TOPO_NLMhwvsy,  TOPO_NLMhwvsz,auxi, auxj, auxk);
		
		for (i=ideb; i<ifin; i++)
		for (j=jdeb; j<jfin; j++)
		for (k=kdeb; k<kfin; k++)
			{
			u.x=u.y=u.z=0; 
			eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, i, j, k, &u);   
    		dx=data[ux][uy][uz].x-u.x;
			dx=dx*dx;
			dy=data[ux][uy][uz].y-u.y;
			dy=dy*dy;
			dz=data[ux][uy][uz].z-u.z;
			dz=dz*dz;
						
			reg+=poids[i-bx0-ux+TOPO_NLMhwnx][j-by0-uy+TOPO_NLMhwny][k-bz0-uz+TOPO_NLMhwnz]*(dx+dy+dz);
			
			}
		
		
		 
	}



free_field3d(champ);
free_dmatrix_3d(poids);
reg=reg*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique
return(reg);
}
