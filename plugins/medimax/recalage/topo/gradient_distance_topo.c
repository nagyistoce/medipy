/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>

#include "recalage/mtch_3d.h"
#include "recalage/chps_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/gradient_distance_topo.h"
#include "math/imx_matrix.h"


/*******************************************************************************
**     gradient_base_quad_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_quad_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
		diff=2.0*(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     gradient_base_quad_sym_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int gradient_base_quad_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 	int 	i,j,k,l,wdth,hght,dpth;
 	double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
 	double gx,gy,gz,diff,diff2,temp,val;
 	double dx,dy,dz;
 	int ux,uy,uz;
 	int bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 	int  D=TOP_D(nb_param);
	int x00,x11,y00,y11,z00,z11;
 	vector3d ***data,***dataG;
 	field3d *champ;
 	vector3d va,vb,vc,vd,ve,vf,vg,vh;
 	double wa,wb,wc,wd,we,wf,wg,wh;
 	double f,px,py,pz;
 	int auxi,auxj,auxk;
 	int i2,j2,k2,i3,j3,k3;
 	double xrel,yrel,zrel;
 	double J=1,dJx,dJy,dJz,xt,yt,zt;
	int aux0,aux1;  								// variables auxiliaires
	double auxdbl[20],aux2,mymono[21];
	double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
	int ideb,ifin,jdeb,jfin,kdeb,kfin;
	double dir_desR3[3];
	int signe=1;
	int resol;
	
 	wdth=I2->width;hght=I2->height;dpth=I2->depth;
 	dataG=gradI2->raw;
 	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 	fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
	resol=BASE3D.resol;
	//-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

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

 // mise a zero du gradient
 grad[0]=grad[1]=grad[2]=0.0;

 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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
 
  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

 // on rempli le vecteur gradient

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	   {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire 
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    	dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    	dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		// Evaluation de Imreca en s=s+u avec interrpolation lineaire 	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
        }

			diff=(double)(val-I2->mri[i][j][k]);
			if (diff!=0)
			{
		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
		signe=1;
		if (J<0) signe=-1;
		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		
		dJx=1.0*dJx/wdth;
		dJy=1.0*dJy/hght;
		dJz=1.0*dJz/dpth;
		
		diff2=1.0*diff*diff; 
      temp=2.0*diff*fx[i-bx0]*fy[j-by0]*fz[k-bz0]*(1.0+fabs(J));
      grad[0]=grad[0]+temp*gx+signe*dJx*diff2;
      grad[1]=grad[1]+temp*gy+signe*dJy*diff2;
      grad[2]=grad[2]+temp*gz+signe*dJz*diff2;
			}
	}
} 
  for (i=0;i<3;i++) grad[i]=grad[i]/(double)(topDx*topDy*topDz);
 
   
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

 
  free_field3d(champ);
  return(1);
}


/*******************************************************************************
**     gradient_base_Lp_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_Lp_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,  double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double p=ORDRE_NORME_LP;
 int signe=1;
 
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
        }
				
	if (((double)I2->mri[i+bx0][j+by0][k+bz0]-val)>0)
		signe=-1;
	else
		signe=1;
		
		diff=1.0*signe*p*pow(fabs(val-(double)I2->mri[i+bx0][j+by0][k+bz0]),(double)(p-1));
		  temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  for (i=0;i<3;i++) grad[i]=grad[i]/(double)(topDx*topDy*topDz);
 
   
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

 
  free_field3d(champ);
  return(1);
}

/*******************************************************************************
**     gradient_base_Lp_sym_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int gradient_base_Lp_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 	int 	i,j,k,l,wdth,hght,dpth;
 	double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
 	double gx,gy,gz,diff,diff2,temp,val;
 	double dx,dy,dz;
 	int ux,uy,uz;
 	int bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 	int  D=TOP_D(nb_param);
	int x00,x11,y00,y11,z00,z11;
 	vector3d ***data,***dataG;
 	field3d *champ;
 	vector3d va,vb,vc,vd,ve,vf,vg,vh;
 	double wa,wb,wc,wd,we,wf,wg,wh;
 	double f,px,py,pz;
 	int auxi,auxj,auxk;
 	int i2,j2,k2,i3,j3,k3;
 	double xrel,yrel,zrel;
 	double J=1,dJx,dJy,dJz,xt,yt,zt;
	int aux0,aux1;  								// variables auxiliaires
	double auxdbl[20],aux2,mymono[21];
	double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
	int ideb,ifin,jdeb,jfin,kdeb,kfin;
	double dir_desR3[3];
	int signeJ=1,signe=1;
	int resol;
	double p=ORDRE_NORME_LP;
	
	
 	wdth=I2->width;hght=I2->height;dpth=I2->depth;
 	dataG=gradI2->raw;
 	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 	fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
	resol=BASE3D.resol;
	//-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

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

 // mise a zero du gradient
 grad[0]=grad[1]=grad[2]=0.0;

 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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
 
  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

 // on rempli le vecteur gradient

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	   {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire 
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    	dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    	dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		// Evaluation de Imreca en s=s+u avec interrpolation lineaire 	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
        }

		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
		signeJ=1;
		if (J<0) signeJ=-1;
		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		dJx=1.0*dJx/wdth;
		dJy=1.0*dJy/hght;
		dJz=1.0*dJz/dpth;
		
		if (((double)I2->mri[i][j][k]-val)>0)
		signe=-1;
	else
		signe=1;
		diff=fabs(val-(double)I2->mri[i][j][k]);
		diff2=pow(diff,(double)p);
		diff=1.0*signe*p*pow(diff,(double)(p-1));
		 
		  temp=diff*fx[i-bx0]*fy[j-by0]*fz[k-bz0]*(1.0+fabs(J));
      grad[0]=grad[0]+temp*gx+signeJ*dJx*diff2;
      grad[1]=grad[1]+temp*gy+signeJ*dJy*diff2;
      grad[2]=grad[2]+temp*gz+signeJ*dJz*diff2;
	}
} 
  for (i=0;i<3;i++) grad[i]=grad[i]/(double)(topDx*topDy*topDz);
 
   
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

  free_field3d(champ);
  return(1);
}

/*******************************************************************************
**     gradient_base_L1L2_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_L1L2_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
			diff=(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
      temp=diff*fx[i]*fy[j]*fz[k];
			temp=temp*pow((diff*diff)+TOPO_EPSILON,-0.5);
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  for (i=0;i<3;i++)
	grad[i]=10000*grad[i]/(double)(topDx*topDy*topDz); /*le facteur 10000 c'est de la bidouille pour que Marquardt marche avec les meme param que Quadratique ... On aurait du normaliser la valeur des gradients et regler les param dans LM en consequent ...*/ 

  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}


  free_field3d(champ);
  return(1);
}

/*******************************************************************************
**     gradient_base_L1L2_sym_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int gradient_base_L1L2_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 	int 	i,j,k,l,wdth,hght,dpth;
 	double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
 	double gx,gy,gz,diff,diff2,temp,val;
 	double dx,dy,dz;
 	int ux,uy,uz;
 	int bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 	int  D=TOP_D(nb_param);
	int x00,x11,y00,y11,z00,z11;
 	vector3d ***data,***dataG;
 	field3d *champ;
 	vector3d va,vb,vc,vd,ve,vf,vg,vh;
 	double wa,wb,wc,wd,we,wf,wg,wh;
 	double f,px,py,pz;
 	int auxi,auxj,auxk;
 	int i2,j2,k2,i3,j3,k3;
 	double xrel,yrel,zrel;
 	double J=1,dJx,dJy,dJz,xt,yt,zt;
	int aux0,aux1;  								// variables auxiliaires
	double auxdbl[20],aux2,mymono[21];
	double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
	int ideb,ifin,jdeb,jfin,kdeb,kfin;
	double dir_desR3[3];
	int signeJ=1;
	int resol;
	
 	wdth=I2->width;hght=I2->height;dpth=I2->depth;
 	dataG=gradI2->raw;
 	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 	fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
	resol=BASE3D.resol;
	//-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

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

 // mise a zero du gradient
 grad[0]=grad[1]=grad[2]=0.0;

 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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
 
  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

 // on rempli le vecteur gradient

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	   {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire 
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i-bx0][j-by0][k-bz0].x;
    	dy=(double)j+data[i-bx0][j-by0][k-bz0].y;
    	dz=(double)k+data[i-bx0][j-by0][k-bz0].z;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		// Evaluation de Imreca en s=s+u avec interrpolation lineaire 	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
        }

		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
		signeJ=1;
		if (J<0) signeJ=-1;
		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		dJx=1.0*dJx/wdth;
		dJy=1.0*dJy/hght;
		dJz=1.0*dJz/dpth;
		
		diff=(double)(val-I2->mri[i][j][k]);
    temp=diff*fx[i]*fy[j]*fz[k];
		temp=temp*pow((diff*diff)+TOPO_EPSILON,-0.5);
    diff2=sqrt(diff*diff+TOPO_EPSILON);
	
		temp=diff*fx[i-bx0]*fy[j-by0]*fz[k-bz0]*(1.0+fabs(J));
    grad[0]=grad[0]+temp*gx+signeJ*dJx*diff2;
    grad[1]=grad[1]+temp*gy+signeJ*dJy*diff2;
    grad[2]=grad[2]+temp*gz+signeJ*dJz*diff2;
	}
} 
  for (i=0;i<3;i++) grad[i]=10000*grad[i]/(double)(topDx*topDy*topDz);
 
   
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

  free_field3d(champ);
  return(1);
}

/*******************************************************************************
**     gradient_base_L1norm_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_L1norm_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk,signe=1;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
			diff=(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
			
			if (diff<0)
				signe = -1;
			else
				signe = 1;
				
			if ((I2->mri[i+bx0][j+by0][k+bz0]+val)!=0)   
		     diff=signe*1.0*I2->mri[i+bx0][j+by0][k+bz0]/((val+I2->mri[i+bx0][j+by0][k+bz0])*(val+I2->mri[i+bx0][j+by0][k+bz0]));
			else
				diff=0;
				
			temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);

  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}


  free_field3d(champ);
  return(1);
}


/*******************************************************************************
**     gradient_base_geman_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_geman_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val,tutu;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
		diff=(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
     tutu=(TOPO_ALPHA_ROBUST+diff*diff);
		 diff=2.0*diff*TOPO_ALPHA_ROBUST/(tutu*tutu);
			temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]*TOPO_ALPHA_ROBUST/(double)(topDx*topDy*topDz);
 
   
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

 
  free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     gradient_regularisation_energie_membrane_local                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_regularisation_energie_membrane_local(int nb_param, double *param, double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
int uxdeb,uydeb,uzdeb;
double Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
double dJxx,dJxy,dJxz,dJyx,dJyy,dJyz,dJzx,dJzy,dJzz;
field3d *champ;
vector3d ***data;
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

l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

grad[0]=0;grad[1]=0;grad[2]=0;


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
				
		 
		 //Jxx=1+	1.0*(data[uxdeb+1][uy][uz].x-data[uxdeb][uy][uz].x);
		 Jxx=		1.0*(data[uxdeb+1][uy][uz].x-data[uxdeb][uy][uz].x);
		 Jxy=		1.0*(data[ux][uydeb+1][uz].x-data[ux][uydeb][uz].x);
		 Jxz=		1.0*(data[ux][uy][uzdeb+1].x-data[ux][uy][uzdeb].x);
		 Jyx=		1.0*(data[uxdeb+1][uy][uz].y-data[uxdeb][uy][uz].y);
		 //Jyy=1+	1.0*(data[ux][uydeb+1][uz].y-data[ux][uydeb][uz].y);
		 Jyy=		1.0*(data[ux][uydeb+1][uz].y-data[ux][uydeb][uz].y);
		 Jyz=		1.0*(data[ux][uy][uzdeb+1].y-data[ux][uy][uzdeb].y);
		 Jzx=		1.0*(data[uxdeb+1][uy][uz].z-data[uxdeb][uy][uz].z);
		 Jzy=		1.0*(data[ux][uydeb+1][uz].z-data[ux][uydeb][uz].z);
		 //Jzz=1+	1.0*(data[ux][uy][uzdeb+1].z-data[ux][uy][uzdeb].z);
		 Jzz=	1.0*(data[ux][uy][uzdeb+1].z-data[ux][uy][uzdeb].z);
		 
		 
		 if (px!=0)
		 {
		 dJxx=dfx[ux]*fy[uy]*fz[uz];
		 dJxy=fx[ux]*dfy[uy]*fz[uz];
		 dJxz=fx[ux]*fy[uy]*dfz[uz];
		 }
		 else
		 {dJxx=dJxy=dJxz=0.0;}
		 
		 
		if (py!=0)
		 {
		 dJyx=dfx[ux]*fy[uy]*fz[uz];
		 dJyy=fx[ux]*dfy[uy]*fz[uz];
		 dJyz=fx[ux]*fy[uy]*dfz[uz];
		 }
		else
		 {dJyx=dJyy=dJyz=0.0;}
		 
		if (pz!=0)
		 {
		 dJzx=dfx[ux]*fy[uy]*fz[uz];
		 dJzy=fx[ux]*dfy[uy]*fz[uz];
		 dJzz=fx[ux]*fy[uy]*dfz[uz];
		 }
		else
		 {dJzx=dJzy=dJzz=0.0;}
		 
		 
		
		 grad[0]+=2.0*(dJxx*Jxx+dJxy*Jxy+dJxz*Jxz);
		 grad[1]+=2.0*(dJyx*Jyx+dJyy*Jyy+dJyz*Jyz);
		 grad[2]+=2.0*(dJzx*Jzx+dJzy*Jzy+dJzz*Jzz);
		 
		 }
		 }

for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


free_field3d(champ);
return(1);
}


/*******************************************************************************
**     gradient_regularisation_energie_membrane_Lp_local                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_regularisation_energie_membrane_Lp_local(int nb_param, double *param,double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
int uxdeb,uydeb,uzdeb;
double Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
double dJxx,dJxy,dJxz,dJyx,dJyy,dJyz,dJzx,dJzy,dJzz;
field3d *champ;
vector3d ***data;
int ideb,jdeb,kdeb,condition;
double ORDRE_REGULARISATION_LP_moins_1=ORDRE_REGULARISATION_LP-1.0,tmp;

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

l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

grad[0]=0;grad[1]=0;grad[2]=0;


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
		 
		 
		 dJxx=dfx[ux]*fy[uy]*fz[uz];
		 dJxy=fx[ux]*dfy[uy]*fz[uz];
		 dJxz=fx[ux]*fy[uy]*dfz[uz];
		 
		 dJyx=dfx[ux]*fy[uy]*fz[uz];
		 dJyy=fx[ux]*dfy[uy]*fz[uz];
		 dJyz=fx[ux]*fy[uy]*dfz[uz];
		
		 dJzx=dfx[ux]*fy[uy]*fz[uz];
		 dJzy=fx[ux]*dfy[uy]*fz[uz];
		 dJzz=fx[ux]*fy[uy]*dfz[uz];
		
		// grad[0]+=2.0*(dJxx*Jxx+dJxy*Jxy+dJxz*Jxz);
		// grad[1]+=2.0*(dJyx*Jyx+dJyy*Jyy+dJyz*Jyz);
		// grad[2]+=2.0*(dJzx*Jzx+dJzy*Jzy+dJzz*Jzz);
		 
		 tmp=dJxx*pow(Jxx,ORDRE_REGULARISATION_LP_moins_1);
		 if (isnan(tmp)) {tmp=0;}
		 grad[0]+=tmp;
	 	 
		  tmp=dJxy*pow(Jxy,ORDRE_REGULARISATION_LP_moins_1);
		 if (isnan(tmp)) {tmp=0;}
		 grad[0]+=tmp;
	 	 
		  tmp=dJxz*pow(Jxz,ORDRE_REGULARISATION_LP_moins_1);
		 if (isnan(tmp)) {tmp=0;}
		 grad[0]+=tmp;
	 	 
		 
		  tmp=dJyx*pow(Jyx,ORDRE_REGULARISATION_LP_moins_1);
			if (isnan(tmp)) {tmp=0;}
			grad[1]+=tmp;
		  
			tmp=dJyy*pow(Jyy,ORDRE_REGULARISATION_LP_moins_1);
			if (isnan(tmp)) {tmp=0;}
			grad[1]+=tmp;
		 
		 tmp=dJyz*pow(Jyz,ORDRE_REGULARISATION_LP_moins_1);
			if (isnan(tmp)) {tmp=0;}
			grad[1]+=tmp;
		 
			
			
			tmp=dJzx*pow(Jzx,ORDRE_REGULARISATION_LP_moins_1);
		 if (isnan(tmp)) {tmp=0;}
		 grad[2]+=tmp;
		 
		 	tmp=dJzy*pow(Jzy,ORDRE_REGULARISATION_LP_moins_1);
		 if (isnan(tmp)) {tmp=0;}
		 grad[2]+=tmp;
		 
		 	tmp=dJzz*pow(Jzz,ORDRE_REGULARISATION_LP_moins_1);
		 if (isnan(tmp)) {tmp=0;}
		 grad[2]+=tmp;
		 }
		 }

for (i=0;i<3;i++)
	grad[i]=1.0*ORDRE_REGULARISATION_LP*grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,ORDRE_REGULARISATION_LP)); // Pour assurer la compatibilite avec le recalage symetrique;


free_field3d(champ);
return(1);
}


/*******************************************************************************
**     gradient_regularisation_log_jacobien_local                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_regularisation_log_jacobien_local(int nb_param, double *param,double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
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
double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
double dir_desR3[3];
double J=1,dJx,dJy,dJz,xt,yt,zt,tmp;

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
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

grad[0]=0;grad[1]=0;grad[2]=0;

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	 {
	xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
	fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
	
	fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
	fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
	fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		
	dJx=1.0*dJx/wdth;
	dJy=1.0*dJy/hght;
	dJz=1.0*dJz/dpth;

    	if (J>0)
	tmp= 2.*log(J)/J;
	else 
	tmp=0;
	
	grad[0]=grad[0]+tmp*dJx;
      	grad[1]=grad[1]+tmp*dJy;
      	grad[2]=grad[2]+tmp*dJz;

	
	 }
 }

for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


return(1);
}

/*******************************************************************************
**     gradient_regularisation_log_jacobien_centre_local                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_regularisation_log_jacobien_centre_local(int nb_param, double *param,double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
{
int resol,wdth,hght,dpth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int i,j,k,auxi,auxj,auxk;
int ideb,ifin,jdeb,jfin,kdeb,kfin;	
int  D=TOP_D(nb_param);
int aux0,aux1;  								// variables auxiliaires
double auxdbl[20],aux2,mymono[21];
double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
double dir_desR3[3];
double J=1,dJx,dJy,dJz,xt,yt,zt,tmp=0,Jmean, dJxmean, dJymean, dJzmean, Jtmp, dJxtmp, dJytmp, dJztmp;

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
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

grad[0]=0;grad[1]=0;grad[2]=0;

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	 {
	xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
	fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
	
	fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
	fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
	fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		
	dJx=1.0*dJx/wdth;
	dJy=1.0*dJy/hght;
	dJz=1.0*dJz/dpth;


	Jmean=0; dJxmean=0; dJymean=0; dJzmean=0;
	
	for (auxi=i-1;auxi<i+2;auxi++)
	for (auxj=j-1;auxj<j+2;auxj++)
	for (auxk=k-1;auxk<k+2;auxk++)
		{
		xt=(double)(auxi)/wdth;yt=(double)(auxj)/hght;zt=(double)(auxk)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &Jtmp, mymono);
		
		fast_eval_fun(pente_x, xt, yt, zt, &dJxtmp, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJytmp, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJztmp, mymono);
	
		dJxtmp=1.0*dJxtmp/wdth;
		dJytmp=1.0*dJytmp/hght;
		dJztmp=1.0*dJztmp/dpth;

		
		Jmean+=log(fabs(Jtmp));

		if (Jtmp>0)
			{
			dJxmean+=dJxtmp/Jtmp;
			dJymean+=dJytmp/Jtmp;
			dJzmean+=dJztmp/Jtmp;
			}
		
		}

	Jmean=Jmean/27.0;
	
	dJxmean=dJxmean/27.0;
	dJymean=dJymean/27.0;
	dJzmean=dJzmean/27.0;
	
	
	
    	
	if (J>0)
	{
	tmp= 2.*(log(J)-Jmean);
	grad[0]=grad[0]+tmp*(dJx/J-dJxmean);
      	grad[1]=grad[1]+tmp*(dJy/J-dJymean);
      	grad[2]=grad[2]+tmp*(dJz/J-dJzmean);
	}
	
	 }
 }

for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


return(1);
}



/*******************************************************************************
**     gradient_regularisation_energie_membrane_jacobien_local                  
**                                                                   
*******************************************************************************/
int gradient_regularisation_energie_membrane_jacobien_local(int nb_param, double *param, double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,f,*dfx,*dfy,*dfz, dfdx, dfdy,dfdz,d2fdxdy, d2fdxdz, d2fdydz, px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
double Jacobien,Jx,Jy,Jz, Jax, Jay, Jaz, Jxax, Jxay, Jxaz, Jyax, Jyay, Jyaz,Jzax, Jzay, Jzaz,J3;
field3d *champ,*champdx, *champdy, *champdz,*champdxy,*champdxz,*champdyz;
mat3d mat;

vector3d ***data,***datadx,***datady,***datadz,***datadxy,***datadxz,***datadyz;
int ideb,jdeb,kdeb,condition,tmp;


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


grad[0]=grad[1]=grad[2]=0.0;

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
		
		
		// calcul de la derivee du jacobien suivant ax
		
		mat.H[0][0]= 	 	dfx[ux]*fy[uy]*fz[uz];
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	fx[ux]*fy[uy]*dfz[uz];
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		
		Jax=determinant_mat3d(&mat);
		
		
		// calcul de la derivee du jacobien suivant ay
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		fx[ux]*fy[uy]*dfz[uz];
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
	
		// calcul de la derivee du jacobien suivant az
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*fy[uy]*fz[uz];
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		fx[ux]*dfy[uy]*fz[uz];
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 		fx[ux]*fy[uy]*dfz[uz];
	
	
			
		// calcul de Jxax
		
		mat.H[0][0]= 	 	dfx[ux]*fy[uy]*fz[uz];
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		dfx[ux]*dfy[uy]*fz[uz];
		mat.H[1][1]= 		datadxy[ux][uy][uz].y;
		mat.H[2][1]= 		datadxy[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	fx[ux]*fy[uy]*dfz[uz];
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		
		Jxax=determinant_mat3d(&mat);
		
		mat.H[0][0]= 	 	dfx[ux]*fy[uy]*fz[uz];
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	dfx[ux]*fy[uy]*dfz[uz];
		mat.H[1][2]= 		datadxz[ux][uy][uz].y;
		mat.H[2][2]= 		datadxz[ux][uy][uz].z;
		
		Jxax=Jxax+determinant_mat3d(&mat);
		
		// calcul de Jxay
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datadxy[ux][uy][uz].x;
		mat.H[1][1]= 		dfx[ux]*dfy[uy]*fz[uz];
		mat.H[2][1]= 		datadxy[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		fx[ux]*fy[uy]*dfz[uz];
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		
		Jxay=determinant_mat3d(&mat);
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadxz[ux][uy][uz].x;
		mat.H[1][2]= 		dfx[ux]*fy[uy]*dfz[uz];
		mat.H[2][2]= 		datadxz[ux][uy][uz].z;
		
		Jxay=Jxay+determinant_mat3d(&mat);
	
		// Jxaz
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*fy[uy]*fz[uz];
			
		mat.H[0][1]= 		datadxy[ux][uy][uz].x;
		mat.H[1][1]= 		datadxy[ux][uy][uz].y;
		mat.H[2][1]= 		dfx[ux]*dfy[uy]*fz[uz];
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 		fx[ux]*fy[uy]*dfz[uz];
		
		
		Jxaz=determinant_mat3d(&mat);
		
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*fy[uy]*fz[uz];
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		fx[ux]*dfy[uy]*fz[uz];
		
		
		mat.H[0][2]= 	 	datadxz[ux][uy][uz].x;
		mat.H[1][2]= 		datadxz[ux][uy][uz].y;
		mat.H[2][2]= 		dfx[ux]*fy[uy]*dfz[uz];
		
		Jxaz=Jxaz+determinant_mat3d(&mat);
		
		
		// Jyax
		
		mat.H[0][0]= 		dfx[ux]*dfy[uy]*fz[uz];
		mat.H[1][0]= 		datadxy[ux][uy][uz].y;
		mat.H[2][0]= 		datadxy[ux][uy][uz].z;
			
		mat.H[0][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	fx[ux]*fy[uy]*dfz[uz];
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		Jyax=determinant_mat3d(&mat);
	
		mat.H[0][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	fx[ux]*dfy[uy]*dfz[uz];
		mat.H[1][2]= 		datadyz[ux][uy][uz].y;
		mat.H[2][2]= 		datadyz[ux][uy][uz].z;
		
		Jyax=Jyax+determinant_mat3d(&mat);
		
		// Jyay
		
		mat.H[0][0]= 		datadxy[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*dfy[uy]*fz[uz];
		mat.H[2][0]= 		datadxy[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		fx[ux]*fy[uy]*dfz[uz];
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		Jyay=determinant_mat3d(&mat);
	
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadyz[ux][uy][uz].x;
		mat.H[1][2]= 		fx[ux]*dfy[uy]*dfz[uz];
		mat.H[2][2]= 		datadyz[ux][uy][uz].z;
		
		Jyay=Jyay+determinant_mat3d(&mat);
		
		
		// Jyaz
		
		mat.H[0][0]= 		datadxy[ux][uy][uz].x;
		mat.H[1][0]= 		datadxy[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*dfy[uy]*fz[uz];
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		fx[ux]*dfy[uy]*fz[uz];
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 		fx[ux]*fy[uy]*dfz[uz];
		
		Jyaz=determinant_mat3d(&mat);
	
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*fy[uy]*fz[uz];
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		fx[ux]*dfy[uy]*fz[uz];
		
		
		mat.H[0][2]= 	 	datadyz[ux][uy][uz].x;
		mat.H[1][2]= 		datadyz[ux][uy][uz].y;
		mat.H[2][2]= 		fx[ux]*dfy[uy]*dfz[uz];
		
		Jyaz=Jyaz+determinant_mat3d(&mat);
		
		
		
		// Jzax
		
		mat.H[0][0]= 	 	dfx[ux]*fy[uy]*dfz[uz];
		mat.H[1][0]= 		datadxz[ux][uy][uz].y;
		mat.H[2][0]= 		datadxz[ux][uy][uz].z;
			
		mat.H[0][1]= 		fx[ux]*dfy[uy]*fz[uz];
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	fx[ux]*fy[uy]*dfz[uz];
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		Jzax=determinant_mat3d(&mat);
	
		mat.H[0][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		fx[ux]*dfy[uy]*dfz[uz];
		mat.H[1][1]= 		datadyz[ux][uy][uz].y;
		mat.H[2][1]= 		datadyz[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	fx[ux]*fy[uy]*dfz[uz];
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;	 
		 
		Jzax=Jzax+determinant_mat3d(&mat);
		
		
		
		// Jzay
		
		mat.H[0][0]= 	 	datadxz[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*fy[uy]*dfz[uz];
		mat.H[2][0]= 		datadxz[ux][uy][uz].z;
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]=		fx[ux]*dfy[uy]*fz[uz];
		mat.H[2][1]= 		datady[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		fx[ux]*fy[uy]*dfz[uz];
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;
		
		Jzay=determinant_mat3d(&mat);
	
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		dfx[ux]*fy[uy]*fz[uz];
		mat.H[2][0]= 		datadx[ux][uy][uz].z;
			
		mat.H[0][1]= 		datadyz[ux][uy][uz].x;
		mat.H[1][1]= 		fx[ux]*dfy[uy]*dfz[uz];
		mat.H[2][1]= 		datadyz[ux][uy][uz].z;
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		fx[ux]*fy[uy]*dfz[uz];
		mat.H[2][2]= 1.0 +	datadz[ux][uy][uz].z;	 
		 
		Jzay=Jzay+determinant_mat3d(&mat);
		
		
		// Jzaz
		
		mat.H[0][0]= 	 	datadxz[ux][uy][uz].x;
		mat.H[1][0]= 		datadxz[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*fy[uy]*dfz[uz];
			
		mat.H[0][1]= 		datady[ux][uy][uz].x;
		mat.H[1][1]= 1.0 +	datady[ux][uy][uz].y;
		mat.H[2][1]= 		fx[ux]*dfy[uy]*fz[uz];
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 		fx[ux]*fy[uy]*dfz[uz];
		
		Jzaz=determinant_mat3d(&mat);
	
		mat.H[0][0]= 1.0 + 	datadx[ux][uy][uz].x;
		mat.H[1][0]= 		datadx[ux][uy][uz].y;
		mat.H[2][0]= 		dfx[ux]*fy[uy]*fz[uz];
			
		mat.H[0][1]= 		datadyz[ux][uy][uz].x;
		mat.H[1][1]= 		datadyz[ux][uy][uz].y;
		mat.H[2][1]= 		fx[ux]*dfy[uy]*dfz[uz];
		
		
		mat.H[0][2]= 	 	datadz[ux][uy][uz].x;
		mat.H[1][2]= 		datadz[ux][uy][uz].y;
		mat.H[2][2]= 		fx[ux]*fy[uy]*dfz[uz];	 
		 
		Jzaz=Jzaz+determinant_mat3d(&mat);
		
		
		
		
		
		if (Jacobien>0)
		{
		tmp=-2.0*(Jx*Jx+Jy*Jy+Jz*Jz);
		J3=Jacobien*Jacobien*Jacobien;
		grad[0]=grad[0]+(2.0*Jacobien*(Jx*Jxax+Jy*Jyax+Jz*Jzax)+tmp*Jax)/J3;
		grad[1]=grad[1]+(2.0*Jacobien*(Jx*Jxay+Jy*Jyay+Jz*Jzay)+tmp*Jay)/J3;
		grad[2]=grad[2]+(2.0*Jacobien*(Jx*Jxaz+Jy*Jyaz+Jz*Jzaz)+tmp*Jaz)/J3;
		
		if (isnan(grad[0]))
			printf("y a un bug grad[0]\n");
		
		if (isnan(grad[1]))
			printf("y a un bug grad[1]\n");
		
		
		if (isnan(grad[2]))
			printf("y a un bug grad[2]\n");
		
		}
		
		
		}

	}










for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


free_field3d(champ);
free_field3d(champdx);
free_field3d(champdy);
free_field3d(champdz);
free_field3d(champdxy);
free_field3d(champdxz);
free_field3d(champdyz);

return(1);
}

/*******************************************************************************
**     gradient_regularisation_dist_identite_local                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_regularisation_dist_identite_local(int nb_param, double *param, double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
field3d *champ;
vector3d ***data;
int ideb,jdeb,kdeb;
double dJxx, dJyy, dJzz;

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

l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

grad[0]=0;grad[1]=0;grad[2]=0;


	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	   {
    	
			
		
		 
		 
		 if (px!=0)
		 {
		 dJxx=dfx[ux]*fy[uy]*fz[uz]*data[ux][uy][uz].x;
		 }
		 else
		 {dJxx=0.0;}
		 
		
		 if (py!=0)
		 {
		 dJyy=fx[ux]*dfy[uy]*fz[uz]*data[ux][uy][uz].y;
		 }
		 else
		 {dJyy=0.0;}
		
		 if (pz!=0)
		 {
		 dJzz=fx[ux]*fy[uy]*dfz[uz]*data[ux][uy][uz].z;
		 }
		 else
		 {dJzz=0.0;}
		 
		 
		
		 grad[0]+=2.0*dJxx;
		 grad[1]+=2.0*dJyy;
		 grad[2]+=2.0*dJzz;
		 		 
		 }

for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


free_field3d(champ);
return(1);
}


/*******************************************************************************
**     gradient_base_quad_locale_pondere_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_quad_locale_pondere_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk,Ntot;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel,Ptot,poids;

Ptot=0;Ntot=0;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			poids=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
	  
		poids=I2->mask->mri[i+bx0][j+by0][k+bz0];
		
		Ptot+=poids;	
		Ntot++;	
        }
				
		diff=2.0*poids*(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  if (Ptot>0)
		{
		for (i=0;i<3;i++)
		grad[i]=grad[i]*Ntot/(double)(topDx*topDy*topDz*Ptot);
  	}
	else
		for (i=0;i<3;i++)
		grad[i]=0;
  	
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
						
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

	
	free_field3d(champ);


  return(1);
}



/*******************************************************************************
**     gradient_base_IM_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_IM_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth,u,v;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,temp,val;
 double dx,dy,dz,x,y;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ,*dpij;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3,i0,j0;
 double xrel,yrel,zrel;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 

dpij=cr_field3d(HISTOJOINT.nbI,HISTOJOINT.nbJ,1);

for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
 	dpij->raw[i][j][0].x=dpij->raw[i][j][0].y=dpij->raw[i][j][0].z=0.0;
	
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
	  x=1.0*val*HISTOJOINT.normI;
		y=1.0*I2->mri[i+bx0][j+by0][k+bz0]*HISTOJOINT.normJ;
		i0=(int)floor(x);
		j0=(int)floor(y);


// calcul de dpij grace � la formule (24) de l'article Thevenaz IEEETMI dec 2000

					if ((x>0)&&(y>0))
					for (u=i0;u<MINI(HISTOJOINT.nbI,i0+2);u++)
					for (v=j0;v<MINI(HISTOJOINT.nbJ,j0+2);v++)
						{
						if (u<x)
							temp=fabs(y-v)*fx[i]*fy[j]*fz[k];
						else
							temp=-1.0*fabs(y-v)*fx[i]*fy[j]*fz[k];
						
						dpij->raw[u][v][0].x+=temp*gx;
						dpij->raw[u][v][0].y+=temp*gy;
						dpij->raw[u][v][0].z+=temp*gz;
						
						}		
				
				}
	}
   
//temp=TOPO_SIZE_HISTOJOINT*TOPO_SIZE_HISTOJOINT*TOPO_SIZE_HISTOJOINT*HISTOJOINT.Ntot;
//temp=HISTOJOINT.Ntot*TOPO_SIZE_HISTOJOINT;
temp=HISTOJOINT.Ntot;
//temp=topDx*topDy*topDz*TOPO_SIZE_HISTOJOINT;
for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
{
dpij->raw[i][j][0].x=dpij->raw[i][j][0].x/temp;
dpij->raw[i][j][0].y=dpij->raw[i][j][0].y/temp;
dpij->raw[i][j][0].z=dpij->raw[i][j][0].z/temp;
}	 
	 
// calcul de dIM grace � la formule (23) de l'article Thevenaz IEEETMI dec 2000

grad[0]=grad[1]=grad[2]=0.0;

for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
	if ((HISTOJOINT.hij[i][j]>0.000001)&&(HISTOJOINT.hi[i]>0.000001))
	{
	
	temp=log(HISTOJOINT.hij[i][j]/HISTOJOINT.hi[i]);
	
	
	grad[0]+=temp*dpij->raw[i][j][0].x;
	grad[1]+=temp*dpij->raw[i][j][0].y;
	grad[2]+=temp*dpij->raw[i][j][0].z;
	
	}
/* grad[0]=grad[0]/(TOPO_SIZE_HISTOJOINT);
 grad[1]=grad[1]/(TOPO_SIZE_HISTOJOINT);
 grad[2]=grad[2]/(TOPO_SIZE_HISTOJOINT);*/

	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
						
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);
free_field3d(dpij);


  return(1);
}


/*******************************************************************************
**     gradient_base_ICP_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_base_ICP_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;

grphic3d *mask=NULL;
mask=I2; 

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
			if (I2->mri[i+bx0][j+by0][k+bz0]>0)	 
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
		diff=2.0*(double)(val);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
								
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     gradient_base_ICP_sym_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int gradient_base_ICP_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double gx,gy,gz,diff,temp,val;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG, ***maskdataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;

grphic3d *mask=NULL;
mask=I2; 

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 maskdataG=maskgradI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
			if (I2->mri[i+bx0][j+by0][k+bz0]>0)	 
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
		diff=2.0*(double)(val);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}


 /*on rempli le vecteur gradient en �chantillonnant sur l'image compl�mentaire */
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
			if (I2->mask->mri[i+bx0][j+by0][k+bz0]>0)	 
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    gx=gy=gz=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dx-(double)i2;yrel=dy-(double)j2;zrel=dz-(double)k2;
	
        va=maskdataG[i2][j2][k2];
        vb=maskdataG[i2][j2][k3];
        vc=maskdataG[i2][j3][k2];
        vd=maskdataG[i2][j3][k3];
        ve=maskdataG[i3][j2][k2];
        vf=maskdataG[i3][j2][k3];
        vg=maskdataG[i3][j3][k2];
        vh=maskdataG[i3][j3][k3];

        gx=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gy=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gz=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mask->mri[i2][j2][k2];
        wb=Im->mask->mri[i2][j2][k3];
        wc=Im->mask->mri[i2][j3][k2];
        wd=Im->mask->mri[i2][j3][k3];
        we=Im->mask->mri[i3][j2][k2];
        wf=Im->mask->mri[i3][j2][k3];
        wg=Im->mask->mri[i3][j3][k2];
        wh=Im->mask->mri[i3][j3][k3];

        val=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
		diff=2.0*(double)(val);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
	}


   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     gradient_atrophie_jacobien_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int gradient_atrophie_jacobien_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 	int 	i,j,k,l,wdth,hght,dpth;
 	double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
 	double gx,gy,gz,diff;
 	int bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 	int  D=TOP_D(nb_param);
 	vector3d ***dataG;
 	double J=1,dJx,dJy,dJz,xt,yt,zt;
	int aux0,aux1;  								// variables auxiliaires
	double auxdbl[20],aux2,mymono[21];
	double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
	int ideb,ifin,jdeb,jfin,kdeb,kfin;
	double dir_desR3[3];
	int resol;
	
 	wdth=I2->width;hght=I2->height;dpth=I2->depth;
 	dataG=gradI2->raw;
 	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 	fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
	resol=BASE3D.resol;
	//-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

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

 // mise a zero du gradient
 grad[0]=grad[1]=grad[2]=0.0;

 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
 bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
 topDx=bx1-bx0;
 topDy=by1-by0;
 topDz=bz1-bz0;


 
  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

 // on rempli le vecteur gradient

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	  	if (I2->mri[i][j][k]>0)
		{
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire 
	    gx=gy=gz=0;

		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);

		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		
		dJx=1.0*dJx/wdth;
		dJy=1.0*dJy/hght;
		dJz=1.0*dJz/dpth;
		
		
		
		diff=2.0*(J-(double)I2->mri[i][j][k]*I2->rcoeff);  
	
      grad[0]=grad[0]+diff*dJx;
      grad[1]=grad[1]+diff*dJy;
      grad[2]=grad[2]+diff*dJz;
			
	}
} 

  for (i=0;i<3;i++) grad[i]=grad[i]/(double)(topDx*topDy*topDz);
 
   
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

 
   return(1);
}

/*******************************************************************************
**     gradient_atrophie_jacobien_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int gradient_atrophie_log_jacobien_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 	int 	i,j,k,l,wdth,hght,dpth;
 	double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 	int	*x0,*y0,*z0,*x1,*y1,*z1;
 	double gx,gy,gz,diff;
 	int bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 	int  D=TOP_D(nb_param);
 	vector3d ***dataG;
 	double J=1,dJx,dJy,dJz,xt,yt,zt;
	int aux0,aux1;  								// variables auxiliaires
	double auxdbl[20],aux2,mymono[21];
	double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
	int ideb,ifin,jdeb,jfin,kdeb,kfin;
	double dir_desR3[3];
	int resol;
	
 	wdth=I2->width;hght=I2->height;dpth=I2->depth;
 	dataG=gradI2->raw;
 	x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 	fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
	resol=BASE3D.resol;
	//-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

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

 // mise a zero du gradient
 grad[0]=grad[1]=grad[2]=0.0;

 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
 bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
 topDx=bx1-bx0;
 topDy=by1-by0;
 topDz=bz1-bz0;


			
 
  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
  
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    aux2 = pow(2.0,-1.0 * resol);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,ordonnee);
    

		//--- initialisation du terme lineaire suivant x ------------------------------ 
			dir_desR3[0]=1.0;dir_desR3[1]=0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_x);
    	
		//--- initialisation du terme lineaire suivant y ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=1.0;dir_desR3[2]=0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_y);
    
		//--- initialisation du terme lineaire suivant z ------------------------------ 
			dir_desR3[0]=0;dir_desR3[1]=0;dir_desR3[2]=1.0;
	    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,resol,auxdbl);
	  	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
			poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,pente_z);
    	
ideb=Slpqr[aux0].b.xm*wdth;
ifin=Slpqr[aux0].b.xM*wdth;
jdeb=Slpqr[aux0].b.ym*hght;
jfin=Slpqr[aux0].b.yM*hght;
kdeb=Slpqr[aux0].b.zm*dpth;
kfin=Slpqr[aux0].b.zM*dpth;

 // on rempli le vecteur gradient

	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
    	if (I2->mri[i][j][k]>0)
		 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire 
	    gx=gy=gz=0;

		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);

		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
		
		dJx=1.0*dJx/wdth;
		dJy=1.0*dJy/hght;
		dJz=1.0*dJz/dpth;
		
			
			diff=2.0*log(J/(double)I2->mri[i][j][k]/I2->rcoeff)/J;
			 
      	if ((!isnan(fabs(diff)))&&(!isinf(fabs(diff))))
					{	grad[0]=grad[0]+diff*dJx;
      			grad[1]=grad[1]+diff*dJy;
      			grad[2]=grad[2]+diff*dJz; }
			
	}
} 

  for (i=0;i<3;i++) grad[i]=grad[i]/(double)(topDx*topDy*topDz);
 
   
	/* int signe=0;
	 for (i=0;i<3;i++)
	 	if (grad[i]>0)
	 		signe++;
			
		if (signe>1)
				{for (i=0;i<3;i++) grad[i]=1.0;}
		else
				{for (i=0;i<3;i++) grad[i]=-1.0;}			*/
	 
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}

 
   return(1);
}


/*******************************************************************************
**     gradient_quad_locale_symetrique_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_quad_locale_symetrique_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp;
 double dxreca,dyreca,dzreca,dxref,dyref,dzref;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG,***dataGref;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double valimreca, valimref;
 double gxreca,gyreca,gzreca,gxref,gyref,gzref;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataGref=maskgradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	    
	dxreca=(double)i+data[i][j][k].x+bx0;
    	dyreca=(double)j+data[i][j][k].y+by0;
    	dzreca=(double)k+data[i][j][k].z+bz0;
    	
	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gxreca=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gyreca=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gzreca=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        valimreca=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
				
		/* Evaluation du gradient de l'image ref en s=s-u avec interrpolation lineaire */
	   
	dxref=(double)i-lambda_ref*data[i][j][k].x+bx0;
    	dyref=(double)j-lambda_ref*data[i][j][k].y+by0;
    	dzref=(double)k-lambda_ref*data[i][j][k].z+bz0;
 
   	i2=(int)dxref;j2=(int)dyref;k2=(int)dzref;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
	
        va=dataGref[i2][j2][k2];
        vb=dataGref[i2][j2][k3];
        vc=dataGref[i2][j3][k2];
        vd=dataGref[i2][j3][k3];
        ve=dataGref[i3][j2][k2];
        vf=dataGref[i3][j2][k3];
        vg=dataGref[i3][j3][k2];
        vh=dataGref[i3][j3][k3];

        gxref=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gyref=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gzref=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imref en s=s-u avec interrpolation lineaire */	
	      wa=I2->mri[i2][j2][k2];
        wb=I2->mri[i2][j2][k3];
        wc=I2->mri[i2][j3][k2];
        wd=I2->mri[i2][j3][k3];
        we=I2->mri[i3][j2][k2];
        wf=I2->mri[i3][j2][k3];
        wg=I2->mri[i3][j3][k2];
        wh=I2->mri[i3][j3][k3];

        valimref=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }		
				
				
				
		diff=2.0*(double)(valimreca-valimref);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*(gxreca+lambda_ref*gxref);
      grad[1]=grad[1]+temp*(gyreca+lambda_ref*gyref);
      grad[2]=grad[2]+temp*(gzreca+lambda_ref*gzref);
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     gradient_quad_locale_symetrique_coupe_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_quad_locale_symetrique_coupe_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp;
 double dxreca,dyreca,dzreca,dxref,dyref,dzref;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG,***dataGref;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double valimreca, valimref;
 double gxreca,gyreca,gzreca,gxref,gyref,gzref;
 int Ntot=0;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;


 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataGref=maskgradI2->raw;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */

	dxreca=(double)i+data[i][j][k].x+bx0;
    	dyreca=(double)j+data[i][j][k].y+by0;
    	dzreca=(double)k+data[i][j][k].z+bz0;

    	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
	
        va=dataG[i2][j2][k2];
        vb=dataG[i2][j2][k3];
        vc=dataG[i2][j3][k2];
        vd=dataG[i2][j3][k3];
        ve=dataG[i3][j2][k2];
        vf=dataG[i3][j2][k3];
        vg=dataG[i3][j3][k2];
        vh=dataG[i3][j3][k3];

        gxreca=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gyreca=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gzreca=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imreca en s=s+u avec interrpolation lineaire */	
	      wa=Im->mri[i2][j2][k2];
        wb=Im->mri[i2][j2][k3];
        wc=Im->mri[i2][j3][k2];
        wd=Im->mri[i2][j3][k3];
        we=Im->mri[i3][j2][k2];
        wf=Im->mri[i3][j2][k3];
        wg=Im->mri[i3][j3][k2];
        wh=Im->mri[i3][j3][k3];

        valimreca=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }
				
				
		/* Evaluation du gradient de l'image ref en s=s-u avec interrpolation lineaire */
	   
	dxref=(double)i-lambda_ref*data[i][j][k].x+bx0;
    	dyref=(double)j-lambda_ref*data[i][j][k].y+by0;
    	dzref=(double)k-lambda_ref*data[i][j][k].z+bz0;

    	i2=(int)dxref;j2=(int)dyref;k2=(int)dzref;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
	
        va=dataGref[i2][j2][k2];
        vb=dataGref[i2][j2][k3];
        vc=dataGref[i2][j3][k2];
        vd=dataGref[i2][j3][k3];
        ve=dataGref[i3][j2][k2];
        vf=dataGref[i3][j2][k3];
        vg=dataGref[i3][j3][k2];
        vh=dataGref[i3][j3][k3];

        gxref=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gyref=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gzref=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Evaluation de Imref en s=s-u avec interrpolation lineaire */	
	      wa=I2->mri[i2][j2][k2];
        wb=I2->mri[i2][j2][k3];
        wc=I2->mri[i2][j3][k2];
        wd=I2->mri[i2][j3][k3];
        we=I2->mri[i3][j2][k2];
        wf=I2->mri[i3][j2][k3];
        wg=I2->mri[i3][j3][k2];
        wh=I2->mri[i3][j3][k3];

        valimref=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		
        }		
				
				
		if ((valimreca>0)&&(valimref>0))		
		{diff=2.0*(double)(valimreca-valimref);
      temp=diff*fx[i]*fy[j]*fz[k];
      grad[0]=grad[0]+temp*(gxreca+lambda_ref*gxref);
      grad[1]=grad[1]+temp*(gyreca+lambda_ref*gyref);
      grad[2]=grad[2]+temp*(gzreca+lambda_ref*gzref);
		Ntot++;	
		}
	}
   
 if (Ntot>0)
 {
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz)/(double) Ntot;
 }
	else
	 for (i=0;i<3;i++)
			grad[i]=0.0;
 
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     gradient_groupwise_variance_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_groupwise_variance_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
int 	i,j,k,l,width,height,depth;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int	*x0,*y0,*z0,*x1,*y1,*z1;
double temp;
double dxreca,dyreca,dzreca,dxref,dyref,dzref;
int ux,uy,uz;
int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
int x00,x11,y00,y11,z00,z11;
vector3d ***data,***Gref,***Gref2, ***Greca;
field3d *champ;
vector3d va,vb,vc,vd,ve,vf,vg,vh,vnull;
double wa,wb,wc,wd,we,wf,wg,wh;
double f,px,py,pz;
int auxi,auxj,auxk;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double valimreca, valimref;
double gxreca,gyreca,gzreca,gxref,gyref,gzref,gxref2,gyref2,gzref2,du;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
	
// vecteur nul
vnull.x=vnull.y=vnull.z=0.0;	
	
width=I2->width;height=I2->height;depth=I2->depth;
 
 
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
Gref=_ParamRecalageBspline.grad_ref_groupwise->raw;
Gref2=_ParamRecalageBspline.grad_ref2_groupwise->raw;
Greca=_ParamRecalageBspline.grad_reca_groupwise->raw;
 
 
 
/*mise a zero du gradient*/
grad[0]=grad[1]=grad[2]=0.0;

D=TOP_D(nb_param);
l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
/* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
/*on rempli le vecteur gradient*/
for (i=0;i<topDx;i++)
for (j=0;j<topDy;j++)
for (k=0;k<topDz;k++)
	 {
	/* Evaluation en s=s+u avec interrpolation lineaire */
	    
	dxreca=(double)i+data[i][j][k].x+bx0;
    dyreca=(double)j+data[i][j][k].y+by0;
    dzreca=(double)k+data[i][j][k].z+bz0;
    	
	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
    if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       	{
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
		
		/* pour interpoler le gradient */
        if (i2==-1 || j2==-1 || k2==-1) va=vnull; else va=Greca[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=vnull; else vb=Greca[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=vnull; else vc=Greca[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=vnull; else vd=Greca[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=vnull; else ve=Greca[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=vnull; else vf=Greca[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=vnull; else vg=Greca[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=vnull; else vh=Greca[i3][j3][k3];
	
        gxreca=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gyreca=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gzreca=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Pour interpoler l'image */	
		if (i2==-1 || j2==-1 || k2==-1) wa=0.0; else wa=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) wb=0.0; else wb=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) wc=0.0; else wc=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) wd=0.0; else wd=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) we=0.0; else we=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) wf=0.0; else wf=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) wg=0.0; else wg=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) wh=0.0; else wh=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k3];
				
        valimreca=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		}
	else
		{
		free_field3d(champ);
		for (i=0;i<3;i++) grad[i]=0.0;
		return(0);
		}
	
				
				
	/* Evaluation en s=s-u avec interrpolation lineaire */
	   
	dxref=(double)i-lambda_ref*data[i][j][k].x+bx0;
    dyref=(double)j-lambda_ref*data[i][j][k].y+by0;
    dzref=(double)k-lambda_ref*data[i][j][k].z+bz0;
 
   	i2=(int)dxref;j2=(int)dyref;k2=(int)dzref;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
    if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
        {
        xrel=dxref-(double)i2;yrel=dyref-(double)j2;zrel=dzref-(double)k2;
	    
		/* pour interpoler le gradient ref*/
        if (i2==-1 || j2==-1 || k2==-1) va=vnull; else va=Gref[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=vnull; else vb=Gref[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=vnull; else vc=Gref[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=vnull; else vd=Gref[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=vnull; else ve=Gref[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=vnull; else vf=Gref[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=vnull; else vg=Gref[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=vnull; else vh=Gref[i3][j3][k3];
	
	    gxref=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gyref=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gzref=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		/* Pour interpoler l'image ref*/	
		if (i2==-1 || j2==-1 || k2==-1) wa=0.0; else wa=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) wb=0.0; else wb=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) wc=0.0; else wc=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) wd=0.0; else wd=(double)_ParamRecalageBspline.ref_groupwise->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) we=0.0; else we=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) wf=0.0; else wf=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) wg=0.0; else wg=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) wh=0.0; else wh=(double)_ParamRecalageBspline.ref_groupwise->mri[i3][j3][k3];
			
        valimref=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		/* pour interpoler le gradient ref2*/
        if (i2==-1 || j2==-1 || k2==-1) va=vnull; else va=Gref2[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=vnull; else vb=Gref2[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=vnull; else vc=Gref2[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=vnull; else vd=Gref2[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=vnull; else ve=Gref2[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=vnull; else vf=Gref2[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=vnull; else vg=Gref2[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=vnull; else vh=Gref2[i3][j3][k3];
	
        gxref2=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
        gyref2=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
        gzref2=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

	   }
	else
		{
		free_field3d(champ);
		for (i=0;i<3;i++) grad[i]=0.0;
		return(0);
		}
	
	valimreca=valimreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
	valimref=valimref*_ParamRecalageBspline.ref_groupwise->rcoeff;
	
				
	du=fx[i]*fy[j]*fz[k]/_ParamRecalageBspline.nb_tot_groupwise;				
	temp=2.0*(valimref+valimreca)/_ParamRecalageBspline.nb_tot_groupwise;
	
	
      grad[0]=grad[0]+du*(2.0*gxreca*valimreca-lambda_ref*gxref2-temp*(gxreca-lambda_ref*gxref));
      grad[1]=grad[1]+du*(2.0*gyreca*valimreca-lambda_ref*gyref2-temp*(gyreca-lambda_ref*gyref));
      grad[2]=grad[2]+du*(2.0*gzreca*valimreca-lambda_ref*gzref2-temp*(gzreca-lambda_ref*gzref));
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);
  	
  
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     gradient_groupwise_variance_locale_nonsym_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_groupwise_variance_locale_nonsym_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
int 	i,j,k,l,width,height,depth;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int	*x0,*y0,*z0,*x1,*y1,*z1;
double temp;
double dxreca,dyreca,dzreca;
int ux,uy,uz;
int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
int x00,x11,y00,y11,z00,z11;
vector3d ***data;
field3d *champ;
vector3d va,vb,vc,vd,ve,vf,vg,vh,vnull;
double wa,wb,wc,wd,we,wf,wg,wh;
double f,px,py,pz;
int auxi,auxj,auxk;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double valimreca, valimref;
double du, gxreca,gyreca,gzreca;
	
// vecteur nul
vnull.x=vnull.y=vnull.z=0.0;	
	
width=I2->width;height=I2->height;depth=I2->depth;
 
 
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 
 
/*mise a zero du gradient*/
grad[0]=grad[1]=grad[2]=0.0;

D=TOP_D(nb_param);
l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
/* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
/*on rempli le vecteur gradient*/
for (i=0;i<topDx;i++)
for (j=0;j<topDy;j++)
for (k=0;k<topDz;k++)
	 {
	/* Evaluation en s=s+u avec interrpolation lineaire */
	    
	dxreca=(double)i+data[i][j][k].x+bx0;
    dyreca=(double)j+data[i][j][k].y+by0;
    dzreca=(double)k+data[i][j][k].z+bz0;
    	
	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    i3=i2+1;j3=j2+1;k3=k2+1;
    
    if (i3>=0 && i2<width && j3>=0 && j2<height && k3>=0 && k2<depth)
       	{
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
		
		/* pour interpoler le gradient */
        if (i2==-1 || j2==-1 || k2==-1) va=vnull; else va=gradI2->raw[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) vb=vnull; else vb=gradI2->raw[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) vc=vnull; else vc=gradI2->raw[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) vd=vnull; else vd=gradI2->raw[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) ve=vnull; else ve=gradI2->raw[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) vf=vnull; else vf=gradI2->raw[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) vg=vnull; else vg=gradI2->raw[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) vh=vnull; else vh=gradI2->raw[i3][j3][k3];
	
        gxreca=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		gxreca=gxreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
		
        gyreca=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		gyreca=gyreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
		
		
        gzreca=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		gzreca=gzreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
		
		/* Pour interpoler l'image */	
		if (i2==-1 || j2==-1 || k2==-1) wa=0.0; else wa=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==depth) wb=0.0; else wb=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j2][k3];
        if (i2==-1 || j3==height || k2==-1) wc=0.0; else wc=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k2];
        if (i2==-1 || j3==height || k3==depth) wd=0.0; else wd=(double)_ParamRecalageBspline.reca_groupwise->mri[i2][j3][k3];
        if (i3==width || j2==-1 || k2==-1) we=0.0; else we=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k2];
        if (i3==width || j2==-1 || k3==depth) wf=0.0; else wf=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j2][k3];
        if (i3==width || j3==height || k2==-1) wg=0.0; else wg=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k2];
        if (i3==width || j3==height || k3==depth) wh=0.0; else wh=(double)_ParamRecalageBspline.reca_groupwise->mri[i3][j3][k3];
				
        valimreca=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		}
	else
		{
		free_field3d(champ);
		for (i=0;i<3;i++) grad[i]=0.0;
		return(0);
		}
	
				
		valimref=_ParamRecalageBspline.ref_groupwise->mri[i+bx0][j+by0][k+bz0];
        
	
		valimreca=valimreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
		valimref=valimref*_ParamRecalageBspline.ref_groupwise->rcoeff;
	
				
		du=2.0*fx[i]*fy[j]*fz[k]/_ParamRecalageBspline.nb_tot_groupwise/_ParamRecalageBspline.nb_tot_groupwise;				
		temp=((_ParamRecalageBspline.nb_tot_groupwise-1.0)*valimreca-valimref)*du;
	
	
      grad[0]=grad[0]+gxreca*temp;
      grad[1]=grad[1]+gyreca*temp;
      grad[2]=grad[2]+gzreca*temp;
	}
   
  for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);
  	
  
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     gradient_quad_sous_champ_locale_3d                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/


int gradient_quad_sous_champ_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, double *grad, int topi, int topj, int topk,TSlpqr *Slpqr, reg_grad_locale_t gradient_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double temp;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data;
 field3d *champ;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 /*mise a zero du gradient*/
 grad[0]=grad[1]=grad[2]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 /* calcul du champ de deformation*/
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

 for (i=MAXI(topi-1,0);i<=MINI(topi+1,D-1);i++)
	for (j=MAXI(topj-1,0);j<=MINI(topj+1,D-1);j++)
		for (k=MAXI(topk-1,0);k<=MINI(topk+1,D-1);k++)
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

 
 /*on rempli le vecteur gradient*/
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   /* Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire */
	
		i2=i+bx0;
		j2=j+by0;
		k2=k+bz0;
		
		dx=data[i][j][k].x-_ParamRecalageBspline.chref->raw[i2][j2][k2].x;
		dy=data[i][j][k].y-_ParamRecalageBspline.chref->raw[i2][j2][k2].y;
		dz=data[i][j][k].z-_ParamRecalageBspline.chref->raw[i2][j2][k2].z;

	
				
	      	temp=2*fx[i]*fy[j]*fz[k];
      		grad[0]=grad[0]+temp*dx;
      		grad[1]=grad[1]+temp*dy;
      		grad[2]=grad[2]+temp*dz;
	}
   
 /* for (i=0;i<3;i++)
	grad[i]=grad[i]/(double)(topDx*topDy*topDz);*/
  
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			double grad_reg[3];
			gradient_reg(nb_param, param, param_norm, grad_reg,  topi,  topj,  topk, NULL,Slpqr);
			
					
		  for (i=0;i<3;i++)
			grad[i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*grad_reg[i];
			}
	
	free_field3d(champ);


  return(1);
}



/*******************************************************************************
**     gradient_regularisation_patch_imagebased_local                  
**                                                                   
**     Calcul du gradient de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int gradient_regularisation_patch_imagebased_local(int nb_param, double *param, double *param_norm, double *grad, int topi, int topj, int topk,grphic3d *mask,TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int x00,x11,y00,y11,z00,z11,bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,f,px,py,pz;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
double tmp;
field3d *champ;
vector3d ***data;
dvector3d u;
int ideb,jdeb,kdeb, ifin, jfin, kfin;
double ***poids, dx,dy,dz, dux, duy, duz , dutx, duty, dutz;
int nhx,nhy,nhz,spx,spy,spz,di,dj,dk;

 	
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
dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

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

l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

grad[0]=0;grad[1]=0;grad[2]=0;


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

	dux=dfx[ux]*fy[uy]*fz[uz];
	duy=fx[ux]*dfy[uy]*fz[uz];
	duz=fx[ux]*fy[uy]*dfz[uz];
	
	
	    	
 	NLMWeightOnePoint3D(_ParamRecalageBspline.ref_groupwise, poids, TOPO_PATCH_COEFF,  
					TOPO_NLMhwnx,  TOPO_NLMhwny,  TOPO_NLMhwnz,  TOPO_NLMhwvsx,  TOPO_NLMhwvsy,  TOPO_NLMhwvsz,auxi, auxj, auxk);
		
		for (i=ideb; i<ifin; i++)
		for (j=jdeb; j<jfin; j++)
		for (k=kdeb; k<kfin; k++)
			{
			u.x=u.y=u.z=0; 
			eval_deplacement_bspline1_3d(param, nb_param, width, height, depth, i, j, k, &u);   
    		dx=data[ux][uy][uz].x-u.x;
			dy=data[ux][uy][uz].y-u.y;
			dz=data[ux][uy][uz].z-u.z;
			
			di = i-bx0;
			dj = j-by0;
			dk = k-bz0;
			
			
						
			tmp=poids[di-ux+TOPO_NLMhwnx][dj-uy+TOPO_NLMhwny][dk-uz+TOPO_NLMhwnz];
			
			if ((di>=0)&&(di<topDx)&&(dj>=0)&&(dj<topDy)&&(dk>=0)&&(dk<topDz))
				{
				dutx=dfx[di]*fy[dj]*fz[dk];
				duty=fx[di]*dfy[dj]*fz[dk];
				dutz=fx[di]*fy[dj]*dfz[dk];
				}
				else
				{
				dutx=duty=dutz=0.0;
				}
				
			grad[0]+=tmp*dx*(dutx-dux);
			grad[1]+=tmp*dy*(duty-duy);
			grad[2]+=tmp*dz*(dutz-duz);
			
			}
		
		
		 
	}



for (i=0;i<3;i++)
	grad[i]=2.0*grad[i]/(double)(topDx*topDy*topDz);


for (i=0;i<3;i++)
	grad[i]=grad[i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


free_dmatrix_3d(poids);
free_field3d(champ);
return(1);
}

