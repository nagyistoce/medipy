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
#include "recalage/topo/hessien_distance_topo.h"

/*******************************************************************************
********************************************************************************
****************** ALLOCATION ET LIBERATION MEMOIRE ****************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**     cr_hessien3d(wdth,hght,dpth)                                   
*/
/*!                                                                   
**     Creation et allocation memoire d'une structure hessien pour les 
**	tailles specifiees en parametres  
**  \param  wdth, hght, dpth : taille du champ
**  \retval hessien3d, le champ alloue, exit(1) si echec
*******************************************************************************/
hessien3d *cr_hessien3d(int wdth, int hght, int dpth)
{ 
  mat3d ***data; 
  mat3d *mt;
  hessien3d  *ch;
  int i,j;
  
  /*allocation memoire de la structure hessien*/
  ch=(hessien3d*)malloc(sizeof(hessien3d));
  if (ch==NULL) {printf("allocation failed in step 1 of cr_hessien3d \n");exit(1);}

  /*allocation memoire du tableau de donnee data*/
  data=(mat3d ***) malloc((size_t)((wdth+1)*sizeof(mat3d **)));
  if (data==NULL) {printf("allocation failed in step 2 of cr_hessien3d \n");exit(1);}

  data[0]=(mat3d **)malloc((size_t)((wdth*hght)*sizeof(mat3d *)));
  if (data[0]==NULL) {printf("allocation failed in step 3 of cr_hessien3d \n");exit(1);}
  
  for(i=1;i<wdth;i++)
     data[i]=data[i-1]+hght; 
 
  data[0][0]=(mat3d *)malloc((size_t)((wdth*hght*dpth)*sizeof(mat3d)));
  if (data[0][0]==NULL) {printf("\n allocation failed in step4  of cr_hessien3d \n"); exit(1);}
  mt=data[0][0];
  
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
   {data[i][j]=mt;mt=mt+dpth;}
  
  ch->raw=data;
  ch->width=wdth;ch->height=hght;ch->depth=dpth;
  
  return (ch);
 }


/*******************************************************************************
**    free_hessien3d(ch)                                   
*/
/*!                                                                   
**     liberation memoire d'une structure hessien:
**	 \param ch: hessien a liberer
**	 \retval 1
*******************************************************************************/
 int free_hessien3d(hessien3d *ch)
{
  mat3d ***data;
  
  data=ch->raw;

  /*liberation memoire des donnees data*/
  free((char *) (data[0][0]));
  free((char*)(data[0]));
  free((char*)data);

  /*liberation memoire de la structure*/
  free(ch);
  
  return (1);
 }

/*******************************************************************************
********************************************************************************
*******************  CALCUL DU HESSIEN D'UNE IMAGE   ***************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
** 	int imx_hessien_3d_p(im,grad,method,t)
*/
/*!	Calcul le hessien d'une image et met le resultat dans un tableau de grphic3d
**		\param im: 	image 
**		\param hessien: 	resultat
**		\param method:	methode de calcul 
**				\li  0: classic
**				\li  1: taylor
**				\li  2: leat square
**		\param t:	taille du filtre si necessaire 
**
**	\attention le champ est suppose alloue a la bonne taille avant l'appel
**	a cette procedure
**
**	\remark  pour plus de details sur les methodes voir le livre:
**	practical handbook on image processing for scientific applications
**	Bernd Jahne, crc press
**
**  hessien.H[i][j]=d./dxi.d./dxj I (gradient suivant x puis suivant y)
**
*******************************************************************************/
int imx_hessien_3d_p(grphic3d *im, hessien3d *hessien, int method, int t)
{
int wdth,hght,dpth;
grphic3d *imtmp;
field3d *grad;
int i,j,k;

wdth=im->width;
hght=im->height;
dpth=im->depth;
imtmp=cr_grphic3d(im);
grad=cr_field3d(wdth,hght,dpth);
//gradtmp=cr_field3d(wdth,hght,dpth);

imx_gradient_3d_p(im,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			imtmp->mri[i][j][k]=(int)grad->raw[i][j][k].x;

	imx_gradient_3d_p(imtmp,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			hessien->raw[i][j][k].H[0][0]=grad->raw[i][j][k].x;
			hessien->raw[i][j][k].H[1][0]=grad->raw[i][j][k].y;
			hessien->raw[i][j][k].H[2][0]=grad->raw[i][j][k].z;
			}

imx_gradient_3d_p(im,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			imtmp->mri[i][j][k]=(int)grad->raw[i][j][k].y;

	imx_gradient_3d_p(imtmp,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			hessien->raw[i][j][k].H[0][1]=grad->raw[i][j][k].x;
			hessien->raw[i][j][k].H[1][1]=grad->raw[i][j][k].y;
			hessien->raw[i][j][k].H[2][1]=grad->raw[i][j][k].z;
			}

imx_gradient_3d_p(im,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			imtmp->mri[i][j][k]=(int)grad->raw[i][j][k].z;

	imx_gradient_3d_p(imtmp,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			hessien->raw[i][j][k].H[0][2]=grad->raw[i][j][k].x;
			hessien->raw[i][j][k].H[1][2]=grad->raw[i][j][k].y;
			hessien->raw[i][j][k].H[2][2]=grad->raw[i][j][k].z;
			}




free_grphic3d(imtmp);free_field3d(grad);//free_field3d(gradtmp);
return(1);
}


/*******************************************************************************
** 	mul_hessien_3d(hessien,hessienres,coeff)                                   
*******************************************************************************/
int	mul_hessien_3d(hessien3d *hessien, hessien3d *hessienres, float coeff)
{
 int	wdth,hght,dpth,i,j,k, u,v;
 mat3d ***data1,***datares;
 
 wdth=hessien->width;hght=hessien->height;dpth=hessien->depth;
 
 
 data1=hessien->raw;datares=hessienres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
   	for (u=0;u<3;u++)
	for (v=0;v<3;v++)
	  (datares[i][j][k]).H[u][v]=(data1[i][j][k]).H[u][v]*coeff;
   }
  
 return(1);
}


/*******************************************************************************
********************************************************************************
************  CALCUL DU HESSIEN RELATIF AU FONCTION DE COUT  *******************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**     hessien_base_quad_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_quad_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
        
		  diff=(double)(I2->mri[i+bx0][j+by0][k+bz0]-val);
      temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*(g[auxi]*g[auxj]-diff*Hval.H[auxi][auxj]);
								
				}	
				
				}		
			}
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
					
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     hessien_base_quad_sym_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_quad_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;
 double J=1,dJx,dJy,dJz,xt,yt,zt,dJ[3];
 int aux0,aux1;  								// variables auxiliaires
 double auxdbl[20],aux2,mymono[21];
 double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
 int ideb,ifin,jdeb,jfin,kdeb,kfin;
 double dir_desR3[3];
 int signeJ=1;
 int resol;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
 resol=BASE3D.resol;
 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
for (i=ideb;i<ifin;i++)
for (j=jdeb;j<jfin;j++)
for (k=kdeb;k<kfin;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
			 
			diff=(double)(I2->mri[i][j][k]-val);
    
		 if (diff!=0)
			{
		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
		signeJ=1;
		if (J<0) signeJ=-1;
		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
	
		dJ[0]=1.0*dJx/wdth;
		dJ[1]=1.0*dJy/hght;
		dJ[2]=1.0*dJz/dpth;
	
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	Hval.H[auxi][auxj]=(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel+(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel+(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel
					+(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel+(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel
					+(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
							  
			temp=fx[i-bx0]*fy[j-by0]*fz[k-bz0];
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{
hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+(1.0+fabs(J))*temp*temp*(g[auxi]*g[auxj]-diff*Hval.H[auxi][auxj])
										-1.0*signeJ*diff*(g[auxi]*temp*dJ[auxj]+g[auxj]*temp*dJ[auxi]);
				}
			}
			else
			{
			xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
		signeJ=1;
						  
		 temp=fx[i-bx0]*fy[j-by0]*fz[k-bz0];	
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{
hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+(1.0+fabs(J))*temp*temp*g[auxi]*g[auxj];
				}
			
			
			}
			
			}
				
			}
   }
tmp=(topDx*topDy*topDz); 
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);

  return(1);
}

/*******************************************************************************
**     hessien_base_L1L2_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_L1L2_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
 hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,val,g[3],diff2;
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	Hval.H[auxi][auxj]=(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel+(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel+(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel
					+(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel+(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel
					+(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];

        
		  diff=(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
			diff2=1.0/sqrt(diff*diff+TOPO_EPSILON);
      temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+temp*((g[auxi]*g[auxj]+diff*Hval.H[auxi][auxj])*diff2-g[auxi]*g[auxj]*diff*diff*diff2*diff2*diff2);
			}
   	}
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);

  return(1);
}

/*******************************************************************************
**     hessien_base_L1L2_sym_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_L1L2_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
 hessien3d *hessienI2,hessien3d *maskhessienI2, mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,diff2,temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;
 double J=1,dJx,dJy,dJz,xt,yt,zt,dJ[3];
 int aux0,aux1;  								// variables auxiliaires
 double auxdbl[20],aux2,mymono[21];
 double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
 int ideb,ifin,jdeb,jfin,kdeb,kfin;
 double dir_desR3[3];
 int signeJ=1;
 int resol;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
 resol=BASE3D.resol;
 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
for (i=ideb;i<ifin;i++)
for (j=jdeb;j<jfin;j++)
for (k=kdeb;k<kfin;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];

		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);
		signeJ=1;
		if (J<0) signeJ=-1;
		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
	
		dJ[0]=1.0*dJx/wdth;
		dJ[1]=1.0*dJy/hght;
		dJ[2]=1.0*dJz/dpth;
	
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	Hval.H[auxi][auxj]=(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel+(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel+(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel
					+(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel+(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel
					+(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];

        


		  diff=(double)(val-I2->mri[i][j][k]);
			diff2=1.0/sqrt(diff*diff+TOPO_EPSILON);
      temp=fx[i-bx0]*fy[j-by0]*fz[k-bz0];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+(1.0+J)*temp*temp*((g[auxi]*g[auxj]+diff*Hval.H[auxi][auxj])*diff2-g[auxi]*g[auxj]*diff*diff*diff2*diff2*diff2)
																											+1.0*signeJ*diff2*temp*(dJ[auxi]*g[auxj]+dJ[auxj]*g[auxi]);
		


/*		  diff=(double)(I2->mri[i][j][k]-val);
      temp=fx[i-bx0]*fy[j-by0]*fz[k-bz0];
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{
hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+(1.0+fabs(J))*temp*temp*(g[auxi]*g[auxj]-diff*Hval.H[auxi][auxj])
										-1.0*signeJ*diff*(g[auxi]*temp*dJ[auxj]+g[auxj]*temp*dJ[auxi]);
				}


diff=(double)(val-I2->mri[i+bx0][j+by0][k+bz0]);
      temp=diff*fx[i]*fy[j]*fz[k];
			temp=temp*pow((diff*diff)+TOPO_EPSILON,-0.5);
      grad[0]=grad[0]+temp*gx;
      grad[1]=grad[1]+temp*gy;
      grad[2]=grad[2]+temp*gz;
			
			*/
			 }
			}
			
   }

tmp=(topDx*topDy*topDz); 
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);

  return(1);
}


/*******************************************************************************
**     hessien_base_geman_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_geman_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, 
hessien3d *hessienI2,hessien3d *maskhessienI2, mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp,tutu;


 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
        
		  diff=(double)(I2->mri[i+bx0][j+by0][k+bz0]-val);
			diff=diff*diff;
			tutu=TOPO_ALPHA_ROBUST+diff;
			diff=(2.0*TOPO_ALPHA_ROBUST*TOPO_ALPHA_ROBUST-6.0*TOPO_ALPHA_ROBUST*diff)/(tutu*tutu*tutu);
      temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*(g[auxi]*g[auxj]-diff*Hval.H[auxi][auxj]);
								
				}	
				
				}		
			}
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp/TOPO_ALPHA_ROBUST;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];

 	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	 
  free_field3d(champ);


  return(1);
}




/*******************************************************************************
**     hessien_base_quad_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_regularisation_energie_membrane_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,px,py,pz;
int i,j;
int ux,uy,uz;
int uxdeb,uydeb,uzdeb;
double dJxx,dJxy,dJxz,dJyx,dJyy,dJyz,dJzx,dJzy,dJzz;
int ideb,jdeb,kdeb,condition;
double tmp;

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


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

for (i=0;i<3;i++)
for (j=0;j<3;j++)
hessien->H[i][j]=0;



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
		 
		 
		 hessien->H[0][0]+=dJxx*dJxx+dJxy*dJxy+dJxz*dJxz;
		 hessien->H[1][1]+=dJyx*dJyx+dJyy*dJyy+dJyz*dJyz;
		 hessien->H[2][2]+=dJzx*dJzx+dJzy*dJzy+dJzz*dJzz;
		 }
		 
		 }

tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
		hessien->H[i][i]=2.0*hessien->H[i][i]/tmp;


 for (i=0;i<3;i++)
		hessien->H[i][i]=hessien->H[i][i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;

return(1);
}

/*******************************************************************************
**     hessien_regularisation_energie_membrane_Lp_local                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_regularisation_energie_membrane_Lp_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
int x00,x11,y00,y11,z00,z11;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,px,py,pz,f;
int i,j,k,auxi,auxj,auxk;
int ux,uy,uz;
int uxdeb,uydeb,uzdeb;
double dJxx,dJxy,dJxz,dJyx,dJyy,dJyz,dJzx,dJzy,dJzz;
int ideb,jdeb,kdeb,condition;
double tmp;
double ORDRE_REGULARISATION_LP_moins_1=ORDRE_REGULARISATION_LP-1.0;
double ORDRE_REGULARISATION_LP_moins_2=ORDRE_REGULARISATION_LP-2.0;
double Jxx,Jxy,Jxz,Jyx,Jyy,Jyz,Jzx,Jzy,Jzz;
field3d *champ;
vector3d ***data;

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

for (i=0;i<3;i++)
for (j=0;j<3;j++)
hessien->H[i][j]=0;



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
	 
		 if (px!=0)
		 {
		 dJxx=dfx[ux]*fy[uy]*fz[uz];
		 dJxy=fx[ux]*dfy[uy]*fz[uz];
		 dJxz=fx[ux]*fy[uy]*dfz[uz];
		 }
		 else
		 {dJxx=dJxy=dJxz=0.0;Jxx=Jxy=Jxz=1;}
		 
		 
		if (py!=0)
		 {
		 dJyx=dfx[ux]*fy[uy]*fz[uz];
		 dJyy=fx[ux]*dfy[uy]*fz[uz];
		 dJyz=fx[ux]*fy[uy]*dfz[uz];
		 }
		else
		 {dJyx=dJyy=dJyz=0.0;Jyx=Jyy=Jyz=1;}
		 
		if (pz!=0)
		 {
		 dJzx=dfx[ux]*fy[uy]*fz[uz];
		 dJzy=fx[ux]*dfy[uy]*fz[uz];
		 dJzz=fx[ux]*fy[uy]*dfz[uz];
		 }
		else
		 {dJzx=dJzy=dJzz=0.0;Jzx=Jzy=Jzz=1;}
		 
		 
		 
		 tmp=dJxx*dJxx*pow(Jxx,ORDRE_REGULARISATION_LP_moins_2);	 
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[0][0]+=tmp;
		 
		 
		 tmp=dJxy*dJxy*pow(Jxy,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[0][0]+=tmp;
		 
		 tmp=dJxz*dJxz*pow(Jxz,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[0][0]+=tmp;
		 
		 
		 
		 tmp=dJyx*dJyx*pow(Jyx,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[1][1]+=tmp;
		 
		 tmp=dJyy*dJyy*pow(Jyy,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[1][1]+=tmp;
		 
		 tmp=dJyz*dJyz*pow(Jyz,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[1][1]+=tmp;
		 
		 
		
		 tmp=dJzx*dJzx*pow(Jzx,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp))
		 	 {tmp=0.0;}
		 hessien->H[2][2]+=tmp;
		 
		 tmp=dJzy*dJzy*pow(Jzy,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[2][2]+=tmp;
		
		 tmp=dJzz*dJzz*pow(Jzz,ORDRE_REGULARISATION_LP_moins_2);
		 if (isnan(tmp)) 
		 	{tmp=0.0;}
		 hessien->H[2][2]+=tmp;
		 }
		 }

tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
		hessien->H[i][i]=1.0*ORDRE_REGULARISATION_LP*ORDRE_REGULARISATION_LP_moins_1*hessien->H[i][i]/tmp;


for (i=0;i<3;i++)
		hessien->H[i][i]=hessien->H[i][i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,ORDRE_REGULARISATION_LP)); // Pour assurer la compatibilite avec le recalage symetrique


free_field3d(champ);

return(1);
}



/*******************************************************************************
**     hessien_regularisation_log_jacobien_local                  
**                                                                   
*******************************************************************************/
int hessien_regularisation_log_jacobien_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int auxi,auxj;
 double tmp;
 double J=1,dJx,dJy,dJz,xt,yt,zt,dJ[3];
 int aux0,aux1;  								// variables auxiliaires
 double auxdbl[20],aux2,mymono[21];
 double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
 int ideb,ifin,jdeb,jfin,kdeb,kfin;
 double dir_desR3[3];
 int resol;


 wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;

 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
 resol=BASE3D.resol;
 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
 bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
 topDx=bx1-bx0;
 topDy=by1-by0;
 topDz=bz1-bz0;


 		
 
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

 
 //on rempli le vecteur gradient
	 
	for (i=ideb;i<ifin;i++)
	for (j=jdeb;j<jfin;j++)
	for (k=kdeb;k<kfin;k++)
	{
 		
		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);

		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
	
		dJ[0]=1.0*dJx/wdth;
		dJ[1]=1.0*dJy/hght;
		dJ[2]=1.0*dJz/dpth;
	
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	if (J>0)
			{
			hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*dJ[auxi]*dJ[auxj]*(1-log(J))/J/J;
			}
	}				
 }
tmp=(topDx*topDy*topDz); 
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=2.0*hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
 for (i=0;i<3;i++)
 for (j=0;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;

return(1);
}



/*******************************************************************************
**     hessien_base_quad_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_quad_locale_pondere_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3,Ntot;
 double xrel,yrel,zrel,Ptot,poids;
 double tmp;

Ptot=0;Ntot=0;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
    	val=0;
			dx=(double)i+data[i][j][k].x+bx0;
    	dy=(double)j+data[i][j][k].y+by0;
    	dz=(double)k+data[i][j][k].z+bz0;
    	i2=(int)dx;j2=(int)dy;k2=(int)dz;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    poids=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
    
			   
		poids=I2->mask->mri[i+bx0][j+by0][k+bz0];
		
		Ptot+=poids;	
		Ntot++;	
		    
		  diff=(double)(I2->mri[i+bx0][j+by0][k+bz0]-val);
      temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp*poids;
    
		
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*(g[auxi]*g[auxj]-diff*Hval.H[auxi][auxj]);
								
				}	
				
				}		
			}
   
if (Ptot>0)
{
tmp=1.0*(topDx*topDy*topDz*Ptot)/Ntot;
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/(tmp);
}
else
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=0;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
 	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
 
  free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     hessien_regularisation_dist_identite_local                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_regularisation_dist_identite_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,px,py,pz;
int i,j;
int ux,uy,uz;
double dJxx,dJyy,dJzz;
int ideb,jdeb,kdeb;
double tmp;

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


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

for (i=0;i<3;i++)
for (j=0;j<3;j++)
hessien->H[i][j]=0;



	for (ux=0;ux<topDx;ux++)
	for (uy=0;uy<topDy;uy++)
	for (uz=0;uz<topDz;uz++)
	   {
    		
		 
		 
		 if (px!=0)
		 {
		 dJxx=dfx[ux]*fy[uy]*fz[uz];
		 }
		 else
		 {dJxx=0.0;}
		 
		 if (py!=0)
		 {
		 dJyy=fx[ux]*dfy[uy]*fz[uz];
		 }
		 else
		 {dJyy=0.0;}
		 
	
		 if (pz!=0)
		 {
		 dJzz=fx[ux]*fy[uy]*dfz[uz];
		 }
		 else
		 {dJzz=0.0;}
	
	
		 
		 hessien->H[0][0]+=dJxx*dJxx;
		 hessien->H[1][1]+=dJyy*dJyy;
		 hessien->H[2][2]+=dJzz*dJzz;
		 
		 
		 }

tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
		hessien->H[i][i]=2.0*hessien->H[i][i]/tmp;


 for (i=0;i<3;i++)
		hessien->H[i][i]=hessien->H[i][i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;

return(1);
}


/*******************************************************************************
**     hessien_base_IM_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_IM_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth,u,v;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ,*dpij;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ***dataH;
 /*mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;*/
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double x,y;
 double dpix,dpiy,dpiz;
 int i0,j0;
/* hessien3d *d2pij; */

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 	dpij=cr_field3d(HISTOJOINT.nbI,HISTOJOINT.nbJ,1);
/*	d2pij=cr_hessien3d(HISTOJOINT.nbI,HISTOJOINT.nbJ,1);*/


for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
 	{
	dpij->raw[i][j][0].x=dpij->raw[i][j][0].y=dpij->raw[i][j][0].z=0.0;
	}
/*	for (u=0;u<3;u++)
	for (v=0;v<3;v++)
		d2pij->raw[i][j][0].H[u][v]=0.0;
		
	}*/
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
	
	/*	// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
    
*/
	
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
						
						dpij->raw[u][v][0].x+=temp*g[0];
						dpij->raw[u][v][0].y+=temp*g[1];
						dpij->raw[u][v][0].z+=temp*g[2];
    
							temp=temp*fx[i]*fy[j]*fz[k];
			/*	// calcul de d2pij 
							for (auxi=0;auxi<3;auxi++)
							for (auxj=auxi;auxj<3;auxj++)
								d2pij->raw[u][v][0].H[auxi][auxj]+=temp*Hval.H[auxi][auxj];*/
		
		
		
						}
						
								

			
   
}}

//temp=TOPO_SIZE_HISTOJOINT*TOPO_SIZE_HISTOJOINT*TOPO_SIZE_HISTOJOINT;
temp=HISTOJOINT.Ntot;
//temp=topDx*topDy*topDz*TOPO_SIZE_HISTOJOINT;
for (i=0;i<HISTOJOINT.nbI;i++)
for (j=0;j<HISTOJOINT.nbJ;j++)
{
dpij->raw[i][j][0].x=dpij->raw[i][j][0].x/temp;
dpij->raw[i][j][0].y=dpij->raw[i][j][0].y/temp;
dpij->raw[i][j][0].z=dpij->raw[i][j][0].z/temp;

/*	for (auxi=0;auxi<3;auxi++)
	for (auxj=auxi;auxj<3;auxj++)
			d2pij->raw[i][j][0].H[auxi][auxj]=d2pij->raw[i][j][0].H[auxi][auxj]/temp;*/
		

}


for (i=0;i<HISTOJOINT.nbI;i++)
if (HISTOJOINT.hi[i]>0.000001)
{
dpix=dpiy=dpiz=0.0;

for (j=0;j<HISTOJOINT.nbJ;j++)
	{
	dpix+=dpij->raw[i][j][0].x;
	dpiy+=dpij->raw[i][j][0].y;
	dpiz+=dpij->raw[i][j][0].z;
	}
	
hessien->H[0][0]+=dpix*dpix/HISTOJOINT.hi[i];
hessien->H[0][1]+=dpix*dpiy/HISTOJOINT.hi[i];
hessien->H[0][2]+=dpix*dpiz/HISTOJOINT.hi[i];
hessien->H[1][1]+=dpiy*dpiy/HISTOJOINT.hi[i];
hessien->H[1][2]+=dpiy*dpiz/HISTOJOINT.hi[i];
hessien->H[2][2]+=dpiz*dpiz/HISTOJOINT.hi[i];

}


for (i=0;i<HISTOJOINT.nbI;i++)
	for (j=0;j<HISTOJOINT.nbJ;j++)
	if (HISTOJOINT.hij[i][j]>0.000001)
	{
	hessien->H[0][0]-=dpij->raw[i][j][0].x*dpij->raw[i][j][0].x/HISTOJOINT.hij[i][j];
	hessien->H[0][1]-=dpij->raw[i][j][0].x*dpij->raw[i][j][0].y/HISTOJOINT.hij[i][j];
	hessien->H[0][2]-=dpij->raw[i][j][0].x*dpij->raw[i][j][0].z/HISTOJOINT.hij[i][j];
	hessien->H[1][1]-=dpij->raw[i][j][0].y*dpij->raw[i][j][0].y/HISTOJOINT.hij[i][j];
	hessien->H[1][2]-=dpij->raw[i][j][0].y*dpij->raw[i][j][0].z/HISTOJOINT.hij[i][j];
	hessien->H[2][2]-=dpij->raw[i][j][0].z*dpij->raw[i][j][0].z/HISTOJOINT.hij[i][j];
	}

 
/*for (i=0;i<HISTOJOINT.nbI;i++)
	for (j=0;j<HISTOJOINT.nbJ;j++)
	if (HISTOJOINT.hij[i][j]>0)
	{
	for (auxi=0;auxi<3;auxi++)
	for (auxj=auxi;auxj<3;auxj++)
			{
			hessien->H[auxi][auxj]-=d2pij->raw[i][j][0].H[auxi][auxj]*log(HISTOJOINT.hij[i][j]/HISTOJOINT.hi[i]);
			}

	}*/
  
 
/* for (i=0;i<3;i++)
	for (j=i;j<3;j++)
 	 hessien->H[i][j]=hessien->H[i][j]/(TOPO_SIZE_HISTOJOINT);*/
 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
free_field3d(champ);
free_field3d(dpij);
/*free_hessien3d(d2pij);*/
  return(1);
}


/*******************************************************************************
**     hessien_base_ICP_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_base_ICP_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, 
hessien3d *hessienI2, hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;

	grphic3d *mask=NULL;
	mask=I2; 

 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 	if(I2->mri[i+bx0][j+by0][k+bz0]>0)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
        
		  temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*(g[auxi]*g[auxj]+val*Hval.H[auxi][auxj]);
								
				}	
				
				}		
			}
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,mask,Slpqr);
					
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     hessien_base_ICP_sym_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/

int hessien_base_ICP_sym_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double temp,val,g[3];
 double dx,dy,dz;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataG,***maskdataG;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hval;
 mat3d ***dataH,***maskdataH;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;

	grphic3d *mask=NULL;
	mask=I2; 


 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 maskdataG=maskgradI2->raw;
 maskdataH=maskhessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient avec l'image compl�mentaire
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 	if(I2->mask->mri[i+bx0][j+by0][k+bz0]>0)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		// Evaluation de Imreca en s=s+u avec interrpolation lineaire	
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=maskdataH[i2][j2][k2];
        hb=maskdataH[i2][j2][k3];
        hc=maskdataH[i2][j3][k2];
        hd=maskdataH[i2][j3][k3];
        he=maskdataH[i3][j2][k2];
        hf=maskdataH[i3][j2][k3];
        hg=maskdataH[i3][j3][k2];
        hh=maskdataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
        
		  temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*(g[auxi]*g[auxj]+val*Hval.H[auxi][auxj]);
								
				}	
				
				}		
			}
 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 	if(I2->mri[i+bx0][j+by0][k+bz0]>0)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
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

        g[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        g[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        g[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
		
		// Evaluation du hessien de l'image en s=s+u avec interpolation lineaire
		    ha=dataH[i2][j2][k2];
        hb=dataH[i2][j2][k3];
        hc=dataH[i2][j3][k2];
        hd=dataH[i2][j3][k3];
        he=dataH[i3][j2][k2];
        hf=dataH[i3][j2][k3];
        hg=dataH[i3][j3][k2];
        hh=dataH[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hval.H[auxi][auxj]=(float)tmp;
					}
        
		  temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*(g[auxi]*g[auxj]+val*Hval.H[auxi][auxj]);
								
				}	
				
				}		
			}
 
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,mask,Slpqr);
					
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}


/*******************************************************************************
**     imx_simul_atrophie_topo_3d_p                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_atrophie_jacobien_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double val,g[3];
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 vector3d ***dataG;
 mat3d ***dataH;
 int auxi,auxj;
 double tmp;
 double J=1,dJx,dJy,dJz,xt,yt,zt,dJ[3];
 int aux0,aux1;  								// variables auxiliaires
 double auxdbl[20],aux2,mymono[21];
 double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
 int ideb,ifin,jdeb,jfin,kdeb,kfin;
 double dir_desR3[3];
 int resol;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
 resol=BASE3D.resol;
 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
 bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
 topDx=bx1-bx0;
 topDy=by1-by0;
 topDz=bz1-bz0;


 		
 
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

 
 //on rempli le vecteur gradient
	 
for (i=ideb;i<ifin;i++)
for (j=jdeb;j<jfin;j++)
for (k=kdeb;k<kfin;k++)
	  if (I2->mri[i][j][k]>0)
		{
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
    	val=0;
		
		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);

		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
	
		dJ[0]=1.0*dJx/wdth;
		dJ[1]=1.0*dJy/hght;
		dJ[2]=1.0*dJz/dpth;
	
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{
hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+dJ[auxi]*dJ[auxj];
				}
			}
			
   }
tmp=(topDx*topDy*topDz); 
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=2.0*hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
 
  return(1);
}

/*******************************************************************************
**     hessien_atrophie_log_jacobien_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_atrophie_log_jacobien_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2,
hessien3d *hessienI2,hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double val,g[3];
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 vector3d ***dataG;
 mat3d ***dataH;
 int auxi,auxj;
 double tmp;
 double J=1,dJx,dJy,dJz,xt,yt,zt,dJ[3];
 int aux0,aux1;  								// variables auxiliaires
 double auxdbl[20],aux2,mymono[21];
 double auxdbl21a[21], ordonnee[21],pente_x[21],pente_y[21],pente_z[21];                      // variables auxiliaire
 int ideb,ifin,jdeb,jfin,kdeb,kfin;
 double dir_desR3[3],J2,auxE;
 int resol;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataG=gradI2->raw;
 dataH=hessienI2->raw;
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
 resol=BASE3D.resol;
 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
 l=TOP_conv_ind(topi,topj,topk,nb_param)/3; 
 
 // calcul du champ de deformation
 bx0=x0[l];bx1=x1[l];by0=y0[l];by1=y1[l];bz0=z0[l];bz1=z1[l];
 topDx=bx1-bx0;
 topDy=by1-by0;
 topDz=bz1-bz0;


 
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

 
 //on rempli le vecteur gradient
	 
for (i=ideb;i<ifin;i++)
for (j=jdeb;j<jfin;j++)
for (k=kdeb;k<kfin;k++)
   if (I2->mri[i][j][k]>0)
	 {
	   // Evaluation du gradient de l'image en s=s+u avec interrpolation lineaire
	    g[0]=g[1]=g[2]=0;
    	val=0;
		
		xt=(double)(i)/wdth;yt=(double)(j)/hght;zt=(double)(k)/dpth;
		fast_eval_fun(ordonnee, xt, yt, zt, &J, mymono);

		fast_eval_fun(pente_x, xt, yt, zt, &dJx, mymono);
		fast_eval_fun(pente_y, xt, yt, zt, &dJy, mymono);
		fast_eval_fun(pente_z, xt, yt, zt, &dJz, mymono);
	
		dJ[0]=1.0*dJx/wdth;
		dJ[1]=1.0*dJy/hght;
		dJ[2]=1.0*dJz/dpth;
	
	
		auxE=log(J/(double)I2->mri[i][j][k]/I2->rcoeff);
		J2=J*J;	
	
		if ((isnan(fabs(auxE)))||(isinf(fabs(auxE))))
				auxE=0;

		auxE=1.0*(1.0-auxE);
	
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{
				hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+dJ[auxi]*dJ[auxj]/J2*auxE;
				if (isnan(hessien->H[auxi][auxj]))
						printf("Oups, not a number en %d %d %d\n",topi,topj,topk);
				
			}
		}
			
   }
tmp=(topDx*topDy*topDz); 
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=2.0*hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
 
  return(1);
}


/*******************************************************************************
**     hessien_quad_locale_symetrique_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_quad_locale_symetrique_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,valimreca,valimref,greca[3],gref[3];
 double dxreca,dyreca,dzreca,dxref,dyref,dzref;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataGreca,***dataGref;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hvalreca,Hvalref;
 mat3d ***dataHreca,***dataHref;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
double lambda_ref2=lambda_ref*lambda_ref;
 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataGreca=gradI2->raw;
 dataHreca=hessienI2->raw;
 
 dataGref=maskgradI2->raw;
 dataHref=maskhessienI2->raw;
 
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   
		 // Evaluation du gradient de l'image Imreca en s=s+u avec interrpolation lineaire

    dxreca=(double)i+data[i][j][k].x+bx0;
    dyreca=(double)j+data[i][j][k].y+by0;
    dzreca=(double)k+data[i][j][k].z+bz0;

    	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
	
        va=dataGreca[i2][j2][k2];
        vb=dataGreca[i2][j2][k3];
        vc=dataGreca[i2][j3][k2];
        vd=dataGreca[i2][j3][k3];
        ve=dataGreca[i3][j2][k2];
        vf=dataGreca[i3][j2][k3];
        vg=dataGreca[i3][j3][k2];
        vh=dataGreca[i3][j3][k3];

        greca[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        greca[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        greca[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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

        valimreca=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		// Evaluation du hessien de l'image imreca en s=s+u avec interpolation lineaire
		    ha=dataHreca[i2][j2][k2];
        hb=dataHreca[i2][j2][k3];
        hc=dataHreca[i2][j3][k2];
        hd=dataHreca[i2][j3][k3];
        he=dataHreca[i3][j2][k2];
        hf=dataHreca[i3][j2][k3];
        hg=dataHreca[i3][j3][k2];
        hh=dataHreca[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hvalreca.H[auxi][auxj]=(float)tmp;
					}
     
		 	
		 
		 // Evaluation du gradient de l'image imref en s=s-u avec interrpolation lineaire

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

        gref[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gref[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gref[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		// Evaluation de Imref en s=s-u avec interrpolation lineaire	
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
		
		// Evaluation du hessien de l'image Imref en s=s-u avec interpolation lineaire
		    ha=dataHref[i2][j2][k2];
        hb=dataHref[i2][j2][k3];
        hc=dataHref[i2][j3][k2];
        hd=dataHref[i2][j3][k3];
        he=dataHref[i3][j2][k2];
        hf=dataHref[i3][j2][k3];
        hg=dataHref[i3][j3][k2];
        hh=dataHref[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hvalref.H[auxi][auxj]=(float)tmp;
					}
     
		    
		  diff=(double)(valimref-valimreca);
      temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*((lambda_ref*gref[auxi]+greca[auxi])*(lambda_ref*gref[auxj]+greca[auxj])
						+diff*(lambda_ref2*Hvalref.H[auxi][auxj]-Hvalreca.H[auxi][auxj]));	
				}	
				
		}}		
	}
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     hessien_quad_locale_symetrique_coupe_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_quad_locale_symetrique_coupe_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double diff,temp,valimreca,valimref,greca[3],gref[3];
 double dxreca,dyreca,dzreca,dxref,dyref,dzref;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data,***dataGreca,***dataGref;
 field3d *champ;
 vector3d va,vb,vc,vd,ve,vf,vg,vh;
 double wa,wb,wc,wd,we,wf,wg,wh;
 mat3d ha,hb,hc,hd,he,hf,hg,hh,Hvalreca,Hvalref;
 mat3d ***dataHreca,***dataHref;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double tmp;
 int Ntot=0;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
double lambda_ref2=lambda_ref*lambda_ref;


 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 dataGreca=gradI2->raw;
 dataHreca=hessienI2->raw;
 
 dataGref=maskgradI2->raw;
 dataHref=maskhessienI2->raw;
 
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	   
		 // Evaluation du gradient de l'image Imreca en s=s+u avec interrpolation lineaire
 	dxreca=(double)i+data[i][j][k].x+bx0;
    	dyreca=(double)j+data[i][j][k].y+by0;
    	dzreca=(double)k+data[i][j][k].z+bz0;
 
     	i2=(int)dxreca;j2=(int)dyreca;k2=(int)dzreca;
    	i3=i2+1;j3=j2+1;k3=k2+1;
    
      	if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
        {
        xrel=dxreca-(double)i2;yrel=dyreca-(double)j2;zrel=dzreca-(double)k2;
	
        va=dataGreca[i2][j2][k2];
        vb=dataGreca[i2][j2][k3];
        vc=dataGreca[i2][j3][k2];
        vd=dataGreca[i2][j3][k3];
        ve=dataGreca[i3][j2][k2];
        vf=dataGreca[i3][j2][k3];
        vg=dataGreca[i3][j3][k2];
        vh=dataGreca[i3][j3][k3];

        greca[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        greca[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        greca[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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

        valimreca=(we-wa)*xrel+(wc-wa)*yrel+(wb-wa)*zrel+(wg-we-wc+wa)*xrel*yrel
		+(wd-wc-wb+wa)*yrel*zrel+(wf-we-wb+wa)*xrel*zrel
		+(wh+we+wc+wb-wg-wf-wd-wa)*xrel*yrel*zrel+wa;
		
		// Evaluation du hessien de l'image imreca en s=s+u avec interpolation lineaire
		    ha=dataHreca[i2][j2][k2];
        hb=dataHreca[i2][j2][k3];
        hc=dataHreca[i2][j3][k2];
        hd=dataHreca[i2][j3][k3];
        he=dataHreca[i3][j2][k2];
        hf=dataHreca[i3][j2][k3];
        hg=dataHreca[i3][j3][k2];
        hh=dataHreca[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hvalreca.H[auxi][auxj]=(float)tmp;
					}
     
		 	
		 
		 // Evaluation du gradient de l'image imref en s=s-u avec interrpolation lineaire
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

        gref[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
		+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
		+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		
        gref[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
		+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
		+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
		
        gref[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
		+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
		+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

		// Evaluation de Imref en s=s-u avec interrpolation lineaire	
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
		
		// Evaluation du hessien de l'image Imref en s=s-u avec interpolation lineaire
		    ha=dataHref[i2][j2][k2];
        hb=dataHref[i2][j2][k3];
        hc=dataHref[i2][j3][k2];
        hd=dataHref[i2][j3][k3];
        he=dataHref[i3][j2][k2];
        hf=dataHref[i3][j2][k3];
        hg=dataHref[i3][j3][k2];
        hh=dataHref[i3][j3][k3];
		
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
					tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
					tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
					tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
					tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
					tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
					tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
					Hvalref.H[auxi][auxj]=(float)tmp;
					}
     
		    
		  if ((valimreca>0)&&(valimref>0))		
			{
			diff=(double)(valimref-valimreca);
      temp=fx[i]*fy[j]*fz[k];
			temp=temp*temp;
    
		for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+2.0*temp*((lambda_ref*gref[auxi]+greca[auxi])*(lambda_ref*gref[auxj]+greca[auxj])
						+diff*(lambda_ref2*Hvalref.H[auxi][auxj]-Hvalreca.H[auxi][auxj]));		
				}	
			Ntot++;
			}
			
			
		}}		
	}
   
if (Ntot>0)
{tmp=1.0*(topDx*topDy*topDz)*Ntot;
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;
}
else
	{
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=0.0;
	}
 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     hessien_quad_sous_champ_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_quad_sous_champ_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
 int 	i,j,k,l,wdth,hght,dpth;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int	*x0,*y0,*z0,*x1,*y1,*z1;
 double temp;
 int ux,uy,uz;
 int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
 int x00,x11,y00,y11,z00,z11;
 vector3d ***data;
 field3d *champ;
 double f,px,py,pz;
 int auxi,auxj,auxk;
 int i2,j2,k2;
 int Ntot=0;


 wdth=I2->width;hght=I2->height;dpth=I2->depth;
 
 
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

 D=TOP_D(nb_param);
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


 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=0.0;

 
 //on rempli le vecteur gradient
	 
 for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 {
	i2=i+bx0;
	j2=j+by0;
	k2=k+bz0;
		
  
      	temp=fx[i]*fy[j]*fz[k];
	temp=2*temp*temp;
	Ntot++;
	
	for (auxi=0;auxi<3;auxi++)
		hessien->H[auxi][auxi]+=temp;
	
	}
	
   
/*if (Ntot>0)
{tmp=1.0*(topDx*topDy*topDz)*Ntot;
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;
}
else
	{
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=0.0;
	}
 
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 */

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}



/*******************************************************************************
**     hessien_groupwise_variance_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_groupwise_variance_locale_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
int i,j,k,l,width,height,depth;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int	*x0,*y0,*z0,*x1,*y1,*z1;
double temp,valimreca,valimref,greca[3],gref[3];
double dxreca,dyreca,dzreca,dxref,dyref,dzref;
int ux,uy,uz;
int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
int x00,x11,y00,y11,z00,z11;
vector3d ***data,***Greca,***Gref;
field3d *champ;
vector3d va,vb,vc,vd,ve,vf,vg,vh,vnull;
double wa,wb,wc,wd,we,wf,wg,wh;
mat3d ha,hb,hc,hd,he,hf,hg,hh,Hvalreca,Hvalref,Hvalref2,hnull;
mat3d ***Hreca,***Href, ***Href2;
double f,px,py,pz;
int auxi,auxj,auxk;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double tmp,du;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
double lambda_ref2=lambda_ref*lambda_ref;
// vecteur nul
vnull.x=vnull.y=vnull.z=0.0;	

// hessien nul
for (auxi=0;auxi<3;auxi++)
for (auxj=0;auxj<3;auxj++)
	hnull.H[auxi][auxj]=0.0;
		
width=I2->width;height=I2->height;depth=I2->depth;
 
Greca=_ParamRecalageBspline.grad_reca_groupwise->raw;
Hreca=_ParamRecalageBspline.hessien_reca_groupwise->raw;
 
Gref=_ParamRecalageBspline.grad_ref_groupwise->raw;
Href=_ParamRecalageBspline.hessien_ref_groupwise->raw;
 
Href2=_ParamRecalageBspline.hessien_ref2_groupwise->raw;
 
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

D=TOP_D(nb_param);
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

 
//on rempli le vecteur gradient
	 
for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 	{
		   
		/* calcul de la valeur interpol�e dans reca */
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
		        
        	greca[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
			+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
			+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
			greca[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
			+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
			+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
        	greca[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
			
			/* Pour interpoler le hessien */	
			if (i2==-1 || j2==-1 || k2==-1) 			ha=hnull; else ha=Hreca[i2][j2][k2];
        	if (i2==-1 || j2==-1 || k3==depth) 			hb=hnull; else hb=Hreca[i2][j2][k3];
        	if (i2==-1 || j3==height || k2==-1) 		hc=hnull; else hc=Hreca[i2][j3][k2];
        	if (i2==-1 || j3==height || k3==depth) 		hd=hnull; else hd=Hreca[i2][j3][k3];
        	if (i3==width || j2==-1 || k2==-1) 			he=hnull; else he=Hreca[i3][j2][k2];
        	if (i3==width || j2==-1 || k3==depth) 		hf=hnull; else hf=Hreca[i3][j2][k3];
        	if (i3==width || j3==height || k2==-1) 		hg=hnull; else hg=Hreca[i3][j3][k2];
        	if (i3==width || j3==height || k3==depth) 	hh=hnull; else hh=Hreca[i3][j3][k3];
			
			for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
				tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
				tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
				tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
				tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
				tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
				tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
				Hvalreca.H[auxi][auxj]=(float)tmp;
				}
			}
		else
			{
		 	free_field3d(champ);
			*hessien=hnull;
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
	
        	gref[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
			+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
			+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
		   	gref[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
			+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
			+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
			gref[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
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
			
		
			/* Pour interpoler le hessien ref */	
			if (i2==-1 || j2==-1 || k2==-1) 			ha=hnull; else ha=Href[i2][j2][k2];
        	if (i2==-1 || j2==-1 || k3==depth) 			hb=hnull; else hb=Href[i2][j2][k3];
        	if (i2==-1 || j3==height || k2==-1) 		hc=hnull; else hc=Href[i2][j3][k2];
        	if (i2==-1 || j3==height || k3==depth) 		hd=hnull; else hd=Href[i2][j3][k3];
        	if (i3==width || j2==-1 || k2==-1) 			he=hnull; else he=Href[i3][j2][k2];
        	if (i3==width || j2==-1 || k3==depth) 		hf=hnull; else hf=Href[i3][j2][k3];
        	if (i3==width || j3==height || k2==-1) 		hg=hnull; else hg=Href[i3][j3][k2];
        	if (i3==width || j3==height || k3==depth) 	hh=hnull; else hh=Href[i3][j3][k3];
			
			for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
				tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
				tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
				tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
				tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
				tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
				tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
				Hvalref.H[auxi][auxj]=(float)tmp;
				}
	
			/* Pour interpoler le hessien ref2 */	
			if (i2==-1 || j2==-1 || k2==-1) 			ha=hnull; else ha=Href2[i2][j2][k2];
        	if (i2==-1 || j2==-1 || k3==depth) 			hb=hnull; else hb=Href2[i2][j2][k3];
        	if (i2==-1 || j3==height || k2==-1) 		hc=hnull; else hc=Href2[i2][j3][k2];
        	if (i2==-1 || j3==height || k3==depth) 		hd=hnull; else hd=Href2[i2][j3][k3];
        	if (i3==width || j2==-1 || k2==-1) 			he=hnull; else he=Href2[i3][j2][k2];
        	if (i3==width || j2==-1 || k3==depth) 		hf=hnull; else hf=Href2[i3][j2][k3];
        	if (i3==width || j3==height || k2==-1) 		hg=hnull; else hg=Href2[i3][j3][k2];
        	if (i3==width || j3==height || k3==depth) 	hh=hnull; else hh=Href2[i3][j3][k3];
			
			for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
				tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
				tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
				tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
				tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
				tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
				tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
				Hvalref2.H[auxi][auxj]=(float)tmp;
				}
			}
		else
			{
			free_field3d(champ);
			*hessien=hnull;
			return(0);
			}
	
	valimreca=valimreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
	valimref=valimref*_ParamRecalageBspline.ref_groupwise->rcoeff;
	du=fx[i]*fy[j]*fz[k];
	du=du*du/_ParamRecalageBspline.nb_tot_groupwise;				
	temp=2.0*(valimref+valimreca)/_ParamRecalageBspline.nb_tot_groupwise;
	 
	for (auxi=0;auxi<3;auxi++)
		for (auxj=auxi;auxj<3;auxj++)
		  	{
			hessien->H[auxi][auxj]=hessien->H[auxi][auxj]
			+du*(
			2.0*(greca[auxi]*greca[auxj]+Hvalreca.H[auxi][auxj]*valimreca)
			+Hvalref2.H[auxi][auxj]*lambda_ref2
			-2.0*(greca[auxi]-gref[auxi]*lambda_ref)*(greca[auxj]-gref[auxj]*lambda_ref)/_ParamRecalageBspline.nb_tot_groupwise
			-1.0*temp*(Hvalreca.H[auxi][auxj]+Hvalref.H[auxi][auxj]*lambda_ref2)
			);
			}						
	}
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

	
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}



/*******************************************************************************
**     hessien_groupwise_variance_locale_nonsym_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_groupwise_variance_locale_nonsym_3d(grphic3d *I2, grphic3d *Im, int nb_param, double *param, double *param_norm, field3d *gradI2, field3d *maskgradI2, hessien3d *hessienI2,
hessien3d *maskhessienI2,mat3d *hessien, int topi, int topj, int topk,TSlpqr *Slpqr, reg_hessien_locale_t hessien_reg)
{
int i,j,k,l,width,height,depth;
double *fx,*fy,*fz,*dfx,*dfy,*dfz;
int	*x0,*y0,*z0,*x1,*y1,*z1;
double temp,valimreca,valimref,greca[3];
double dxreca,dyreca,dzreca;
int ux,uy,uz;
int D,bx0,bx1,by0,by1,bz0,bz1,topDx,topDy,topDz;
int x00,x11,y00,y11,z00,z11;
vector3d ***data;
field3d *champ;
vector3d va,vb,vc,vd,ve,vf,vg,vh,vnull;
double wa,wb,wc,wd,we,wf,wg,wh;
mat3d ha,hb,hc,hd,he,hf,hg,hh,Hvalreca,hnull;
double f,px,py,pz;
int auxi,auxj,auxk;
int i2,j2,k2,i3,j3,k3;
double xrel,yrel,zrel;
double tmp,du;
double Nmoins1=_ParamRecalageBspline.nb_tot_groupwise-1.0;

// vecteur nul
vnull.x=vnull.y=vnull.z=0.0;	

// hessien nul
for (auxi=0;auxi<3;auxi++)
for (auxj=0;auxj<3;auxj++)
	hnull.H[auxi][auxj]=0.0;
		
width=I2->width;height=I2->height;depth=I2->depth;
 
 
x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;

 
 
 //mise a zero du hessien
for (i=0;i<3;i++)
for (j=0;j<3;j++)
	hessien->H[i][j]=0.0;

D=TOP_D(nb_param);
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

 
//on rempli le vecteur gradient
	 
for (i=0;i<topDx;i++)
	for (j=0;j<topDy;j++)
		for (k=0;k<topDz;k++)
	 	{
		   
		/* calcul de la valeur interpol�e dans reca */
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
		        
        	greca[0]=(ve.x-va.x)*xrel+(vc.x-va.x)*yrel+(vb.x-va.x)*zrel+(vg.x-ve.x-vc.x+va.x)*xrel*yrel
			+(vd.x-vc.x-vb.x+va.x)*yrel*zrel+(vf.x-ve.x-vb.x+va.x)*xrel*zrel
			+(vh.x+ve.x+vc.x+vb.x-vg.x-vf.x-vd.x-va.x)*xrel*yrel*zrel+va.x;
		
			greca[0]=greca[0]*_ParamRecalageBspline.reca_groupwise->rcoeff;
		
			greca[1]=(ve.y-va.y)*xrel+(vc.y-va.y)*yrel+(vb.y-va.y)*zrel+(vg.y-ve.y-vc.y+va.y)*xrel*yrel
			+(vd.y-vc.y-vb.y+va.y)*yrel*zrel+(vf.y-ve.y-vb.y+va.y)*xrel*zrel
			+(vh.y+ve.y+vc.y+vb.y-vg.y-vf.y-vd.y-va.y)*xrel*yrel*zrel+va.y;
		
			greca[1]=greca[1]*_ParamRecalageBspline.reca_groupwise->rcoeff;
		
        	greca[2]=(ve.z-va.z)*xrel+(vc.z-va.z)*yrel+(vb.z-va.z)*zrel+(vg.z-ve.z-vc.z+va.z)*xrel*yrel
			+(vd.z-vc.z-vb.z+va.z)*yrel*zrel+(vf.z-ve.z-vb.z+va.z)*xrel*zrel
			+(vh.z+ve.z+vc.z+vb.z-vg.z-vf.z-vd.z-va.z)*xrel*yrel*zrel+va.z;

			greca[2]=greca[2]*_ParamRecalageBspline.reca_groupwise->rcoeff;
		
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
			
			/* Pour interpoler le hessien */	
			if (i2==-1 || j2==-1 || k2==-1) 			ha=hnull; else ha=hessienI2->raw[i2][j2][k2];
        	if (i2==-1 || j2==-1 || k3==depth) 			hb=hnull; else hb=hessienI2->raw[i2][j2][k3];
        	if (i2==-1 || j3==height || k2==-1) 		hc=hnull; else hc=hessienI2->raw[i2][j3][k2];
        	if (i2==-1 || j3==height || k3==depth) 		hd=hnull; else hd=hessienI2->raw[i2][j3][k3];
        	if (i3==width || j2==-1 || k2==-1) 			he=hnull; else he=hessienI2->raw[i3][j2][k2];
        	if (i3==width || j2==-1 || k3==depth) 		hf=hnull; else hf=hessienI2->raw[i3][j2][k3];
        	if (i3==width || j3==height || k2==-1) 		hg=hnull; else hg=hessienI2->raw[i3][j3][k2];
        	if (i3==width || j3==height || k3==depth) 	hh=hnull; else hh=hessienI2->raw[i3][j3][k3];
			
			for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
			  	{
				tmp=(double)(he.H[auxi][auxj]-ha.H[auxi][auxj])*xrel;
				tmp=tmp+(double)(hc.H[auxi][auxj]-ha.H[auxi][auxj])*yrel+(hb.H[auxi][auxj]-ha.H[auxi][auxj])*zrel;
				tmp=tmp+(double)(hg.H[auxi][auxj]-he.H[auxi][auxj]-hc.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*yrel;
				tmp=tmp+(double)(hd.H[auxi][auxj]-hc.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*yrel*zrel;
				tmp=tmp+(double)(hf.H[auxi][auxj]-he.H[auxi][auxj]-hb.H[auxi][auxj]+ha.H[auxi][auxj])*xrel*zrel;
				tmp=tmp+(double)(hh.H[auxi][auxj]+he.H[auxi][auxj]+hc.H[auxi][auxj]+hb.H[auxi][auxj]-hg.H[auxi][auxj]-hf.H[auxi][auxj]-hd.H[auxi][auxj]-ha.H[auxi][auxj])*xrel*yrel*zrel+ha.H[auxi][auxj];
				Hvalreca.H[auxi][auxj]=(float)tmp*_ParamRecalageBspline.reca_groupwise->rcoeff;
				}
			}
		else
			{
		 	free_field3d(champ);
			*hessien=hnull;
			return(0);
			}
	
	
	valimref=_ParamRecalageBspline.ref_groupwise->mri[i+bx0][j+by0][k+bz0];
  
	valimreca=valimreca*_ParamRecalageBspline.reca_groupwise->rcoeff;
	valimref=valimref*_ParamRecalageBspline.ref_groupwise->rcoeff;
	du=fx[i]*fy[j]*fz[k];
	du=2.0*du*du/_ParamRecalageBspline.nb_tot_groupwise/_ParamRecalageBspline.nb_tot_groupwise;				
	
	temp=(Nmoins1*valimreca-valimref);
	 
	for (auxi=0;auxi<3;auxi++)
		for (auxj=auxi;auxj<3;auxj++)
		  	{
			hessien->H[auxi][auxj]=hessien->H[auxi][auxj]
									+du*(Hvalreca.H[auxi][auxj]*temp+(Nmoins1*greca[auxi]*greca[auxj]));
			}						
	}
   
tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
	for (j=i;j<3;j++)
		hessien->H[i][j]=hessien->H[i][j]/tmp;

	
 hessien->H[1][0]=hessien->H[0][1];
 hessien->H[2][0]=hessien->H[0][2];
 hessien->H[2][1]=hessien->H[1][2];
 

 
  	
	if (TOPO_REGULARISATION_MEMBRANE_ELASTIQUE!=0)
			{
			mat3d H_reg;
			hessien_reg(nb_param, param, param_norm, &H_reg,  topi,  topj,  topk,NULL,Slpqr);
							
		  for (i=0;i<3;i++)
			hessien->H[i][i]+=TOPO_REGULARISATION_MEMBRANE_ELASTIQUE_NORMALISEE*TOPO_REGULARISATION_MEMBRANE_ELASTIQUE*H_reg.H[i][i];
			}
	
  free_field3d(champ);


  return(1);
}

/*******************************************************************************
**     hessien_base_quad_locale_3d                  
**                                                                   
**     Calcul du hessien de la transfo sur la base de fonctions                     
**     I2: image final = imref
**     Im: image transforme de I1 (Imreca) a l'iteration k
**     gradI2= gradient de l'image I2	
**     grad: vecteur 1D contenant le gradient de l'energie  
*******************************************************************************/
int hessien_regularisation_patch_imagebased_local(int nb_param, double *param, double *param_norm, mat3d *hessien, int topi, int topj, int topk,grphic3d *mask, TSlpqr *Slpqr)
{
int resol,width,height,depth,topD,l,topDx,topDy,topDz;
int *x0,*x1,*y0,*y1,*z0,*z1;
int bx0,bx1,by0,by1,bz0,bz1;
double *fx,*fy,*fz,*dfx,*dfy,*dfz,px,py,pz;
int i,j,k,di,dj,dk;
int ux,uy,uz;
int ideb,jdeb,kdeb, ifin,jfin,kfin;
double ***poids, tmp;
int nhx,nhy,nhz,spx,spy,spz,auxi,auxj,auxk;
double du[3], dut[3];

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


l=TOP_conv_ind(topi,topj,topk,nb_param)/3;
px=param[3*l];py=param[3*l+1];pz=param[3*l+2];
   		

ideb=topDx/2;
jdeb=topDy/2;
kdeb=topDz/2;

for (i=0;i<3;i++)
for (j=0;j<3;j++)
hessien->H[i][j]=0;




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

	du[0]=dfx[ux]*fy[uy]*fz[uz];
	du[1]=fx[ux]*dfy[uy]*fz[uz];
	du[2]=fx[ux]*fy[uy]*dfz[uz];
	
	
	    	
 	NLMWeightOnePoint3D(_ParamRecalageBspline.ref_groupwise, poids, TOPO_PATCH_COEFF,  
					TOPO_NLMhwnx,  TOPO_NLMhwny,  TOPO_NLMhwnz,  TOPO_NLMhwvsx,  TOPO_NLMhwvsy,  TOPO_NLMhwvsz,auxi, auxj, auxk);
		
		for (i=ideb; i<ifin; i++)
		for (j=jdeb; j<jfin; j++)
		for (k=kdeb; k<kfin; k++)
			{
			
			di = i-bx0;
			dj = j-by0;
			dk = k-bz0;
			
			
						
			tmp=poids[di-ux+TOPO_NLMhwnx][dj-uy+TOPO_NLMhwny][dk-uz+TOPO_NLMhwnz];
			
			if ((di>=0)&&(di<topDx)&&(dj>=0)&&(dj<topDy)&&(dk>=0)&&(dk<topDz))
				{
				dut[0]=dfx[di]*fy[dj]*fz[dk];
				dut[1]=fx[di]*dfy[dj]*fz[dk];
				dut[2]=fx[di]*fy[dj]*dfz[dk];
				}
				else
				{
				dut[0]=dut[1]=dut[2]=0.0;
				}
				
			for (auxi=0;auxi<3;auxi++)
			for (auxj=auxi;auxj<3;auxj++)
		  	{
			hessien->H[auxi][auxj]=hessien->H[auxi][auxj]+tmp*(dut[auxi]-du[auxi])*(dut[auxj]-du[auxj]);
			}						

			}
		
		
		 
	}


tmp=1.0*(topDx*topDy*topDz);
 for (i=0;i<3;i++)
		hessien->H[i][i]=2.0*hessien->H[i][i]/tmp;


 for (i=0;i<3;i++)
		hessien->H[i][i]=hessien->H[i][i]*(pow(1.0*TOPO_PONDERATION_RECALAGE_SYMETRIQUE,2.0)); // Pour assurer la compatibilite avec le recalage symetrique;


free_dmatrix_3d(poids);

return(1);
}

