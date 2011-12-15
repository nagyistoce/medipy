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
#include "noyau/io/imx_log.h"
#include "noyau/gui/imx_picture3d.h"

#ifdef __GTK_GUI
#include "noyau/gui/imx_3dmo.h"
#endif /* __GTK_GUI */

#include "math/imx_matrix.h"
#include "math/imx_bspline.h"	/* pour les fonctions d'interpolation */

#include "noyau/imx_lang.h"

#include "outils/imx_misc.h"

#include "traitement/trai_3d.h"

#include "recalage/distance_3d.h"
#include "recalage/mtch_robust_3d.h"
#include "recalage/chps_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/validation_topo_3d.h"
#include "recalage/topo/divers_topo.h"
#include "recalage/topo/distance_topo.h"




/*******************************************************************************
**       Superposition_damier_3d                                               
**                                                                        
**       g�n�re une superposition d'image en damier                          
*******************************************************************************/
void superposition_damier_3d(void)
{
 int im_ref,im_reca,im_res,step,err;
 
 /*question sur les images a recaler*/
 im_reca=GET_PLACE3D(TEXT0030);
 im_ref= GET_PLACE3D(TEXT0031);
 im_res= GET_PLACE3D(TEXT0006);

 step = GET_INT("pas du damier",5,&err);

superposition_damier_3d_p(im_reca,im_ref,im_res,step);

  
 show_picture_3d(im_res);
}


/*******************************************************************************
**       superposition_damier_3d_p                                               
**                                                                        
**       g�n�re une superposition d'image en damier                          
*******************************************************************************/
void superposition_damier_3d_p(int im_reca,int im_ref,int im_res,int step)
{
grphic3d *imref,*imreca,*imres;
int wdth,hght,dpth,i,j,k,tmp;

imref=ptr_img_3d(im_ref);
imreca=ptr_img_3d(im_reca);
imres=ptr_img_3d(im_res);

wdth=MINI(imref->width,imreca->width);
hght=MINI(imref->height,imreca->height);
dpth=MINI(imref->depth,imreca->depth);

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	{
	tmp=i/step+j/step+k/step;
	if ((tmp%2)==0)
		imres->mri[i][j][k]=imref->mri[i][j][k];
	else
		imres->mri[i][j][k]=imreca->mri[i][j][k];
	}

imx_copie_param_3d_p(imref,imres); 
}



/*******************************************************************************
**        synthese_champ_topo_3d                                               
**                                                                        
**       genere un fichier chp contenant un champ de synthese                          
*******************************************************************************/
void synthese_champ_topo_3d(void)
{
 char nomfichres[500],str[500],*ext;
 int err,width,height,depth,e,id, size_min;
 float dx,dy,dz;
 
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 
 /*test de l'extension*/
  ext=strstr(nomfichres,".chp");
  if (ext!=NULL) *ext=0;
 
 width= GET_INT("Width", 0, &e);
 height= GET_INT("Height", 0, &e);
 depth= GET_INT("Depth", 0, &e);
 
 dx= GET_FLOAT("dx", 1, &e);
 dy= GET_FLOAT("dy", 1, &e);
 dz= GET_FLOAT("dz", 1, &e);
 
 
 id = GET_INT("Identite = 0", 1, &e);
 
 if (id==0)
 id_champ(nomfichres,width,height,depth,dx,dy,dz);
 else
 {
 size_min= GET_INT("Taille mini d'un bloc", 10, &e);
 rand_champ_topo(nomfichres,width,height,depth,dx,dy,dz,size_min);
 }
 return;
}


void rand_champ_topo(char *nomfichres,int width, int height, int depth,float dx, float dy, float dz, int size_min)
{
int x0,y0,z0,i,j,k,ix[3],iy[3],iz[3];
double x,y,z,phi,dim,dimprec;
double xmid,ymid,zmid,i0,j0,k0;
double xm,xM,ym,yM,zm,zM,maxX,minX,maxY,minY,maxZ,minZ,tmp1,tmp2,res,minres;
int err=1;
TOPTliste myL,myLfinal1,myLfinal2;
TOPTbox myB,myBaux;
field3d  *ch;
transf3d *transfo;
vector3d ***data;
char nomfichinv[500];
//srand(1); // initialise le generateur de nombre aleatoire


//----------initialisation de la liste chainee-------------
myL.nelem = 1;
myL.tete = (TOPTbox*) malloc(sizeof(TOPTbox));
myL.tete->next = NULL;
myL.tete->xm = 0;   	myL.tete->xM = width;
myL.tete->ym = 0;   	myL.tete->yM = height;
myL.tete->zm = 0;   	myL.tete->zM = depth;

myLfinal1.nelem = 0; myLfinal1.tete = NULL;


//-------------Partition 1 de l'image---------------
while (myL.nelem>0)
	{
	myB = depile_end(&myL);
	
	if (dim_min_box(myB)>size_min)
		{
		//------on casse en 8 sous-boite---------------
		x0=0.5*(random()*(myB.xM-myB.xm))/RAND_MAX+0.75*myB.xm+0.25*myB.xM;
		y0=0.5*(random()*(myB.yM-myB.ym))/RAND_MAX+0.75*myB.ym+0.25*myB.yM;
		z0=0.5*(random()*(myB.zM-myB.zm))/RAND_MAX+0.75*myB.zm+0.25*myB.zM;
		
		ix[0]=myB.xm; ix[1]=x0; ix[2]=myB.xM;
		iy[0]=myB.ym; iy[1]=y0; iy[2]=myB.yM;
		iz[0]=myB.zm; iz[1]=z0; iz[2]=myB.zM;
		
		for (i=0;i<2;i++)
			for (j=0;j<2;j++)
				for (k=0;k<2;k++)
				{
				myBaux.xm=ix[i]; myBaux.xM=ix[i+1];
				myBaux.ym=iy[j]; myBaux.yM=iy[j+1];
				myBaux.zm=iz[k]; myBaux.zM=iz[k+1];
				empile(&myL,myBaux);
				}				
		}
	else 
		{
		//-----------Tirage al�toire des coeff du champ------------
		/*myB.ax=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.xM-myB.xm)/(3.0*PI));
		myB.ay=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.yM-myB.ym)/(3.0*PI));
		myB.az=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.zM-myB.zm)/(3.0*PI));*/
		res=10;
		while (res>0.9/PI)
			{
				myB.ax=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.xM-myB.xm)/(1.0*PI));
				myB.ay=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.yM-myB.ym)/(1.0*PI));
				myB.az=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.zM-myB.zm)/(1.0*PI));
				res=1.0*fabs(myB.ax)/(myB.xM-myB.xm)+1.0*fabs(myB.ay)/(myB.yM-myB.ym)+1.0*fabs(myB.az)/(myB.zM-myB.zm);
			}
		
		empile(&myLfinal1,myB);
		}
	}

//----------initialisation de la liste chainee-------------
myL.nelem = 1;
myL.tete = (TOPTbox*) malloc(sizeof(TOPTbox));
myL.tete->next = NULL;
myL.tete->xm = 0;   	myL.tete->xM = width;
myL.tete->ym = 0;   	myL.tete->yM = height;
myL.tete->zm = 0;   	myL.tete->zM = depth;

myLfinal2.nelem = 0; myLfinal2.tete = NULL;


//-------------Partition 2 de l'image---------------
while (myL.nelem>0)
	{
	myB = depile_end(&myL);
	
	if (dim_min_box(myB)>size_min)
		{
		//------on casse en 8 sous-boite---------------
		x0=0.5*(random()*(myB.xM-myB.xm))/RAND_MAX+0.75*myB.xm+0.25*myB.xM;
		y0=0.5*(random()*(myB.yM-myB.ym))/RAND_MAX+0.75*myB.ym+0.25*myB.yM;
		z0=0.5*(random()*(myB.zM-myB.zm))/RAND_MAX+0.75*myB.zm+0.25*myB.zM;
		
		ix[0]=myB.xm; ix[1]=x0; ix[2]=myB.xM;
		iy[0]=myB.ym; iy[1]=y0; iy[2]=myB.yM;
		iz[0]=myB.zm; iz[1]=z0; iz[2]=myB.zM;
		
		for (i=0;i<2;i++)
			for (j=0;j<2;j++)
				for (k=0;k<2;k++)
				{
				myBaux.xm=ix[i]; myBaux.xM=ix[i+1];
				myBaux.ym=iy[j]; myBaux.yM=iy[j+1];
				myBaux.zm=iz[k]; myBaux.zM=iz[k+1];
				empile(&myL,myBaux);
				}				
		}
	else 
		{
		//-----------Tirage al�toire des coeff du champ------------
		/*myB.ax=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.xM-myB.xm)/(3.0*PI));
		myB.ay=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.yM-myB.ym)/(3.0*PI));
		myB.az=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.zM-myB.zm)/(3.0*PI));*/
	
		res=10;
		while (res>0.9/PI)
			{
				myB.ax=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.xM-myB.xm)/(1.0*PI));
				myB.ay=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.yM-myB.ym)/(1.0*PI));
				myB.az=(2.0*rand()/RAND_MAX-1.0)*(1.0*(myB.zM-myB.zm)/(1.0*PI));
				res=1.0*fabs(myB.ax)/(myB.xM-myB.xm)+1.0*fabs(myB.ay)/(myB.yM-myB.ym)+1.0*fabs(myB.az)/(myB.zM-myB.zm);
			}
		empile(&myLfinal2,myB);
		}
	}


//-------------generation du champ de deformation par combinaison des 2 transfos----------------

ch=cr_field3d(width,height,depth);
data=ch->raw;
 
for (i=0;i<width;i++)
	for (j=0;j<height;j++)
		for (k=0;k<depth;k++)
		{
		myBaux=find_box((double)i,(double)j,(double)k,myLfinal1,&err);
		
	
		
		if (err==0) { printf("Y a pas de boite qui marche !!!\n"); myBaux.ax=myBaux.ay=myBaux.az=0;}
		phi=1.0*sin(PI*1.0*(i-myBaux.xm)/(myBaux.xM-myBaux.xm))*sin(PI*1.0*(j-myBaux.ym)/(myBaux.yM-myBaux.ym))*sin(PI*1.0*(k-myBaux.zm)/(myBaux.zM-myBaux.zm));
		x=i+myBaux.ax*phi;
		y=j+myBaux.ay*phi;
		z=k+myBaux.az*phi;
		
		myBaux=find_box(x,y,z,myLfinal2,&err);
		
		if (err==0) { printf("Y a pas de boite qui marche !!!\n"); myBaux.ax=myBaux.ay=myBaux.az=0;}
		phi=1.0*sin(PI*1.0*(x-myBaux.xm)/(myBaux.xM-myBaux.xm))*sin(PI*1.0*(y-myBaux.ym)/(myBaux.yM-myBaux.ym))*sin(PI*1.0*(z-myBaux.zm)/(myBaux.zM-myBaux.zm));
		
		
		data[i][j][k].x=1.0*(x+myBaux.ax*phi-i)*dx;
		data[i][j][k].y=1.0*(y+myBaux.ay*phi-j)*dy;
		data[i][j][k].z=1.0*(z+myBaux.az*phi-k)*dz;
		}


transfo=field_to_transf_3d(ch,NULL,NULL);

/* C'est vrai pour la plupart des images IRM ... */
transfo->dx=dx;
transfo->dy=dy;
transfo->dz=dz;

save_transf_3d(transfo,nomfichres);

//--------------------- generation du champ de transformation inverse -----------------
minres=0;
for (i=0;i<width;i++)
	for (j=0;j<height;j++)
		for (k=0;k<depth;k++)
			{
			data[i][j][k].x=0;
			data[i][j][k].y=0;
			data[i][j][k].z=0;


			myBaux=find_box((double)i,(double)j,(double)k,myLfinal2,&err);
			xm=myBaux.xm;xM=myBaux.xM;
			ym=myBaux.ym;yM=myBaux.yM;
			zm=myBaux.zm;zM=myBaux.zM;
			xmid=1.0*(myBaux.xM+myBaux.xm)/2;
			ymid=1.0*(myBaux.yM+myBaux.ym)/2;
			zmid=1.0*(myBaux.zM+myBaux.zm)/2;
		
			
			if (err==0) { printf("Y a pas de boite qui marche pour la transfo inverse !!!\n");}
			else
				{dimprec=0;
					while	(((dim=dim_max_box(myBaux))>0.01)&&(fabs(dim-dimprec)>0.0001))
						{
							dimprec=dim;
							//---- Suivant x ---------
							if ((xmid>myBaux.xm)&&(xmid<myBaux.xM))
									{
									maxX=1;
									tmp1=fabs(1.0*sin(PI*(1.0*(myBaux.xm-xm)/(xM-xm))));
									tmp2=fabs(1.0*sin(PI*(1.0*(myBaux.xM-xm)/(xM-xm))));
									if (tmp1<tmp2)
										minX=tmp1;
									else minX=tmp2;
									}
							else
									{
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.xm-xm)/(xM-xm)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.xM-xm)/(xM-xm)));
									if (tmp1<tmp2)
										{minX=tmp1;maxX=tmp2;}
									else {minX=tmp2;maxX=tmp1;}
									}
							
							//---- Suivant y ---------
							if ((ymid>myBaux.ym)&&(ymid<myBaux.yM))
									{
									maxY=1;
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.ym-ym)/(yM-ym)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.yM-ym)/(yM-ym)));
									if (tmp1<tmp2)
										minY=tmp1;
									else minY=tmp2;
									}
							else
									{
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.ym-ym)/(yM-ym)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.yM-ym)/(yM-ym)));
									if (tmp1<tmp2)
										{minY=tmp1;maxY=tmp2;}
									else {minY=tmp2;maxY=tmp1;}
									}
							
							//---- Suivant z ---------
							if ((zmid>myBaux.zm)&&(zmid<myBaux.zM))
									{
									maxZ=1;
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.zm-zm)/(zM-zm)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.zM-zm)/(zM-zm)));
									if (tmp1<tmp2)
										minZ=tmp1;
									else minZ=tmp2;
									}
							else
									{
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.zm-zm)/(zM-zm)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.zM-zm)/(zM-zm)));
									if (tmp1<tmp2)
										{minZ=tmp1;maxZ=tmp2;}
									else {minZ=tmp2;maxZ=tmp1;}
									}
							
							//------Mise �jour de la boite -------
							
							if (myBaux.ax>0) {myBaux.xm=MAXI(i-myBaux.ax*maxX*maxY*maxZ,xm);myBaux.xM=MINI(i-myBaux.ax*minX*minY*minZ,xM);}
							else {myBaux.xM=MINI(i-myBaux.ax*maxX*maxY*maxZ,xM);myBaux.xm=MAXI(i-myBaux.ax*minX*minY*minZ,xm);}
							
							if (myBaux.ay>0) {myBaux.ym=MAXI(j-myBaux.ay*maxX*maxY*maxZ,ym);myBaux.yM=MINI(j-myBaux.ay*minX*minY*minZ,yM);}
							else {myBaux.yM=MINI(j-myBaux.ay*maxX*maxY*maxZ,yM);myBaux.ym=MAXI(j-myBaux.ay*minX*minY*minZ,ym);}
							
							if (myBaux.az>0) {myBaux.zm=MAXI(k-myBaux.az*maxX*maxY*maxZ,zm);myBaux.zM=MINI(k-myBaux.az*minX*minY*minZ,zM);}
							else {myBaux.zM=MINI(k-myBaux.az*maxX*maxY*maxZ,zM);myBaux.zm=MAXI(k-myBaux.az*minX*minY*minZ,zm);}
									
							

						}
								
				}
			
				i0=1.0*(myBaux.xM+myBaux.xm)/2;
				j0=1.0*(myBaux.yM+myBaux.ym)/2;
				k0=1.0*(myBaux.zM+myBaux.zm)/2;
	
			
			myBaux=find_box(i0,j0,k0,myLfinal1,&err);
			xm=myBaux.xm;xM=myBaux.xM;
			ym=myBaux.ym;yM=myBaux.yM;
			zm=myBaux.zm;zM=myBaux.zM;
			xmid=1.0*(myBaux.xM+myBaux.xm)/2;
			ymid=1.0*(myBaux.yM+myBaux.ym)/2;
			zmid=1.0*(myBaux.zM+myBaux.zm)/2;
			
			if (err==0) { printf("Y a pas de boite qui marche pour la transfo inverse !!!\n");}
			else
				{
					dimprec=0;
					while	(((dim=dim_max_box(myBaux))>0.01)&&(fabs(dim-dimprec)>0.0001))
						{
							dimprec=dim;
							//---- Suivant x ---------
							if ((xmid>myBaux.xm)&&(xmid<myBaux.xM))
									{
									maxX=1;
									tmp1=fabs(1.0*sin(PI*(1.0*(myBaux.xm-xm)/(xM-xm))));
									tmp2=fabs(1.0*sin(PI*(1.0*(myBaux.xM-xm)/(xM-xm))));
									if (tmp1<tmp2)
										minX=tmp1;
									else minX=tmp2;
									}
							else
									{
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.xm-xm)/(xM-xm)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.xM-xm)/(xM-xm)));
									if (tmp1<tmp2)
										{minX=tmp1;maxX=tmp2;}
									else {minX=tmp2;maxX=tmp1;}
									}
							
							//---- Suivant y ---------
							if ((ymid>myBaux.ym)&&(ymid<myBaux.yM))
									{
									maxY=1;
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.ym-ym)/(yM-ym)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.yM-ym)/(yM-ym)));
									if (tmp1<tmp2)
										minY=tmp1;
									else minY=tmp2;
									}
							else
									{
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.ym-ym)/(yM-ym)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.yM-ym)/(yM-ym)));
									if (tmp1<tmp2)
										{minY=tmp1;maxY=tmp2;}
									else {minY=tmp2;maxY=tmp1;}
									}
							
							//---- Suivant z ---------
							if ((zmid>myBaux.zm)&&(zmid<myBaux.zM))
									{
									maxZ=1;
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.zm-zm)/(zM-zm)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.zM-zm)/(zM-zm)));
									if (tmp1<tmp2)
										minZ=tmp1;
									else minZ=tmp2;
									}
							else
									{
									tmp1=fabs(1.0*sin(PI*1.0*(myBaux.zm-zm)/(zM-zm)));
									tmp2=fabs(1.0*sin(PI*1.0*(myBaux.zM-zm)/(zM-zm)));
									if (tmp1<tmp2)
										{minZ=tmp1;maxZ=tmp2;}
									else {minZ=tmp2;maxZ=tmp1;}
									}
							
							//------Mise �jour de la boite -------
							
							if (myBaux.ax>0) {myBaux.xm=MAXI(i0-myBaux.ax*maxX*maxY*maxZ,xm);myBaux.xM=MINI(i0-myBaux.ax*minX*minY*minZ,xM);}
							else {myBaux.xM=MINI(i0-myBaux.ax*maxX*maxY*maxZ,xM);myBaux.xm=MAXI(i0-myBaux.ax*minX*minY*minZ,xm);}
							
							if (myBaux.ay>0) {myBaux.ym=MAXI(j0-myBaux.ay*maxX*maxY*maxZ,ym);myBaux.yM=MINI(j0-myBaux.ay*minX*minY*minZ,yM);}
							else {myBaux.yM=MINI(j0-myBaux.ay*maxX*maxY*maxZ,yM);myBaux.ym=MAXI(j0-myBaux.ay*minX*minY*minZ,ym);}
							
							if (myBaux.az>0) {myBaux.zm=MAXI(k0-myBaux.az*maxX*maxY*maxZ,zm);myBaux.zM=MINI(k0-myBaux.az*minX*minY*minZ,zM);}
							else {myBaux.zM=MINI(k0-myBaux.az*maxX*maxY*maxZ,zM);myBaux.zm=MAXI(k0-myBaux.az*minX*minY*minZ,zm);}
									
							

						}		
				}
				
			data[i][j][k].x=(1.0*(myBaux.xM+myBaux.xm)/2-i)*dx;
			data[i][j][k].y=(1.0*(myBaux.yM+myBaux.ym)/2-j)*dy;
			data[i][j][k].z=(1.0*(myBaux.zM+myBaux.zm)/2-k)*dz;
			printf ("%d %d %d \r",i,j,k);
			res=my_funct(1.0*(myBaux.xM+myBaux.xm)/2, 1.0*(myBaux.yM+myBaux.ym)/2, 1.0*(myBaux.zM+myBaux.zm)/2, myLfinal1, myLfinal2, i, j, k  );
			if (res>minres)
			minres=res;
	    
			}


printf("l'erreur maximale commise en module pour l'inversion du champ est %f \n",minres); 		
			
strcpy(nomfichinv,nomfichres);
strcat(nomfichinv,"_inv");

transfo=field_to_transf_3d(ch,NULL,NULL);
transfo->dx=dx;
transfo->dy=dy;
transfo->dz=dz;
save_transf_3d(transfo,nomfichinv);

free_transf3d(transfo);
free_field3d(ch);


}


void id_champ(char *nomfichres,int width, int height, int depth, float dx, float dy, float dz)
{
int i,j,k;
field3d  *ch;
transf3d *transfo;
vector3d ***data;

ch=cr_field3d(width,height,depth);
data=ch->raw;
 
for (i=0;i<width;i++)
	for (j=0;j<height;j++)
		for (k=0;k<depth;k++)
		data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0.0;

transfo=field_to_transf_3d(ch,NULL,NULL);
transfo->dx=dx;
transfo->dy=dy;
transfo->dz=dz;
save_transf_3d(transfo,nomfichres);

free_transf3d(transfo);
free_field3d(ch);


}


/*******************************************************************************
**        distance_landmark                                                
**                                                                                                   
*******************************************************************************/

void distance_landmark(void)
{
  int im_1, im_2,im_res, err=0;
  char nomfichier[PATH_LEN],nomfichres[500],str[500];
	grphic3d *im1,*im2,*imres;
	field3d *champ;
	transf3d *transfo;
  /*question sur les images*/
  im_1 = GET_PLACE3D(TEXT0030);
  im_2 = GET_PLACE3D(TEXT0031);
  im_res=GET_PLACE3D(TEXT0006);

  /*fichier contenant la transformation*/
  {
    sprintf(nomfichier,"%s",GET_FILE("*.trf",&err));
	
    if (err)
      return ;
	  
	put_file_extension(nomfichier,".trf",nomfichier);

    if (test_fichier_transf_3d(nomfichier) != 1)
    {
      PUT_ERROR("This file contains no 3D transformation");
      return ;
    }
  }

 	sprintf(str,"File ? in %s",_CurrentPath);
 	strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 
	im1=ptr_img_3d(im_1);	
	im2=ptr_img_3d(im_2);
	imres=ptr_img_3d(im_res);
	
	transfo=load_transf_3d(nomfichier);
	
	champ=transf_to_field_3d(transfo,im1,im2);

if ((im1->width!=im2->width)||(im1->height!=im2->height)||(im1->depth!=im2->depth))
	{printf("Pas les meme tailles d'images !!!\n");
 	free_field3d(champ); 
  free_transf3d(transfo);
	return;
	}

	distance_landmark_p(im2,im1,imres,champ,nomfichres);
	show_picture_3d(im_res);
	free_field3d(champ); 
  free_transf3d(transfo);
	
}

/*******************************************************************************
**        distance_landmark_p                                                
**                                                                                                   
*******************************************************************************/

void distance_landmark_p(grphic3d *im2,grphic3d *im1,grphic3d *imres,field3d *champ,char* nomfichres)
{
int ind,nb_landmark,i,j,k,l;
int wdth,hght,dpth;
int ix=-1,iy=-1,iz=-1,xref=-1,yref=-1,zref=-1,Ntot;
double x,y,z,*dist,mean,std;

ind=init_flog(nomfichres);

wdth=im2->width;
hght=im2->height;
dpth=im2->depth;

imx_copie_param_3d_p(im1,imres);
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		{
		imres->mri[i][j][k]=0;
		}

	
// recherche du nombre de landmark
nb_landmark = 0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im2->mri[i][j][k] > nb_landmark) nb_landmark = im2->mri[i][j][k];

dist=CALLOC(nb_landmark,double);

		
for (l=1;l<=nb_landmark;l++)
	{
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		{
		if (im2->mri[i][j][k] == l) {ix=i; iy=j; iz=k;}
		if (im1->mri[i][j][k] == l) {xref=i; yref=j; zref=k;}
		}
		
	if (((ix != 0)||(iy != 0)||(iz != 0))&&((xref != 0)||(yref != 0)||(zref != 0)))
	{
	x=ix+champ->raw[ix][iy][iz].x;
	y=iy+champ->raw[ix][iy][iz].y;
	z=iz+champ->raw[ix][iy][iz].z;
		
	dist[l-1]=sqrt(pow(x-xref,2)+pow(y-yref,2)+pow(z-zref,2));
	aff_flog(ind,"%d \t %f\n",l,dist[l-1]);
	}
	else
		{dist[l-1]=-1;
		aff_flog(ind,"%d \t %f\n",l,dist[l-1]);	
		}
	/*x2=floor(x+0.5);
	y2=floor(y+0.5);
	z2=floor(z+0.5);*/
	
	/** remplissage de l'image resultat illustrant la distance entre 2 landmarks **/
	/*mmax=fabs(x2-xref)+fabs(y2-yref)+fabs(z2-zref);
	dx=x2-xref;
	dy=y2-yref;
	dz=z2-zref;*/
	/*for (m=0;m<=mmax;m++)
		{
		imres->mri[(int)floor(xref+1.0*m*dx/mmax+0.5)][(int)floor(yref+1.0*m*dy/mmax+0.5)][(int)floor(zref+1.0*m*dz/mmax+0.5)]=l;
		}*/
	
	}

//calcul moyenne ecart-type
mean=0;
Ntot=0;
for (l=0;l<nb_landmark;l++)
	if (dist[l]>-1)
	{
	mean=mean+dist[l];
	Ntot++;
	}
mean=mean/Ntot;

std=0;
for (l=0;l<nb_landmark;l++)
	if (dist[l]>-1)
	{
	std=std+(dist[l]-mean)*(dist[l]-mean);
	}
std=sqrt(std/(Ntot-1));

aff_flog(ind,"\n moyenne : %f \n ecart-type : %f\n",mean,std);

free(dist);
end_flog(ind);
}



/*******************************************************************************
**        info_landmark                                                
**                                                                                                   
*******************************************************************************/

void info_landmark(void)
{
  int im_1,err=0,ans;
  grphic3d *im1;
	char nomfichier[PATH_LEN],s[256];
	field3d *champ=NULL;
	transf3d *transfo;
 
	im_1 = GET_PLACE3D(TEXT0030);
 
	im1=ptr_img_3d(im_1);	
	
	
  sprintf(s,"Utilsation d'une transformation?\n");
  ans=GET_EXIST(s,&err);
  if (ans==1)
	 /*fichier contenant la transformation*/
  {
    sprintf(nomfichier,"%s",GET_FILE("*.trf",&err));
	
    if (err)
      return ;
	  
	put_file_extension(nomfichier,".trf",nomfichier);

    if (test_fichier_transf_3d(nomfichier) != 1)
    {
      PUT_ERROR("This file contains no 3D transformation");
      return ;
    }
	transfo=load_transf_3d(nomfichier);
	champ=transf_to_field_3d(transfo,im1,im1);
 	free_transf3d(transfo);
  }
	
	info_landmark_p(im1,champ);
	
	if (champ != NULL)
		free_field3d(champ); 
 
	
}

/*******************************************************************************
**        info_landmark                                                
**                                                                                                   
*******************************************************************************/

void info_landmark_p(grphic3d *im1, field3d *champ)
{
int nb_landmark,wdth,hght,dpth,i,j,k,nb,l,ix,iy,iz;
float itmp,jtmp,ktmp;
wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

// recherche du nombre de landmark
nb_landmark = 0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im1->mri[i][j][k] > nb_landmark) nb_landmark = im1->mri[i][j][k];

printf("Nb landmarks : %d \n",nb_landmark);
	printf("label \t x \t y \t z  \n");

for (l=1;l<=nb_landmark;l++)
	{nb=0;
	ix=iy=iz=-1;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		{
		if (im1->mri[i][j][k] == l) {ix=i; iy=j; iz=k;nb++;}
		}
if ((ix==-1)||(iy==-1)||(iz==-1))
	printf("%d \t -1 \t -1 \t -1  \n",l);
else
	{
	if (champ != NULL)
		{itmp=ix+champ->raw[ix][iy][iz].x;
		jtmp=iy+champ->raw[ix][iy][iz].y;
		ktmp=iz+champ->raw[ix][iy][iz].z;
		}
	else
		{itmp=ix; jtmp=iy; ktmp=iz;}
	printf("%d \t %2f \t %2f \t %2f  \n",l,itmp,jtmp,ktmp);
	}
if (nb>1)
	printf("Attention ce landmark est declare %d fois ...\n",nb);		

	}
	
}


/*******************************************************************************
**        LVV_operator                                                
**                                                                                                   
*******************************************************************************/

void LVV_operator(void)
{
	int im_1,im_res,e;
	grphic3d *im1,*imres;
	float sigma;
	
	im_1 = GET_PLACE3D(TEXT0231);
	im_res=GET_PLACE3D(TEXT0006);
	
	
	im1=ptr_img_3d(im_1);	
	imres=ptr_img_3d(im_res);
	
	sigma = GET_FLOAT("sigma", 1, &e);
 
	LVV_operator_p(im1,imres,sigma);
	
	show_picture_3d(im_res);

}

/*******************************************************************************
**        LVV_operator_p                                                
**                                                                                                   
*******************************************************************************/

void LVV_operator_p(grphic3d *im1, grphic3d *imres,float sigma)
{
int wdth,hght,dpth;
grphic3d *imx,*imy,*imz,*imxx,*imyy,*imzz,*imxy,*imyz,*imxz;
field3d *grad;
int i,j,k;
int method=2,t=5;
int Ix,Iy,Iz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz;
double w2,max;
double tmp,***LVV;
int wnd=3;

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;
imx=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imx);
imy=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imy);
imz=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imz);
imxx=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imxx);
imyy=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imyy);
imzz=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imzz);
imxy=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imxy);
imxz=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imxz);
imyz=cr_grphic3d(im1);
imx_copie_param_3d_p(im1,imyz);


LVV=alloc_dmatrix_3d(wdth,hght,dpth);

grad=cr_field3d(wdth,hght,dpth);
imx_copie_param_3d_p(im1,imres);
imx_gaussian_filter_3d_p(im1, imres, wnd, sigma);

imx_gradient_3d_p(imres,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			imx->mri[i][j][k]=(int)grad->raw[i][j][k].x;
			imy->mri[i][j][k]=(int)grad->raw[i][j][k].y;
			imz->mri[i][j][k]=(int)grad->raw[i][j][k].z;
			}
	

imx_gaussian_filter_3d_p(imx, imx, wnd, sigma);
imx_gaussian_filter_3d_p(imy, imy, wnd, sigma);
imx_gaussian_filter_3d_p(imz, imz, wnd, sigma);
	
imx_gradient_3d_p(imx,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			imxx->mri[i][j][k]=(int)grad->raw[i][j][k].x;
			imxy->mri[i][j][k]=(int)grad->raw[i][j][k].y;
			imxz->mri[i][j][k]=(int)grad->raw[i][j][k].z;
			}
			
imx_gradient_3d_p(imy,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			imyy->mri[i][j][k]=(int)grad->raw[i][j][k].y;
			imyz->mri[i][j][k]=(int)grad->raw[i][j][k].z;
			}

			
imx_gradient_3d_p(imz,grad,method,t);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			imzz->mri[i][j][k]=(int)grad->raw[i][j][k].z;
			}


for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			{
			Ix=imx->mri[i][j][k];
			Iy=imy->mri[i][j][k];
			Iz=imz->mri[i][j][k];
			Ixx=imxx->mri[i][j][k];
			Ixy=imxy->mri[i][j][k];
			Ixz=imxz->mri[i][j][k];
			Iyy=imyy->mri[i][j][k];
			Iyz=imyz->mri[i][j][k];
			Izz=imzz->mri[i][j][k];
			
			w2=2.0*sqrt(Ix*Ix+Iy*Iy+Iz*Iz);
			
			//imres->mri[i][j][k]=(int)(-1.0/w2*((Ix*Ix*(Iyy+Izz)-2.0*Iy*Iz*Iyz)+(Iy*Iy*(Ixx+Izz)-2.0*Ix*Iz*Ixz)+(Iz*Iz*(Iyy+Ixx)-2.0*Ix*Iy*Ixy)));
			tmp=(-1.0/w2*((Ix*Ix*(Iyy+Izz)-2.0*Iy*Iz*Iyz)+(Iy*Iy*(Ixx+Izz)-2.0*Ix*Iz*Ixz)+(Iz*Iz*(Iyy+Ixx)-2.0*Ix*Iy*Ixy)));
			if (isnan(tmp))
				tmp=0;
				
			LVV[i][j][k]=tmp;
			//imres->mri[i][j][k]=(int)tmp;
			//imres->mri[i][j][k]=(int)w2;
			}
			
max = 0;
for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			max=MAXI(max,fabs(LVV[i][j][k]));

printf ("max =%f \n",max);

for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
			imres->mri[i][j][k]=(int)(LVV[i][j][k]*1000.0/max);

imx_inimaxminpixel_3d_p(imres);
			
free_grphic3d(imx);
free_grphic3d(imy);
free_grphic3d(imz);
free_grphic3d(imxx);
free_grphic3d(imyy);
free_grphic3d(imzz);
free_grphic3d(imxy);
free_grphic3d(imxz);
free_grphic3d(imyz);
free_field3d(grad);
free_dmatrix_3d(LVV);

}



/*******************************************************************************
**     verif_topo()                  
**                                                                   
**     A partir des param d'une deformation Bspline verifie si la transformation
**		viole la topologie et pr�cise dans quel bloc 
*******************************************************************************/
void verif_topo(void)
{

  int err;
  char nomfichier[256];
  transf3d *transfo;
  int wdth, hght, dpth, i, j, k, nb_param,topi,topj,topk;
  double Jm=0,JM=HUGE_VAL,minJ;
  double dir_desR3[3];
  scal_func_t scal_func;
  double *param;
  int aux0,aux1;  								// variables auxiliaires
  double auxdbl[20],aux2;
  double auxdbl4[4];
  double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaires
  TSlpqr Slpqr[8];
  int signe;  			      // D : valeur max des indices i, j, k
  int	*x0,*x1,*y0,*y1,*z0,*z1;
  double xm,xM,ym,yM,zm,zM;
  int nb_func,nbfx,nbfy,nbfz;
  int D, echelle;
  int	im_res,max_pix;
  grphic3d *imres;
 
  
  //fichier contenant la transformation
  {
    char *nomoriginal = GET_FILE("*.trf",&err);

    if (err)
      return ;

    strcpy(nomfichier, nomoriginal);
    if (test_fichier_transf_3d(nomfichier) != 1)
    {
      PUT_ERROR("This file contains no 3D transformation");
      return ;
    }
  }

//question sur les images
 im_res=GET_PLACE3D(TEXT0031);
 imres=ptr_img_3d(im_res);

//question sur l'encadrement du Jacobien
Jm=(double) GET_FLOAT("Jmin<1", 0, &err);
JM=(double) GET_FLOAT("Jmax>1", 0, &err);

  //chargement de la transformation
  transfo = load_transf_3d(nomfichier);
  if (! transfo)
  {
    err = 1;
	PUT_ERROR("Probleme lors du chargement de la transformation");
	return;
  }
 
 wdth=transfo->width;hght=transfo->height;dpth=transfo->depth;
 imres->width=wdth;imres->height=hght;imres->depth=dpth;
 imres->rcoeff=1.0;imres->icomp=0.0;

	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		imres->mri[i][j][k]=0;

 
 if (transfo->typetrans != BSPLINE3D)
  {
  	PUT_ERROR("Ce n'est pas une transformation Bspline");
    return ;
  }
  
   switch (transfo->degre)
    {
    case 0: scal_func=Bsplined0;break;
    case 1: scal_func=Bsplined1;break;
    case 2: scal_func=Bsplined2;break;
    default: scal_func=Bsplined0;break;
    }

	echelle=transfo->resol;
	D=(int)pow(2.0,1.0*echelle)-1;
    param = CALLOC(transfo->nb_param,double);  
    
	nb_param = init_base_3d(wdth, hght, dpth, scal_func);
    for (i = BASE3D.resol; i<transfo->resol; i++)
      nb_param = base_resol_up_3d(param,nb_param);
    if (nb_param != transfo->nb_param)
      printf("ERREUR dans transf_to_field_3d\n");

for(i=0;i<nb_param;i++) param[i]=transfo->param[i];


 for (i=0;i<nb_param/3;i++) 
 	{
	param[3*i]=param[3*i]/wdth;
	param[3*i+1]=param[3*i+1]/hght;
	param[3*i+2]=param[3*i+2]/dpth;
	}

 
 	D= TOP_D(nb_param);	
  	nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;echelle=BASE3D.resol;
	x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  
 	dir_desR3[0]=dir_desR3[1]=dir_desR3[2]=0.0;
	printf("Debut de la verification ..... \n");
	
	for (i=0; i<D; i++)
  		for (j=0; j<D; j++)
  			for (k=0; k<D; k++)
  		
	{
  	 	//-------------------------------------
  		// initialisation des boites  ---------
  		//-------------------------------------
       
  		aux0 = TOP_conv_ind(i,j,k,nb_param);
  		xm = 1.0*x0[aux0/3]/wdth; xM = 1.0*x1[aux0/3]/wdth; ym = 1.0*y0[aux0/3]/hght; 
  		yM = 1.0*y1[aux0/3]/hght; zm = 1.0*z0[aux0/3]/dpth; zM =1.0*z1[aux0/3]/dpth;

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

  		TOP_init_slpqr(Slpqr  ,i-1,j-1,k-1,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+1,i-1,j-1,k  ,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+2,i-1,j  ,k-1,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+3,i-1,j  ,k  ,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+4,i  ,j-1,k-1,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+5,i  ,j-1,k  ,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+6,i  ,j  ,k-1,D,nb_param,param);
  		TOP_init_slpqr(Slpqr+7,i  ,j  ,k  ,D,nb_param,param);


	
  		//-------------------------------------
  		// initialisation ind, bdelta, Jm, JM -
  		//-------------------------------------
  		for (aux0=0; aux0<8; aux0++)
  		{ 	Slpqr[aux0].ind = 7 - aux0;
    		Slpqr[aux0].bdeltam = HUGE_VAL;
    		Slpqr[aux0].bdeltaM = HUGE_VAL;
    		Slpqr[aux0].JM = JM;
    		Slpqr[aux0].Jm = Jm;
  		}

  		//-------------------------------------
  		// initialisation des alpha  ----------
  		//-------------------------------------
  		for (aux0=0; aux0<8; aux0++)
  		{ //--- initialisation du terme constant ------------------------------
    	TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,echelle,auxdbl);
  
    	//--- changement de variables pour le polynome ----------------------
    	for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
    	aux2 = pow(2.0,-1.0 * echelle);
		poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,auxdbl21b);
    	for (aux1=0; aux1<20; aux1++) Slpqr[aux0].alpha[aux1].ordorg = auxdbl21b[aux1+1];

		}


  		//-------------------------------------
  		// verification step  -----------------
  		//-------------------------------------
		signe=0;
  		minJ=0;
		for (aux0=0; aux0<8; aux0++) 
  		{ TOP_verification_delta(Slpqr+aux0,0.0,auxdbl4);
    		if (auxdbl4[1]<=0) 
				{signe=signe-1;
				if (auxdbl4[1]<=minJ) minJ=auxdbl4[1];
				}
			else signe=signe+1;
		
  		}
		aux0 = TOP_conv_ind(i,j,k,nb_param);
	if ((signe) != 8)
		{ 
		//printf("Violation de la topologie dans la region x [%d %d] ; y [%d %d] ; z[%d %d] \n",x0[aux0/3],x1[aux0/3],y0[aux0/3],y1[aux0/3],z0[aux0/3],z1[aux0/3]);
		for(topi=x0[aux0/3];topi<x1[aux0/3];topi++)
		for(topj=y0[aux0/3];topj<y1[aux0/3];topj++)
		for(topk=z0[aux0/3];topk<z1[aux0/3];topk++)
			imres->mri[topi][topj][topk]=imres->mri[topi][topj][topk]+floor(-100.0*minJ);
		
		}
	}
 

 
 
	printf("Fin de la verification ..... \n");
	max_pix=0;	
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
		if (imres->mri[i][j][k]>max_pix) max_pix=imres->mri[i][j][k];
		
	imres->max_pixel=imres->cutoff_max=max_pix;
	show_picture_3d(im_res);
  
    end_base_3d();
    free(param);
 
}

/*******************************************************************************
**     verif_poly()                  
**                                                                   
**     A partir des param d'une deformation Bspline verifie si la transformation
**		viole la topologie et pr�cise dans quel bloc 
*******************************************************************************/


void verif_poly(int i,int j,int k,int ind,int echelle,int nb_param,double *param,TSlpqr *Slpqr)
{
  int l;
  int aux0,aux1;  								// variables auxiliaires
  double auxdbl[20],aux2;
  double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaires
  int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
  double xm,xM,ym,yM,zm,zM;
  int width,height,depth;
  
 
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
 
 
  for (l=0;l<nb_param/3;l++)							// on se place sur [0,1]^3
  { param[3*l]   = 1.0*param[3*l]   / width;
	param[3*l+1] = 1.0*param[3*l+1] / height;
	param[3*l+2] = 1.0*param[3*l+2] / depth;
  }
  
  
  aux0 = TOP_conv_ind(i,j,k,nb_param);
  xm=1.0*(double)i/(1.0*pow(2.0,(double)echelle));
  xM=1.0*((double)i+2.0)/(1.0*pow(2.0,(double)echelle));
  ym=1.0*(double)j/(1.0*pow(2.0,(double)echelle));
  yM=1.0*((double)j+2.0)/(1.0*pow(2.0,(double)echelle));
  zm=1.0*(double)k/(1.0*pow(2.0,(double)echelle));
  zM=1.0*((double)k+2.0)/(1.0*pow(2.0,(double)echelle));

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

  TOP_init_slpqr(Slpqr  ,i-1,j-1,k-1,D,nb_param,param);
  TOP_init_slpqr(Slpqr+1,i-1,j-1,k  ,D,nb_param,param);
  TOP_init_slpqr(Slpqr+2,i-1,j  ,k-1,D,nb_param,param);
  TOP_init_slpqr(Slpqr+3,i-1,j  ,k  ,D,nb_param,param);
  TOP_init_slpqr(Slpqr+4,i  ,j-1,k-1,D,nb_param,param);
  TOP_init_slpqr(Slpqr+5,i  ,j-1,k  ,D,nb_param,param);
  TOP_init_slpqr(Slpqr+6,i  ,j  ,k-1,D,nb_param,param);
  TOP_init_slpqr(Slpqr+7,i  ,j  ,k  ,D,nb_param,param) ;


	
  //-------------------------------------
  // initialisation des alpha  ----------
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { //--- initialisation du terme constant ------------------------------
    TOP_poly_init(Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,echelle,auxdbl);
  	
    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
	
    aux2 = pow(2.0,-1.0 * echelle);
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,auxdbl21b);
    
	for (aux1=0; aux1<20; aux1++) Slpqr[aux0].alpha[aux1].ordorg = auxdbl21b[aux1+1];

	for (aux1=0; aux1<20; aux1++) Slpqr[aux0].alpha[aux1].pente = 0;

  }

  for (l=0;l<nb_param/3;l++)							// on se place sur [0,1]^3
  { param[3*l]   = 1.0*param[3*l]   * width;
	param[3*l+1] = 1.0*param[3*l+1] * height;
	param[3*l+2] = 1.0*param[3*l+2] * depth;
  }

return;
}

/*******************************************************************************
**        topo_synth_3d                                                
**                                                                        
**       genere un fichier chp contenant un champ de synthese                          
*******************************************************************************/
void topo_synth_3d(void)
{
 char nomfichres[500],str[500],*ext;
 int err,resol,taille,e;
 double Jmin,Jmax;
 
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 
 /*test de l'extension*/
  ext=strstr(nomfichres,".chp");
  if (ext!=NULL) *ext=0;
 
 taille= GET_INT("Taille", 0, &e);
 
 resol= GET_INT("resolution", 0, &e);
 
 Jmin=(double) GET_FLOAT("Jmin", 0, &err);
 Jmax=(double) GET_FLOAT("Jmax", 0, &err);

 
 rand_topo_field(nomfichres,taille,resol,Jmin,Jmax);
 
 return;
}

/*******************************************************************************
**        					rand_topo_field()                                              
**                                                                        
**       Cree le champs de deformation aleatoire conservant la topologie
**			   
*******************************************************************************/

void rand_topo_field(char *nomfichres,int taille, int resolf, double Jmin, double Jmax)                                              
{

int nb_param,resol;
double *param,*p,*grad,cfmax,cf,*param_norm;
int nbmax_param=760000;
int wdth,hght,dpth,topD;
int r,l,i,j,k,t;
transf3d *transfo;      
int *x00,*x11,*y00,*y11,*z00,*z11;
double xm,xM,ym,yM,zm,zM;
TSlpqr Slpqr[8];  
wdth=hght=dpth=taille;
 if(resolf==1)nbmax_param=3+10;
 if(resolf==2)nbmax_param=81+10;
 if(resolf==3)nbmax_param=1029+10;
 if(resolf==4)nbmax_param=10125+10;
 if(resolf==5)nbmax_param=89373+10;
 if(resolf==6)nbmax_param=750141+10;

 if((param = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return; }
 if((p = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return; }
 if((grad = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return; }
 if((param_norm = calloc(nbmax_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return; }
  
  
  nb_param=init_base_3d(wdth,hght,dpth,Bsplined1);
  resol=BASE3D.resol;

 init_log(NULL);

//srand(1); // initialise le generateur de nombre aleatoire

  /* Initialisation tableau */
  for (l=0;l<nb_param;l++)
   {param[l]=0.0;}
      
  /* Debut traitement */
  for (r=resol;r<=resolf;r++)
  {
    aff_log("\nRESOLUTION %d    nb param %d \n",r,nb_param);
	TOP_verification_all(nb_param,param,NULL,Jmin,Jmax); 
	
	for (l=0;l<nb_param;l++)
 	  {p[l]=param[l];}
  
	for (l=0;l<nb_param/3;l++)
		{
		param_norm[3*l]  =param[3*l]  /wdth;
		param_norm[3*l+1]=param[3*l+1]/hght;
		param_norm[3*l+2]=param[3*l+2]/dpth;
		}
		
    topD=(int)pow(2.0,1.0*BASE3D.resol)-1;
 	x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
  
	for (i=0;i<topD;i++)
	for (j=0;j<topD;j++)
	for (k=0;k<topD;k++)
	{
	
	 l = TOP_conv_ind(i,j,k,nb_param);
    	
  //---------------------------------------------------------------
  //----- initialisation Slpqr            -------------------------
  //---------------------------------------------------------------

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
  
	/* tirage aleatoire d'un gradient */
	do {
	grad[l]	 =((double)(2*rand()-RAND_MAX)/RAND_MAX);
	grad[l+1]=((double)(2*rand()-RAND_MAX)/RAND_MAX);
	grad[l+2]=((double)(2*rand()-RAND_MAX)/RAND_MAX);
	} while ((grad[l]==0) && (grad[l+1]==0) && (grad[l+2]==0));

	for (t=0; t<nb_param; t++)
   		{param[t]=1.0*param[t]/taille;}
	
	
	cfmax = TOP_linemin_maxpas_3d(nb_param,param,grad,i,j,k,Jmin,Jmax,Slpqr);
	cfmax=0.95*cfmax;
 	for (t=0; t<nb_param; t++)
   		{param[t]=1.0*param[t]*taille;}


	cf = ((double)rand()/RAND_MAX)*cfmax*taille;
	
	for (t=l; t<l+3; t++) p[t] = param[t] - cf*grad[t];
   
	param_norm[l]  =p[l]  /wdth;
	param_norm[l+1]=p[l+1]/hght;
	param_norm[l+2]=p[l+2]/dpth;
	printf("%d  %d  %d \r",i,j,k);  	 
	while ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,i,j,k,Slpqr)==0))
  		{ 
		cf = ((double)rand()/RAND_MAX)*cfmax*taille;
		for (t=l; t<l+3; t++) p[t] = param[t] - cf*grad[t];
    	}
		
 	for (t=l; t<l+3; t++) param[t] = p[t];
    
	}
 
   if (r!=resolf)
    {/*on passe a la resolution superieure*/
     nb_param=base_resol_up_3d(param,nb_param);
    }


  }

  
  /*enregistrement du champ resultat dans un fichier*/
   if (nomfichres!=NULL)
   { 
     transfo=cr_transf3d(wdth,hght,dpth,NULL);
     transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
     transfo->typetrans=BSPLINE3D;
     transfo->resol=resolf;transfo->degre=1;
     for (l=0;l<nb_param;l++) transfo->param[l]=param[l]; 
  printf("param : [%f  %f   %f]\n",param[0],param[1],param[2]);  
	save_transf_3d(transfo,nomfichres);
    free_transf3d(transfo);
   } 


  /*liberation de la base*/
   end_base_3d();

 free(param);

}


/*******************************************************************************
**        					chps_test_3d()                                              
**                                                                        
**       Cree le champs de deformation engendre par ma transformation
**			x -> x^2
**			y -> y^3
**			z -> z^4
**                        
*******************************************************************************/

void chps_test_3d(char *nomfichres)                                              
{
 field3d  *ch;
 transf3d *transfo;
 vector3d ***data;
 int wdth,hght,dpth,i,j,k;
 wdth=128;
 hght=128;
 dpth=128;
 ch=cr_field3d(wdth,hght,dpth);
 data=ch->raw;
 
for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
		for (k=0;k<dpth;k++)
		{
		data[i][j][k].x=1.0*wdth*pow((1.0*i/(double)wdth),0.8)-i;
		data[i][j][k].y=1.0*hght*pow((1.0*j/(double)hght),1)-j;
		data[i][j][k].z=1.0*dpth*pow((1.0*k/(double)dpth),1.2)-k;
		}
transfo=field_to_transf_3d(ch,NULL,NULL);
transfo->dx = 1.0;
transfo->dy = 1.0;
transfo->dz = 1.0;
save_transf_3d(transfo,nomfichres/*"/home/noblet/chps_test.trf"*/);
  
free_transf3d(transfo);
free_field3d(ch);

}


/*******************************************************************************
** 	void compare_critere(grphic3d *imref, grphic3d *imreca, grphic3d *imres,
**		field3d *champ_fin,double *param, double
**		param_norm, int nb_param, int ***masque_param, int nb_iter)
**
** 				Calcul tout une s�rie de crit�re 
*******************************************************************************/

void compare_critere(grphic3d *imref, grphic3d *imreca, grphic3d *imres,field3d *champ_fin,double *param, double
*param_norm, int nb_param, int ***masque_param, int nb_iter)
{
int resol,max_ref,max_reca;
double L2,L2sym,Lp,Lpsym,L1L2,quad_robust,woods,woods_robust,IM,IM_simple,IM_studholme,IM_MAES_1,IM_MAES_2,ec,cr,woods_vrai;
dist_func_locale_t dist;
ptr_distance_param dist_param=CALLOC(1, distance_param);
 
resol=BASE3D.resol;
base_to_field_3d(nb_param,param,champ_fin,imref,imreca);
inter_qsinc3_3d(imreca,champ_fin,imres); 


//------------ Quadratique ---------------
dist=Energie_quad_locale_3d;
L2=Energie_globale_3d(imref,imreca,nb_param,param,param_norm,dist,NULL,masque_param);

//---------- Quadratique symetrique ----------
dist=Energie_quad_sym_locale_3d;
L2sym=Energie_globale_3d(imref,imreca,nb_param,param,param_norm,dist,NULL,masque_param);

//------------- Lp ---------------
dist=Energie_Lp_locale_3d;
Lp=Energie_globale_3d(imref,imreca,nb_param,param,param_norm,dist,NULL,masque_param);

//-------------- Lp symetrique -------
dist=Energie_Lp_sym_locale_3d;
Lpsym=Energie_globale_3d(imref,imreca,nb_param,param,param_norm,dist,NULL,masque_param);

//-------------- L1L2 -------
dist=Energie_L1L2_locale_3d;
L1L2=Energie_globale_3d(imref,imreca,nb_param,param,param_norm,dist,NULL,masque_param);

calc_robust_stat_3d(imref,imreca,0,0);

dist_param->imreca=imref;
dist_param->imref=imres;

max_ref=imref->max_pixel;
max_reca=imreca->max_pixel;

dist_param->maxRef=max_ref;
dist_param->maxReca=max_reca;

dist_param->nbclRef=imx_calculer_nombre_de_classes_3d_p(imref,&max_ref);
dist_param->nbclReca=imx_calculer_nombre_de_classes_3d_p(imreca,&max_reca);

//------- QUAD ROBUST ---------------
quad_robust=erreur_quad_robust_geman_3d(dist_param);

//--------- Woods ----------------
woods=erreur_woods_3d(dist_param);

//------ Woods robust ------------
woods_robust=erreur_woods_robust_3d(dist_param);

//---------- IM ----------
IM=100.0/erreur_IM_3d(dist_param);

//---------- IM simple ---------
IM_simple=erreur_IM_simple_3d(dist_param);
IM_simple=-1.0/IM_simple;

//-------- IM Normalisee studholme ---------
IM_studholme=erreur_IM_Normalisee_Studholme_simple_3d(dist_param);
IM_studholme=-1.0/IM_studholme;

//--------- IM Normalisee maes1 -----------
IM_MAES_1=erreur_IM_Normalisee_Maes1_simple_3d(dist_param);
IM_MAES_1=-1.0/IM_MAES_1;

//--------- IM Normalisee maes2 -----------
IM_MAES_2=erreur_IM_Normalisee_Maes2_simple_3d(dist_param);
IM_MAES_2=IM_MAES_2;

//--------- Entropie conjointe -----------
ec=erreur_entropie_conjointe_simple_3d(dist_param);

//---------- Ratio de correlation --------
cr=erreur_CR_3d(dist_param);

//---------- Vrai critere de woods ---------
woods_vrai=erreur_Woods_vraie_3d(dist_param);


aff_flog(LOG_ENERGY,"%d \t %d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n",resol,nb_iter,L2,L2sym,Lp,Lpsym,L1L2,quad_robust,woods,woods_robust,IM,IM_simple,IM_studholme,IM_MAES_1,IM_MAES_2,ec,cr,woods_vrai);
printf("L2:%f \t Lp:%f \t L1L2:%f \t L2robuste:%f \t woods:%f \t IM:%f \t entropie_conjointe:%f \t ratio_correlation:%f \t woods_vrai:%f \n",L2,Lp,L1L2,quad_robust,woods,IM,ec,cr,woods_vrai);
//printf("IM:%f \n",IM);

FREE(dist_param);
return;
}

/******************************************************************************/
/*      distances_inter_images_3d() ----------------------------------------------
****************************************************************************/

void distances_inter_images_3d(void)
{
  int im_1,im_2;

  im_1=GET_PLACE3D("Image 1 ?");
  im_2=GET_PLACE3D("Image 2 ?");

  imx_distances_inter_images_3d(im_1,im_2);

 }

/****************************************************************************/

void imx_distances_inter_images_3d(int im_1, int im_2)
{
  grphic3d *im1,*im2,*imres;
  field3d *champ;
	double param[3];
	scal_func_t scal_func;
	int nb_param;
	int ***mask;
	
	mask=alloc_imatrix_3d(1,1,1);
	mask[0][0][0]=1;
	
  im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
	
	scal_func=Bsplined1;
	param[0]=param[1]=param[2]=0;
	nb_param=init_base_3d(im1->width,im1->height,im1->depth,scal_func);	

imres=cr_grphic3d(im2);imx_copie_param_3d_p(im2,imres);
champ=cr_field3d(im1->width,im1->height,im1->depth);
compare_critere(im1, im2, imres,champ,param,param, nb_param, mask,0);

end_base_3d();
free_imatrix_3d(mask);
free_grphic3d(imres);
free_field3d(champ);

}


/******************************************************************************/
/*      distances_inter_images_3d() ----------------------------------------------
****************************************************************************/

void correlation_inter_images_3d(void)
{
  int im_1,im_2;

  im_1=GET_PLACE3D("Image 1 ?");
  im_2=GET_PLACE3D("Image 2 ?");

  imx_correlation_inter_images_3d(im_1,im_2);

 }

/****************************************************************************/

void imx_correlation_inter_images_3d(int im_1, int im_2)
{
  grphic3d *im1,*im2;
	int i,j,k,l,m,Ntot;
	int height,width,depth;
	double moy1,moy2,var1,var2,sygma1,sygma2,cov;
	double rcoeff1,rcoeff2;

	im1=ptr_img_3d(im_1);
  im2=ptr_img_3d(im_2);
	



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

/*cacul des moyennes des images*/	
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
     {
      if (im1->mri[i][j][k])
      {
        l++;
        moy1 += im1->mri[i][j][k];
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
        moy2 += im2->mri[i][j][k];
      }
     }

moy1= rcoeff1*moy1/l;
moy2= rcoeff2*moy2/m;

/*calcul des variances des images*/
height=im1->height;
width=im1->width;
depth=im1->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im1->mri[i][j][k])
          var1 +=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff1*im1->mri[i][j][k]-moy1);
      }

height=im2->height;
width=im2->width;
depth=im2->depth;
for(i=0;i<width;i++)
  for(j=0;j<height;j++)
     for(k=0;k<depth;k++)
      {
        if (im2->mri[i][j][k])
          var2 +=(rcoeff2*im2->mri[i][j][k]-moy2)*(rcoeff2*im2->mri[i][j][k]-moy2);
      }

var1= var1/l;
var2=var2/m;
sygma1=sqrt(var1);
sygma2=sqrt(var2);

/*calcul de la covariance*/
cov=0;
Ntot=0;
for(i=0;i<width;i++)
	for(j=0;j<height;j++)
		for(k=0;k<depth;k++)
		 	 if ((im2->mri[i][j][k])||(im1->mri[i][j][k]))
			{
			cov+=(rcoeff1*im1->mri[i][j][k]-moy1)*(rcoeff2*im2->mri[i][j][k]-moy2); 
			Ntot++;
			}

cov=1.0*cov/(Ntot*sygma1*sygma2);
cov=cov*cov;			
printf("Correlation : %f \n",cov);				
}
