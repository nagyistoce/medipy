/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>

#include "noyau/io/imx_log.h"

#include "recalage/topo/analyse_intervalle.h"
#include "recalage/topo/mtch_topo_3d.h"

/*******************************************************************************
********************************************************************************
************** CONSERVATION TOPOLOGIE (ANALYSE PAR INTERVALLES) ****************
********************************************************************************
*******************************************************************************/

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
int TOP_verification_all(int nb_param,double *param,int ***masque_param,double Jm,double JM)
{ int i,j,k,l,wdth,dpth,hght;
  int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
  int aux,status=1;
  double *param_norm;
	TSlpqr Slpqr[8];
	int *x00,*x11,*y00,*y11,*z00,*z11;
	double xm,xM,ym,yM,zm,zM;
	int topD1;

	
	wdth=BASE3D.width;hght=BASE3D.height;dpth=BASE3D.depth;
  x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
	RENVERSEMENT_SEUIL=1;
	if((param_norm = calloc(nb_param, sizeof(double))) == NULL) {
   PUT_ERROR("[imx_matching_Bspline_topo_3d_p] memory allocation error !\n"); return(-1); }
 
	for (l=0;l<nb_param/3;l++)
		{
		param_norm[3*l]  =param[3*l]  /wdth;
		param_norm[3*l+1]=param[3*l+1]/hght;
		param_norm[3*l+2]=param[3*l+2]/dpth;
		}
 
 					
if (masque_param != NULL)
	{  
  	for (i=0; i<D; i++)
  	for (j=0; j<D; j++)
  	for (k=0; k<D; k++)
  	if(masque_param[i][j][k]!=0)
  		{
  		//---------------------------------------------------------------
  		//----- initialisation Slpqr            -------------------------
  		//---------------------------------------------------------------
			
			topD1 = TOP_conv_ind(i,j,k,nb_param);

  		xm = 1.0*x00[topD1/3]/wdth; xM = 1.0*x11[topD1/3]/wdth; ym = 1.0*y00[topD1/3]/hght; 
  		yM = 1.0*y11[topD1/3]/hght; zm = 1.0*z00[topD1/3]/dpth; zM =1.0*z11[topD1/3]/dpth;
	
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

			aux=TOP_verification_locale(nb_param,param,param_norm,Jm,JM,i,j,k,Slpqr);
  		if (aux==0) status=0; 
  		}
  	}
else
	{  
  	for (i=0; i<D; i++)
  	for (j=0; j<D; j++)
  	for (k=0; k<D; k++)
  		{

  		topD1 = TOP_conv_ind(i,j,k,nb_param);

			//---------------------------------------------------------------
  		//----- initialisation Slpqr            -------------------------
  		//---------------------------------------------------------------
				
  		xm = 1.0*x00[topD1/3]/wdth; xM = 1.0*x11[topD1/3]/wdth; ym = 1.0*y00[topD1/3]/hght; 
  		yM = 1.0*y11[topD1/3]/hght; zm = 1.0*z00[topD1/3]/dpth; zM =1.0*z11[topD1/3]/dpth;
	
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

			aux=TOP_verification_locale(nb_param,param,param_norm,Jm,JM,i,j,k,Slpqr);
  		if (aux==0) status=0; 
  		}
  	}
RENVERSEMENT_SEUIL=0;
return(status);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int TOP_verification_locale(int nb_param,double *param,double *param_norm,double Jm,double JM,int i,int j,int k,TSlpqr *Slpqr)
{ 
  int aux0,aux1;  								// variables auxiliaires
  double auxdbl[20],aux2;
  double auxdbl4[4];
  double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaires
  int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
  int	*x0,*x1,*y0,*y1,*z0,*z1;
  int nb_func,nbfx,nbfy,nbfz,echelle;
  int width,height,depth;
  int status=1;
  
	
if (param_norm != NULL)
	{ 
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;echelle=BASE3D.resol;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
 

  aux0 = TOP_conv_ind(i,j,k,nb_param);
 
  //-------------------------------------
  // initialisation des alx,aly,alz  ----
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { for (aux1=0; aux1<8; aux1++)
    { Slpqr[aux0].alx[aux1] = 0; Slpqr[aux0].aly[aux1] = 0; Slpqr[aux0].alz[aux1] = 0;
    }
  }

  TOP_init_slpqr(Slpqr  ,i-1,j-1,k-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+1,i-1,j-1,k  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+2,i-1,j  ,k-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+3,i-1,j  ,k  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+4,i  ,j-1,k-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+5,i  ,j-1,k  ,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+6,i  ,j  ,k-1,D,nb_param,param_norm);
  TOP_init_slpqr(Slpqr+7,i  ,j  ,k  ,D,nb_param,param_norm);


	
  //-------------------------------------
  // initialisation ind, bdelta, Jm, JM -
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++)
  { Slpqr[aux0].ind = 7 - aux0;
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
    for (aux1=0; aux1<20; aux1++) 
		{
		   Slpqr[aux0].alpha[aux1].ordorg = auxdbl21b[aux1+1];
		   Slpqr[aux0].alpha[aux1].pente  = 0.0;
	  }

  /*if ((Slpqr[aux0].b.xm==0.125)&&(Slpqr[aux0].b.ym==0.875)&&(Slpqr[aux0].b.zm==0.125))
		printf("coucou \n");

  if ((Slpqr[aux0].b.xm==0.1875)&&(Slpqr[aux0].b.ym==0.875)&&(Slpqr[aux0].b.zm==0.1875))
		printf("coucou \n");*/
		/*if ((Slpqr[aux0].b.xm==0.75)&&(Slpqr[aux0].b.ym==0.5)&&(Slpqr[aux0].b.zm==0.125))
		printf("coucou \n");*/

	}


  //-------------------------------------
  // verification step  -----------------
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++) 
  { TOP_verification_delta(&(Slpqr[aux0]),0.0,auxdbl4);
    if (auxdbl4[0]<=Slpqr[aux0].Jm)
		{ 
  	#ifndef TOPO_VERBOSE
		aff_log("Probleme conservation topologie!! General en l'origine (%d) : %e %e %e  %d %d %d\n",aux0,auxdbl4[0],Slpqr[aux0].Jm,auxdbl4[1],i,j,k);
	#endif
		status=0;
		}
    if (auxdbl4[3]>=Slpqr[aux0].JM)
		{ 
  	#ifndef TOPO_VERBOSE
		aff_log("Probleme conservation topologie!! General en l'origine (%d) : %e %e %e  %d %d %d\n",aux0,auxdbl4[2],Slpqr[aux0].JM,auxdbl4[3],i,j,k);
	#endif
		status=0;
		}
	if ((auxdbl4[1]<Jm)||(auxdbl4[2]>JM))
	  	{
	#ifndef TOPO_VERBOSE
		aff_log("Probleme  %e  Jm=%e  %e   :  %e  JM=%e   %e\n",auxdbl4[0],Jm,auxdbl4[1],auxdbl4[2],JM,auxdbl4[3]);
	#endif
  		status=0;
		}

  }

	}

return(status);
}


/*---------------------------------------------------------------------------*/
/*----- retourne la valeur max autorisée du pas -----------------------------*/
/*---------------------------------------------------------------------------*/
double TOP_linemin_maxpas_3d(int nb_param,double *param,double *moins_direc,int i,int j,int k, double Jm, double JM, TSlpqr *Slpqr)
{ int aux0,aux1;  								// variables auxiliaires
  double auxdbl[20],aux2;
  double auxdbl4[4];
  double auxdbl21a[21], auxdbl21b[21];                      // variables auxiliaires
  int D = TOP_D(nb_param); 			      // D : valeur max des indices i, j, k
  double dir_desR3[3];          				// direction de descente dans R^3
  double delta_i, delta_s;
  int	*x0,*x1,*y0,*y1,*z0,*z1;
  int nb_func,nbfx,nbfy,nbfz,echelle;
  int width,height,depth;
  nb_func=BASE3D.nb_func;nbfx=BASE3D.nbfx;nbfy=BASE3D.nbfy;nbfz=BASE3D.nbfz;echelle=BASE3D.resol;
  x0=BASE3D.x0;x1=BASE3D.x1;y0=BASE3D.y0;y1=BASE3D.y1;z0=BASE3D.z0;z1=BASE3D.z1;
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  
  //-------------------------------------
  // initialisation des boites  ---------
  //-------------------------------------
  
       
  aux0 = TOP_conv_ind(i,j,k,nb_param);
  dir_desR3[0] = - moins_direc[aux0]; dir_desR3[1] = - moins_direc[aux0+1]; dir_desR3[2] = - moins_direc[aux0+2];
  
  if ((dir_desR3[0]==0)&&(dir_desR3[1]==0)&&(dir_desR3[2]==0)) return(HUGE_VAL);
  
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
  { Slpqr[aux0].ind = 7 - aux0;
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

    //--- initialisation du terme lineaire ------------------------------ 
    TOP_poly_pente_coeff(dir_desR3,Slpqr[aux0].ind+1,Slpqr[aux0].alx,Slpqr[aux0].aly,Slpqr[aux0].alz,echelle,auxdbl);

    //--- changement de variables pour le polynome ----------------------
    for (aux1=1; aux1<=20; aux1++) auxdbl21a[aux1] = auxdbl[aux1-1];
	poly_rescale(auxdbl21a,aux2,Slpqr[aux0].b.xm,aux2,Slpqr[aux0].b.ym,aux2,Slpqr[aux0].b.zm,auxdbl21b);
    for (aux1=0; aux1<20; aux1++) Slpqr[aux0].alpha[aux1].pente = auxdbl21b[aux1+1];
  }


  //-------------------------------------
  // verification step  -----------------
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++) 
  { TOP_verification_delta(Slpqr+aux0,0.0,auxdbl4);
    if (auxdbl4[0]<=Slpqr[aux0].Jm){ 
	#ifndef TOPO_VERBOSE
	aff_log("Probleme conservation topologie!!! En l'origine (%d) : %e %e %e\n",aux0,auxdbl4[0],Slpqr[aux0].Jm,auxdbl4[1]);
	#endif
	RENVERSEMENT_SEUIL=1;	/* C'est pour eviter les petits pb numeriques aux changements de resolution ... */

	}
	if (auxdbl4[1]<Slpqr[aux0].Jm)
	  auxdbl4[3] = 1.0;
		
   if (auxdbl4[3]>=Slpqr[aux0].JM){ 
	#ifndef TOPO_VERBOSE
	aff_log("Probleme conservation topologie!!  En l'origine (%d) : %e %e %e\n",aux0,auxdbl4[2],Slpqr[aux0].JM,auxdbl4[3]);
	#endif
	RENVERSEMENT_SEUIL=1;	/* C'est pour eviter les petits pb numeriques aux changements de resolution ... */

	}
  }
	
	if (RENVERSEMENT_SEUIL==1)
	for (aux0=0; aux0<8; aux0++) 
  { TOP_verification_delta(Slpqr+aux0,0.0,auxdbl4);
    if (auxdbl4[0]<=Slpqr[aux0].Jm){ 
	//#ifndef TOPO_VERBOSE
	aff_log("Ca marche toujours pas !!!! En l'origine (%d) : %e %e %e\n",aux0,auxdbl4[0],Slpqr[aux0].Jm,auxdbl4[1]);
	//#endif
	}
	if (auxdbl4[1]<Slpqr[aux0].Jm)
	  auxdbl4[3] = 1.0;
		
   if (auxdbl4[3]>=Slpqr[aux0].JM){ 
	//#ifndef TOPO_VERBOSE
	aff_log("Ca marche toujours pas !!!!  En l'origine (%d) : %e %e %e\n",aux0,auxdbl4[2],Slpqr[aux0].JM,auxdbl4[3]);
	//#endif
	/*if ((auxdbl4[3]<=Slpqr[aux0].JM)&&(auxdbl4[0]>=Slpqr[aux0].Jm))
		aff_log("\n\n\n Cette fois ci ca a marche !!!!\n\n\n");*/
	}
  }
	
  //-------------------------------------
  // initialization step  ---------------
  //-------------------------------------
  for (aux0=0; aux0<8; aux0++) TOP_initialization_delta(Slpqr+aux0);

  //-------------------------------------
  // refinement step --------------------
  //-------------------------------------
  TOP_refinement(Slpqr,&delta_i,&delta_s);
 
  //-------------------------------------
  // verification step  -----------------
  //-------------------------------------
/*  for (aux0=0; aux0<8; aux0++) 
  { TOP_verification_delta(Slpqr+aux0,delta_i,auxdbl4);
    if (auxdbl4[0]<=Slpqr[aux0].Jm) 
	aff_log("En delta_i (%d) : %e %e %e\n",aux0,auxdbl4[0],Slpqr[aux0].Jm,auxdbl4[1]);
    if (auxdbl4[3]>=Slpqr[aux0].JM) 
	aff_log("En delta_i (%d) : %e %e %e\n",aux0,auxdbl4[2],Slpqr[aux0].JM,auxdbl4[3]);
  }*/

  //-------------------------------------
  // calcul final -----------------------
  //-------------------------------------
  RENVERSEMENT_SEUIL=0;
	return (delta_i);
}


/*---------------------------------------------------------------------------*/
/*----- refinement step de l'article ----------------------------------------*/
/*---------------------------------------------------------------------------*/
void TOP_refinement(TSlpqr *Slpqr, double *delta_i, double *delta_s)
{ double aux, aux_i, aux_s;
  int i, arg_inf=0, status; int arg_inf_optim=1;

  //-------------------------------------
  // initialisation ---------------------
  //-------------------------------------
  aux = HUGE_VAL;
  for (i=0; i<8; i++)
  { if (Slpqr[i].bdeltam<aux) { aux = Slpqr[i].bdeltam; arg_inf = i; arg_inf_optim =  1; }
    if (Slpqr[i].bdeltaM<aux) { aux = Slpqr[i].bdeltaM; arg_inf = i; arg_inf_optim = -1; }
  }
  *delta_i = HUGE_VAL;
  *delta_s = aux;

  //-------------------------------------
  // iterations -------------------------
  //-------------------------------------
  make_propo(Slpqr+arg_inf,*delta_s,arg_inf_optim,&aux_i,&aux_s);
  *delta_s = MIN2(*delta_s,aux_s);
  *delta_i = MIN2(*delta_i,aux_i);

  for (i=0; i<8; i++)
  { status = acknowledge(Slpqr+i,*delta_s,1);
    if (status==0)
    { make_propo(Slpqr+i,*delta_s,1,&aux_i,&aux_s);
      *delta_s = MIN2(*delta_s,aux_s);
      *delta_i = MIN2(*delta_i,aux_i);
    }

    status = acknowledge(Slpqr+i,*delta_s,-1);
    if (status==0)
    { make_propo(Slpqr+i,*delta_s,-1,&aux_i,&aux_s);
      *delta_s = MIN2(*delta_s,aux_s);
      *delta_i = MIN2(*delta_i,aux_i);
    }
  }

// ATTENTION :  on fait maintenant la acknoledge sur delta_i et non sur delta_s
/*stop = 1;
while (stop == 1)
	{
	stop=0;
	 for (i=0; i<8; i++)
  	{ status = acknowledge(Slpqr+i,*delta_i,1);
    	if (status==0)
    	{ make_propo(Slpqr+i,*delta_s,1,&aux_i,&aux_s);
      	*delta_s = MIN2(*delta_s,aux_s);
      	*delta_i = MIN2(*delta_i,aux_i);
		stop=1;
    	}

    	status = acknowledge(Slpqr+i,*delta_i,-1);
    	if (status==0)
    	{ make_propo(Slpqr+i,*delta_s,-1,&aux_i,&aux_s);
    	  *delta_s = MIN2(*delta_s,aux_s);
    	  *delta_i = MIN2(*delta_i,aux_i);
		  stop=1;
    	}
  	}
	
	}  
*/
}

/*---------------------------------------------------------------------------*/
/*----- make proposition ----------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void make_propo(TSlpqr *Slpqr, double delta_s, int minimiser, double *aux_i, double *aux_s)
{ double a_i[21], a_s[21];  double aux_i_old, aux_s_old, pen, seuil;
  double binf_i, bsup_i, bsupx_i, bsupy_i, bsupz_i; double binf_s, bsup_s, bsupx_s, bsupy_s, bsupz_s;
  int i,j; double gamma_i, gamma_s; int fin;
  double diff;
  int nb_iter=0;
  //-------------------------------------
  // initialisations --------------------
  //-------------------------------------
  if (minimiser==1) seuil =   Slpqr->Jm;
  else              seuil = - Slpqr->JM;

  *aux_i = 0; *aux_s = delta_s;

  //-------------------------------------
  // iterations -------------------------
  //-------------------------------------
do
{ aux_i_old = *aux_i; aux_s_old = *aux_s;

  for (i=1; i<=20; i++)
  { a_i[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + aux_i_old * Slpqr->alpha[i-1].pente);
    a_s[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + aux_s_old * Slpqr->alpha[i-1].pente);
  }

  master_minim_main (seuil,a_i,Slpqr->b,&binf_i,&bsup_i,&bsupx_i,&bsupy_i,&bsupz_i);
  master_minim_main (seuil,a_s,Slpqr->b,&binf_s,&bsup_s,&bsupx_s,&bsupy_s,&bsupz_s);

  //----- mise a jour de aux_s ---------------------------------------------
  if (bsup_s<seuil)
  { pen = TOP_eval_poly(bsupx_s,bsupy_s,bsupz_s,Slpqr->alpha,1.0,0.0);
    if (pen<0)  																// pbs numériques possibles (bsup_s<seuil et pen>0, d'où pb si pas d'accolades)
	{ gamma_s = (seuil-bsup_s) / pen + aux_s_old;
      *aux_s = gamma_s;
	}  
  }

  if (bsup_i>seuil)
  { pen = TOP_eval_poly(bsupx_i,bsupy_i,bsupz_i,Slpqr->alpha,1.0,0.0);
    if (pen<0)  gamma_i = (seuil-bsup_i) / pen + aux_i_old;
    else        gamma_i = aux_s_old;
    *aux_s = MIN2(*aux_s,gamma_i);
  }

  //----- mise a jour de aux_i ---------------------------------------------
  if ( (binf_i > seuil) && (binf_s < seuil) )
  { *aux_i = aux_i_old - (aux_i_old - aux_s_old) / (binf_i - binf_s) * (binf_i - seuil);
    if (*aux_i>=*aux_s) *aux_i = aux_i_old;										// pb numerique possible
  }

  //----- verification de aux_i ---------------------------------------------
	for (i=1; i<=20; i++)
		a_i[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + *aux_i * Slpqr->alpha[i-1].pente);
  
	master_minim_main (seuil,a_i,Slpqr->b,&binf_i,&bsup_i,&bsupx_i,&bsupy_i,&bsupz_i);
  	if (binf_i<=seuil)
		*aux_i = aux_i_old;

  //----- pour être certain de se premunir contre d'improbables problemes numeriques...----
  *aux_s = MAX2(*aux_s,0.0);
  *aux_i = MAX2(*aux_i,0.0);
  binf_i = MIN2(*aux_i,*aux_s);
  bsup_i = MAX2(*aux_i,*aux_s);
  *aux_i = binf_i;
  *aux_s = bsup_i;

if(*aux_s<=0) 
{
#ifndef TOPO_VERBOSE
printf("ATTENTION !!!!!!!!!!!! Probable pb numérique !!!!!!!!!!!! *aux_s = %f    *aux_i= %f \n ",*aux_s,*aux_i);
#endif
}
 nb_iter++;
  fin = ( fabs((*aux_i - *aux_s)) / (*aux_i + *aux_s) < 1e-3) || ((aux_i_old==*aux_i)&&(aux_s_old==*aux_s))||(nb_iter>100);
 
//  printf("%d %e  \n",fin, fabs((*aux_i - *aux_s)) / (*aux_i + *aux_s));
 
} while (fin==0);

/* } while ( (aux_i_old!=*aux_i) || (aux_s_old!=*aux_s) ); */


//----- Affinement lineaire de la valeur de aux_i (resoud le pb si aux_i =0)-------------
if ( fabs((*aux_i - *aux_s)) / (*aux_i + *aux_s) > 1e-3)
 { 
  	aux_i_old = *aux_i;

  	for (i=1; i<=20; i++)
  	{ a_i[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + aux_i_old * Slpqr->alpha[i-1].pente);
    	a_s[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + *aux_s * Slpqr->alpha[i-1].pente);
  	}

  	master_minim_main (seuil,a_i,Slpqr->b,&binf_i,&bsup_i,&bsupx_i,&bsupy_i,&bsupz_i);
  	master_minim_main (seuil,a_s,Slpqr->b,&binf_s,&bsup_s,&bsupx_s,&bsupy_s,&bsupz_s);

 	*aux_i = aux_i_old - (aux_i_old - *aux_s) / (binf_i - binf_s) * (binf_i - seuil);

	diff=*aux_i-aux_i_old;
	j=9;
	do
	{
		*aux_i=aux_i_old+0.1*j*diff;
		for (i=1; i<=20; i++)
  			 a_i[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + *aux_i * Slpqr->alpha[i-1].pente);
  		

  	master_minim_main (seuil,a_i,Slpqr->b,&binf_i,&bsup_i,&bsupx_i,&bsupy_i,&bsupz_i);
  	j--;
	} while ((binf_i<=seuil)&&(j>=0));
 }

}

/*---------------------------------------------------------------------------*/
/*----- acknowledgement -----------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// minimiser =  1 pour minimisation
// minimiser = -1 pour maximisation

int acknowledge(TSlpqr *Slpqr, double delta_s, int minimiser)
{ double a[21]; int i, status;
  double binf, bsup, bsupx, bsupy, bsupz;

  for (i=1; i<=20; i++)  a[i] = 1.0 * minimiser * (Slpqr->alpha[i-1].ordorg + delta_s * Slpqr->alpha[i-1].pente);

  if (minimiser==1)
    master_minim_main(Slpqr->Jm,a,Slpqr->b,&binf,&bsup,&bsupx,&bsupy,&bsupz);
  else
    master_minim_main(-Slpqr->JM,a,Slpqr->b,&binf,&bsup,&bsupx,&bsupy,&bsupz);

  if (minimiser==1) status = ( binf >=   Slpqr->Jm );
  else              status = ( binf >= - Slpqr->JM );

  return status;
}


/*---------------------------------------------------------------------------*/
/*----- verification step ---------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void TOP_verification_delta(TSlpqr *Slpqr, double coeff, double *result)
{ 
	int i; double a[21]; 
	double binf, bsup, bsupx, bsupy, bsupz; 
	
	a[0]=0; // Pour eviter un warning dans valgrind
  for (i=1; i<=20; i++) 
  	a[i] = Slpqr->alpha[i-1].ordorg + coeff * Slpqr->alpha[i-1].pente;

  master_minim_main (Slpqr->Jm,a,Slpqr->b,&binf,&bsup,&bsupx,&bsupy,&bsupz);
  result[0] = binf; 
  result[1] = bsup;

  for (i=1; i<=20; i++) 
  	a[i] = - a[i]; 
  master_minim_main (- Slpqr->JM,a,Slpqr->b,&binf,&bsup,&bsupx,&bsupy,&bsupz);
  result[2] = - bsup; result[3] = - binf;

}

/*---------------------------------------------------------------------------*/
/*----- initialization step de l'article ------------------------------------*/
/*---------------------------------------------------------------------------*/
void TOP_initialization_delta(TSlpqr *Slpqr)
{ int i,j,k;  double pente, ordorg; double aux0;
  double x[2], y[2], z[2];

  x[0] = Slpqr->b.xm; x[1] = Slpqr->b.xM;
  y[0] = Slpqr->b.ym; y[1] = Slpqr->b.yM;
  z[0] = Slpqr->b.zm; z[1] = Slpqr->b.zM;

/*  dx = 0.05*(x[1]-x[0]);   dy = 0.05*(y[1]-y[0]);   dz = 0.05*(z[1]-z[0]);
  x[0] = x[0] + dx; x[1] = x[1] - dx;
  y[0] = y[0] + dy; y[1] = y[1] - dy;
  z[0] = z[0] + dz; z[1] = z[1] - dz;*/

  for (i=0; i<2; i++)
  { for (j=0; j<2; j++)
    { for (k=0; k<2; k++)
      { pente  = TOP_eval_poly(x[i],y[j],z[k],Slpqr->alpha,1,0);
        ordorg = TOP_eval_poly(x[i],y[j],z[k],Slpqr->alpha,0,1);
        if ((ordorg<=Slpqr->Jm)||(ordorg>=Slpqr->JM))    		// possible suite à pb numérique
        { Slpqr->bdeltam = 0 ; Slpqr->bdeltaM = 0;
        }
        else if (pente>0)
        { aux0 = (Slpqr->JM-ordorg) / pente;
          if (aux0<Slpqr->bdeltaM) Slpqr->bdeltaM = aux0;
        }
        else if (pente<0)
        { aux0 = (Slpqr->Jm-ordorg) / pente;
          if (aux0<Slpqr->bdeltam) Slpqr->bdeltam = aux0;
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
double TOP_eval_poly(double xt, double yt, double zt, Tdroite *droite, double cpente, double cordorg)
{ double valeur = 0,newpol[20]; int i;

  for (i=0; i<20; i++)
    newpol[i] = droite[i].pente * cpente + droite[i].ordorg * cordorg;

  valeur = newpol[0] + newpol[1] * xt + newpol[2] * yt + newpol[3] * zt + newpol[4] * xt*xt \
          + newpol[5] * xt*yt + newpol[6] * xt*zt + newpol[7] * yt*yt + newpol[8] * yt*zt + newpol[9] * zt*zt \
          + newpol[10] * xt*xt*yt + newpol[11] * xt*xt*zt + newpol[12] * xt*yt*yt + newpol[13] * xt*yt*zt \
          + newpol[14] * xt*zt*zt + newpol[15] * yt*yt*zt + newpol[16] * yt*zt*zt +  newpol[17] * xt*xt*yt*zt \
          + newpol[18] * xt*yt*yt*zt +  newpol[19] * xt*yt*zt*zt;

  return(valeur);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void TOP_init_slpqr(TSlpqr *Slpqr,int p,int q, int r, int D, int nb_param, double *param)
{ int aux0, auxp, auxq, auxr, auxcoeff = 0;

  for (auxp=p; auxp<p+2; auxp++)
  { for (auxq=q; auxq<q+2; auxq++)
    { for (auxr=r; auxr<r+2; auxr++)
      { if ((auxp>=0)&&(auxq>=0)&&(auxr>=0)&&(auxp<=D-1)&&(auxq<=D-1)&&(auxr<=D-1))
        { aux0 = TOP_conv_ind(auxp,auxq,auxr,nb_param);
		 Slpqr->alx[auxcoeff] = param[aux0]; Slpqr->aly[auxcoeff] = param[aux0+1]; Slpqr->alz[auxcoeff] = param[aux0+2];
       
        }
        auxcoeff++;
      }
    }
  }
}



/*----------------------------------------------------------------------*/
/*----- permet d'adapter l'interface d'appel par rapport a minim_main --*/
/*----------------------------------------------------------------------*/
void master_minim_main(double seuil_arret, double *a, Tparall b, double *binf, double *bsup, double *bsupx, double *bsupy, double *bsupz)
{ double precis,J;

  precis = dim_max_parall(b) / 100.0;
	J=seuil_arret;
  seuil_arret = HUGE_VAL;
  minim_main(seuil_arret,precis, a, b.xm, b.xM, b.ym, b.yM, b.zm, b.zM, binf, bsup, bsupx, bsupy, bsupz,J);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
%-----------------------------------------------------------
%----- fichier : minim_main.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : procédure de minimisation
%-----           routine principale
%-----
%-----           retourne bornes inf et sup du minimum
%-----                    argument de la borne sup
%-----
%-----           recherche dans le pavé xinf,xsup,yinf,ysup,zinf,zsup
%-----
%-----           precis : précision de la recherche
%-----           binf et bsup sont calculées sur des pavés dont les dim < precis
%-----------------------------------------------------------


function [binf,bsup,bsupx,bsupy,bsupz] = minim_main(precis,a,xinf,xsup,yinf,ysup,zinf,zsup);*/

void minim_main(double seuil_arret, double precis, double *a, double xinf, double xsup, double yinf, double ysup, double zinf, double zsup, double *binf, double *bsup, double *bsupx, double *bsupy, double *bsupz, double J)
{ int GRIMAX, i, status, niter = 0; double c_bar;
  int aux_coinx1, aux_coinx2, aux_coiny1, aux_coiny2, aux_coinz1, aux_coinz2, aux_coinx, aux_coiny, aux_coinz;
  double c_tildex, c_tildey, c_tildez;
  double xcass[2], ycass[2], zcass[2];
  double myX[2], myY[2], myZ[2];   int ix, iy, iz;
  double auxdbl, auxdbl2[2], auxdbl21[21], seuil;
  TOPTliste myQ, myL, myLaux;
  TOPTbox myB, myBsmall, *myptrBsmall;

  /*%-----------------------------------------------------------
  %----- initialisation --------------------------------------
  %-----------------------------------------------------------*/
  GRIMAX = 15;           // % nombre de boucles sur pave_grignote
  myQ.nelem = 1;
  myQ.tete = (TOPTbox*) malloc(sizeof(TOPTbox));
  myQ.tete->next = NULL;
  myQ.tete->xm = xinf;   myQ.tete->xM = xsup;
  myQ.tete->ym = yinf;   myQ.tete->yM = ysup;
  myQ.tete->zm = zinf;   myQ.tete->zM = zsup;

  for (i=0; i<GRIMAX; i++)
  { pave_grignote(precis,a,myQ.tete,&status);
    if (status==0) break;
    if (dim_max_box(*myQ.tete)<0.9*precis) break;
  }

  myB = *(myQ.tete);
  if_natural(a,myB.xm,myB.xM,myB.ym,myB.yM,myB.zm,myB.zM,auxdbl2);
  myQ.tete->fm = MIN2(auxdbl2[0],auxdbl2[1]);
  myQ.tete->fM = MAX2(auxdbl2[0],auxdbl2[1]);

  if (myQ.tete->fm == myQ.tete->fM)
  { *binf = myQ.tete->fm; *bsup = *binf;
    *bsupx = (myB.xm+myB.xM)/2.0; *bsupy = (myB.ym+myB.yM)/2.0; *bsupz = (myB.zm+myB.zM)/2.0;
	free(myQ.tete); 
	return;
  }
  
  if (myQ.tete->fm > seuil_arret)
  { *binf = myQ.tete->fm; 
    *bsupx = (myB.xm+myB.xM)/2.0; *bsupy = (myB.ym+myB.yM)/2.0; *bsupz = (myB.zm+myB.zM)/2.0;
    eval_fun(a,*bsupx,*bsupy,*bsupz,bsup,auxdbl21);
    free(myQ.tete);
	return;
  }

  c_bar = MAX2(auxdbl2[0],auxdbl2[1]);                              // c_bar : borne sup du minimum global
  myL.nelem = 0;     myL.tete = NULL;
  myLaux.nelem = 0;  myLaux.tete = NULL;

  //-----------------------------------------------------------
  //----- boucle principale - scan du domaine -----------------
  //-----------------------------------------------------------
  while (myQ.nelem>0)
  {
    niter++;

//    printf("%4d %4d %4d \n",niter,myQ.nelem,myL.nelem);

    myB = depile_end(&myQ);                     //myB = myQ(:,length(myQ(1,:)));    myQ = myQ(:,1:length(myQ(1,:))-1);         % on traite les petits pavés d'abord pour éviter d'avoir une longue liste chainée (casser des gros pavés ne donne pas d'indication précise sur l'encadrement du minimum)

    if (myB.fm<=c_bar)										// myBox = (*myListe).tete.b;
    {
      if (dim_max_box(myB)<=precis)
      { empile(&myL,myB);
        seuil = myB.fm - (myB.fM - myB.fm);
				//seuil = myB.fm;
				
				if ((seuil<J)&&(myB.fm>J)&&(RENVERSEMENT_SEUIL==1))
					seuil= 1.0*(myB.fm+J)/2;
				
				
				c_bar = MIN2(seuil,c_bar);
      }
      else
      {
        minim_iter(precis,a,myB,&c_tildex,&c_tildey,&c_tildez);

        aux_coinx1 = (c_tildex-myB.xm<=precis); aux_coinx2 = (myB.xM-c_tildex<=precis);
        aux_coiny1 = (c_tildey-myB.ym<=precis); aux_coiny2 = (myB.yM-c_tildey<=precis);
        aux_coinz1 = (c_tildez-myB.zm<=precis); aux_coinz2 = (myB.zM-c_tildez<=precis);

        aux_coinx = (aux_coinx1 || aux_coinx2); aux_coiny = (aux_coiny1 || aux_coiny2); aux_coinz = (aux_coinz1 || aux_coinz2);

        if (aux_coinx + aux_coiny + aux_coinz == 3)

        { xcass[0] = (myB.xm + myB.xM) / 2.0;
          if      (aux_coinx1==1) xcass[1] = myB.xm + precis;
          else if (aux_coinx2==1) xcass[1] = myB.xM - precis;

          ycass[0] = (myB.ym + myB.yM) / 2.0;
          if      (aux_coiny1==1) ycass[1] = myB.ym + precis;
          else if (aux_coiny2==1) ycass[1] = myB.yM - precis;

          zcass[0] = (myB.zm + myB.zM) / 2.0;
          if      (aux_coinz1==1) zcass[1] = myB.zm + precis;
          else if (aux_coinz2==1) zcass[1] = myB.zM - precis;

          pave_casser(precis,myB,xcass,ycass,zcass,&myLaux);
        }

        else
        { xcass[0] = c_tildex - precis / 3.0;  xcass[1] = c_tildex + precis / 3.0;
          ycass[0] = c_tildey - precis / 3.0;  ycass[1] = c_tildey + precis / 3.0;
          zcass[0] = c_tildez - precis / 3.0;  zcass[1] = c_tildez + precis / 3.0;
          pave_casser(precis,myB,xcass,ycass,zcass,&myLaux);
          eval_fun(a,c_tildex,c_tildey,c_tildez,&seuil,auxdbl21);
          c_bar = MIN2(c_bar,seuil);
        }

        myptrBsmall = NULL;
        while (myLaux.nelem>0)
        { myB = depile_end(&myLaux);
          if_natural(a,myB.xm,myB.xM,myB.ym,myB.yM,myB.zm,myB.zM,auxdbl2);

          if (MIN2(auxdbl2[0],auxdbl2[1])<=c_bar)
          { for (i=0; i<GRIMAX; i++)
            { pave_grignote(precis,a,&myB,&status);
              if (status==0) break;
              if_natural(a,myB.xm,myB.xM,myB.ym,myB.yM,myB.zm,myB.zM,auxdbl2);
              if (MIN2(auxdbl2[0],auxdbl2[1])>c_bar) break;
              if (dim_max_box(myB)<precis) break;
            }

            myB.fm = MIN2(auxdbl2[0],auxdbl2[1]);
            myB.fM = MAX2(auxdbl2[0],auxdbl2[1]);
            c_bar = MIN2(c_bar,myB.fM);
            if (myB.fm<=c_bar)
            { empile(&myQ,myB);
              if (dim_max_box(myB)<=precis)
              { myBsmall = myB; myptrBsmall = &myBsmall; }
            }
          }    // fin if (MIN2() <= c_bar)
        }      // fin while (myLaux.nelem>0)

        if (myptrBsmall!=NULL) empile(&myQ,myBsmall);
      }
    }       //% fin "if (myB(7)<=c_bar)"
  }         // fin tant que myQ nonj vide

  //%-----------------------------------------------------------
  //%----- fin scan du domaine ---------------------------------
  //%-----------------------------------------------------------

  //%-----------------------------------------------------------
  //%----- extraction des bornes de l'optimum ------------------
  //%-----------------------------------------------------------
  //% borne sup retenue : min sur les sommets (+ l'argument)

  //binf = min(myL(7,:));                   //% et pas max(myL(7,:)) car ce n'est pas forcément le même minimum qui est encadré dans chaque cas

  *binf = c_bar;
  *bsup = HUGE_VAL;

  /*while (myL.nelem>0)
  { myB = depile_end(&myL);

    eval_fun(a,(myB.xm+myB.xM)/2.0,(myB.ym+myB.yM)/2.0,(myB.zm+myB.zM)/2.0,&auxdbl,auxdbl21);
    if (auxdbl<*bsup)
      { *bsup = auxdbl; *bsupx = (myB.xm+myB.xM)/2.0 ; *bsupy = (myB.ym+myB.yM)/2.0; *bsupz = (myB.zm+myB.zM)/2.0; }
  }*/


  while (myL.nelem>0)
  { myB = depile_end(&myL);

    myX[0] = myB.xm; myX[1] = myB.xM;
    myY[0] = myB.ym; myY[1] = myB.yM;
    myZ[0] = myB.zm; myZ[1] = myB.zM;

    for (ix=0; ix<2; ix++)
    for (iy=0; iy<2; iy++)
    for (iz=0; iz<2; iz++)
    { eval_fun(a,myX[ix],myY[iy],myZ[iz],&auxdbl,auxdbl21);
      if (auxdbl<*bsup)
      { *bsup = auxdbl; *bsupx = myX[ix]; *bsupy = myY[iy]; *bsupz = myZ[iz]; }
    }
  }  
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double dim_max_parall(Tparall b)
{ double wx, wy, wz;
  wx = b.xM - b.xm;  wy = b.yM - b.ym; wz = b.zM - b.zm;
  if ((wx>=wy)&&(wx>=wz)) return wx;
  else if (wy>=wz) return wy;
  else return wz;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double dim_max_box(TOPTbox b)
{ double wx, wy, wz;
  wx = b.xM - b.xm;  wy = b.yM - b.ym; wz = b.zM - b.zm;
  if ((wx>=wy)&&(wx>=wz)) return wx;
  else if (wy>=wz) return wy;
  else return wz;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : pave_grignote.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : contraction d'un pavé par encradrement du gradient
%-----           Voir classeur OPT 4 et 5
%-----------------------------------------------------------


function [xinf,xsup,yinf,ysup,zinf,zsup] = pave_grignote(precis,a,xinf,xsup,yinf,ysup,zinf,zsup);*/

void pave_grignote(double precis,double *a,TOPTbox *b, int *status)
{ double mygrad[21][4]; double big_mat[21][10];
  double alpha[21], beta[21];
  double alpha_min, beta_min, alpha_max, beta_max;
  double xinf, xsup, yinf, ysup, zinf, zsup;
  double myinf[4], mysup[4];
  double xi, xs, xih, xsh, xil, xsl;
  int modxi, modxs;
  double auxdbl2[2];
  int i, j;

  /*-----------------------------------------------------------*/
  /*-----------------------------------------------------------*/
  /*-----------------------------------------------------------*/
  *status = 0;
  xinf = (*b).xm; xsup = (*b).xM;
  yinf = (*b).ym; ysup = (*b).yM;
  zinf = (*b).zm; zsup = (*b).zM;

  if (dim_max_box(*b)<0.9*precis) return;

  poly_grad(a,mygrad);
  myinf[1] = xinf; myinf[2] = yinf; myinf[3] = zinf;
  mysup[1] = xsup; mysup[2] = ysup; mysup[3] = zsup;

  //-----------------------------------------------------------
  for (i=1; i<=3; i++)                           // on traite successivement chacune des trois directions
  {
    //-----------------------------------------------------------
    xi = myinf[i]; xs = mysup[i];
    modxi = 0; modxs = 0;                        // booléens indiquant la modification des bornes

    poly_marg_poly_fast(i,mygrad,big_mat);		 // ATTENTION : modification par rapport à la version matlab
    //big_mat = poly_marg_poly_fast(i,mygrad(:,i));

    for (j=1; j<=20; j++)
    { alpha[j] = big_mat[j][3*(i-1)+2];
      beta[j]  = big_mat[j][3*(i-1)+3];
    }

    //alpha = big_mat(:,3*(i-1)+2);
    //beta  = big_mat(:,3*(i-1)+3);

    if_natural(alpha,myinf[1],mysup[1],myinf[2],mysup[2],myinf[3],mysup[3],auxdbl2);
    alpha_min = MIN2(auxdbl2[0],auxdbl2[1]);
    alpha_max = MAX2(auxdbl2[0],auxdbl2[1]);
    if_natural(beta,myinf[1],mysup[1],myinf[2],mysup[2],myinf[3],mysup[3],auxdbl2);
    beta_min = MIN2(auxdbl2[0],auxdbl2[1]);
    beta_max = MAX2(auxdbl2[0],auxdbl2[1]);

    /*aux = if_natural(alpha,myinf(1),mysup(1),myinf(2),mysup(2),myinf(3),mysup(3));
    alpha_min = min(aux);
    alpha_max = max(aux);
    aux = if_natural(beta,myinf(1),mysup(1),myinf(2),mysup(2),myinf(3),mysup(3));
    beta_min = min(aux);
    beta_max = max(aux);*/

    //-----------------------------------------------------------
    if (xi<0)
    { xih = alpha_min * xi + beta_max;
      if (xih<0)
      { xsh = beta_max;
        if (xsh<=0) { xi = 0;                       modxi = 1; }
        else        { xi = - beta_max / alpha_min;  modxi = 1; }
      }
    }

    //-----------------------------------------------------------
    if (xi>=0)
    { xih = alpha_max * xi + beta_max;
      if (xih<0)
      { xsh = alpha_max * xs + beta_max;
        if (xsh<=0) { xi = xs;                      modxi = 1; }
        else        { xi = - beta_max / alpha_max;  modxi = 1; }
      }
    }

    //-----------------------------------------------------------
    if (xs>0)
    { xsl = alpha_min * xs + beta_min;
      if (xsl>0)
      { xil = beta_min;
        if (xil>=0) { xs = 0;                       modxs = 1; }
        else        { xs = - beta_min / alpha_min;  modxs = 1; }
      }
    }

    //-----------------------------------------------------------
    if (xs<=0)
    { xsl = alpha_max * xs + beta_min;
      if (xsl>0)
      { xil = alpha_max * xi + beta_min;
        if (xil>=0) { xs = xi;                            modxs = 1; }
        else        { xs = - beta_min / alpha_max;        modxs = 1; }
      }
    }

    //-----------------------------------------------------------
    if (modxi==1)
    { if (xi-precis/3 > myinf[i])
      { myinf[i] = xi - precis / 3;
        myinf[i] = MIN2(mysup[i], myinf[i]);           // min([mysup(i) , myinf(i)]);   //% au cas où on pourrait avoir myinf(i)>mysup(i) avec la ligne précédente, ce qui est possible
      }
    }

    if (modxs==1)
    { if (xs+precis/3 < mysup[i])
      { mysup[i] = xs + precis / 3;
        mysup[i] = MAX2(mysup[i], myinf[i]);                             // max([mysup(i) , myinf(i)]);                  % au cas où on pourrait avoir mysup(i)<myinf(i) avec la ligne précédente, ce qui est possible
      }
    }

  }    // fin for (i=1; i<=3; i++)


   //-----------------------------------------------------------
   //-----------------------------------------------------------
   //-----------------------------------------------------------

   //xinf = myinf(1); yinf = myinf(2); zinf = myinf(3);
   //xsup = mysup(1); ysup = mysup(2); zsup = mysup(3);

   if ( ((*b).xm!=myinf[1]) || ((*b).xM!=mysup[1]) || ((*b).ym!=myinf[2]) || ((*b).yM!=mysup[2]) || ((*b).zm!=myinf[3]) || ((*b).zM!=mysup[3]) )
     *status = 1;

  (*b).xm = myinf[1]; (*b).xM = mysup[1];
  (*b).ym = myinf[2]; (*b).yM = mysup[2];
  (*b).zm = myinf[3]; (*b).zM = mysup[3];

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : poly_marg.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : marginales du polynome
%-----           Valeurs des coefficients en un point donné
%-----------------------------------------------------------
function [marg] = poly_marg(indice,a,xt,yt,zt);*/

void poly_marg(int indice, double *a, double xt, double yt, double zt, double *marg)
{

  switch (indice)

  {  case 1 :                                       // marginale en x

     marg[1] = a[5] + yt*a[11] + zt*a[12] + yt*zt*a[18];
     marg[2] = a[2] + yt*a[6] + zt*a[7] + yt*zt*a[14] + yt*yt*a[13] + zt*zt*a[15] + yt*yt*zt*a[19] + yt*zt*zt*a[20];
     marg[3] = a[1] + yt*a[3] + zt*a[4] + yt*zt*a[9] + yt*yt*a[8] + zt*zt*a[10] + yt*yt*zt*a[16] + yt*zt*zt*a[17];
     break;

     case 2 :                                       // marginale en y

     marg[1] = a[8] + xt*a[13] + zt*a[16] + xt*zt*a[19];
     marg[2] = a[3] + xt*a[6] + zt*a[9] + xt*zt*a[14] + xt*xt*a[11] + zt*zt*a[17] + xt*xt*zt*a[18] + xt*zt*zt*a[20];
     marg[3] = a[1] + xt*a[2] + zt*a[4] + xt*zt*a[7] + xt*xt*a[5] + zt*zt*a[10] + xt*xt*zt*a[12] + xt*zt*zt*a[15];
     break;

     case 3 :                                       // marginale en z

     marg[1] = a[10] + xt*a[15] + yt*a[17] + xt*yt*a[20];
     marg[2] = a[4] + xt*a[7] + yt*a[9] + xt*yt*a[14] + xt*xt*a[12] + yt*yt*a[16] + xt*xt*yt*a[18] + xt*yt*yt*a[19];
     marg[3] = a[1] + xt*a[2] + yt*a[3] + xt*yt*a[6] + xt*xt*a[5] + yt*yt*a[8] + xt*xt*yt*a[11] + xt*yt*yt*a[13];

  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
%-----------------------------------------------------------
%----- fichier : poly_marg_poly_fast.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : marginales du polynome sous forme de polynome
%-----           voir poly_marg.txt (MuPad)
%-----           matrice big_mat remplie seulement partiellement
%-----           Remplissage complet : voir poly_marg_poly.m
%-----------------------------------------------------------

%-----------------------------------------------------------
% Sortie de la procédure : matrice
%   [ x2 x1 x0 y2 y1 y0 z2 z1 z0 ]
% chaque colonne contient les coefficients du polynome.
% Ex : colonne 4 : marginale en y, coefficient de y2 sous
% forme de vecteur de coordonnées dans la base des polynomes
% considérés.
%-----------------------------------------------------------

function [big_mat] = poly_marg_poly_fast(param,a); */

//
// ATTENTION : modification des paramètres d'appel par rapport à la version matlab
//
// a ci-dessus correspond dans la fonction appelante à my_grad(:,param)
//

void poly_marg_poly_fast(int param, double my_grad[][4], double big_mat[][10])
{ int i, j;
  double r[21];
  double a[21];

  for (i=1; i<=20; i++) a[i] = my_grad[i][param];

  //big_mat = zeros(length(a),9);
  //r = zeros(size(a));                                 //% présence de r car commodité MuPad

  for (i=1; i<=20; i++)
  { for (j=1; j<=9; j++)
    { big_mat[i][j] = 0.0;
    }
  }

  for (i=1; i<=20; i++) r[i] = 0.0;

//%-----------------------------------------------------------
//%-----------------------------------------------------------
if (param==1)
{
big_mat[1 ][ 2] = a[2] + r[1];
big_mat[2 ][ 2] = r[2];
big_mat[3 ][ 2] = a[6] + r[3];
big_mat[4 ][ 2] = a[7] + r[4];
big_mat[5 ][ 2] = r[5];
big_mat[6 ][ 2] = r[6];
big_mat[7 ][ 2] = r[7];
big_mat[8 ][ 2] = a[13] + r[8];
big_mat[9 ][ 2] = a[14] + r[9];
big_mat[10 ][ 2] = a[15] + r[10];
big_mat[11 ][ 2] = r[11];
big_mat[12 ][ 2] = r[12];
big_mat[13 ][ 2] = r[13];
big_mat[14 ][ 2] = r[14];
big_mat[15 ][ 2] = r[15];
big_mat[16 ][ 2] = a[19] + r[16];
big_mat[17 ][ 2] = a[20] + r[17];
big_mat[18 ][ 2] = r[18];
big_mat[19 ][ 2] = r[19];
big_mat[20 ][ 2] = r[20];

big_mat[1 ][ 3] = a[1] + r[1];
big_mat[2 ][ 3] = r[2];
big_mat[3 ][ 3] = a[3] + r[3];
big_mat[4 ][ 3] = a[4] + r[4];
big_mat[5 ][ 3] = r[5];
big_mat[6 ][ 3] = r[6];
big_mat[7 ][ 3] = r[7];
big_mat[8 ][ 3] = a[8] + r[8];
big_mat[9 ][ 3] = a[9] + r[9];
big_mat[10 ][ 3] = a[10] + r[10];
big_mat[11 ][ 3] = r[11];
big_mat[12 ][ 3] = r[12];
big_mat[13 ][ 3] = r[13];
big_mat[14 ][ 3] = r[14];
big_mat[15 ][ 3] = r[15];
big_mat[16 ][ 3] = a[16] + r[16];
big_mat[17 ][ 3] = a[17] + r[17];
big_mat[18 ][ 3] = r[18];
big_mat[19 ][ 3] = r[19];
big_mat[20 ][ 3] = r[20];
}
//%-----------------------------------------------------------
//%-----------------------------------------------------------
else if (param==2)
{
big_mat[1 ][ 5] = a[3] + r[1];
big_mat[2 ][ 5] = a[6] + r[2];
big_mat[3 ][ 5] = r[3];
big_mat[4 ][ 5] = a[9] + r[4];
big_mat[5 ][ 5] = a[11] + r[5];
big_mat[6 ][ 5] = r[6];
big_mat[7 ][ 5] = a[14] + r[7];
big_mat[8 ][ 5] = r[8];
big_mat[9 ][ 5] = r[9];
big_mat[10 ][ 5] = a[17] + r[10];
big_mat[11 ][ 5] = r[11];
big_mat[12 ][ 5] = a[18] + r[12];
big_mat[13 ][ 5] = r[13];
big_mat[14 ][ 5] = r[14];
big_mat[15 ][ 5] = a[20] + r[15];
big_mat[16 ][ 5] = r[16];
big_mat[17 ][ 5] = r[17];
big_mat[18 ][ 5] = r[18];
big_mat[19 ][ 5] = r[19];
big_mat[20 ][ 5] = r[20];

big_mat[1 ][ 6] = a[1] + r[1];
big_mat[2 ][ 6] = a[2] + r[2];
big_mat[3 ][ 6] = r[3];
big_mat[4 ][ 6] = a[4] + r[4];
big_mat[5 ][ 6] = a[5] + r[5];
big_mat[6 ][ 6] = r[6];
big_mat[7 ][ 6] = a[7] + r[7];
big_mat[8 ][ 6] = r[8];
big_mat[9 ][ 6] = r[9];
big_mat[10 ][ 6] = a[10] + r[10];
big_mat[11 ][ 6] = r[11];
big_mat[12 ][ 6] = a[12] + r[12];
big_mat[13 ][ 6] = r[13];
big_mat[14 ][ 6] = r[14];
big_mat[15 ][ 6] = a[15] + r[15];
big_mat[16 ][ 6] = r[16];
big_mat[17 ][ 6] = r[17];
big_mat[18 ][ 6] = r[18];
big_mat[19 ][ 6] = r[19];
big_mat[20 ][ 6] = r[20];
}
//%-----------------------------------------------------------
//%-----------------------------------------------------------
else
{
big_mat[1 ][ 8] = a[4] + r[1];
big_mat[2 ][ 8] = a[7] + r[2];
big_mat[3 ][ 8] = a[9] + r[3];
big_mat[4 ][ 8] = r[4];
big_mat[5 ][ 8] = a[12] + r[5];
big_mat[6 ][ 8] = a[14] + r[6];
big_mat[7 ][ 8] = r[7];
big_mat[8 ][ 8] = a[16] + r[8];
big_mat[9 ][ 8] = r[9];
big_mat[10 ][ 8] = r[10];
big_mat[11 ][ 8] = a[18] + r[11];
big_mat[12 ][ 8] = r[12];
big_mat[13 ][ 8] = a[19] + r[13];
big_mat[14 ][ 8] = r[14];
big_mat[15 ][ 8] = r[15];
big_mat[16 ][ 8] = r[16];
big_mat[17 ][ 8] = r[17];
big_mat[18 ][ 8] = r[18];
big_mat[19 ][ 8] = r[19];
big_mat[20 ][ 8] = r[20];
big_mat[1 ][ 9] = a[1] + r[1];
big_mat[2 ][ 9] = a[2] + r[2];
big_mat[3 ][ 9] = a[3] + r[3];
big_mat[4 ][ 9] = r[4];
big_mat[5 ][ 9] = a[5] + r[5];
big_mat[6 ][ 9] = a[6] + r[6];
big_mat[7 ][ 9] = r[7];
big_mat[8 ][ 9] = a[8] + r[8];
big_mat[9 ][ 9] = r[9];
big_mat[10 ][ 9] = r[10];
big_mat[11 ][ 9] = a[11] + r[11];
big_mat[12 ][ 9] = r[12];
big_mat[13 ][ 9] = a[13] + r[13];
big_mat[14 ][ 9] = r[14];
big_mat[15 ][ 9] = r[15];
big_mat[16 ][ 9] = r[16];
big_mat[17 ][ 9] = r[17];
big_mat[18 ][ 9] = r[18];
big_mat[19 ][ 9] = r[19];
big_mat[20 ][ 9] = r[20];
}
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : poly_grad.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : gradient du polynome
%-----------------------------------------------------------


function [mygrad] = poly_grad(a);*/

void poly_grad(double *a, double mygrad[][4])
{ int i,j;

  //mygrad = zeros(length(a),3);

  for (i=1; i<=20; i++)
  { for (j=1; j<=3; j++)
      mygrad[i][j] = 0.0;
  }


  mygrad[1][1] = a[2];
  mygrad[2][1] = 2*a[5];
  mygrad[3][1] = a[6];
  mygrad[4][1] = a[7];
  mygrad[6][1] = 2*a[11];
  mygrad[7][1] = 2*a[12];
  mygrad[8][1] = a[13];
  mygrad[9][1] = a[14];
  mygrad[10][1] = a[15];
  mygrad[14][1] = 2*a[18];
  mygrad[16][1] = a[19];
  mygrad[17][1] = a[20];


  mygrad[1][2] = a[3];
  mygrad[2][2] = a[6];
  mygrad[3][2] = 2*a[8];
  mygrad[4][2] = a[9];
  mygrad[5][2] = a[11];
  mygrad[6][2] = 2*a[13];
  mygrad[7][2] = a[14];
  mygrad[9][2] = 2*a[16];
  mygrad[10][2] = a[17];
  mygrad[12][2] = a[18];
  mygrad[14][2] = 2*a[19];
  mygrad[15][2] = a[20];

  mygrad[1][3] = a[4];
  mygrad[2][3] = a[7];
  mygrad[3][3] = a[9];
  mygrad[4][3] = 2*a[10];
  mygrad[5][3] = a[12];
  mygrad[6][3] = a[14];
  mygrad[7][3] = 2*a[15];
  mygrad[8][3] = a[16];
  mygrad[9][3] = 2*a[17];
  mygrad[11][3] = a[18];
  mygrad[13][3] = a[19];
  mygrad[14][3] = 2*a[20];
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : if_natural.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : natural inclusion function
%-----
%-----           A implanter : forcer les arrondis - IEEE 754
%-----------------------------------------------------------

function [bornes] = if_natural(a,xinf,xsup,yinf,ysup,zinf,zsup);
%-----------------------------------------------------------
%-----------------------------------------------------------
%-----------------------------------------------------------
bornes = zeros(2,1);

%-----------------------------------------------------------
pol1 = poly_transl(a,-xinf,-yinf,-zinf);

[aux , monomes] = eval_fun(pol1,xsup-xinf,ysup-yinf,zsup-zinf);

aux = monomes(2:length(monomes));

bornes(1) = monomes(1) + sum( aux .* (aux<0) );
bornes(2) = monomes(1) + sum( aux .* (aux>0) );
*/


void if_natural(double *a, double xinf, double xsup, double yinf, double ysup, double zinf, double zsup, double *auxdbl2)
{ double pol1[21]; int i;
  double aux, monomes[21];

  poly_transl(a,-xinf,-yinf,-zinf,pol1);

  eval_fun(pol1,xsup-xinf,ysup-yinf,zsup-zinf,&aux,monomes);

  auxdbl2[0] = monomes[1];
  auxdbl2[1] = monomes[1];
  for (i=2; i<=20; i++)
  { if (monomes[i]<0)
      auxdbl2[0] += monomes[i];
    else
      auxdbl2[1] += monomes[i];
  }

}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : poly_transl.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : changement de variable
%-----
%-----           nxt = xt + b1
%-----           nyt = yt + b2
%-----           nzt = zt + b3
%-----------------------------------------------------------


function [alpha] = poly_transl(a,b1,b2,b3);*/

void poly_transl(double *a, double b1, double b2, double b3, double *alpha)
{
  //alpha = zeros(20,1);

  alpha[1] = a[1] - b1*(a[2] - b1*a[5]) - b2*(a[3] - b1*(a[6] - b1*a[11]) - b2*(a[8] - b1*a[13])) - b3*(a[4] - b1*(a[7] - b1*a[12]) - b3*(a[10] - b1*a[15] - b2*(a[17] - b1*a[20])) - b2*(a[9] - b1*(a[14] - b1*a[18]) - b2*(a[16] - b1*a[19])));
  alpha[2] = a[2] - 2*b1*a[5] - b2*(a[6] - 2*b1*a[11] - b2*a[13]) - b3*(a[7] - 2*b1*a[12] - b3*(a[15] - b2*a[20]) - b2*(a[14] - 2*b1*a[18] - b2*a[19]));
  alpha[3] = a[3] - b1*(a[6] - b1*a[11]) - 2*b2*(a[8] - b1*a[13]) - b3*(a[9] - b1*(a[14] - b1*a[18]) - 2*b2*(a[16] - b1*a[19]) - b3*(a[17] - b1*a[20]));
  alpha[4] = a[4] - b1*(a[7] - b1*a[12]) - 2*b3*(a[10] - b1*a[15] - b2*(a[17] - b1*a[20])) - b2*(a[9] - b1*(a[14] - b1*a[18]) - b2*(a[16] - b1*a[19]));
  alpha[5] = a[5] - b2*a[11] - b3*(a[12] - b2*a[18]);
  alpha[6] = a[6] - 2*b1*a[11] - 2*b2*a[13] - b3*(a[14] - 2*b1*a[18] - 2*b2*a[19] - b3*a[20]);
  alpha[7] = a[7] - 2*b1*a[12] - 2*b3*(a[15] - b2*a[20]) - b2*(a[14] - 2*b1*a[18] - b2*a[19]);
  alpha[8] = a[8] - b1*a[13] - b3*(a[16] - b1*a[19]);
  alpha[9] = a[9] - b1*(a[14] - b1*a[18]) - 2*b2*(a[16] - b1*a[19]) - 2*b3*(a[17] - b1*a[20]);
  alpha[10] = a[10] - b1*a[15] - b2*(a[17] - b1*a[20]);
  alpha[11] = a[11] - b3*a[18];
  alpha[12] = a[12] - b2*a[18];
  alpha[13] = a[13] - b3*a[19];
  alpha[14] = a[14] - 2*b1*a[18] - 2*b2*a[19] - 2*b3*a[20];
  alpha[15] = a[15] - b2*a[20];
  alpha[16] = a[16] - b1*a[19];
  alpha[17] = a[17] - b1*a[20];
  alpha[18] = a[18];
  alpha[19] = a[19];
  alpha[20] = a[20];
}

/*%-----------------------------------------------------------
%----- fichier : poly_ch_var.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : changement de variable
%-----
%-----           nxt = a1 xt + b1
%-----           nyt = a2 yt + b2
%-----           nzt = a3 zt + b3
%-----------------------------------------------------------


function [alpha] = poly_rescale(a,a1,b1,a2,b2,a3,b3);*/

void poly_rescale(double *a, double a1, double b1, double a2, double b2, double a3, double b3, double *alpha)
{
//alpha = zeros(20,1);

alpha[1] = a[1] - 1.0/a1*b1*(a[2] - 1.0/a1*b1*a[5]) - 1.0/a2*b2*(a[3] - 1.0/a1*b1*(a[6] - 1.0/a1*b1*a[11]) - 1.0/a2*b2*(a[8] - 1.0/a1*b1*a[13])) - 1.0/a3*b3*(a[4] - 1.0/a1*b1*(a[7] - 1.0/a1*b1*a[12]) - 1.0/a3*b3*(a[10] - 1.0/a1*b1*a[15] - 1.0/a2*b2*(a[17] - 1.0/a1*b1*a[20])) - 1.0/a2*b2*(a[9] - 1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) - 1.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])));
alpha[2] = 1.0/a1*(a[2] - 2.0/a1*b1*a[5]) - 1.0/a2*b2*(1.0/a1*(a[6] - 2.0/a1*b1*a[11]) - 1.0/a1/a2*b2*a[13]) - 1.0/a3*b3*(1/a1*(a[7] - 2.0/a1*b1*a[12]) - 1.0/a3*b3*(1.0/a1*a[15] - 1.0/a1/a2*b2*a[20]) - 1.0/a2*b2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) - 1.0/a1/a2*b2*a[19]));
alpha[3] = 1.0/a2*(a[3] - 1.0/a1*b1*(a[6] - 1.0/a1*b1*a[11]) - 2.0/a2*b2*(a[8] - 1.0/a1*b1*a[13])) - 1.0/a3*b3*(1.0/a2*(a[9] - 1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) - 2.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])) - 1.0/a2/a3*b3*(a[17] - 1.0/a1*b1*a[20]));
alpha[4] = 1.0/a3*(a[4] - 1.0/a1*b1*(a[7] - 1.0/a1*b1*a[12]) - 2.0/a3*b3*(a[10] - 1.0/a1*b1*a[15] - 1.0/a2*b2*(a[17] - 1.0/a1*b1*a[20])) - 1.0/a2*b2*(a[9] - 1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) - 1.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])));
alpha[5] = 1.0/a1/a1*a[5] - 1.0/a1/a1/a2*b2*a[11] - 1.0/a3*b3*(1.0/a1/a1*a[12] - 1.0/a1/a1/a2*b2*a[18]);
alpha[6] = 1.0/a2*(1.0/a1*(a[6] - 2.0/a1*b1*a[11]) - 2.0/a1/a2*b2*a[13]) - 1.0/a3*b3*(1.0/a2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) - 2.0/a1/a2*b2*a[19]) - 1.0/a1/a2/a3*b3*a[20]);
alpha[7] = 1.0/a3*(1.0/a1*(a[7] - 2.0/a1*b1*a[12]) - 2.0/a3*b3*(1.0/a1*a[15] - 1.0/a1/a2*b2*a[20]) - 1.0/a2*b2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) - 1.0/a1/a2*b2*a[19]));
alpha[8] = 1.0/a2/a2*(a[8] - 1.0/a1*b1*a[13]) - 1.0/a2/a2/a3*b3*(a[16] - 1.0/a1*b1*a[19]);
alpha[9] = 1.0/a3*(1.0/a2*(a[9] - 1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) - 2.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])) - 2.0/a2/a3*b3*(a[17] - 1.0/a1*b1*a[20]));
alpha[10] = 1.0/a3/a3*(a[10] - 1.0/a1*b1*a[15] - 1.0/a2*b2*(a[17] - 1.0/a1*b1*a[20]));
alpha[11] = 1.0/a1/a1/a2*a[11] - 1.0/a1/a1/a2/a3*b3*a[18];
alpha[12] = 1.0/a3*(1.0/a1/a1*a[12] - 1.0/a1/a1/a2*b2*a[18]);
alpha[13] = 1.0/a1/a2/a2*a[13] - 1.0/a1/a2/a2/a3*b3*a[19];
alpha[14] = 1.0/a3*(1.0/a2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) - 2.0/a1/a2*b2*a[19]) - 2.0/a1/a2/a3*b3*a[20]);
alpha[15] = 1.0/a3/a3*(1.0/a1*a[15] - 1.0/a1/a2*b2*a[20]);
alpha[16] = 1.0/a2/a2/a3*(a[16] - 1.0/a1*b1*a[19]);
alpha[17] = 1.0/a2/a3/a3*(a[17] - 1.0/a1*b1*a[20]);
alpha[18] = 1.0/(a1*a1)/a2/a3*a[18];
alpha[19] = 1.0/a1/(a2*a2)/a3*a[19];
alpha[20] = 1.0/a1/a2/(a3*a3)*a[20];
/*alpha[1] = a[1] - 1.0/a1*b1*(a[2] - 1.0/a1*b1*a[5]) - 1.0/a2*b2*(a[3] -
1.0/a1*b1*(a[6] - 1.0/a1*b1*a[11]) - 1.0/a2*b2*(a[8] - 1.0/a1*b1*a[13])) -
1.0/a3*b3*(a[4] - 1.0/a1*b1*(a[7] - 1.0/a1*b1*a[12]) - 1.0/a3*b3*(a[10] -
1.0/a1*b1*a[15] - 1.0/a2*b2*(a[17] - 1.0/a1*b1*a[20])) - 1.0/a2*b2*(a[9] -
1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) - 1.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])));
alpha[2] = 1.0/a1*(a[2] - 2.0/a1*b1*a[5]) - 1.0/a2*b2*(1.0/a1*(a[6] -
2.0/a1*b1*a[11]) - 1.0/a1/a2*b2*a[13]) - 1.0/a3*b3*(1.0/a1*(a[7] -
2.0/a1*b1*a[12]) - 1.0/a3*b3*(1.0/a1*a[15] - 1.0/a1/a2*b2*a[20]) -
1.0/a2*b2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) - 1.0/a1/a2*b2*a[19]));
alpha[3] = 1.0/a2*(a[3] - 1.0/a1*b1*(a[6] - 1.0/a1*b1*a[11]) - 2.0/a2*b2*(a[8]
- 1.0/a1*b1*a[13])) - 1.0/a3*b3*(1.0/a2*(a[9] - 1.0/a1*b1*(a[14] -
1.0/a1*b1*a[18]) - 2.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])) - 1.0/a2/a3*b3*(a[17] -
1.0/a1*b1*a[20]));
alpha[4] = 1.0/a3*(a[4] - 1.0/a1*b1*(a[7] - 1.0/a1*b1*a[12]) - 2.0/a3*b3*(a[10]
- 1.0/a1*b1*a[15] - 1.0/a2*b2*(a[17] - 1.0/a1*b1*a[20])) - 1.0/a2*b2*(a[9] -
1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) - 1.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])));
alpha[5] = 1.0/(a1*a1)*a[5] - 1.0/(a1*a1)/a2*b2*a[11] -
1.0/a3*b3*(1.0/(a1*a1)*a[12] - 1.0/(a1*a1)/a2*b2*a[18]);
alpha[6] = 1.0/a2*(1.0/a1*(a[6] - 2.0/a1*b1*a[11]) - 2.0/a1/a2*b2*a[13]) -
1.0/a3*b3*(1.0/a2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) - 2.0/a1/a2*b2*a[19]) -
1.0/a1/a2/a3*b3*a[20]);
alpha[7] = 1.0/a3*(1.0/a1*(a[7] - 2.0/a1*b1*a[12]) - 2.0/a3*b3*(1.0/a1*a[15] -
1.0/a1/a2*b2*a[20]) - 1.0/a2*b2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) -
1.0/a1/a2*b2*a[19]));
alpha[8] = 1.0/(a2*a2)*(a[8] - 1.0/a1*b1*a[13]) - 1.0/(a2*a2)/a3*b3*(a[16] -
1.0/a1*b1*a[19]);
alpha[9] = 1.0/a3*(1.0/a2*(a[9] - 1.0/a1*b1*(a[14] - 1.0/a1*b1*a[18]) -
2.0/a2*b2*(a[16] - 1.0/a1*b1*a[19])) - 2.0/a2/a3*b3*(a[17] - 1.0/a1*b1*a[20]));
alpha[10] = 1.0/(a3*a3)*(a[10] - 1.0/a1*b1*a[15] - 1.0/a2*b2*(a[17] -
1.0/a1*b1*a[20]));
alpha[11] = 1.0/(a1*a1)/a2*a[11] - 1.0/(a1*a1)/a2/a3*b3*a[18];
alpha[12] = 1.0/a3*(1.0/(a1*a1)*a[12] - 1.0/(a1*a1)/a2*b2*a[18]);
alpha[13] = 1.0/a1/(a2*a2)*a[13] - 1.0/a1/(a2*a2)/a3*b3*a[19];
alpha[14] = 1.0/a3*(1.0/a2*(1.0/a1*(a[14] - 2.0/a1*b1*a[18]) -
2.0/a1/a2*b2*a[19]) - 2.0/a1/a2/a3*b3*a[20]);
alpha[15] = 1.0/(a3*a3)*(1.0/a1*a[15] - 1.0/a1/a2*b2*a[20]);
alpha[16] = 1.0/(a2*a2)/a3*(a[16] - 1.0/a1*b1*a[19]);
alpha[17] = 1.0/a2/(a3*a3)*(a[17] - 1.0/a1*b1*a[20]);
alpha[18] = 1.0/(a1*a1)/a2/a3*a[18];
alpha[19] = 1.0/a1/(a2*a2)/a3*a[19];
alpha[20] = 1.0/a1/a2/(a3*a3)*a[20];*/

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : eval_fun.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : évaluation de la fonction associée au
%-----           polynome
%-----------------------------------------------------------


function [myfunc, mymono] = eval_fun(a,xt,yt,zt);*/


void eval_fun(double *a, double xt, double yt, double zt, double *aux, double *mymono)
{ int i;
  //mymono = zeros(size(a));

mymono[1]  = a[1];
mymono[2]  = a[2]*xt;
mymono[3]  = a[3]*yt;
mymono[4]  = a[4]*zt;
mymono[5]  = a[5]*xt*xt;
mymono[6]  = a[6]*xt*yt;
mymono[7]  = a[7]*xt*zt;
mymono[8]  = a[8]*yt*yt;
mymono[9]  = a[9]*yt*zt;
mymono[10] = a[10]*zt*zt;
mymono[11] = a[11]*xt*xt*yt;
mymono[12] = a[12]*xt*xt*zt;
mymono[13] = a[13]*xt*yt*yt;
mymono[14] = a[14]*xt*yt*zt;
mymono[15] = a[15]*xt*zt*zt;
mymono[16] = a[16]*yt*yt*zt;
mymono[17] = a[17]*yt*zt*zt;
mymono[18] = a[18]*xt*xt*yt*zt;
mymono[19] = a[19]*xt*yt*yt*zt;
mymono[20] = a[20]*xt*yt*zt*zt;

*aux = 0;
for (i=1; i<=20; i++)
  *aux += mymono[i];


//myfunc = sum(mymono);
}


void fast_eval_fun(double *a, double xt, double yt, double zt, double *aux, double *mymono)
{ 
*aux=a[1]+zt*(a[4]+a[10]*zt)
		+yt*(a[3]+a[8]*yt+zt*(a[9]+a[16]*yt+a[17]*zt))
		+xt*(a[2]+a[5]*xt+zt*(a[7]+a[12]*xt+a[15]*zt)
		+yt*(a[6]+a[11]*xt+a[13]*yt+zt*(a[14]+a[18]*xt+a[19]*yt+a[20]*zt)));
}


/*%-----------------------------------------------------------
%----- fichier : eval_integrale_fun
%----- objet   : évaluation de l'integrale du jacobien
%-----          
%-----------------------------------------------------------*/

void eval_integrale_fun(double *a, double x0, double y0, double z0, double x1, double y1, double z1, double *aux)
{

double x02=x0*x0;
double y02=y0*y0;
double z02=z0*z0;
double x12=x1*x1;
double y12=y1*y1;
double z12=z1*z1;


double x03=x0*x0*x0;
double y03=y0*y0*y0;
double z03=z0*z0*z0;
double x13=x1*x1*x1;
double y13=y1*y1*y1;
double z13=z1*z1*z1;

double demi=1.0/2.0;
double sixieme=1.0/6.0;
double tiers=1.0/3.0;
double quart=1.0/4.0;
double huitieme=1.0/8.0;
double douzieme=1.0/12.0;

double a1=a[1];
double a2=a[2];
double a3=a[3];
double a4=a[4];
double a5=a[5];
double a6=a[6];
double a7=a[7];
double a8=a[8];
double a9=a[9];
double a10=a[10];
double a11=a[11];
double a12=a[12];
double a13=a[13];
double a14=a[14];
double a15=a[15];
double a16=a[16];
double a17=a[17];
double a18=a[18];
double a19=a[19];
double a20=a[20];


*aux=-sixieme*a13*x12*y13*z0-demi*a4*z02*x0*y0+demi*a4*z02*x1*y0-demi*a4*z02*x1*y1+demi*a4*z02*x0*y1+demi*a4*z12*x0*y0-demi*a4*z12*x0*y1+demi*a4*z12*x1*y1-tiers*a8*y13*x1*z0-sixieme*a11*x13*y12*z0-quart*a6*x12*y12*z0+quart*a6*x02*y12*z0-a1*x1*y1*z0+huitieme*a14*x02*y12*z02+demi*a2*x02*y1*z0-sixieme*a13*x12*y03*z1-tiers*a8*y03*x1*z1-demi*a2*x12*y0*z1+sixieme*a12*x03*z12*y0-sixieme*a12*x13*z12*y0-huitieme*a14*x12*y02*z12-douzieme*a20*x12*y02*z13-sixieme*a16*y03*z12*x1-sixieme*a15*x12*z13*y0+sixieme*a15*x02*z13*y0-douzieme*a19*x12*y03*z12+demi*a2*x12*y1*z1-a1*x0*y1*z1-douzieme*a18*x13*y02*z12+sixieme*a17*y02*z13*x0+tiers*a8*y13*x1*z1-douzieme*a18*x03*y12*z12+douzieme*a19*x12*y13*z12-tiers*a8*y13*x0*z1-quart*a6*x02*y12*z1+quart*a9*y02*z12*x0+huitieme*a14*x12*y12*z12+a1*x1*y1*z1-sixieme*a11*x03*y12*z1+sixieme*a12*x13*z12*y1-quart*a7*x02*z12*y1+tiers*a10*z13*x1*y1-demi*a2*x02*y1*z1+sixieme*a11*x13*y12*z1+quart*a6*x12*y12*z1+douzieme*a20*x02*y02*z13+quart*a7*x02*z12*y0+huitieme*a14*x02*y02*z12+sixieme*a16*y03*z12*x0+douzieme*a20*x12*y12*z13+sixieme*a16*y13*z12*x1+sixieme*a17*y12*z13*x1+quart*a7*x12*z12*y1+sixieme*a15*x12*z13*y1+tiers*a10*z13*x0*y0-sixieme*a11*x13*y02*z1+tiers*a5*x03*y0*z1-a1*x1*y0*z1+sixieme*a11*x03*y02*z1-quart*a6*x12*y02*z1-demi*a3*y02*x1*z1+sixieme*a13*x02*y03*z1+a1*x0*y0*z1+tiers*a8*y03*x0*z1-quart*a7*x12*z12*y0+douzieme*a18*x03*y02*z12-tiers*a10*z13*x0*y1-sixieme*a12*x03*z12*y1-sixieme*a15*x02*z13*y1-sixieme*a13*x02*y13*z1+demi*a3*y12*x1*z1+sixieme*a13*x12*y13*z1-tiers*a5*x13*y0*z1-demi*a3*y12*x0*z1+tiers*a5*x13*y1*z1-tiers*a5*x03*y1*z1+quart*a6*x02*y02*z1-quart*a9*y12*z12*x0+quart*a9*y12*z12*x1+douzieme*a18*x13*y12*z12-douzieme*a19*x02*y13*z12+demi*a2*x02*y0*z1-tiers*a10*z13*x1*y0+douzieme*a19*x02*y03*z12+demi*a3*y02*x0*z1-sixieme*a17*y02*z13*x1-quart*a9*y02*z12*x1-sixieme*a16*y13*z12*x0-sixieme*a17*y12*z13*x0-huitieme*a14*x02*y12*z12-douzieme*a20*x02*y12*z13+douzieme*a19*x02*y13*z02-huitieme*a14*x12*y12*z02+douzieme*a18*x03*y12*z02-douzieme*a19*x12*y13*z02+douzieme*a20*x02*y12*z03+quart*a9*y12*z02*x0-quart*a9*y12*z02*x1+huitieme*a14*x12*y02*z02+douzieme*a18*x13*y02*z02+sixieme*a17*y02*z03*x1-quart*a9*y02*z02*x0+quart*a7*x12*z02*y0+tiers*a10*z03*x1*y0-douzieme*a19*x02*y03*z02-demi*a3*y02*x0*z0+sixieme*a11*x03*y12*z0-douzieme*a20*x12*y12*z03-sixieme*a12*x03*z02*y0-sixieme*a12*x13*z02*y1+quart*a7*x02*z02*y1-tiers*a10*z03*x1*y1+quart*a6*x12*y02*z0+demi*a3*y02*x1*z0-sixieme*a13*x02*y03*z0-a1*x0*y0*z0+tiers*a8*y13*x0*z0+demi*a3*y12*x0*z0-tiers*a5*x03*y0*z0-sixieme*a16*y13*z02*x1-sixieme*a17*y12*z03*x1+sixieme*a15*x12*z03*y0-sixieme*a15*x02*z03*y0+douzieme*a19*x12*y03*z02-tiers*a5*x13*y1*z0+tiers*a5*x03*y1*z0-quart*a6*x02*y02*z0+sixieme*a11*x13*y02*z0-sixieme*a11*x03*y02*z0+sixieme*a17*y12*z03*x0-douzieme*a18*x13*y12*z02+tiers*a10*z03*x0*y1+sixieme*a12*x03*z02*y1+sixieme*a15*x02*z03*y1+sixieme*a13*x02*y13*z0-demi*a3*y12*x1*z0+quart*a9*y02*z02*x1+sixieme*a16*y13*z02*x0+sixieme*a12*x13*z02*y0+tiers*a5*x13*y0*z0+a1*x1*y0*z0-tiers*a8*y03*x0*z0+a1*x0*y1*z0-douzieme*a18*x03*y02*z02+demi*a2*x12*y0*z0-demi*a2*x02*y0*z0+sixieme*a13*x12*y03*z0+tiers*a8*y03*x1*z0-quart*a7*x12*z02*y1-sixieme*a15*x12*z03*y1-douzieme*a20*x02*y02*z03-tiers*a10*z03*x0*y0-quart*a7*x02*z02*y0+sixieme*a16*y03*z02*x1-huitieme*a14*x02*y02*z02-sixieme*a16*y03*z02*x0-demi*a2*x12*y1*z0+douzieme*a20*x12*y02*z03-sixieme*a17*y02*z03*x0-demi*a4*z12*x1*y0;


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TOPTbox depile_end(TOPTliste *myListe)
{ TOPTbox myBox, *myPtr; unsigned int i;

  if (myListe->nelem == 1)
  { myBox = *(myListe->tete);
    free(myListe->tete);
    myListe->tete = NULL;
  }
  else
  { myPtr = myListe->tete;
    for (i=1; i<=myListe->nelem-2; i++) myPtr = myPtr->next;
    myBox = *(myPtr->next);
    free(myPtr->next);
    myPtr->next = NULL;
  }

  myListe->nelem--;
  return myBox;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void empile(TOPTliste *myListe, TOPTbox myB)
{ TOPTbox *ptr;

  myB.next = NULL;

  if (myListe->nelem ==0)
  { myListe->tete = (TOPTbox*) malloc(sizeof(TOPTbox));
    *myListe->tete = myB;
  }
  else
  { ptr = myListe->tete;
    while (ptr->next != NULL) ptr = ptr->next;
    ptr->next = (TOPTbox*) malloc(sizeof(TOPTbox));
    ptr = ptr->next;
    *ptr = myB;
  }
  myListe->nelem++;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void empile_front(TOPTliste *myListe, TOPTbox myB)
{
  myB.next = myListe->tete;

  if (myListe->nelem ==0)
  { myListe->tete = (TOPTbox*) malloc(sizeof(TOPTbox));
		 myB.next=NULL;
    *myListe->tete = myB;
		
  }
  else
  {  myListe->tete = (TOPTbox*) malloc(sizeof(TOPTbox));
     *myListe->tete  = myB;
  }
  myListe->nelem++;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : minim_iter.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : procédure de minimisation
%-----           minimisation itérative
%-----
%-----           Gauss-Seidel car quadratique en les
%-----           marginales
%-----
%-----           coin : indique si le minimum est atteint en un coin
%-----
%-----           Voir les commentaires en fin de procédure
%-----------------------------------------------------------


function [x_hat,y_hat,z_hat,coin] = minim_iter(precis,a,x0,y0,z0,xinf,xsup,yinf,ysup,zinf,zsup); */

void minim_iter(double precis, double *a, TOPTbox myB, double *x_hat, double *y_hat, double *z_hat)
{ double poly[4];    			// c'est un poly_marg
  double xinf, yinf, zinf, xsup, ysup, zsup;
  double x0, y0, z0;
  double x1, y1, z1;
  int continuer;

  xinf = myB.xm; yinf = myB.ym; zinf = myB.zm;
  xsup = myB.xM; ysup = myB.yM; zsup = myB.zM;
  x0 = (xinf+xsup) / 2.0;
  y0 = (yinf+ysup) / 2.0;
  z0 = (zinf+zsup) / 2.0;

  //-----------------------------------------------------------
  //-----------------------------------------------------------
  //-----------------------------------------------------------
  continuer = 1;

  while (continuer==1)
  {
    if (xsup-xinf>=precis)
    { poly_marg(1,a,x0,y0,z0,poly);
      x1 = minim_min_ligne(poly,xinf,xsup);
    }
    else
      x1 = x0;

    if (ysup-yinf>=precis)
    { poly_marg(2,a,x1,y0,z0,poly);
      y1 = minim_min_ligne(poly,yinf,ysup);
    }
    else
      y1 = y0;

    if (zsup-zinf>=precis)
    { poly_marg(3,a,x1,y1,z0,poly);
      z1 = minim_min_ligne(poly,zinf,zsup);
    }
    else
      z1 = z0;

    if ( (fabs(x1-x0)<precis/10.0) && (fabs(y1-y0)<precis/10.0) && (fabs(z1-z0)<precis/10.0) )
      continuer = 0;

    x0 = x1; y0 = y1; z0 = z1;
  }

  //----------------------------------------------------------
  *x_hat = x0; *y_hat = y0; *z_hat = z0;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : minim_min_ligne.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : procédure de minimisation
%-----           minimisation de ligne
%-----------------------------------------------------------


function [myarg , extr] = minim_min_ligne(poly_marg,vinf,vsup); */

double minim_min_ligne(double *poly_marg,double vinf, double vsup)
{ double a, b, c;
  double finf, fsup, faux=0, vaux=0;
  int fauxinit;

  //-----------------------------------------------------------
  //-----------------------------------------------------------
  //-----------------------------------------------------------
  a = poly_marg[1];
  b = poly_marg[2];
  c = poly_marg[3];

/*  finf = a * vinf * vinf + b * vinf + c;
  fsup = a * vsup * vsup + b * vsup + c;*/
  finf = a * vinf * vinf + b * vinf;
  fsup = a * vsup * vsup + b * vsup;

  //  faux = nan;
  fauxinit = 0;

  if (a!=0.0)
  {  vaux = - b / 2.0 / a;
     if ((vaux>vinf)&&(vaux<vsup))
     { 
		 	/*faux = a * vaux * vaux + b * vaux + c;*/
			faux = a * vaux * vaux + b * vaux;
       fauxinit = 1;
     }
  }

  if (fauxinit==0)
  { 
	
	if (finf<=fsup)
         return vinf;
    else return vsup;
  }
  else
  { if ((finf<=fsup)&&(finf<=faux))
      return vinf;
    else if (faux<=fsup)
      return vaux;
    else
      return vsup;
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/


//function [lout] = vecteur_insere(lin,val,precis);

void vecteur_insere(int *ncomp, double *lin, double val, double precis)
{ int i; double aux;

  if ( (val>lin[0]) && (val<lin[*ncomp-1]) && (lin[*ncomp-1]-lin[0]>=precis) )
  { *ncomp = *ncomp + 1;
    lin[*ncomp-1] = val;
    for (i=*ncomp-2; (i>=0)&&(lin[i]>lin[i+1]); i--)
    { aux = lin[i];
      lin[i] = lin[i+1];
      lin[i+1] = aux;
    }
  }

}
/*lout = lin;
dim = length(lin);

if ( (val>lin(1)) & (val<lin(dim)) & (lin(dim)-lin(1)>=precis) )
  lout = sort([lin,val]);
end
*/


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : pave_casser.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : retourne une liste de pavés resultant de
%-----           casser myB en xcass, ycass, zcass
%-----           bornes non initialisées
%-----------------------------------------------------------


function [myL] = pave_casser(precis,myB,xcass,ycass,zcass);*/

// ATTENTION : xcass, ycass, zcass de dimension 2 nécessairement


void pave_casser(double precis, TOPTbox myB, double *xcass, double *ycass, double *zcass, TOPTliste *myL)
{ double ix[4], iy[4], iz[4];
  int nix, niy, niz;
  int i, j, k;
  TOPTbox myBaux;

  //--------------------------------------------------------
  nix = 2; ix[0] = myB.xm; ix[1] = myB.xM;
  vecteur_insere(&nix,ix,xcass[0],precis);
  vecteur_insere(&nix,ix,xcass[1],precis);

  //--------------------------------------------------------
  niy = 2; iy[0] = myB.ym; iy[1] = myB.yM;
  vecteur_insere(&niy,iy,ycass[0],precis);
  vecteur_insere(&niy,iy,ycass[1],precis);

  //--------------------------------------------------------
  niz = 2; iz[0] = myB.zm; iz[1] = myB.zM;
  vecteur_insere(&niz,iz,zcass[0],precis);
  vecteur_insere(&niz,iz,zcass[1],precis);

  //--------------------------------------------------------
  for (i=0; i<=nix-2; i++)
  { for (j=0; j<=niy-2; j++)
    { for (k=0; k<=niz-2; k++)
      { myBaux.xm = ix[i]; myBaux.xM = ix[i+1];
        myBaux.ym = iy[j]; myBaux.yM = iy[j+1];
        myBaux.zm = iz[k]; myBaux.zM = iz[k+1];
        myBaux.fm = EDOM;  myBaux.fM = EDOM;
        empile(myL,myBaux);
      }
    }
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : poly_pente_coeff.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : donne la pente (en fonction du pas) de
%-----           chaque coefficient
%-----
%-----           dir_des : vecteur colonne de dimension 3, normé
%-----           indice  : entier, entre 1 et 8 (i de alxi)
%-----------------------------------------------------------*/

// function [pente] = poly_pente_coeff(dir_des,indice,alx,aly,alz,l);

void TOP_poly_pente_coeff(double *dir_des, int indice, double *alx, double *aly, double *alz, int l, double *pente)
{
double big_jac[21][25];
double alx1,alx2,alx3,alx4,alx5,alx6,alx7,alx8;
double aly1,aly2,aly3,aly4,aly5,aly6,aly7,aly8;
double alz1,alz2,alz3,alz4,alz5,alz6,alz7,alz8;
int i,j;


//myco = 2^(-5*l/2);
// double myco = pow(2.0,-5.0*l/2);
/* Convention O.Musse : facteur de la fonction d'echelle = 1 */
//l=0;
double myco=pow(2.0,-1.0*l);

alx1 = alx[0];
alx2 = alx[1];
alx3 = alx[2];
alx4 = alx[3];
alx5 = alx[4];
alx6 = alx[5];
alx7 = alx[6];
alx8 = alx[7];

aly1 = aly[0];
aly2 = aly[1];
aly3 = aly[2];
aly4 = aly[3];
aly5 = aly[4];
aly6 = aly[5];
aly7 = aly[6];
aly8 = aly[7];

alz1 = alz[0];
alz2 = alz[1];
alz3 = alz[2];
alz4 = alz[3];
alz5 = alz[4];
alz6 = alz[5];
alz7 = alz[6];
alz8 = alz[7];

switch(indice)
{ case 1:

big_jac[1 ][ 1] = aly2*alz3 - aly3*alz2 - aly2*alz5 + alz2*aly5 + aly3*alz5 - alz3*aly5 - aly3*myco - alz2*myco + aly5*myco + alz5*myco - myco*myco;
big_jac[2 ][ 1] = 2*aly3*alz2 - 2*aly2*alz3 + aly2*alz5 - alz2*aly5 - aly3*alz5 + alz3*aly5 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 + aly5*alz6 - aly6*alz5 - aly5*alz7 + alz5*aly7 + aly3*myco + alz2*myco - aly7*myco - alz6*myco;
big_jac[3 ][ 1] = aly3*alz2 - aly2*alz3 + 2*aly2*alz5 - aly3*alz4 - 2*alz2*aly5 + aly4*alz3 - aly3*alz5 + alz3*aly5 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 + aly3*alz7 - alz3*aly7 + 2*alz2*myco - aly5*myco - alz4*myco - 2*alz5*myco + aly7*myco + alz7*myco + myco*myco;
big_jac[4 ][ 1] = aly3*alz2 - aly2*alz3 + aly2*alz4 - alz2*aly4 + aly2*alz5 - alz2*aly5 - aly2*alz6 - 2*aly3*alz5 + alz2*aly6 + 2*alz3*aly5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 + 2*aly3*myco - aly4*myco - 2*aly5*myco + aly6*myco - alz5*myco + alz6*myco + myco*myco;
big_jac[5 ][ 1] = aly2*alz3 - aly3*alz2 - aly2*alz7 + aly3*alz6 + alz2*aly7 - alz3*aly6 + aly6*alz7 - aly7*alz6;
big_jac[6 ][ 1] = 2*aly2*alz3 - 2*aly3*alz2 - 2*aly2*alz5 + 2*aly3*alz4 + 2*alz2*aly5 - 2*aly4*alz3 + aly3*alz5 - alz3*aly5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 - aly3*alz7 + alz3*aly7 - aly3*alz8 + aly4*alz7 + alz3*aly8 - 2*aly5*alz6 - alz4*aly7 + 2*aly6*alz5 + aly5*alz7 - alz5*aly7 + aly5*alz8 - aly6*alz7 - alz5*aly8 + aly7*alz6 - 2*alz2*myco + alz4*myco + 2*alz6*myco - alz8*myco;
big_jac[7 ][ 1] = 2*aly2*alz3 - 2*aly3*alz2 - 2*aly2*alz4 + 2*alz2*aly4 - aly2*alz5 + alz2*aly5 + aly2*alz6 + 2*aly3*alz5 - alz2*aly6 - 2*alz3*aly5 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 + aly2*alz8 - alz2*aly8 - aly4*alz6 + alz4*aly6 - aly5*alz6 + aly6*alz5 + 2*aly5*alz7 - 2*alz5*aly7 - aly5*alz8 - aly6*alz7 + alz5*aly8 + aly7*alz6 - 2*aly3*myco + aly4*myco + 2*aly7*myco - aly8*myco;
big_jac[8 ][ 1] = alz2*aly5 - aly2*alz5 + aly2*alz7 - alz2*aly7 + aly4*alz5 - aly5*alz4 - aly4*alz7 + alz4*aly7 - alz2*myco + alz4*myco + alz5*myco - alz7*myco;
big_jac[9 ][ 1] = aly2*alz3 - aly3*alz2 - aly2*alz4 + alz2*aly4 - 2*aly2*alz5 + aly3*alz4 + 2*alz2*aly5 - aly4*alz3 + 2*aly2*alz6 + 2*aly3*alz5 - 2*alz2*aly6 - 2*alz3*aly5 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 - aly2*alz8 - 2*aly3*alz7 + alz2*aly8 - aly4*alz6 + 2*alz3*aly7 + alz4*aly6 + aly3*alz8 + aly4*alz7 - alz3*aly8 - alz4*aly7 + 2*aly5*myco - aly6*myco + 2*alz5*myco - 2*aly7*myco - 2*alz6*myco + aly8*myco - alz7*myco + alz8*myco - myco*myco;
big_jac[10 ][ 1] = aly3*alz5 - alz3*aly5 - aly3*alz6 - aly4*alz5 + alz3*aly6 + aly5*alz4 + aly4*alz6 - alz4*aly6 - aly3*myco + aly4*myco + aly5*myco - aly6*myco;
big_jac[11 ][ 1] = aly3*alz2 - aly2*alz3 - aly3*alz4 + aly4*alz3 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 + aly3*alz8 - aly4*alz7 - alz3*aly8 + alz4*aly7 - aly6*alz7 + aly7*alz6 - aly7*alz8 + aly8*alz7;
big_jac[12 ][ 1] = aly3*alz2 - aly2*alz3 + aly2*alz4 - alz2*aly4 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 - aly2*alz8 + alz2*aly8 + aly4*alz6 - alz4*aly6 - aly6*alz7 + aly7*alz6 + aly6*alz8 - alz6*aly8;
big_jac[13 ][ 1] = aly2*alz5 - alz2*aly5 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 + aly4*alz7 + aly5*alz6 - alz4*aly7 - aly6*alz5 - aly5*alz8 + aly6*alz7 + alz5*aly8 - aly7*alz6 + aly7*alz8 - aly8*alz7 + alz2*myco - alz4*myco - alz6*myco + alz8*myco;
big_jac[14 ][ 1] = 2*aly3*alz2 - 2*aly2*alz3 + 2*aly2*alz4 - 2*alz2*aly4 + 2*aly2*alz5 - 2*aly3*alz4 - 2*alz2*aly5 + 2*aly4*alz3 - 2*aly2*alz6 - 2*aly3*alz5 + 2*alz2*aly6 + 2*alz3*aly5 + 2*aly3*alz7 + 2*aly4*alz6 - 2*alz3*aly7 - 2*alz4*aly6 - 2*aly4*alz7 + 2*aly5*alz6 + 2*alz4*aly7 - 2*aly6*alz5 - 2*aly5*alz7 + 2*alz5*aly7 + 2*aly6*alz7 - 2*aly7*alz6;
big_jac[15 ][ 1] = alz3*aly5 - aly3*alz5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 - aly4*alz6 + alz4*aly6 - aly5*alz7 + alz5*aly7 + aly5*alz8 + aly6*alz7 - alz5*aly8 - aly7*alz6 - aly6*alz8 + alz6*aly8 + aly3*myco - aly4*myco - aly7*myco + aly8*myco;
big_jac[16 ][ 1] = aly2*alz5 - alz2*aly5 - aly2*alz6 + alz2*aly6 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 + aly2*alz8 - alz2*aly8 + aly4*alz6 - alz4*aly6 + aly4*alz7 - alz4*aly7 - aly4*alz8 + alz4*aly8 - alz5*myco + alz6*myco + alz7*myco - alz8*myco;
big_jac[17 ][ 1] = alz3*aly5 - aly3*alz5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 + aly3*alz7 - aly4*alz6 - alz3*aly7 + alz4*aly6 - aly3*alz8 - aly4*alz7 + alz3*aly8 + alz4*aly7 + aly4*alz8 - alz4*aly8 - aly5*myco + aly6*myco + aly7*myco - aly8*myco;
big_jac[18 ][ 1] = aly2*alz3 - aly3*alz2 - aly2*alz4 + alz2*aly4 + aly3*alz4 - aly4*alz3 - aly2*alz7 + aly3*alz6 + alz2*aly7 - alz3*aly6 + aly2*alz8 - alz2*aly8 - aly4*alz6 + alz4*aly6 - aly3*alz8 + aly4*alz7 + alz3*aly8 - alz4*aly7 + aly6*alz7 - aly7*alz6 - aly6*alz8 + alz6*aly8 + aly7*alz8 - aly8*alz7;
big_jac[19 ][ 1] = alz2*aly5 - aly2*alz5 + aly2*alz6 - alz2*aly6 + aly2*alz7 - alz2*aly7 + aly4*alz5 - aly5*alz4 - aly2*alz8 + alz2*aly8 - aly4*alz6 + alz4*aly6 - aly4*alz7 - aly5*alz6 + alz4*aly7 + aly6*alz5 + aly4*alz8 - alz4*aly8 + aly5*alz8 - aly6*alz7 - alz5*aly8 + aly7*alz6 - aly7*alz8 + aly8*alz7;
big_jac[20 ][ 1] = aly3*alz5 - alz3*aly5 - aly3*alz6 - aly4*alz5 + alz3*aly6 + aly5*alz4 - aly3*alz7 + aly4*alz6 + alz3*aly7 - alz4*aly6 + aly3*alz8 + aly4*alz7 - alz3*aly8 - alz4*aly7 - aly4*alz8 + aly5*alz7 + alz4*aly8 - alz5*aly7 - aly5*alz8 - aly6*alz7 + alz5*aly8 + aly7*alz6 + aly6*alz8 - alz6*aly8;
big_jac[1 ][ 2] = alx3*alz2 - alx2*alz3 + alx2*alz5 - alz2*alx5 - alx3*alz5 + alx5*alz3 + alx3*myco - alz2*myco - alx5*myco + alz3*myco - myco*myco;
big_jac[2 ][ 2] = 2*alx2*alz3 - 2*alx3*alz2 - alx2*alz5 + alz2*alx5 + alx3*alz5 - alx5*alz3 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 - alx5*alz6 + alx6*alz5 + alx5*alz7 - alx7*alz5 - alx3*myco + 2*alz2*myco - 2*alz3*myco + alx7*myco - alz6*myco + alz7*myco + myco*myco;
big_jac[3 ][ 2] = alx2*alz3 - alx3*alz2 - 2*alx2*alz5 + alx3*alz4 - alx4*alz3 + 2*alz2*alx5 + alx3*alz5 - alx5*alz3 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 - alx3*alz7 + alz3*alx7 + alz2*myco + alx5*myco - alz4*myco - alx7*myco;
big_jac[4 ][ 2] = alx2*alz3 - alx3*alz2 - alx2*alz4 + alx4*alz2 - alx2*alz5 + alz2*alx5 + alx2*alz6 + 2*alx3*alz5 - alz2*alx6 - 2*alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 - 2*alx3*myco + alx4*myco + 2*alx5*myco - alz3*myco - alx6*myco + alz4*myco + myco*myco;
big_jac[5 ][ 2] = alx3*alz2 - alx2*alz3 + alx2*alz7 - alx3*alz6 - alz2*alx7 + alz3*alx6 - alx6*alz7 + alx7*alz6 - alz2*myco + alz3*myco + alz6*myco - alz7*myco;
big_jac[6 ][ 2] = 2*alx3*alz2 - 2*alx2*alz3 + 2*alx2*alz5 - 2*alx3*alz4 + 2*alx4*alz3 - 2*alz2*alx5 - alx3*alz5 + alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 + alx3*alz7 - alz3*alx7 + alx3*alz8 - alx4*alz7 + 2*alx5*alz6 - alz3*alx8 - 2*alx6*alz5 + alz4*alx7 - alx5*alz7 + alx7*alz5 - alx5*alz8 + alx6*alz7 - alx7*alz6 + alz5*alx8 - 2*alz2*myco + 2*alz4*myco + alz6*myco - alz8*myco;
big_jac[7 ][ 2] = 2*alx3*alz2 - 2*alx2*alz3 + 2*alx2*alz4 - 2*alx4*alz2 + alx2*alz5 - alz2*alx5 - alx2*alz6 - 2*alx3*alz5 + alz2*alx6 + 2*alx5*alz3 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 - alx2*alz8 + alx4*alz6 + alz2*alx8 - alx6*alz4 + alx5*alz6 - alx6*alz5 - 2*alx5*alz7 + 2*alx7*alz5 + alx5*alz8 + alx6*alz7 - alx7*alz6 - alz5*alx8 + 2*alx3*myco - alx4*myco + 2*alz3*myco - 2*alz4*myco - 2*alx7*myco + alx8*myco - alz7*myco + alz8*myco - myco*myco;
big_jac[8 ][ 2] = alx2*alz5 - alz2*alx5 - alx2*alz7 - alx4*alz5 + alz2*alx7 + alx5*alz4 + alx4*alz7 - alz4*alx7;
big_jac[9 ][ 2] = alx3*alz2 - alx2*alz3 + alx2*alz4 - alx4*alz2 + 2*alx2*alz5 - alx3*alz4 + alx4*alz3 - 2*alz2*alx5 - 2*alx2*alz6 - 2*alx3*alz5 + 2*alz2*alx6 + 2*alx5*alz3 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 + alx2*alz8 + 2*alx3*alz7 + alx4*alz6 - alz2*alx8 - 2*alz3*alx7 - alx6*alz4 - alx3*alz8 - alx4*alz7 + alz3*alx8 + alz4*alx7 - 2*alx5*myco + alx6*myco + 2*alx7*myco - alx8*myco;
big_jac[10 ][ 2] = alx5*alz3 - alx3*alz5 + alx3*alz6 + alx4*alz5 - alx5*alz4 - alz3*alx6 - alx4*alz6 + alx6*alz4 + alx3*myco - alx4*myco - alx5*myco + alx6*myco;
big_jac[11 ][ 2] = alx2*alz3 - alx3*alz2 + alx3*alz4 - alx4*alz3 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 - alx3*alz8 + alx4*alz7 + alz3*alx8 - alz4*alx7 + alx6*alz7 - alx7*alz6 + alx7*alz8 - alx8*alz7 + alz2*myco - alz4*myco - alz6*myco + alz8*myco;
big_jac[12 ][ 2] = alx2*alz3 - alx3*alz2 - alx2*alz4 + alx4*alz2 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 + alx2*alz8 - alx4*alz6 - alz2*alx8 + alx6*alz4 + alx6*alz7 - alx7*alz6 - alx6*alz8 + alx8*alz6 - alz3*myco + alz4*myco + alz7*myco - alz8*myco;
big_jac[13 ][ 2] = alz2*alx5 - alx2*alz5 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 - alx4*alz7 - alx5*alz6 + alx6*alz5 + alz4*alx7 + alx5*alz8 - alx6*alz7 + alx7*alz6 - alz5*alx8 - alx7*alz8 + alx8*alz7;
big_jac[14 ][ 2] = 2*alx2*alz3 - 2*alx3*alz2 - 2*alx2*alz4 + 2*alx4*alz2 - 2*alx2*alz5 + 2*alx3*alz4 - 2*alx4*alz3 + 2*alz2*alx5 + 2*alx2*alz6 + 2*alx3*alz5 - 2*alz2*alx6 - 2*alx5*alz3 - 2*alx3*alz7 - 2*alx4*alz6 + 2*alz3*alx7 + 2*alx6*alz4 + 2*alx4*alz7 - 2*alx5*alz6 + 2*alx6*alz5 - 2*alz4*alx7 + 2*alx5*alz7 - 2*alx7*alz5 - 2*alx6*alz7 + 2*alx7*alz6;
big_jac[15 ][ 2] = alx3*alz5 - alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 + alx4*alz6 - alx6*alz4 + alx5*alz7 - alx7*alz5 - alx5*alz8 - alx6*alz7 + alx7*alz6 + alz5*alx8 + alx6*alz8 - alx8*alz6 - alx3*myco + alx4*myco + alx7*myco - alx8*myco;
big_jac[16 ][ 2] = alz2*alx5 - alx2*alz5 + alx2*alz6 - alz2*alx6 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 - alx2*alz8 - alx4*alz6 + alz2*alx8 + alx6*alz4 - alx4*alz7 + alz4*alx7 + alx4*alz8 - alz4*alx8;
big_jac[17 ][ 2] = alx3*alz5 - alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 - alx3*alz7 + alx4*alz6 + alz3*alx7 - alx6*alz4 + alx3*alz8 + alx4*alz7 - alz3*alx8 - alz4*alx7 - alx4*alz8 + alz4*alx8 + alx5*myco - alx6*myco - alx7*myco + alx8*myco;
big_jac[18 ][ 2] = alx3*alz2 - alx2*alz3 + alx2*alz4 - alx4*alz2 - alx3*alz4 + alx4*alz3 + alx2*alz7 - alx3*alz6 - alz2*alx7 + alz3*alx6 - alx2*alz8 + alx4*alz6 + alz2*alx8 - alx6*alz4 + alx3*alz8 - alx4*alz7 - alz3*alx8 + alz4*alx7 - alx6*alz7 + alx7*alz6 + alx6*alz8 - alx8*alz6 - alx7*alz8 + alx8*alz7;
big_jac[19 ][ 2] = alx2*alz5 - alz2*alx5 - alx2*alz6 + alz2*alx6 - alx2*alz7 - alx4*alz5 + alz2*alx7 + alx5*alz4 + alx2*alz8 + alx4*alz6 - alz2*alx8 - alx6*alz4 + alx4*alz7 + alx5*alz6 - alx6*alz5 - alz4*alx7 - alx4*alz8 + alz4*alx8 - alx5*alz8 + alx6*alz7 - alx7*alz6 + alz5*alx8 + alx7*alz8 - alx8*alz7;
big_jac[20 ][ 2] = alx5*alz3 - alx3*alz5 + alx3*alz6 + alx4*alz5 - alx5*alz4 - alz3*alx6 + alx3*alz7 - alx4*alz6 - alz3*alx7 + alx6*alz4 - alx3*alz8 - alx4*alz7 + alz3*alx8 + alz4*alx7 + alx4*alz8 - alx5*alz7 - alz4*alx8 + alx7*alz5 + alx5*alz8 + alx6*alz7 - alx7*alz6 - alz5*alx8 - alx6*alz8 + alx8*alz6;
big_jac[1 ][ 3] = alx2*aly3 - alx3*aly2 - alx2*aly5 + aly2*alx5 + alx3*aly5 - aly3*alx5 + alx2*myco + aly2*myco - aly3*myco - alx5*myco - myco*myco;
big_jac[2 ][ 3] = 2*alx3*aly2 - 2*alx2*aly3 + alx2*aly5 - aly2*alx5 - alx3*aly5 + aly3*alx5 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 + alx5*aly6 - alx6*aly5 - alx5*aly7 + aly5*alx7 - alx2*myco - 2*aly2*myco + 2*aly3*myco + alx6*myco + aly6*myco - aly7*myco + myco*myco;
big_jac[3 ][ 3] = alx3*aly2 - alx2*aly3 + 2*alx2*aly5 - alx3*aly4 - 2*aly2*alx5 + alx4*aly3 - alx3*aly5 + aly3*alx5 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 + alx3*aly7 - aly3*alx7 - 2*alx2*myco - aly2*myco + alx4*myco + 2*alx5*myco + aly4*myco - alx7*myco + myco*myco;
big_jac[4 ][ 3] = alx3*aly2 - alx2*aly3 + alx2*aly4 - aly2*alx4 + alx2*aly5 - aly2*alx5 - alx2*aly6 - 2*alx3*aly5 + aly2*alx6 + 2*aly3*alx5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 + aly3*myco + alx5*myco - aly4*myco - alx6*myco;
big_jac[5 ][ 3] = alx2*aly3 - alx3*aly2 - alx2*aly7 + alx3*aly6 + aly2*alx7 - aly3*alx6 + alx6*aly7 - alx7*aly6 + aly2*myco - aly3*myco - aly6*myco + aly7*myco;
big_jac[6 ][ 3] = 2*alx2*aly3 - 2*alx3*aly2 - 2*alx2*aly5 + 2*alx3*aly4 + 2*aly2*alx5 - 2*alx4*aly3 + alx3*aly5 - aly3*alx5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 - alx3*aly7 + aly3*alx7 - alx3*aly8 + alx4*aly7 + aly3*alx8 - 2*alx5*aly6 - aly4*alx7 + 2*alx6*aly5 + alx5*aly7 - aly5*alx7 + alx5*aly8 - alx6*aly7 - aly5*alx8 + alx7*aly6 + 2*alx2*myco + 2*aly2*myco - alx4*myco - 2*aly4*myco - 2*alx6*myco - aly6*myco + alx8*myco + aly8*myco - myco*myco;
big_jac[7 ][ 3] = 2*alx2*aly3 - 2*alx3*aly2 - 2*alx2*aly4 + 2*aly2*alx4 - alx2*aly5 + aly2*alx5 + alx2*aly6 + 2*alx3*aly5 - aly2*alx6 - 2*aly3*alx5 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 + alx2*aly8 - aly2*alx8 - alx4*aly6 + aly4*alx6 - alx5*aly6 + alx6*aly5 + 2*alx5*aly7 - 2*aly5*alx7 - alx5*aly8 - alx6*aly7 + aly5*alx8 + alx7*aly6 - 2*aly3*myco + 2*aly4*myco + aly7*myco - aly8*myco;
big_jac[8 ][ 3] = aly2*alx5 - alx2*aly5 + alx2*aly7 - aly2*alx7 + alx4*aly5 - alx5*aly4 - alx4*aly7 + aly4*alx7 + alx2*myco - alx4*myco - alx5*myco + alx7*myco;
big_jac[9 ][ 3] = alx2*aly3 - alx3*aly2 - alx2*aly4 + aly2*alx4 - 2*alx2*aly5 + alx3*aly4 + 2*aly2*alx5 - alx4*aly3 + 2*alx2*aly6 + 2*alx3*aly5 - 2*aly2*alx6 - 2*aly3*alx5 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 - alx2*aly8 - 2*alx3*aly7 + aly2*alx8 - alx4*aly6 + 2*aly3*alx7 + aly4*alx6 + alx3*aly8 + alx4*aly7 - aly3*alx8 - aly4*alx7 - 2*alx5*myco + 2*alx6*myco + alx7*myco - alx8*myco;
big_jac[10 ][ 3] = alx3*aly5 - aly3*alx5 - alx3*aly6 - alx4*aly5 + aly3*alx6 + alx5*aly4 + alx4*aly6 - aly4*alx6;
big_jac[11 ][ 3] = alx3*aly2 - alx2*aly3 - alx3*aly4 + alx4*aly3 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 + alx3*aly8 - alx4*aly7 - aly3*alx8 + aly4*alx7 - alx6*aly7 + alx7*aly6 - alx7*aly8 + alx8*aly7 - aly2*myco + aly4*myco + aly6*myco - aly8*myco;
big_jac[12 ][ 3] = alx3*aly2 - alx2*aly3 + alx2*aly4 - aly2*alx4 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 - alx2*aly8 + aly2*alx8 + alx4*aly6 - aly4*alx6 - alx6*aly7 + alx7*aly6 + alx6*aly8 - aly6*alx8 + aly3*myco - aly4*myco - aly7*myco + aly8*myco;
big_jac[13 ][ 3] = alx2*aly5 - aly2*alx5 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 + alx4*aly7 + alx5*aly6 - aly4*alx7 - alx6*aly5 - alx5*aly8 + alx6*aly7 + aly5*alx8 - alx7*aly6 + alx7*aly8 - alx8*aly7 - alx2*myco + alx4*myco + alx6*myco - alx8*myco;
big_jac[14 ][ 3] = 2*alx3*aly2 - 2*alx2*aly3 + 2*alx2*aly4 - 2*aly2*alx4 + 2*alx2*aly5 - 2*alx3*aly4 - 2*aly2*alx5 + 2*alx4*aly3 - 2*alx2*aly6 - 2*alx3*aly5 + 2*aly2*alx6 + 2*aly3*alx5 + 2*alx3*aly7 + 2*alx4*aly6 - 2*aly3*alx7 - 2*aly4*alx6 - 2*alx4*aly7 + 2*alx5*aly6 + 2*aly4*alx7 - 2*alx6*aly5 - 2*alx5*aly7 + 2*aly5*alx7 + 2*alx6*aly7 - 2*alx7*aly6;
big_jac[15 ][ 3] = aly3*alx5 - alx3*aly5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 - alx4*aly6 + aly4*alx6 - alx5*aly7 + aly5*alx7 + alx5*aly8 + alx6*aly7 - aly5*alx8 - alx7*aly6 - alx6*aly8 + aly6*alx8;
big_jac[16 ][ 3] = alx2*aly5 - aly2*alx5 - alx2*aly6 + aly2*alx6 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 + alx2*aly8 - aly2*alx8 + alx4*aly6 - aly4*alx6 + alx4*aly7 - aly4*alx7 - alx4*aly8 + aly4*alx8 + alx5*myco - alx6*myco - alx7*myco + alx8*myco;
big_jac[17 ][ 3] = aly3*alx5 - alx3*aly5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 + alx3*aly7 - alx4*aly6 - aly3*alx7 + aly4*alx6 - alx3*aly8 - alx4*aly7 + aly3*alx8 + aly4*alx7 + alx4*aly8 - aly4*alx8;
big_jac[18 ][ 3] = alx2*aly3 - alx3*aly2 - alx2*aly4 + aly2*alx4 + alx3*aly4 - alx4*aly3 - alx2*aly7 + alx3*aly6 + aly2*alx7 - aly3*alx6 + alx2*aly8 - aly2*alx8 - alx4*aly6 + aly4*alx6 - alx3*aly8 + alx4*aly7 + aly3*alx8 - aly4*alx7 + alx6*aly7 - alx7*aly6 - alx6*aly8 + aly6*alx8 + alx7*aly8 - alx8*aly7;
big_jac[19 ][ 3] = aly2*alx5 - alx2*aly5 + alx2*aly6 - aly2*alx6 + alx2*aly7 - aly2*alx7 + alx4*aly5 - alx5*aly4 - alx2*aly8 + aly2*alx8 - alx4*aly6 + aly4*alx6 - alx4*aly7 - alx5*aly6 + aly4*alx7 + alx6*aly5 + alx4*aly8 - aly4*alx8 + alx5*aly8 - alx6*aly7 - aly5*alx8 + alx7*aly6 - alx7*aly8 + alx8*aly7;
big_jac[20 ][ 3] = alx3*aly5 - aly3*alx5 - alx3*aly6 - alx4*aly5 + aly3*alx6 + alx5*aly4 - alx3*aly7 + alx4*aly6 + aly3*alx7 - aly4*alx6 + alx3*aly8 + alx4*aly7 - aly3*alx8 - aly4*alx7 - alx4*aly8 + alx5*aly7 + aly4*alx8 - aly5*alx7 - alx5*aly8 - alx6*aly7 + aly5*alx8 + alx7*aly6 + alx6*aly8 - aly6*alx8;

break;
case 2:

big_jac[1 ][ 4] = alz1*aly3 - aly1*alz3 + aly1*alz5 - alz1*aly5 - aly3*alz5 + alz3*aly5 + alz1*myco - alz5*myco;
big_jac[2 ][ 4] = 2*aly1*alz3 - 2*alz1*aly3 - aly1*alz5 + alz1*aly5 - aly1*alz7 + alz1*aly7 + 2*aly3*alz5 - 2*alz3*aly5 + aly5*alz7 - alz5*aly7 - alz1*myco + alz5*myco;
big_jac[3 ][ 4] = aly1*alz3 - alz1*aly3 - 2*aly1*alz5 + 2*alz1*aly5 + aly1*alz7 - alz1*aly7 + 2*aly3*alz5 - 2*alz3*aly5 - aly3*alz7 + alz3*aly7 - 2*alz1*myco + alz3*myco + 2*alz5*myco - alz7*myco;
big_jac[4 ][ 4] = aly1*alz3 - alz1*aly3 - aly1*alz4 + alz1*aly4 - aly1*alz5 + alz1*aly5 + aly1*alz6 - alz1*aly6 + 2*aly3*alz5 - 2*alz3*aly5 - aly3*alz6 - aly4*alz5 + alz3*aly6 + aly5*alz4 - aly3*myco + aly5*myco + alz5*myco - alz6*myco - myco*myco;
big_jac[5 ][ 4] = alz1*aly3 - aly1*alz3 + aly1*alz7 - alz1*aly7 - aly3*alz5 + alz3*aly5 - aly5*alz7 + alz5*aly7;
big_jac[6 ][ 4] = 2*alz1*aly3 - 2*aly1*alz3 + 2*aly1*alz5 - 2*alz1*aly5 - 3*aly3*alz5 + 3*alz3*aly5 + aly3*alz7 - alz3*aly7 - aly5*alz7 + alz5*aly7 + 2*alz1*myco - alz3*myco - 2*alz5*myco + alz7*myco;
big_jac[7 ][ 4] = 2*alz1*aly3 - 2*aly1*alz3 + 2*aly1*alz4 - 2*alz1*aly4 + aly1*alz5 - alz1*aly5 - aly1*alz6 + alz1*aly6 + aly1*alz7 - alz1*aly7 - 3*aly3*alz5 + 3*alz3*aly5 - aly1*alz8 + alz1*aly8 + aly3*alz6 + 2*aly4*alz5 - alz3*aly6 - 2*aly5*alz4 + aly5*alz6 - aly6*alz5 - 2*aly5*alz7 + 2*alz5*aly7 + aly5*alz8 + aly6*alz7 - alz5*aly8 - aly7*alz6 + aly3*myco - aly7*myco;
big_jac[8 ][ 4] = aly1*alz5 - alz1*aly5 - aly1*alz7 + alz1*aly7 - aly3*alz5 + alz3*aly5 + aly3*alz7 - alz3*aly7 + alz1*myco - alz3*myco - alz5*myco + alz7*myco;
big_jac[9 ][ 4] = alz1*aly3 - aly1*alz3 + aly1*alz4 - alz1*aly4 + 2*aly1*alz5 - 2*alz1*aly5 - 2*aly1*alz6 + 2*alz1*aly6 - aly3*alz4 + aly4*alz3 - aly1*alz7 + alz1*aly7 - 3*aly3*alz5 + 3*alz3*aly5 + aly1*alz8 - alz1*aly8 + 2*aly3*alz6 + aly4*alz5 - 2*alz3*aly6 - aly5*alz4 + 2*aly3*alz7 - 2*alz3*aly7 - aly3*alz8 - aly4*alz7 + alz3*aly8 + alz4*aly7 - aly5*myco - 2*alz5*myco + aly7*myco + 2*alz6*myco + alz7*myco - alz8*myco + myco*myco;
big_jac[10 ][ 4] = alz3*aly5 - aly3*alz5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 - aly4*alz6 + alz4*aly6 + aly3*myco - aly4*myco - aly5*myco + aly6*myco;
big_jac[11 ][ 4] = aly1*alz3 - alz1*aly3 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly5*alz7 - alz5*aly7;
big_jac[12 ][ 4] = aly1*alz3 - alz1*aly3 - aly1*alz4 + alz1*aly4 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly1*alz8 - alz1*aly8 - aly4*alz5 + aly5*alz4 + aly5*alz7 - alz5*aly7 - aly5*alz8 + alz5*aly8;
big_jac[13 ][ 4] = alz1*aly5 - aly1*alz5 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly3*alz7 + alz3*aly7 - alz1*myco + alz3*myco + alz5*myco - alz7*myco;
big_jac[14 ][ 4] = 2*aly1*alz3 - 2*alz1*aly3 - 2*aly1*alz4 + 2*alz1*aly4 - 2*aly1*alz5 + 2*alz1*aly5 + 2*aly1*alz6 - 2*alz1*aly6 + 2*aly3*alz4 - 2*aly4*alz3 + 4*aly3*alz5 - 4*alz3*aly5 - 2*aly3*alz6 - 2*aly4*alz5 + 2*alz3*aly6 + 2*aly5*alz4 - 2*aly3*alz7 + 2*alz3*aly7 + 2*aly4*alz7 - 2*aly5*alz6 - 2*alz4*aly7 + 2*aly6*alz5 + 2*aly5*alz7 - 2*alz5*aly7 - 2*aly6*alz7 + 2*aly7*alz6;
big_jac[15 ][ 4] = aly3*alz5 - alz3*aly5 - aly3*alz6 - aly4*alz5 + alz3*aly6 + aly5*alz4 + aly4*alz6 - alz4*aly6 + aly5*alz7 - alz5*aly7 - aly5*alz8 - aly6*alz7 + alz5*aly8 + aly7*alz6 + aly6*alz8 - alz6*aly8 - aly3*myco + aly4*myco + aly7*myco - aly8*myco;
big_jac[16 ][ 4] = alz1*aly5 - aly1*alz5 + aly1*alz6 - alz1*aly6 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly1*alz8 + alz1*aly8 - aly3*alz6 + alz3*aly6 - aly3*alz7 + alz3*aly7 + aly3*alz8 - alz3*aly8 + alz5*myco - alz6*myco - alz7*myco + alz8*myco;
big_jac[17 ][ 4] = aly3*alz5 - alz3*aly5 - aly3*alz6 - aly4*alz5 + alz3*aly6 + aly5*alz4 - aly3*alz7 + aly4*alz6 + alz3*aly7 - alz4*aly6 + aly3*alz8 + aly4*alz7 - alz3*aly8 - alz4*aly7 - aly4*alz8 + alz4*aly8 + aly5*myco - aly6*myco - aly7*myco + aly8*myco;
big_jac[18 ][ 4] = alz1*aly3 - aly1*alz3 + aly1*alz4 - alz1*aly4 - aly3*alz4 + aly4*alz3 + aly1*alz7 - alz1*aly7 - aly3*alz5 + alz3*aly5 - aly1*alz8 + alz1*aly8 + aly4*alz5 - aly5*alz4 + aly3*alz8 - aly4*alz7 - alz3*aly8 + alz4*aly7 - aly5*alz7 + alz5*aly7 + aly5*alz8 - alz5*aly8 - aly7*alz8 + aly8*alz7;
big_jac[19 ][ 4] = aly1*alz5 - alz1*aly5 - aly1*alz6 + alz1*aly6 - aly1*alz7 + alz1*aly7 - aly3*alz5 + alz3*aly5 + aly1*alz8 - alz1*aly8 + aly3*alz6 - alz3*aly6 + aly3*alz7 - alz3*aly7 - aly3*alz8 + alz3*aly8 + aly5*alz6 - aly6*alz5 - aly5*alz8 + aly6*alz7 + alz5*aly8 - aly7*alz6 + aly7*alz8 - aly8*alz7;
big_jac[20 ][ 4] = alz3*aly5 - aly3*alz5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 + aly3*alz7 - aly4*alz6 - alz3*aly7 + alz4*aly6 - aly3*alz8 - aly4*alz7 + alz3*aly8 + alz4*aly7 + aly4*alz8 - aly5*alz7 - alz4*aly8 + alz5*aly7 + aly5*alz8 + aly6*alz7 - alz5*aly8 - aly7*alz6 - aly6*alz8 + alz6*aly8;
big_jac[1 ][ 5] = alx1*alz3 - alx3*alz1 - alx1*alz5 + alz1*alx5 + alx3*alz5 - alx5*alz3 + alz1*myco - alz3*myco;
big_jac[2 ][ 5] = 2*alx3*alz1 - 2*alx1*alz3 + alx1*alz5 - alz1*alx5 + alx1*alz7 - 2*alx3*alz5 - alz1*alx7 + 2*alx5*alz3 - alx5*alz7 + alx7*alz5 - 2*alz1*myco + 2*alz3*myco + alz5*myco - alz7*myco;
big_jac[3 ][ 5] = alx3*alz1 - alx1*alz3 + 2*alx1*alz5 - 2*alz1*alx5 - alx1*alz7 - 2*alx3*alz5 + alz1*alx7 + 2*alx5*alz3 + alx3*alz7 - alz3*alx7 - alz1*myco + alz3*myco;
big_jac[4 ][ 5] = alx3*alz1 - alx1*alz3 + alx1*alz4 - alz1*alx4 + alx1*alz5 - alz1*alx5 - alx1*alz6 + alz1*alx6 - 2*alx3*alz5 + 2*alx5*alz3 + alx3*alz6 + alx4*alz5 - alx5*alz4 - alz3*alx6 + alx3*myco - alx5*myco + alz3*myco - alz4*myco - myco*myco;
big_jac[5 ][ 5] = alx1*alz3 - alx3*alz1 - alx1*alz7 + alx3*alz5 + alz1*alx7 - alx5*alz3 + alx5*alz7 - alx7*alz5 + alz1*myco - alz3*myco - alz5*myco + alz7*myco;
big_jac[6 ][ 5] = 2*alx1*alz3 - 2*alx3*alz1 - 2*alx1*alz5 + 2*alz1*alx5 + 3*alx3*alz5 - 3*alx5*alz3 - alx3*alz7 + alz3*alx7 + alx5*alz7 - alx7*alz5 + 2*alz1*myco - 2*alz3*myco - alz5*myco + alz7*myco;
big_jac[7 ][ 5] = 2*alx1*alz3 - 2*alx3*alz1 - 2*alx1*alz4 + 2*alz1*alx4 - alx1*alz5 + alz1*alx5 + alx1*alz6 - alz1*alx6 - alx1*alz7 + 3*alx3*alz5 + alz1*alx7 - 3*alx5*alz3 + alx1*alz8 - alx3*alz6 - alz1*alx8 - 2*alx4*alz5 + 2*alx5*alz4 + alz3*alx6 - alx5*alz6 + alx6*alz5 + 2*alx5*alz7 - 2*alx7*alz5 - alx5*alz8 - alx6*alz7 + alx7*alz6 + alz5*alx8 - alx3*myco - 2*alz3*myco + 2*alz4*myco + alx7*myco + alz7*myco - alz8*myco + myco*myco;
big_jac[8 ][ 5] = alz1*alx5 - alx1*alz5 + alx1*alz7 + alx3*alz5 - alz1*alx7 - alx5*alz3 - alx3*alz7 + alz3*alx7;
big_jac[9 ][ 5] = alx1*alz3 - alx3*alz1 - alx1*alz4 + alz1*alx4 - 2*alx1*alz5 + 2*alz1*alx5 + 2*alx1*alz6 + alx3*alz4 - 2*alz1*alx6 - alx4*alz3 + alx1*alz7 + 3*alx3*alz5 - alz1*alx7 - 3*alx5*alz3 - alx1*alz8 - 2*alx3*alz6 + alz1*alx8 - alx4*alz5 + alx5*alz4 + 2*alz3*alx6 - 2*alx3*alz7 + 2*alz3*alx7 + alx3*alz8 + alx4*alz7 - alz3*alx8 - alz4*alx7 + alx5*myco - alx7*myco;
big_jac[10 ][ 5] = alx3*alz5 - alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 + alx4*alz6 - alx6*alz4 - alx3*myco + alx4*myco + alx5*myco - alx6*myco;
big_jac[11 ][ 5] = alx3*alz1 - alx1*alz3 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx5*alz7 + alx7*alz5 - alz1*myco + alz3*myco + alz5*myco - alz7*myco;
big_jac[12 ][ 5] = alx3*alz1 - alx1*alz3 + alx1*alz4 - alz1*alx4 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx1*alz8 + alz1*alx8 + alx4*alz5 - alx5*alz4 - alx5*alz7 + alx7*alz5 + alx5*alz8 - alz5*alx8 + alz3*myco - alz4*myco - alz7*myco + alz8*myco;
big_jac[13 ][ 5] = alx1*alz5 - alz1*alx5 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx3*alz7 - alz3*alx7;
big_jac[14 ][ 5] = 2*alx3*alz1 - 2*alx1*alz3 + 2*alx1*alz4 - 2*alz1*alx4 + 2*alx1*alz5 - 2*alz1*alx5 - 2*alx1*alz6 - 2*alx3*alz4 + 2*alz1*alx6 + 2*alx4*alz3 - 4*alx3*alz5 + 4*alx5*alz3 + 2*alx3*alz6 + 2*alx4*alz5 - 2*alx5*alz4 - 2*alz3*alx6 + 2*alx3*alz7 - 2*alz3*alx7 - 2*alx4*alz7 + 2*alx5*alz6 - 2*alx6*alz5 + 2*alz4*alx7 - 2*alx5*alz7 + 2*alx7*alz5 + 2*alx6*alz7 - 2*alx7*alz6;
big_jac[15 ][ 5] = alx5*alz3 - alx3*alz5 + alx3*alz6 + alx4*alz5 - alx5*alz4 - alz3*alx6 - alx4*alz6 + alx6*alz4 - alx5*alz7 + alx7*alz5 + alx5*alz8 + alx6*alz7 - alx7*alz6 - alz5*alx8 - alx6*alz8 + alx8*alz6 + alx3*myco - alx4*myco - alx7*myco + alx8*myco;
big_jac[16 ][ 5] = alx1*alz5 - alz1*alx5 - alx1*alz6 + alz1*alx6 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx1*alz8 + alx3*alz6 - alz1*alx8 - alz3*alx6 + alx3*alz7 - alz3*alx7 - alx3*alz8 + alz3*alx8;
big_jac[17 ][ 5] = alx5*alz3 - alx3*alz5 + alx3*alz6 + alx4*alz5 - alx5*alz4 - alz3*alx6 + alx3*alz7 - alx4*alz6 - alz3*alx7 + alx6*alz4 - alx3*alz8 - alx4*alz7 + alz3*alx8 + alz4*alx7 + alx4*alz8 - alz4*alx8 - alx5*myco + alx6*myco + alx7*myco - alx8*myco;
big_jac[18 ][ 5] = alx1*alz3 - alx3*alz1 - alx1*alz4 + alz1*alx4 + alx3*alz4 - alx4*alz3 - alx1*alz7 + alx3*alz5 + alz1*alx7 - alx5*alz3 + alx1*alz8 - alz1*alx8 - alx4*alz5 + alx5*alz4 - alx3*alz8 + alx4*alz7 + alz3*alx8 - alz4*alx7 + alx5*alz7 - alx7*alz5 - alx5*alz8 + alz5*alx8 + alx7*alz8 - alx8*alz7;
big_jac[19 ][ 5] = alz1*alx5 - alx1*alz5 + alx1*alz6 - alz1*alx6 + alx1*alz7 + alx3*alz5 - alz1*alx7 - alx5*alz3 - alx1*alz8 - alx3*alz6 + alz1*alx8 + alz3*alx6 - alx3*alz7 + alz3*alx7 + alx3*alz8 - alx5*alz6 - alz3*alx8 + alx6*alz5 + alx5*alz8 - alx6*alz7 + alx7*alz6 - alz5*alx8 - alx7*alz8 + alx8*alz7;
big_jac[20 ][ 5] = alx3*alz5 - alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 - alx3*alz7 + alx4*alz6 + alz3*alx7 - alx6*alz4 + alx3*alz8 + alx4*alz7 - alz3*alx8 - alz4*alx7 - alx4*alz8 + alx5*alz7 + alz4*alx8 - alx7*alz5 - alx5*alz8 - alx6*alz7 + alx7*alz6 + alz5*alx8 + alx6*alz8 - alx8*alz6;
big_jac[1 ][ 6] = aly1*alx3 - alx1*aly3 + alx1*aly5 - aly1*alx5 - alx3*aly5 + aly3*alx5 - alx1*myco - aly1*myco + aly3*myco + alx5*myco + myco*myco;
big_jac[2 ][ 6] = 2*alx1*aly3 - 2*aly1*alx3 - alx1*aly5 + aly1*alx5 - alx1*aly7 + aly1*alx7 + 2*alx3*aly5 - 2*aly3*alx5 + alx5*aly7 - aly5*alx7 + alx1*myco + 2*aly1*myco - 2*aly3*myco - alx5*myco - aly5*myco + aly7*myco - myco*myco;
big_jac[3 ][ 6] = alx1*aly3 - aly1*alx3 - 2*alx1*aly5 + 2*aly1*alx5 + alx1*aly7 - aly1*alx7 + 2*alx3*aly5 - 2*aly3*alx5 - alx3*aly7 + aly3*alx7 + 2*alx1*myco + aly1*myco - alx3*myco - aly3*myco - 2*alx5*myco + alx7*myco - myco*myco;
big_jac[4 ][ 6] = alx1*aly3 - aly1*alx3 - alx1*aly4 + aly1*alx4 - alx1*aly5 + aly1*alx5 + alx1*aly6 - aly1*alx6 + 2*alx3*aly5 - 2*aly3*alx5 - alx3*aly6 - alx4*aly5 + aly3*alx6 + alx5*aly4 - aly3*myco - alx5*myco + aly4*myco + alx6*myco;
big_jac[5 ][ 6] = aly1*alx3 - alx1*aly3 + alx1*aly7 - aly1*alx7 - alx3*aly5 + aly3*alx5 - alx5*aly7 + aly5*alx7 - aly1*myco + aly3*myco + aly5*myco - aly7*myco;
big_jac[6 ][ 6] = 2*aly1*alx3 - 2*alx1*aly3 + 2*alx1*aly5 - 2*aly1*alx5 - 3*alx3*aly5 + 3*aly3*alx5 + alx3*aly7 - aly3*alx7 - alx5*aly7 + aly5*alx7 - 2*alx1*myco - 2*aly1*myco + alx3*myco + 2*aly3*myco + 2*alx5*myco + aly5*myco - alx7*myco - aly7*myco + myco*myco;
big_jac[7 ][ 6] = 2*aly1*alx3 - 2*alx1*aly3 + 2*alx1*aly4 - 2*aly1*alx4 + alx1*aly5 - aly1*alx5 - alx1*aly6 + aly1*alx6 + alx1*aly7 - aly1*alx7 - 3*alx3*aly5 + 3*aly3*alx5 - alx1*aly8 + aly1*alx8 + alx3*aly6 + 2*alx4*aly5 - aly3*alx6 - 2*alx5*aly4 + alx5*aly6 - alx6*aly5 - 2*alx5*aly7 + 2*aly5*alx7 + alx5*aly8 + alx6*aly7 - aly5*alx8 - alx7*aly6 + 2*aly3*myco - 2*aly4*myco - aly7*myco + aly8*myco;
big_jac[8 ][ 6] = alx1*aly5 - aly1*alx5 - alx1*aly7 + aly1*alx7 - alx3*aly5 + aly3*alx5 + alx3*aly7 - aly3*alx7 - alx1*myco + alx3*myco + alx5*myco - alx7*myco;
big_jac[9 ][ 6] = aly1*alx3 - alx1*aly3 + alx1*aly4 - aly1*alx4 + 2*alx1*aly5 - 2*aly1*alx5 - 2*alx1*aly6 + 2*aly1*alx6 - alx3*aly4 + alx4*aly3 - alx1*aly7 + aly1*alx7 - 3*alx3*aly5 + 3*aly3*alx5 + alx1*aly8 - aly1*alx8 + 2*alx3*aly6 + alx4*aly5 - 2*aly3*alx6 - alx5*aly4 + 2*alx3*aly7 - 2*aly3*alx7 - alx3*aly8 - alx4*aly7 + aly3*alx8 + aly4*alx7 + 2*alx5*myco - 2*alx6*myco - alx7*myco + alx8*myco;
big_jac[10 ][ 6] = aly3*alx5 - alx3*aly5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 - alx4*aly6 + aly4*alx6;
big_jac[11 ][ 6] = alx1*aly3 - aly1*alx3 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx5*aly7 - aly5*alx7 + aly1*myco - aly3*myco - aly5*myco + aly7*myco;
big_jac[12 ][ 6] = alx1*aly3 - aly1*alx3 - alx1*aly4 + aly1*alx4 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx1*aly8 - aly1*alx8 - alx4*aly5 + alx5*aly4 + alx5*aly7 - aly5*alx7 - alx5*aly8 + aly5*alx8 - aly3*myco + aly4*myco + aly7*myco - aly8*myco;
big_jac[13 ][ 6] = aly1*alx5 - alx1*aly5 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx3*aly7 + aly3*alx7 + alx1*myco - alx3*myco - alx5*myco + alx7*myco;
big_jac[14 ][ 6] = 2*alx1*aly3 - 2*aly1*alx3 - 2*alx1*aly4 + 2*aly1*alx4 - 2*alx1*aly5 + 2*aly1*alx5 + 2*alx1*aly6 - 2*aly1*alx6 + 2*alx3*aly4 - 2*alx4*aly3 + 4*alx3*aly5 - 4*aly3*alx5 - 2*alx3*aly6 - 2*alx4*aly5 + 2*aly3*alx6 + 2*alx5*aly4 - 2*alx3*aly7 + 2*aly3*alx7 + 2*alx4*aly7 - 2*alx5*aly6 - 2*aly4*alx7 + 2*alx6*aly5 + 2*alx5*aly7 - 2*aly5*alx7 - 2*alx6*aly7 + 2*alx7*aly6;
big_jac[15 ][ 6] = alx3*aly5 - aly3*alx5 - alx3*aly6 - alx4*aly5 + aly3*alx6 + alx5*aly4 + alx4*aly6 - aly4*alx6 + alx5*aly7 - aly5*alx7 - alx5*aly8 - alx6*aly7 + aly5*alx8 + alx7*aly6 + alx6*aly8 - aly6*alx8;
big_jac[16 ][ 6] = aly1*alx5 - alx1*aly5 + alx1*aly6 - aly1*alx6 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx1*aly8 + aly1*alx8 - alx3*aly6 + aly3*alx6 - alx3*aly7 + aly3*alx7 + alx3*aly8 - aly3*alx8 - alx5*myco + alx6*myco + alx7*myco - alx8*myco;
big_jac[17 ][ 6] = alx3*aly5 - aly3*alx5 - alx3*aly6 - alx4*aly5 + aly3*alx6 + alx5*aly4 - alx3*aly7 + alx4*aly6 + aly3*alx7 - aly4*alx6 + alx3*aly8 + alx4*aly7 - aly3*alx8 - aly4*alx7 - alx4*aly8 + aly4*alx8;
big_jac[18 ][ 6] = aly1*alx3 - alx1*aly3 + alx1*aly4 - aly1*alx4 - alx3*aly4 + alx4*aly3 + alx1*aly7 - aly1*alx7 - alx3*aly5 + aly3*alx5 - alx1*aly8 + aly1*alx8 + alx4*aly5 - alx5*aly4 + alx3*aly8 - alx4*aly7 - aly3*alx8 + aly4*alx7 - alx5*aly7 + aly5*alx7 + alx5*aly8 - aly5*alx8 - alx7*aly8 + alx8*aly7;
big_jac[19 ][ 6] = alx1*aly5 - aly1*alx5 - alx1*aly6 + aly1*alx6 - alx1*aly7 + aly1*alx7 - alx3*aly5 + aly3*alx5 + alx1*aly8 - aly1*alx8 + alx3*aly6 - aly3*alx6 + alx3*aly7 - aly3*alx7 - alx3*aly8 + aly3*alx8 + alx5*aly6 - alx6*aly5 - alx5*aly8 + alx6*aly7 + aly5*alx8 - alx7*aly6 + alx7*aly8 - alx8*aly7;
big_jac[20 ][ 6] = aly3*alx5 - alx3*aly5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 + alx3*aly7 - alx4*aly6 - aly3*alx7 + aly4*alx6 - alx3*aly8 - alx4*aly7 + aly3*alx8 + aly4*alx7 + alx4*aly8 - alx5*aly7 - aly4*alx8 + aly5*alx7 + alx5*aly8 + alx6*aly7 - aly5*alx8 - alx7*aly6 - alx6*aly8 + aly6*alx8;

break;
case 3:

big_jac[1 ][ 7] = aly1*alz2 - aly2*alz1 - aly1*alz5 + alz1*aly5 + aly2*alz5 - alz2*aly5 + aly1*myco - aly5*myco;
big_jac[2 ][ 7] = 2*aly2*alz1 - 2*aly1*alz2 + aly1*alz5 - alz1*aly5 + aly1*alz6 - 2*aly2*alz5 - alz1*aly6 + 2*alz2*aly5 - aly5*alz6 + aly6*alz5 - aly1*myco + aly5*myco;
big_jac[3 ][ 7] = aly2*alz1 - aly1*alz2 + aly1*alz4 - alz1*aly4 + aly1*alz5 - alz1*aly5 - 2*aly2*alz5 + 2*alz2*aly5 - aly1*alz7 + alz1*aly7 + aly2*alz7 - alz2*aly7 + aly4*alz5 - aly5*alz4 - alz2*myco + aly5*myco + alz5*myco - aly7*myco - myco*myco;
big_jac[4 ][ 7] = aly2*alz1 - aly1*alz2 + 2*aly1*alz5 - 2*alz1*aly5 - aly1*alz6 - 2*aly2*alz5 + alz1*aly6 + 2*alz2*aly5 + aly2*alz6 - alz2*aly6 - 2*aly1*myco + aly2*myco + 2*aly5*myco - aly6*myco;
big_jac[5 ][ 7] = aly1*alz2 - aly2*alz1 - aly1*alz6 + aly2*alz5 + alz1*aly6 - alz2*aly5 + aly5*alz6 - aly6*alz5;
big_jac[6 ][ 7] = 2*aly1*alz2 - 2*aly2*alz1 - 2*aly1*alz4 + 2*alz1*aly4 - aly1*alz5 + alz1*aly5 - aly1*alz6 + 3*aly2*alz5 + alz1*aly6 - 3*alz2*aly5 + aly1*alz7 - alz1*aly7 + aly1*alz8 - aly2*alz7 - alz1*aly8 + alz2*aly7 - 2*aly4*alz5 + 2*aly5*alz4 + 2*aly5*alz6 - 2*aly6*alz5 - aly5*alz7 + alz5*aly7 - aly5*alz8 + aly6*alz7 + alz5*aly8 - aly7*alz6 + alz2*myco - alz6*myco;
big_jac[7 ][ 7] = 2*aly1*alz2 - 2*aly2*alz1 - 2*aly1*alz5 + 2*alz1*aly5 + 3*aly2*alz5 - 3*alz2*aly5 - aly2*alz6 + alz2*aly6 + aly5*alz6 - aly6*alz5 + 2*aly1*myco - aly2*myco - 2*aly5*myco + aly6*myco;
big_jac[8 ][ 7] = aly2*alz5 - alz2*aly5 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 + aly4*alz7 - alz4*aly7 + alz2*myco - alz4*myco - alz5*myco + alz7*myco;
big_jac[9 ][ 7] = aly1*alz2 - aly2*alz1 - aly1*alz4 + alz1*aly4 - 2*aly1*alz5 + aly2*alz4 + 2*alz1*aly5 - alz2*aly4 + aly1*alz6 + 3*aly2*alz5 - alz1*aly6 - 3*alz2*aly5 + 2*aly1*alz7 - 2*aly2*alz6 - 2*alz1*aly7 + 2*alz2*aly6 - aly1*alz8 - 2*aly2*alz7 + alz1*aly8 + 2*alz2*aly7 - aly4*alz5 + aly5*alz4 + aly2*alz8 - alz2*aly8 + aly4*alz6 - alz4*aly6 - 2*aly5*myco + aly6*myco - alz5*myco + 2*aly7*myco + alz6*myco - aly8*myco + myco*myco;
big_jac[10 ][ 7] = alz1*aly5 - aly1*alz5 + aly1*alz6 + aly2*alz5 - alz1*aly6 - alz2*aly5 - aly2*alz6 + alz2*aly6 + aly1*myco - aly2*myco - aly5*myco + aly6*myco;
big_jac[11 ][ 7] = aly2*alz1 - aly1*alz2 + aly1*alz4 - alz1*aly4 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly1*alz8 + alz1*aly8 + aly4*alz5 - aly5*alz4 - aly5*alz6 + aly6*alz5 + aly5*alz8 - alz5*aly8;
big_jac[12 ][ 7] = aly2*alz1 - aly1*alz2 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly5*alz6 + aly6*alz5;
big_jac[13 ][ 7] = alz2*aly5 - aly2*alz5 + aly2*alz7 - alz2*aly7 + aly4*alz5 - aly5*alz4 - aly4*alz7 - aly5*alz6 + alz4*aly7 + aly6*alz5 + aly5*alz8 - aly6*alz7 - alz5*aly8 + aly7*alz6 - aly7*alz8 + aly8*alz7 - alz2*myco + alz4*myco + alz6*myco - alz8*myco;
big_jac[14 ][ 7] = 2*aly2*alz1 - 2*aly1*alz2 + 2*aly1*alz4 - 2*alz1*aly4 + 2*aly1*alz5 - 2*aly2*alz4 - 2*alz1*aly5 + 2*alz2*aly4 - 4*aly2*alz5 + 4*alz2*aly5 - 2*aly1*alz7 + 2*aly2*alz6 + 2*alz1*aly7 - 2*alz2*aly6 + 2*aly2*alz7 - 2*alz2*aly7 + 2*aly4*alz5 - 2*aly5*alz4 - 2*aly4*alz6 + 2*alz4*aly6 - 2*aly5*alz6 + 2*aly6*alz5 + 2*aly5*alz7 - 2*alz5*aly7 - 2*aly6*alz7 + 2*aly7*alz6;
big_jac[15 ][ 7] = aly1*alz5 - alz1*aly5 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 + aly2*alz6 - alz2*aly6 - aly1*myco + aly2*myco + aly5*myco - aly6*myco;
big_jac[16 ][ 7] = alz2*aly5 - aly2*alz5 + aly2*alz6 - alz2*aly6 + aly2*alz7 - alz2*aly7 + aly4*alz5 - aly5*alz4 - aly2*alz8 + alz2*aly8 - aly4*alz6 + alz4*aly6 - aly4*alz7 + alz4*aly7 + aly4*alz8 - alz4*aly8 + alz5*myco - alz6*myco - alz7*myco + alz8*myco;
big_jac[17 ][ 7] = aly1*alz5 - alz1*aly5 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 - aly1*alz7 + aly2*alz6 + alz1*aly7 - alz2*aly6 + aly1*alz8 + aly2*alz7 - alz1*aly8 - alz2*aly7 - aly2*alz8 + alz2*aly8 + aly5*myco - aly6*myco - aly7*myco + aly8*myco;
big_jac[18 ][ 7] = aly1*alz2 - aly2*alz1 - aly1*alz4 + alz1*aly4 + aly2*alz4 - alz2*aly4 - aly1*alz6 + aly2*alz5 + alz1*aly6 - alz2*aly5 + aly1*alz8 - alz1*aly8 - aly4*alz5 + aly5*alz4 - aly2*alz8 + alz2*aly8 + aly4*alz6 - alz4*aly6 + aly5*alz6 - aly6*alz5 - aly5*alz8 + alz5*aly8 + aly6*alz8 - alz6*aly8;
big_jac[19 ][ 7] = aly2*alz5 - alz2*aly5 - aly2*alz6 + alz2*aly6 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 + aly2*alz8 - alz2*aly8 + aly4*alz6 - alz4*aly6 + aly4*alz7 + aly5*alz6 - alz4*aly7 - aly6*alz5 - aly4*alz8 + alz4*aly8 - aly5*alz8 + aly6*alz7 + alz5*aly8 - aly7*alz6 + aly7*alz8 - aly8*alz7;
big_jac[20 ][ 7] = alz1*aly5 - aly1*alz5 + aly1*alz6 + aly2*alz5 - alz1*aly6 - alz2*aly5 + aly1*alz7 - aly2*alz6 - alz1*aly7 + alz2*aly6 - aly1*alz8 - aly2*alz7 + alz1*aly8 + alz2*aly7 + aly2*alz8 - alz2*aly8 - aly5*alz7 + alz5*aly7 + aly5*alz8 + aly6*alz7 - alz5*aly8 - aly7*alz6 - aly6*alz8 + alz6*aly8;
big_jac[1 ][ 8] = alx2*alz1 - alx1*alz2 + alx1*alz5 - alz1*alx5 - alx2*alz5 + alz2*alx5 - alx1*myco - alz1*myco + alz2*myco + alx5*myco + myco*myco;
big_jac[2 ][ 8] = 2*alx1*alz2 - 2*alx2*alz1 - alx1*alz5 + alz1*alx5 - alx1*alz6 + 2*alx2*alz5 + alz1*alx6 - 2*alz2*alx5 + alx5*alz6 - alx6*alz5 + alx1*myco + 2*alz1*myco - 2*alz2*myco - alx5*myco - alz5*myco + alz6*myco - myco*myco;
big_jac[3 ][ 8] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alz1*alx4 - alx1*alz5 + alz1*alx5 + 2*alx2*alz5 - 2*alz2*alx5 + alx1*alz7 - alz1*alx7 - alx2*alz7 - alx4*alz5 + alz2*alx7 + alx5*alz4 - alz2*myco - alx5*myco + alz4*myco + alx7*myco;
big_jac[4 ][ 8] = alx1*alz2 - alx2*alz1 - 2*alx1*alz5 + 2*alz1*alx5 + alx1*alz6 + 2*alx2*alz5 - alz1*alx6 - 2*alz2*alx5 - alx2*alz6 + alz2*alx6 + 2*alx1*myco - alx2*myco + alz1*myco - alz2*myco - 2*alx5*myco + alx6*myco - myco*myco;
big_jac[5 ][ 8] = alx2*alz1 - alx1*alz2 + alx1*alz6 - alx2*alz5 - alz1*alx6 + alz2*alx5 - alx5*alz6 + alx6*alz5 - alz1*myco + alz2*myco + alz5*myco - alz6*myco;
big_jac[6 ][ 8] = 2*alx2*alz1 - 2*alx1*alz2 + 2*alx1*alz4 - 2*alz1*alx4 + alx1*alz5 - alz1*alx5 + alx1*alz6 - 3*alx2*alz5 - alz1*alx6 + 3*alz2*alx5 - alx1*alz7 + alz1*alx7 - alx1*alz8 + alx2*alz7 + alz1*alx8 + 2*alx4*alz5 - alz2*alx7 - 2*alx5*alz4 - 2*alx5*alz6 + 2*alx6*alz5 + alx5*alz7 - alx7*alz5 + alx5*alz8 - alx6*alz7 + alx7*alz6 - alz5*alx8 + 2*alz2*myco - 2*alz4*myco - alz6*myco + alz8*myco;
big_jac[7 ][ 8] = 2*alx2*alz1 - 2*alx1*alz2 + 2*alx1*alz5 - 2*alz1*alx5 - 3*alx2*alz5 + 3*alz2*alx5 + alx2*alz6 - alz2*alx6 - alx5*alz6 + alx6*alz5 - 2*alx1*myco + alx2*myco - 2*alz1*myco + 2*alz2*myco + 2*alx5*myco - alx6*myco + alz5*myco - alz6*myco + myco*myco;
big_jac[8 ][ 8] = alz2*alx5 - alx2*alz5 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 - alx4*alz7 + alz4*alx7;
big_jac[9 ][ 8] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alz1*alx4 + 2*alx1*alz5 - alx2*alz4 - 2*alz1*alx5 + alx4*alz2 - alx1*alz6 - 3*alx2*alz5 + alz1*alx6 + 3*alz2*alx5 - 2*alx1*alz7 + 2*alx2*alz6 + 2*alz1*alx7 - 2*alz2*alx6 + alx1*alz8 + 2*alx2*alz7 - alz1*alx8 + alx4*alz5 - 2*alz2*alx7 - alx5*alz4 - alx2*alz8 - alx4*alz6 + alz2*alx8 + alx6*alz4 + 2*alx5*myco - alx6*myco - 2*alx7*myco + alx8*myco;
big_jac[10 ][ 8] = alx1*alz5 - alz1*alx5 - alx1*alz6 - alx2*alz5 + alz1*alx6 + alz2*alx5 + alx2*alz6 - alz2*alx6 - alx1*myco + alx2*myco + alx5*myco - alx6*myco;
big_jac[11 ][ 8] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alz1*alx4 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx1*alz8 - alz1*alx8 - alx4*alz5 + alx5*alz4 + alx5*alz6 - alx6*alz5 - alx5*alz8 + alz5*alx8 - alz2*myco + alz4*myco + alz6*myco - alz8*myco;
big_jac[12 ][ 8] = alx1*alz2 - alx2*alz1 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx5*alz6 - alx6*alz5 + alz1*myco - alz2*myco - alz5*myco + alz6*myco;
big_jac[13 ][ 8] = alx2*alz5 - alz2*alx5 - alx2*alz7 - alx4*alz5 + alz2*alx7 + alx5*alz4 + alx4*alz7 + alx5*alz6 - alx6*alz5 - alz4*alx7 - alx5*alz8 + alx6*alz7 - alx7*alz6 + alz5*alx8 + alx7*alz8 - alx8*alz7;
big_jac[14 ][ 8] = 2*alx1*alz2 - 2*alx2*alz1 - 2*alx1*alz4 + 2*alz1*alx4 - 2*alx1*alz5 + 2*alx2*alz4 + 2*alz1*alx5 - 2*alx4*alz2 + 4*alx2*alz5 - 4*alz2*alx5 + 2*alx1*alz7 - 2*alx2*alz6 - 2*alz1*alx7 + 2*alz2*alx6 - 2*alx2*alz7 - 2*alx4*alz5 + 2*alz2*alx7 + 2*alx5*alz4 + 2*alx4*alz6 - 2*alx6*alz4 + 2*alx5*alz6 - 2*alx6*alz5 - 2*alx5*alz7 + 2*alx7*alz5 + 2*alx6*alz7 - 2*alx7*alz6;
big_jac[15 ][ 8] = alz1*alx5 - alx1*alz5 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 - alx2*alz6 + alz2*alx6 + alx1*myco - alx2*myco - alx5*myco + alx6*myco;
big_jac[16 ][ 8] = alx2*alz5 - alz2*alx5 - alx2*alz6 + alz2*alx6 - alx2*alz7 - alx4*alz5 + alz2*alx7 + alx5*alz4 + alx2*alz8 + alx4*alz6 - alz2*alx8 - alx6*alz4 + alx4*alz7 - alz4*alx7 - alx4*alz8 + alz4*alx8;
big_jac[17 ][ 8] = alz1*alx5 - alx1*alz5 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 + alx1*alz7 - alx2*alz6 - alz1*alx7 + alz2*alx6 - alx1*alz8 - alx2*alz7 + alz1*alx8 + alz2*alx7 + alx2*alz8 - alz2*alx8 - alx5*myco + alx6*myco + alx7*myco - alx8*myco;
big_jac[18 ][ 8] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alz1*alx4 - alx2*alz4 + alx4*alz2 + alx1*alz6 - alx2*alz5 - alz1*alx6 + alz2*alx5 - alx1*alz8 + alz1*alx8 + alx4*alz5 - alx5*alz4 + alx2*alz8 - alx4*alz6 - alz2*alx8 + alx6*alz4 - alx5*alz6 + alx6*alz5 + alx5*alz8 - alz5*alx8 - alx6*alz8 + alx8*alz6;
big_jac[19 ][ 8] = alz2*alx5 - alx2*alz5 + alx2*alz6 - alz2*alx6 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 - alx2*alz8 - alx4*alz6 + alz2*alx8 + alx6*alz4 - alx4*alz7 - alx5*alz6 + alx6*alz5 + alz4*alx7 + alx4*alz8 - alz4*alx8 + alx5*alz8 - alx6*alz7 + alx7*alz6 - alz5*alx8 - alx7*alz8 + alx8*alz7;
big_jac[20 ][ 8] = alx1*alz5 - alz1*alx5 - alx1*alz6 - alx2*alz5 + alz1*alx6 + alz2*alx5 - alx1*alz7 + alx2*alz6 + alz1*alx7 - alz2*alx6 + alx1*alz8 + alx2*alz7 - alz1*alx8 - alz2*alx7 - alx2*alz8 + alz2*alx8 + alx5*alz7 - alx7*alz5 - alx5*alz8 - alx6*alz7 + alx7*alz6 + alz5*alx8 + alx6*alz8 - alx8*alz6;
big_jac[1 ][ 9] = alx1*aly2 - alx2*aly1 - alx1*aly5 + aly1*alx5 + alx2*aly5 - aly2*alx5 + aly1*myco - aly2*myco;
big_jac[2 ][ 9] = 2*alx2*aly1 - 2*alx1*aly2 + alx1*aly5 - aly1*alx5 + alx1*aly6 - 2*alx2*aly5 - aly1*alx6 + 2*aly2*alx5 - alx5*aly6 + alx6*aly5 - 2*aly1*myco + 2*aly2*myco + aly5*myco - aly6*myco;
big_jac[3 ][ 9] = alx2*aly1 - alx1*aly2 + alx1*aly4 - aly1*alx4 + alx1*aly5 - aly1*alx5 - 2*alx2*aly5 + 2*aly2*alx5 - alx1*aly7 + aly1*alx7 + alx2*aly7 - aly2*alx7 + alx4*aly5 - alx5*aly4 + alx2*myco + aly2*myco - alx5*myco - aly4*myco - myco*myco;
big_jac[4 ][ 9] = alx2*aly1 - alx1*aly2 + 2*alx1*aly5 - 2*aly1*alx5 - alx1*aly6 - 2*alx2*aly5 + aly1*alx6 + 2*aly2*alx5 + alx2*aly6 - aly2*alx6 - aly1*myco + aly2*myco;
big_jac[5 ][ 9] = alx1*aly2 - alx2*aly1 - alx1*aly6 + alx2*aly5 + aly1*alx6 - aly2*alx5 + alx5*aly6 - alx6*aly5 + aly1*myco - aly2*myco - aly5*myco + aly6*myco;
big_jac[6 ][ 9] = 2*alx1*aly2 - 2*alx2*aly1 - 2*alx1*aly4 + 2*aly1*alx4 - alx1*aly5 + aly1*alx5 - alx1*aly6 + 3*alx2*aly5 + aly1*alx6 - 3*aly2*alx5 + alx1*aly7 - aly1*alx7 + alx1*aly8 - alx2*aly7 - aly1*alx8 + aly2*alx7 - 2*alx4*aly5 + 2*alx5*aly4 + 2*alx5*aly6 - 2*alx6*aly5 - alx5*aly7 + aly5*alx7 - alx5*aly8 + alx6*aly7 + aly5*alx8 - alx7*aly6 - alx2*myco - 2*aly2*myco + 2*aly4*myco + alx6*myco + aly6*myco - aly8*myco + myco*myco;
big_jac[7 ][ 9] = 2*alx1*aly2 - 2*alx2*aly1 - 2*alx1*aly5 + 2*aly1*alx5 + 3*alx2*aly5 - 3*aly2*alx5 - alx2*aly6 + aly2*alx6 + alx5*aly6 - alx6*aly5 + 2*aly1*myco - 2*aly2*myco - aly5*myco + aly6*myco;
big_jac[8 ][ 9] = alx2*aly5 - aly2*alx5 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 + alx4*aly7 - aly4*alx7 - alx2*myco + alx4*myco + alx5*myco - alx7*myco;
big_jac[9 ][ 9] = alx1*aly2 - alx2*aly1 - alx1*aly4 + aly1*alx4 - 2*alx1*aly5 + alx2*aly4 + 2*aly1*alx5 - aly2*alx4 + alx1*aly6 + 3*alx2*aly5 - aly1*alx6 - 3*aly2*alx5 + 2*alx1*aly7 - 2*alx2*aly6 - 2*aly1*alx7 + 2*aly2*alx6 - alx1*aly8 - 2*alx2*aly7 + aly1*alx8 + 2*aly2*alx7 - alx4*aly5 + alx5*aly4 + alx2*aly8 - aly2*alx8 + alx4*aly6 - aly4*alx6 + alx5*myco - alx6*myco;
big_jac[10 ][ 9] = aly1*alx5 - alx1*aly5 + alx1*aly6 + alx2*aly5 - aly1*alx6 - aly2*alx5 - alx2*aly6 + aly2*alx6;
big_jac[11 ][ 9] = alx2*aly1 - alx1*aly2 + alx1*aly4 - aly1*alx4 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx1*aly8 + aly1*alx8 + alx4*aly5 - alx5*aly4 - alx5*aly6 + alx6*aly5 + alx5*aly8 - aly5*alx8 + aly2*myco - aly4*myco - aly6*myco + aly8*myco;
big_jac[12 ][ 9] = alx2*aly1 - alx1*aly2 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx5*aly6 + alx6*aly5 - aly1*myco + aly2*myco + aly5*myco - aly6*myco;
big_jac[13 ][ 9] = aly2*alx5 - alx2*aly5 + alx2*aly7 - aly2*alx7 + alx4*aly5 - alx5*aly4 - alx4*aly7 - alx5*aly6 + aly4*alx7 + alx6*aly5 + alx5*aly8 - alx6*aly7 - aly5*alx8 + alx7*aly6 - alx7*aly8 + alx8*aly7 + alx2*myco - alx4*myco - alx6*myco + alx8*myco;
big_jac[14 ][ 9] = 2*alx2*aly1 - 2*alx1*aly2 + 2*alx1*aly4 - 2*aly1*alx4 + 2*alx1*aly5 - 2*alx2*aly4 - 2*aly1*alx5 + 2*aly2*alx4 - 4*alx2*aly5 + 4*aly2*alx5 - 2*alx1*aly7 + 2*alx2*aly6 + 2*aly1*alx7 - 2*aly2*alx6 + 2*alx2*aly7 - 2*aly2*alx7 + 2*alx4*aly5 - 2*alx5*aly4 - 2*alx4*aly6 + 2*aly4*alx6 - 2*alx5*aly6 + 2*alx6*aly5 + 2*alx5*aly7 - 2*aly5*alx7 - 2*alx6*aly7 + 2*alx7*aly6;
big_jac[15 ][ 9] = alx1*aly5 - aly1*alx5 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 + alx2*aly6 - aly2*alx6;
big_jac[16 ][ 9] = aly2*alx5 - alx2*aly5 + alx2*aly6 - aly2*alx6 + alx2*aly7 - aly2*alx7 + alx4*aly5 - alx5*aly4 - alx2*aly8 + aly2*alx8 - alx4*aly6 + aly4*alx6 - alx4*aly7 + aly4*alx7 + alx4*aly8 - aly4*alx8 - alx5*myco + alx6*myco + alx7*myco - alx8*myco;
big_jac[17 ][ 9] = alx1*aly5 - aly1*alx5 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 - alx1*aly7 + alx2*aly6 + aly1*alx7 - aly2*alx6 + alx1*aly8 + alx2*aly7 - aly1*alx8 - aly2*alx7 - alx2*aly8 + aly2*alx8;
big_jac[18 ][ 9] = alx1*aly2 - alx2*aly1 - alx1*aly4 + aly1*alx4 + alx2*aly4 - aly2*alx4 - alx1*aly6 + alx2*aly5 + aly1*alx6 - aly2*alx5 + alx1*aly8 - aly1*alx8 - alx4*aly5 + alx5*aly4 - alx2*aly8 + aly2*alx8 + alx4*aly6 - aly4*alx6 + alx5*aly6 - alx6*aly5 - alx5*aly8 + aly5*alx8 + alx6*aly8 - aly6*alx8;
big_jac[19 ][ 9] = alx2*aly5 - aly2*alx5 - alx2*aly6 + aly2*alx6 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 + alx2*aly8 - aly2*alx8 + alx4*aly6 - aly4*alx6 + alx4*aly7 + alx5*aly6 - aly4*alx7 - alx6*aly5 - alx4*aly8 + aly4*alx8 - alx5*aly8 + alx6*aly7 + aly5*alx8 - alx7*aly6 + alx7*aly8 - alx8*aly7;
big_jac[20 ][ 9] = aly1*alx5 - alx1*aly5 + alx1*aly6 + alx2*aly5 - aly1*alx6 - aly2*alx5 + alx1*aly7 - alx2*aly6 - aly1*alx7 + aly2*alx6 - alx1*aly8 - alx2*aly7 + aly1*alx8 + aly2*alx7 + alx2*aly8 - aly2*alx8 - alx5*aly7 + aly5*alx7 + alx5*aly8 + alx6*aly7 - aly5*alx8 - alx7*aly6 - alx6*aly8 + aly6*alx8;

break;
case 4:

big_jac[1 ][ 10] = 0;
big_jac[2 ][ 10] = 0;
big_jac[3 ][ 10] = alz1*aly3 - aly1*alz3 + aly1*alz5 - alz1*aly5 - aly3*alz5 + alz3*aly5 + alz1*myco - alz5*myco;
big_jac[4 ][ 10] = aly1*alz2 - aly2*alz1 - aly1*alz5 + alz1*aly5 + aly2*alz5 - alz2*aly5 + aly1*myco - aly5*myco;
big_jac[5 ][ 10] = 0;
big_jac[6 ][ 10] = 2*aly1*alz3 - 2*alz1*aly3 - aly1*alz5 + alz1*aly5 - aly1*alz7 + alz1*aly7 + 2*aly3*alz5 - 2*alz3*aly5 + aly5*alz7 - alz5*aly7 - alz1*myco + alz5*myco;
big_jac[7 ][ 10] = 2*aly2*alz1 - 2*aly1*alz2 + aly1*alz5 - alz1*aly5 + aly1*alz6 - 2*aly2*alz5 - alz1*aly6 + 2*alz2*aly5 - aly5*alz6 + aly6*alz5 - aly1*myco + aly5*myco;
big_jac[8 ][ 10] = alz1*aly5 - aly1*alz5 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly3*alz7 + alz3*aly7 - alz1*myco + alz3*myco + alz5*myco - alz7*myco;
big_jac[9 ][ 10] = aly2*alz1 - aly1*alz2 + aly1*alz3 - alz1*aly3 - aly2*alz3 + aly3*alz2 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 + aly5*myco + alz5*myco - aly7*myco - alz6*myco - myco*myco;
big_jac[10 ][ 10] = aly1*alz5 - alz1*aly5 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 + aly2*alz6 - alz2*aly6 - aly1*myco + aly2*myco + aly5*myco - aly6*myco;
big_jac[11 ][ 10] = alz1*aly3 - aly1*alz3 + aly1*alz7 - alz1*aly7 - aly3*alz5 + alz3*aly5 - aly5*alz7 + alz5*aly7;
big_jac[12 ][ 10] = aly1*alz2 - aly2*alz1 - aly1*alz6 + aly2*alz5 + alz1*aly6 - alz2*aly5 + aly5*alz6 - aly6*alz5;
big_jac[13 ][ 10] = aly1*alz5 - alz1*aly5 - aly1*alz7 + alz1*aly7 - aly3*alz5 + alz3*aly5 + aly3*alz7 - alz3*aly7 + alz1*myco - alz3*myco - alz5*myco + alz7*myco;
big_jac[14 ][ 10] = 2*aly1*alz2 - 2*aly2*alz1 - 2*aly1*alz3 + 2*alz1*aly3 + 2*aly2*alz3 - 2*aly3*alz2 - 2*aly1*alz6 + 2*aly2*alz5 + 2*alz1*aly6 - 2*alz2*aly5 + 2*aly1*alz7 - 2*alz1*aly7 - 2*aly3*alz5 + 2*alz3*aly5 - 2*aly2*alz7 + 2*aly3*alz6 + 2*alz2*aly7 - 2*alz3*aly6 + 2*aly5*alz6 - 2*aly6*alz5 - 2*aly5*alz7 + 2*alz5*aly7 + 2*aly6*alz7 - 2*aly7*alz6;
big_jac[15 ][ 10] = alz1*aly5 - aly1*alz5 + aly1*alz6 + aly2*alz5 - alz1*aly6 - alz2*aly5 - aly2*alz6 + alz2*aly6 + aly1*myco - aly2*myco - aly5*myco + aly6*myco;
big_jac[16 ][ 10] = aly1*alz5 - alz1*aly5 - aly1*alz6 + alz1*aly6 - aly1*alz7 + alz1*aly7 - aly3*alz5 + alz3*aly5 + aly1*alz8 - alz1*aly8 + aly3*alz6 - alz3*aly6 + aly3*alz7 - alz3*aly7 - aly3*alz8 + alz3*aly8 - alz5*myco + alz6*myco + alz7*myco - alz8*myco;
big_jac[17 ][ 10] = alz1*aly5 - aly1*alz5 + aly1*alz6 + aly2*alz5 - alz1*aly6 - alz2*aly5 + aly1*alz7 - aly2*alz6 - alz1*aly7 + alz2*aly6 - aly1*alz8 - aly2*alz7 + alz1*aly8 + alz2*aly7 + aly2*alz8 - alz2*aly8 - aly5*myco + aly6*myco + aly7*myco - aly8*myco;
big_jac[18 ][ 10] = aly2*alz1 - aly1*alz2 + aly1*alz3 - alz1*aly3 - aly2*alz3 + aly3*alz2 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 - aly5*alz6 + aly6*alz5 + aly5*alz7 - alz5*aly7 - aly6*alz7 + aly7*alz6;
big_jac[19 ][ 10] = alz1*aly5 - aly1*alz5 + aly1*alz6 - alz1*aly6 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly1*alz8 + alz1*aly8 - aly3*alz6 + alz3*aly6 - aly3*alz7 + alz3*aly7 + aly3*alz8 - alz3*aly8 - aly5*alz6 + aly6*alz5 + aly5*alz8 - aly6*alz7 - alz5*aly8 + aly7*alz6 - aly7*alz8 + aly8*alz7;
big_jac[20 ][ 10] = aly1*alz5 - alz1*aly5 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 - aly1*alz7 + aly2*alz6 + alz1*aly7 - alz2*aly6 + aly1*alz8 + aly2*alz7 - alz1*aly8 - alz2*aly7 - aly2*alz8 + alz2*aly8 + aly5*alz7 - alz5*aly7 - aly5*alz8 - aly6*alz7 + alz5*aly8 + aly7*alz6 + aly6*alz8 - alz6*aly8;
big_jac[1 ][ 11] = 0;
big_jac[2 ][ 11] = 0;
big_jac[3 ][ 11] = alx1*alz3 - alx3*alz1 - alx1*alz5 + alz1*alx5 + alx3*alz5 - alx5*alz3 + alz1*myco - alz3*myco;
big_jac[4 ][ 11] = alx2*alz1 - alx1*alz2 + alx1*alz5 - alz1*alx5 - alx2*alz5 + alz2*alx5 - alx1*myco - alz1*myco + alz2*myco + alx5*myco + myco*myco;
big_jac[5 ][ 11] = 0;
big_jac[6 ][ 11] = 2*alx3*alz1 - 2*alx1*alz3 + alx1*alz5 - alz1*alx5 + alx1*alz7 - 2*alx3*alz5 - alz1*alx7 + 2*alx5*alz3 - alx5*alz7 + alx7*alz5 - 2*alz1*myco + 2*alz3*myco + alz5*myco - alz7*myco;
big_jac[7 ][ 11] = 2*alx1*alz2 - 2*alx2*alz1 - alx1*alz5 + alz1*alx5 - alx1*alz6 + 2*alx2*alz5 + alz1*alx6 - 2*alz2*alx5 + alx5*alz6 - alx6*alz5 + alx1*myco + 2*alz1*myco - 2*alz2*myco - alx5*myco - alz5*myco + alz6*myco - myco*myco;
big_jac[8 ][ 11] = alx1*alz5 - alz1*alx5 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx3*alz7 - alz3*alx7;
big_jac[9 ][ 11] = alx1*alz2 - alx2*alz1 - alx1*alz3 + alx3*alz1 + alx2*alz3 - alx3*alz2 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 - alx5*myco + alx7*myco;
big_jac[10 ][ 11] = alz1*alx5 - alx1*alz5 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 - alx2*alz6 + alz2*alx6 + alx1*myco - alx2*myco - alx5*myco + alx6*myco;
big_jac[11 ][ 11] = alx1*alz3 - alx3*alz1 - alx1*alz7 + alx3*alz5 + alz1*alx7 - alx5*alz3 + alx5*alz7 - alx7*alz5 + alz1*myco - alz3*myco - alz5*myco + alz7*myco;
big_jac[12 ][ 11] = alx2*alz1 - alx1*alz2 + alx1*alz6 - alx2*alz5 - alz1*alx6 + alz2*alx5 - alx5*alz6 + alx6*alz5 - alz1*myco + alz2*myco + alz5*myco - alz6*myco;
big_jac[13 ][ 11] = alz1*alx5 - alx1*alz5 + alx1*alz7 + alx3*alz5 - alz1*alx7 - alx5*alz3 - alx3*alz7 + alz3*alx7;
big_jac[14 ][ 11] = 2*alx2*alz1 - 2*alx1*alz2 + 2*alx1*alz3 - 2*alx3*alz1 - 2*alx2*alz3 + 2*alx3*alz2 + 2*alx1*alz6 - 2*alx2*alz5 - 2*alz1*alx6 + 2*alz2*alx5 - 2*alx1*alz7 + 2*alx3*alz5 + 2*alz1*alx7 - 2*alx5*alz3 + 2*alx2*alz7 - 2*alx3*alz6 - 2*alz2*alx7 + 2*alz3*alx6 - 2*alx5*alz6 + 2*alx6*alz5 + 2*alx5*alz7 - 2*alx7*alz5 - 2*alx6*alz7 + 2*alx7*alz6;
big_jac[15 ][ 11] = alx1*alz5 - alz1*alx5 - alx1*alz6 - alx2*alz5 + alz1*alx6 + alz2*alx5 + alx2*alz6 - alz2*alx6 - alx1*myco + alx2*myco + alx5*myco - alx6*myco;
big_jac[16 ][ 11] = alz1*alx5 - alx1*alz5 + alx1*alz6 - alz1*alx6 + alx1*alz7 + alx3*alz5 - alz1*alx7 - alx5*alz3 - alx1*alz8 - alx3*alz6 + alz1*alx8 + alz3*alx6 - alx3*alz7 + alz3*alx7 + alx3*alz8 - alz3*alx8;
big_jac[17 ][ 11] = alx1*alz5 - alz1*alx5 - alx1*alz6 - alx2*alz5 + alz1*alx6 + alz2*alx5 - alx1*alz7 + alx2*alz6 + alz1*alx7 - alz2*alx6 + alx1*alz8 + alx2*alz7 - alz1*alx8 - alz2*alx7 - alx2*alz8 + alz2*alx8 + alx5*myco - alx6*myco - alx7*myco + alx8*myco;
big_jac[18 ][ 11] = alx1*alz2 - alx2*alz1 - alx1*alz3 + alx3*alz1 + alx2*alz3 - alx3*alz2 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 + alx5*alz6 - alx6*alz5 - alx5*alz7 + alx7*alz5 + alx6*alz7 - alx7*alz6;
big_jac[19 ][ 11] = alx1*alz5 - alz1*alx5 - alx1*alz6 + alz1*alx6 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx1*alz8 + alx3*alz6 - alz1*alx8 - alz3*alx6 + alx3*alz7 - alz3*alx7 - alx3*alz8 + alx5*alz6 + alz3*alx8 - alx6*alz5 - alx5*alz8 + alx6*alz7 - alx7*alz6 + alz5*alx8 + alx7*alz8 - alx8*alz7;
big_jac[20 ][ 11] = alz1*alx5 - alx1*alz5 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 + alx1*alz7 - alx2*alz6 - alz1*alx7 + alz2*alx6 - alx1*alz8 - alx2*alz7 + alz1*alx8 + alz2*alx7 + alx2*alz8 - alz2*alx8 - alx5*alz7 + alx7*alz5 + alx5*alz8 + alx6*alz7 - alx7*alz6 - alz5*alx8 - alx6*alz8 + alx8*alz6;
big_jac[1 ][ 12] = 0;
big_jac[2 ][ 12] = 0;
big_jac[3 ][ 12] = aly1*alx3 - alx1*aly3 + alx1*aly5 - aly1*alx5 - alx3*aly5 + aly3*alx5 - alx1*myco - aly1*myco + aly3*myco + alx5*myco + myco*myco;
big_jac[4 ][ 12] = alx1*aly2 - alx2*aly1 - alx1*aly5 + aly1*alx5 + alx2*aly5 - aly2*alx5 + aly1*myco - aly2*myco;
big_jac[5 ][ 12] = 0;
big_jac[6 ][ 12] = 2*alx1*aly3 - 2*aly1*alx3 - alx1*aly5 + aly1*alx5 - alx1*aly7 + aly1*alx7 + 2*alx3*aly5 - 2*aly3*alx5 + alx5*aly7 - aly5*alx7 + alx1*myco + 2*aly1*myco - 2*aly3*myco - alx5*myco - aly5*myco + aly7*myco - myco*myco;
big_jac[7 ][ 12] = 2*alx2*aly1 - 2*alx1*aly2 + alx1*aly5 - aly1*alx5 + alx1*aly6 - 2*alx2*aly5 - aly1*alx6 + 2*aly2*alx5 - alx5*aly6 + alx6*aly5 - 2*aly1*myco + 2*aly2*myco + aly5*myco - aly6*myco;
big_jac[8 ][ 12] = aly1*alx5 - alx1*aly5 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx3*aly7 + aly3*alx7 + alx1*myco - alx3*myco - alx5*myco + alx7*myco;
big_jac[9 ][ 12] = alx2*aly1 - alx1*aly2 + alx1*aly3 - aly1*alx3 - alx2*aly3 + alx3*aly2 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 - alx5*myco + alx6*myco;
big_jac[10 ][ 12] = alx1*aly5 - aly1*alx5 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 + alx2*aly6 - aly2*alx6;
big_jac[11 ][ 12] = aly1*alx3 - alx1*aly3 + alx1*aly7 - aly1*alx7 - alx3*aly5 + aly3*alx5 - alx5*aly7 + aly5*alx7 - aly1*myco + aly3*myco + aly5*myco - aly7*myco;
big_jac[12 ][ 12] = alx1*aly2 - alx2*aly1 - alx1*aly6 + alx2*aly5 + aly1*alx6 - aly2*alx5 + alx5*aly6 - alx6*aly5 + aly1*myco - aly2*myco - aly5*myco + aly6*myco;
big_jac[13 ][ 12] = alx1*aly5 - aly1*alx5 - alx1*aly7 + aly1*alx7 - alx3*aly5 + aly3*alx5 + alx3*aly7 - aly3*alx7 - alx1*myco + alx3*myco + alx5*myco - alx7*myco;
big_jac[14 ][ 12] = 2*alx1*aly2 - 2*alx2*aly1 - 2*alx1*aly3 + 2*aly1*alx3 + 2*alx2*aly3 - 2*alx3*aly2 - 2*alx1*aly6 + 2*alx2*aly5 + 2*aly1*alx6 - 2*aly2*alx5 + 2*alx1*aly7 - 2*aly1*alx7 - 2*alx3*aly5 + 2*aly3*alx5 - 2*alx2*aly7 + 2*alx3*aly6 + 2*aly2*alx7 - 2*aly3*alx6 + 2*alx5*aly6 - 2*alx6*aly5 - 2*alx5*aly7 + 2*aly5*alx7 + 2*alx6*aly7 - 2*alx7*aly6;
big_jac[15 ][ 12] = aly1*alx5 - alx1*aly5 + alx1*aly6 + alx2*aly5 - aly1*alx6 - aly2*alx5 - alx2*aly6 + aly2*alx6;
big_jac[16 ][ 12] = alx1*aly5 - aly1*alx5 - alx1*aly6 + aly1*alx6 - alx1*aly7 + aly1*alx7 - alx3*aly5 + aly3*alx5 + alx1*aly8 - aly1*alx8 + alx3*aly6 - aly3*alx6 + alx3*aly7 - aly3*alx7 - alx3*aly8 + aly3*alx8 + alx5*myco - alx6*myco - alx7*myco + alx8*myco;
big_jac[17 ][ 12] = aly1*alx5 - alx1*aly5 + alx1*aly6 + alx2*aly5 - aly1*alx6 - aly2*alx5 + alx1*aly7 - alx2*aly6 - aly1*alx7 + aly2*alx6 - alx1*aly8 - alx2*aly7 + aly1*alx8 + aly2*alx7 + alx2*aly8 - aly2*alx8;
big_jac[18 ][ 12] = alx2*aly1 - alx1*aly2 + alx1*aly3 - aly1*alx3 - alx2*aly3 + alx3*aly2 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 - alx5*aly6 + alx6*aly5 + alx5*aly7 - aly5*alx7 - alx6*aly7 + alx7*aly6;
big_jac[19 ][ 12] = aly1*alx5 - alx1*aly5 + alx1*aly6 - aly1*alx6 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx1*aly8 + aly1*alx8 - alx3*aly6 + aly3*alx6 - alx3*aly7 + aly3*alx7 + alx3*aly8 - aly3*alx8 - alx5*aly6 + alx6*aly5 + alx5*aly8 - alx6*aly7 - aly5*alx8 + alx7*aly6 - alx7*aly8 + alx8*aly7;
big_jac[20 ][ 12] = alx1*aly5 - aly1*alx5 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 - alx1*aly7 + alx2*aly6 + aly1*alx7 - aly2*alx6 + alx1*aly8 + alx2*aly7 - aly1*alx8 - aly2*alx7 - alx2*aly8 + aly2*alx8 + alx5*aly7 - aly5*alx7 - alx5*aly8 - alx6*aly7 + aly5*alx8 + alx7*aly6 + alx6*aly8 - aly6*alx8;

break;
case 5:

big_jac[1 ][ 13] = aly2*alz1 - aly1*alz2 + aly1*alz3 - alz1*aly3 - aly2*alz3 + aly3*alz2 - aly1*myco - alz1*myco + aly3*myco + alz2*myco + myco*myco;
big_jac[2 ][ 13] = aly1*alz2 - aly2*alz1 - aly1*alz3 + alz1*aly3 + 2*aly2*alz3 - 2*aly3*alz2 - aly1*alz6 + alz1*aly6 + aly1*alz7 - alz1*aly7 - aly2*alz7 + aly3*alz6 + alz2*aly7 - alz3*aly6 - aly3*myco - alz2*myco + aly7*myco + alz6*myco;
big_jac[3 ][ 13] = 2*aly1*alz2 - 2*aly2*alz1 - aly1*alz3 + alz1*aly3 - aly1*alz4 + 2*aly2*alz3 + alz1*aly4 - 2*aly3*alz2 + aly3*alz4 - aly4*alz3 + aly1*myco + 2*alz1*myco - aly3*myco - 2*alz2*myco - alz3*myco + alz4*myco - myco*myco;
big_jac[4 ][ 13] = aly1*alz2 - aly2*alz1 - 2*aly1*alz3 + 2*alz1*aly3 + aly1*alz4 + 2*aly2*alz3 - alz1*aly4 - 2*aly3*alz2 - aly2*alz4 + alz2*aly4 + 2*aly1*myco - aly2*myco + alz1*myco - 2*aly3*myco - alz2*myco + aly4*myco - myco*myco;
big_jac[5 ][ 13] = aly3*alz2 - aly2*alz3 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 - aly6*alz7 + aly7*alz6;
big_jac[6 ][ 13] = 2*aly2*alz1 - 2*aly1*alz2 + aly1*alz3 - alz1*aly3 + aly1*alz4 - 3*aly2*alz3 - alz1*aly4 + 3*aly3*alz2 + 2*aly1*alz6 - 2*alz1*aly6 - 2*aly3*alz4 + 2*aly4*alz3 - aly1*alz7 + alz1*aly7 - aly1*alz8 + aly2*alz7 + alz1*aly8 - 2*aly3*alz6 - alz2*aly7 + 2*alz3*aly6 + aly3*alz7 - alz3*aly7 + aly3*alz8 - aly4*alz7 - alz3*aly8 + alz4*aly7 + 2*alz2*myco - alz4*myco - 2*alz6*myco + alz8*myco;
big_jac[7 ][ 13] = aly2*alz1 - aly1*alz2 + 2*aly1*alz3 - 2*alz1*aly3 - aly1*alz4 - 3*aly2*alz3 + alz1*aly4 + 3*aly3*alz2 + 2*aly2*alz4 - 2*alz2*aly4 + aly1*alz6 - alz1*aly6 - 2*aly1*alz7 - aly2*alz6 + 2*alz1*aly7 + alz2*aly6 + aly1*alz8 + 2*aly2*alz7 - alz1*aly8 - aly3*alz6 - 2*alz2*aly7 + alz3*aly6 - aly2*alz8 + alz2*aly8 + aly4*alz6 - alz4*aly6 + 2*aly3*myco - aly4*myco - 2*aly7*myco + aly8*myco;
big_jac[8 ][ 13] = aly2*alz1 - aly1*alz2 + aly1*alz4 - aly2*alz3 - alz1*aly4 + aly3*alz2 - aly3*alz4 + aly4*alz3 - alz1*myco + alz2*myco + alz3*myco - alz4*myco;
big_jac[9 ][ 13] = 2*aly2*alz1 - 2*aly1*alz2 + 2*aly1*alz3 - 2*alz1*aly3 - 3*aly2*alz3 + 3*aly3*alz2 + aly2*alz4 - alz2*aly4 - aly3*alz4 + aly4*alz3 - 2*aly1*myco + aly2*myco - 2*alz1*myco + 2*aly3*myco + 2*alz2*myco - aly4*myco + alz3*myco - alz4*myco + myco*myco;
big_jac[10 ][ 13] = aly1*alz3 - alz1*aly3 - aly1*alz4 - aly2*alz3 + alz1*aly4 + aly3*alz2 + aly2*alz4 - alz2*aly4 - aly1*myco + aly2*myco + aly3*myco - aly4*myco;
big_jac[11 ][ 13] = aly2*alz3 - aly3*alz2 + aly3*alz4 - aly4*alz3 - aly2*alz7 + aly3*alz6 + alz2*aly7 - alz3*aly6 - aly3*alz8 + aly4*alz7 + alz3*aly8 - alz4*aly7 + aly6*alz7 - aly7*alz6 + aly7*alz8 - aly8*alz7;
big_jac[12 ][ 13] = aly2*alz3 - aly3*alz2 - aly2*alz4 + alz2*aly4 - aly2*alz7 + aly3*alz6 + alz2*aly7 - alz3*aly6 + aly2*alz8 - alz2*aly8 - aly4*alz6 + alz4*aly6 + aly6*alz7 - aly7*alz6 - aly6*alz8 + alz6*aly8;
big_jac[13 ][ 13] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 - aly1*alz6 + alz1*aly6 + aly3*alz4 - aly4*alz3 + aly1*alz8 - alz1*aly8 + aly3*alz6 - alz3*aly6 - aly3*alz8 + alz3*aly8 - alz2*myco + alz4*myco + alz6*myco - alz8*myco;
big_jac[14 ][ 13] = 2*aly1*alz2 - 2*aly2*alz1 - 2*aly1*alz3 + 2*alz1*aly3 + 4*aly2*alz3 - 4*aly3*alz2 - 2*aly2*alz4 + 2*alz2*aly4 - 2*aly1*alz6 + 2*alz1*aly6 + 2*aly3*alz4 - 2*aly4*alz3 + 2*aly1*alz7 + 2*aly2*alz6 - 2*alz1*aly7 - 2*alz2*aly6 - 2*aly2*alz7 + 2*aly3*alz6 + 2*alz2*aly7 - 2*alz3*aly6 - 2*aly3*alz7 - 2*aly4*alz6 + 2*alz3*aly7 + 2*alz4*aly6 + 2*aly4*alz7 - 2*alz4*aly7;
big_jac[15 ][ 13] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 - aly2*alz4 + alz2*aly4 + aly1*alz7 - alz1*aly7 - aly1*alz8 - aly2*alz7 + alz1*aly8 + alz2*aly7 + aly2*alz8 - alz2*aly8 - aly3*myco + aly4*myco + aly7*myco - aly8*myco;
big_jac[16 ][ 13] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 + aly3*alz4 - aly4*alz3 + alz1*myco - alz2*myco - alz3*myco + alz4*myco;
big_jac[17 ][ 13] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 - aly2*alz4 + alz2*aly4 + aly1*myco - aly2*myco - aly3*myco + aly4*myco;
big_jac[18 ][ 13] = aly3*alz2 - aly2*alz3 + aly2*alz4 - alz2*aly4 - aly3*alz4 + aly4*alz3 + aly2*alz7 - aly3*alz6 - alz2*aly7 + alz3*aly6 - aly2*alz8 + alz2*aly8 + aly4*alz6 - alz4*aly6 + aly3*alz8 - aly4*alz7 - alz3*aly8 + alz4*aly7 - aly6*alz7 + aly7*alz6 + aly6*alz8 - alz6*aly8 - aly7*alz8 + aly8*alz7;
big_jac[19 ][ 13] = aly2*alz1 - aly1*alz2 + aly1*alz4 - aly2*alz3 - alz1*aly4 + aly3*alz2 + aly1*alz6 - alz1*aly6 - aly3*alz4 + aly4*alz3 - aly2*alz6 + alz2*aly6 - aly1*alz8 + alz1*aly8 - aly3*alz6 + alz3*aly6 + aly2*alz8 - alz2*aly8 + aly4*alz6 - alz4*aly6 + aly3*alz8 - alz3*aly8 - aly4*alz8 + alz4*aly8;
big_jac[20 ][ 13] = aly1*alz3 - alz1*aly3 - aly1*alz4 - aly2*alz3 + alz1*aly4 + aly3*alz2 + aly2*alz4 - alz2*aly4 - aly1*alz7 + alz1*aly7 + aly1*alz8 + aly2*alz7 - alz1*aly8 - alz2*aly7 - aly2*alz8 + aly3*alz7 + alz2*aly8 - alz3*aly7 - aly3*alz8 - aly4*alz7 + alz3*aly8 + alz4*aly7 + aly4*alz8 - alz4*aly8;
big_jac[1 ][ 14] = alx1*alz2 - alx2*alz1 - alx1*alz3 + alx3*alz1 + alx2*alz3 - alx3*alz2 + alx1*myco - alx3*myco;
big_jac[2 ][ 14] = alx2*alz1 - alx1*alz2 + alx1*alz3 - alx3*alz1 - 2*alx2*alz3 + 2*alx3*alz2 + alx1*alz6 - alz1*alx6 - alx1*alz7 + alz1*alx7 + alx2*alz7 - alx3*alz6 - alz2*alx7 + alz3*alx6 + alx3*myco - alz2*myco + alz3*myco - alx7*myco - myco*myco;
big_jac[3 ][ 14] = 2*alx2*alz1 - 2*alx1*alz2 + alx1*alz3 - alx3*alz1 + alx1*alz4 - 2*alx2*alz3 + 2*alx3*alz2 - alz1*alx4 - alx3*alz4 + alx4*alz3 - alx1*myco + alx3*myco;
big_jac[4 ][ 14] = alx2*alz1 - alx1*alz2 + 2*alx1*alz3 - 2*alx3*alz1 - alx1*alz4 - 2*alx2*alz3 + 2*alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - 2*alx1*myco + alx2*myco + 2*alx3*myco - alx4*myco;
big_jac[5 ][ 14] = alx2*alz3 - alx3*alz2 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 + alx6*alz7 - alx7*alz6 + alz2*myco - alz3*myco - alz6*myco + alz7*myco;
big_jac[6 ][ 14] = 2*alx1*alz2 - 2*alx2*alz1 - alx1*alz3 + alx3*alz1 - alx1*alz4 + 3*alx2*alz3 - 3*alx3*alz2 + alz1*alx4 - 2*alx1*alz6 + 2*alx3*alz4 + 2*alz1*alx6 - 2*alx4*alz3 + alx1*alz7 - alz1*alx7 + alx1*alz8 - alx2*alz7 + 2*alx3*alz6 - alz1*alx8 + alz2*alx7 - 2*alz3*alx6 - alx3*alz7 + alz3*alx7 - alx3*alz8 + alx4*alz7 + alz3*alx8 - alz4*alx7 + alz2*myco - alz4*myco;
big_jac[7 ][ 14] = alx1*alz2 - alx2*alz1 - 2*alx1*alz3 + 2*alx3*alz1 + alx1*alz4 + 3*alx2*alz3 - 3*alx3*alz2 - alz1*alx4 - 2*alx2*alz4 + 2*alx4*alz2 - alx1*alz6 + alz1*alx6 + 2*alx1*alz7 + alx2*alz6 - 2*alz1*alx7 - alz2*alx6 - alx1*alz8 - 2*alx2*alz7 + alx3*alz6 + alz1*alx8 + 2*alz2*alx7 - alz3*alx6 + alx2*alz8 - alx4*alz6 - alz2*alx8 + alx6*alz4 - 2*alx3*myco + alx4*myco - alz3*myco + alz4*myco + 2*alx7*myco - alx8*myco + myco*myco;
big_jac[8 ][ 14] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alx2*alz3 - alx3*alz2 + alz1*alx4 + alx3*alz4 - alx4*alz3;
big_jac[9 ][ 14] = 2*alx1*alz2 - 2*alx2*alz1 - 2*alx1*alz3 + 2*alx3*alz1 + 3*alx2*alz3 - 3*alx3*alz2 - alx2*alz4 + alx4*alz2 + alx3*alz4 - alx4*alz3 + 2*alx1*myco - alx2*myco - 2*alx3*myco + alx4*myco;
big_jac[10 ][ 14] = alx3*alz1 - alx1*alz3 + alx1*alz4 + alx2*alz3 - alx3*alz2 - alz1*alx4 - alx2*alz4 + alx4*alz2 + alx1*myco - alx2*myco - alx3*myco + alx4*myco;
big_jac[11 ][ 14] = alx3*alz2 - alx2*alz3 - alx3*alz4 + alx4*alz3 + alx2*alz7 - alx3*alz6 - alz2*alx7 + alz3*alx6 + alx3*alz8 - alx4*alz7 - alz3*alx8 + alz4*alx7 - alx6*alz7 + alx7*alz6 - alx7*alz8 + alx8*alz7 - alz2*myco + alz4*myco + alz6*myco - alz8*myco;
big_jac[12 ][ 14] = alx3*alz2 - alx2*alz3 + alx2*alz4 - alx4*alz2 + alx2*alz7 - alx3*alz6 - alz2*alx7 + alz3*alx6 - alx2*alz8 + alx4*alz6 + alz2*alx8 - alx6*alz4 - alx6*alz7 + alx7*alz6 + alx6*alz8 - alx8*alz6 + alz3*myco - alz4*myco - alz7*myco + alz8*myco;
big_jac[13 ][ 14] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 + alx1*alz6 - alx3*alz4 - alz1*alx6 + alx4*alz3 - alx1*alz8 - alx3*alz6 + alz1*alx8 + alz3*alx6 + alx3*alz8 - alz3*alx8;
big_jac[14 ][ 14] = 2*alx2*alz1 - 2*alx1*alz2 + 2*alx1*alz3 - 2*alx3*alz1 - 4*alx2*alz3 + 4*alx3*alz2 + 2*alx2*alz4 - 2*alx4*alz2 + 2*alx1*alz6 - 2*alx3*alz4 - 2*alz1*alx6 + 2*alx4*alz3 - 2*alx1*alz7 - 2*alx2*alz6 + 2*alz1*alx7 + 2*alz2*alx6 + 2*alx2*alz7 - 2*alx3*alz6 - 2*alz2*alx7 + 2*alz3*alx6 + 2*alx3*alz7 + 2*alx4*alz6 - 2*alz3*alx7 - 2*alx6*alz4 - 2*alx4*alz7 + 2*alz4*alx7;
big_jac[15 ][ 14] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - alx1*alz7 + alz1*alx7 + alx1*alz8 + alx2*alz7 - alz1*alx8 - alz2*alx7 - alx2*alz8 + alz2*alx8 + alx3*myco - alx4*myco - alx7*myco + alx8*myco;
big_jac[16 ][ 14] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 - alx3*alz4 + alx4*alz3;
big_jac[17 ][ 14] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - alx1*myco + alx2*myco + alx3*myco - alx4*myco;
big_jac[18 ][ 14] = alx2*alz3 - alx3*alz2 - alx2*alz4 + alx4*alz2 + alx3*alz4 - alx4*alz3 - alx2*alz7 + alx3*alz6 + alz2*alx7 - alz3*alx6 + alx2*alz8 - alx4*alz6 - alz2*alx8 + alx6*alz4 - alx3*alz8 + alx4*alz7 + alz3*alx8 - alz4*alx7 + alx6*alz7 - alx7*alz6 - alx6*alz8 + alx8*alz6 + alx7*alz8 - alx8*alz7;
big_jac[19 ][ 14] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alx2*alz3 - alx3*alz2 + alz1*alx4 - alx1*alz6 + alx3*alz4 + alz1*alx6 - alx4*alz3 + alx2*alz6 - alz2*alx6 + alx1*alz8 + alx3*alz6 - alz1*alx8 - alz3*alx6 - alx2*alz8 - alx4*alz6 + alz2*alx8 + alx6*alz4 - alx3*alz8 + alz3*alx8 + alx4*alz8 - alz4*alx8;
big_jac[20 ][ 14] = alx3*alz1 - alx1*alz3 + alx1*alz4 + alx2*alz3 - alx3*alz2 - alz1*alx4 - alx2*alz4 + alx4*alz2 + alx1*alz7 - alz1*alx7 - alx1*alz8 - alx2*alz7 + alz1*alx8 + alz2*alx7 + alx2*alz8 - alx3*alz7 - alz2*alx8 + alz3*alx7 + alx3*alz8 + alx4*alz7 - alz3*alx8 - alz4*alx7 - alx4*alz8 + alz4*alx8;
big_jac[1 ][ 15] = alx2*aly1 - alx1*aly2 + alx1*aly3 - aly1*alx3 - alx2*aly3 + alx3*aly2 + alx1*myco - alx2*myco;
big_jac[2 ][ 15] = alx1*aly2 - alx2*aly1 - alx1*aly3 + aly1*alx3 + 2*alx2*aly3 - 2*alx3*aly2 - alx1*aly6 + aly1*alx6 + alx1*aly7 - aly1*alx7 - alx2*aly7 + alx3*aly6 + aly2*alx7 - aly3*alx6 + alx2*myco + aly2*myco - aly3*myco - alx6*myco - myco*myco;
big_jac[3 ][ 15] = 2*alx1*aly2 - 2*alx2*aly1 - alx1*aly3 + aly1*alx3 - alx1*aly4 + 2*alx2*aly3 + aly1*alx4 - 2*alx3*aly2 + alx3*aly4 - alx4*aly3 - 2*alx1*myco + 2*alx2*myco + alx3*myco - alx4*myco;
big_jac[4 ][ 15] = alx1*aly2 - alx2*aly1 - 2*alx1*aly3 + 2*aly1*alx3 + alx1*aly4 + 2*alx2*aly3 - aly1*alx4 - 2*alx3*aly2 - alx2*aly4 + aly2*alx4 - alx1*myco + alx2*myco;
big_jac[5 ][ 15] = alx3*aly2 - alx2*aly3 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 - alx6*aly7 + alx7*aly6 - aly2*myco + aly3*myco + aly6*myco - aly7*myco;
big_jac[6 ][ 15] = 2*alx2*aly1 - 2*alx1*aly2 + alx1*aly3 - aly1*alx3 + alx1*aly4 - 3*alx2*aly3 - aly1*alx4 + 3*alx3*aly2 + 2*alx1*aly6 - 2*aly1*alx6 - 2*alx3*aly4 + 2*alx4*aly3 - alx1*aly7 + aly1*alx7 - alx1*aly8 + alx2*aly7 + aly1*alx8 - 2*alx3*aly6 - aly2*alx7 + 2*aly3*alx6 + alx3*aly7 - aly3*alx7 + alx3*aly8 - alx4*aly7 - aly3*alx8 + aly4*alx7 - 2*alx2*myco - aly2*myco + alx4*myco + aly4*myco + 2*alx6*myco - alx8*myco + myco*myco;
big_jac[7 ][ 15] = alx2*aly1 - alx1*aly2 + 2*alx1*aly3 - 2*aly1*alx3 - alx1*aly4 - 3*alx2*aly3 + aly1*alx4 + 3*alx3*aly2 + 2*alx2*aly4 - 2*aly2*alx4 + alx1*aly6 - aly1*alx6 - 2*alx1*aly7 - alx2*aly6 + 2*aly1*alx7 + aly2*alx6 + alx1*aly8 + 2*alx2*aly7 - aly1*alx8 - alx3*aly6 - 2*aly2*alx7 + aly3*alx6 - alx2*aly8 + aly2*alx8 + alx4*aly6 - aly4*alx6 + aly3*myco - aly4*myco;
big_jac[8 ][ 15] = alx2*aly1 - alx1*aly2 + alx1*aly4 - alx2*aly3 - aly1*alx4 + alx3*aly2 - alx3*aly4 + alx4*aly3 + alx1*myco - alx2*myco - alx3*myco + alx4*myco;
big_jac[9 ][ 15] = 2*alx2*aly1 - 2*alx1*aly2 + 2*alx1*aly3 - 2*aly1*alx3 - 3*alx2*aly3 + 3*alx3*aly2 + alx2*aly4 - aly2*alx4 - alx3*aly4 + alx4*aly3 + 2*alx1*myco - 2*alx2*myco - alx3*myco + alx4*myco;
big_jac[10 ][ 15] = alx1*aly3 - aly1*alx3 - alx1*aly4 - alx2*aly3 + aly1*alx4 + alx3*aly2 + alx2*aly4 - aly2*alx4;
big_jac[11 ][ 15] = alx2*aly3 - alx3*aly2 + alx3*aly4 - alx4*aly3 - alx2*aly7 + alx3*aly6 + aly2*alx7 - aly3*alx6 - alx3*aly8 + alx4*aly7 + aly3*alx8 - aly4*alx7 + alx6*aly7 - alx7*aly6 + alx7*aly8 - alx8*aly7 + aly2*myco - aly4*myco - aly6*myco + aly8*myco;
big_jac[12 ][ 15] = alx2*aly3 - alx3*aly2 - alx2*aly4 + aly2*alx4 - alx2*aly7 + alx3*aly6 + aly2*alx7 - aly3*alx6 + alx2*aly8 - aly2*alx8 - alx4*aly6 + aly4*alx6 + alx6*aly7 - alx7*aly6 - alx6*aly8 + aly6*alx8 - aly3*myco + aly4*myco + aly7*myco - aly8*myco;
big_jac[13 ][ 15] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 - alx1*aly6 + aly1*alx6 + alx3*aly4 - alx4*aly3 + alx1*aly8 - aly1*alx8 + alx3*aly6 - aly3*alx6 - alx3*aly8 + aly3*alx8 + alx2*myco - alx4*myco - alx6*myco + alx8*myco;
big_jac[14 ][ 15] = 2*alx1*aly2 - 2*alx2*aly1 - 2*alx1*aly3 + 2*aly1*alx3 + 4*alx2*aly3 - 4*alx3*aly2 - 2*alx2*aly4 + 2*aly2*alx4 - 2*alx1*aly6 + 2*aly1*alx6 + 2*alx3*aly4 - 2*alx4*aly3 + 2*alx1*aly7 + 2*alx2*aly6 - 2*aly1*alx7 - 2*aly2*alx6 - 2*alx2*aly7 + 2*alx3*aly6 + 2*aly2*alx7 - 2*aly3*alx6 - 2*alx3*aly7 - 2*alx4*aly6 + 2*aly3*alx7 + 2*aly4*alx6 + 2*alx4*aly7 - 2*aly4*alx7;
big_jac[15 ][ 15] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 - alx2*aly4 + aly2*alx4 + alx1*aly7 - aly1*alx7 - alx1*aly8 - alx2*aly7 + aly1*alx8 + aly2*alx7 + alx2*aly8 - aly2*alx8;
big_jac[16 ][ 15] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 + alx3*aly4 - alx4*aly3 - alx1*myco + alx2*myco + alx3*myco - alx4*myco;
big_jac[17 ][ 15] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 - alx2*aly4 + aly2*alx4;
big_jac[18 ][ 15] = alx3*aly2 - alx2*aly3 + alx2*aly4 - aly2*alx4 - alx3*aly4 + alx4*aly3 + alx2*aly7 - alx3*aly6 - aly2*alx7 + aly3*alx6 - alx2*aly8 + aly2*alx8 + alx4*aly6 - aly4*alx6 + alx3*aly8 - alx4*aly7 - aly3*alx8 + aly4*alx7 - alx6*aly7 + alx7*aly6 + alx6*aly8 - aly6*alx8 - alx7*aly8 + alx8*aly7;
big_jac[19 ][ 15] = alx2*aly1 - alx1*aly2 + alx1*aly4 - alx2*aly3 - aly1*alx4 + alx3*aly2 + alx1*aly6 - aly1*alx6 - alx3*aly4 + alx4*aly3 - alx2*aly6 + aly2*alx6 - alx1*aly8 + aly1*alx8 - alx3*aly6 + aly3*alx6 + alx2*aly8 - aly2*alx8 + alx4*aly6 - aly4*alx6 + alx3*aly8 - aly3*alx8 - alx4*aly8 + aly4*alx8;
big_jac[20 ][ 15] = alx1*aly3 - aly1*alx3 - alx1*aly4 - alx2*aly3 + aly1*alx4 + alx3*aly2 + alx2*aly4 - aly2*alx4 - alx1*aly7 + aly1*alx7 + alx1*aly8 + alx2*aly7 - aly1*alx8 - aly2*alx7 - alx2*aly8 + alx3*aly7 + aly2*alx8 - aly3*alx7 - alx3*aly8 - alx4*aly7 + aly3*alx8 + aly4*alx7 + alx4*aly8 - aly4*alx8;

break;
case 6:

big_jac[1 ][ 16] = 0;
big_jac[2 ][ 16] = alz1*aly3 - aly1*alz3 + aly1*alz5 - alz1*aly5 - aly3*alz5 + alz3*aly5 + alz1*myco - alz5*myco;
big_jac[3 ][ 16] = 0;
big_jac[4 ][ 16] = aly2*alz1 - aly1*alz2 + aly1*alz3 - alz1*aly3 - aly2*alz3 + aly3*alz2 - aly1*myco - alz1*myco + aly3*myco + alz2*myco + myco*myco;
big_jac[5 ][ 16] = aly1*alz3 - alz1*aly3 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly5*alz7 - alz5*aly7;
big_jac[6 ][ 16] = aly1*alz3 - alz1*aly3 - 2*aly1*alz5 + 2*alz1*aly5 + aly1*alz7 - alz1*aly7 + 2*aly3*alz5 - 2*alz3*aly5 - aly3*alz7 + alz3*aly7 - 2*alz1*myco + alz3*myco + 2*alz5*myco - alz7*myco;
big_jac[7 ][ 16] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 - aly1*alz5 + alz1*aly5 + aly2*alz5 - alz2*aly5 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 - aly3*myco + aly7*myco;
big_jac[8 ][ 16] = 0;
big_jac[9 ][ 16] = 2*aly1*alz2 - 2*aly2*alz1 - aly1*alz3 + alz1*aly3 - aly1*alz4 + 2*aly2*alz3 + alz1*aly4 - 2*aly3*alz2 + aly3*alz4 - aly4*alz3 + aly1*myco + 2*alz1*myco - aly3*myco - 2*alz2*myco - alz3*myco + alz4*myco - myco*myco;
big_jac[10 ][ 16] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 - aly2*alz4 + alz2*aly4 + aly1*myco - aly2*myco - aly3*myco + aly4*myco;
big_jac[11 ][ 16] = alz1*aly3 - aly1*alz3 + aly1*alz7 - alz1*aly7 - aly3*alz5 + alz3*aly5 - aly5*alz7 + alz5*aly7;
big_jac[12 ][ 16] = alz1*aly3 - aly1*alz3 + aly1*alz4 - alz1*aly4 + aly1*alz7 - alz1*aly7 - aly3*alz5 + alz3*aly5 - aly1*alz8 + alz1*aly8 + aly4*alz5 - aly5*alz4 - aly5*alz7 + alz5*aly7 + aly5*alz8 - alz5*aly8;
big_jac[13 ][ 16] = aly1*alz5 - alz1*aly5 - aly1*alz7 + alz1*aly7 - aly3*alz5 + alz3*aly5 + aly3*alz7 - alz3*aly7 + alz1*myco - alz3*myco - alz5*myco + alz7*myco;
big_jac[14 ][ 16] = 2*aly2*alz1 - 2*aly1*alz2 + 2*aly1*alz4 - 2*aly2*alz3 - 2*alz1*aly4 + 2*aly3*alz2 + 2*aly1*alz5 - 2*alz1*aly5 - 2*aly2*alz5 - 2*aly3*alz4 + 2*alz2*aly5 + 2*aly4*alz3 - 2*aly1*alz7 + 2*alz1*aly7 - 2*aly3*alz5 + 2*alz3*aly5 + 2*aly2*alz7 - 2*alz2*aly7 + 2*aly4*alz5 - 2*aly5*alz4 + 2*aly3*alz7 - 2*alz3*aly7 - 2*aly4*alz7 + 2*alz4*aly7;
big_jac[15 ][ 16] = aly1*alz3 - alz1*aly3 - aly1*alz4 - aly2*alz3 + alz1*aly4 + aly3*alz2 + aly2*alz4 - alz2*aly4 - aly1*alz7 + alz1*aly7 + aly1*alz8 + aly2*alz7 - alz1*aly8 - alz2*aly7 - aly2*alz8 + alz2*aly8 + aly3*myco - aly4*myco - aly7*myco + aly8*myco;
big_jac[16 ][ 16] = aly2*alz1 - aly1*alz2 + aly1*alz4 - aly2*alz3 - alz1*aly4 + aly3*alz2 - aly3*alz4 + aly4*alz3 - alz1*myco + alz2*myco + alz3*myco - alz4*myco;
big_jac[17 ][ 16] = aly1*alz3 - alz1*aly3 - aly1*alz4 - aly2*alz3 + alz1*aly4 + aly3*alz2 + aly2*alz4 - alz2*aly4 - aly1*myco + aly2*myco + aly3*myco - aly4*myco;
big_jac[18 ][ 16] = aly1*alz3 - alz1*aly3 - aly1*alz4 + alz1*aly4 + aly3*alz4 - aly4*alz3 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly1*alz8 - alz1*aly8 - aly4*alz5 + aly5*alz4 - aly3*alz8 + aly4*alz7 + alz3*aly8 - alz4*aly7 + aly5*alz7 - alz5*aly7 - aly5*alz8 + alz5*aly8 + aly7*alz8 - aly8*alz7;
big_jac[19 ][ 16] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 - aly1*alz5 + alz1*aly5 + aly2*alz5 + aly3*alz4 - alz2*aly5 - aly4*alz3 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly2*alz7 + alz2*aly7 - aly4*alz5 + aly5*alz4 - aly3*alz7 + alz3*aly7 + aly4*alz7 - alz4*aly7;
big_jac[20 ][ 16] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 - aly2*alz4 + alz2*aly4 + aly1*alz7 - alz1*aly7 - aly1*alz8 - aly2*alz7 + alz1*aly8 + alz2*aly7 + aly2*alz8 - aly3*alz7 - alz2*aly8 + alz3*aly7 + aly3*alz8 + aly4*alz7 - alz3*aly8 - alz4*aly7 - aly4*alz8 + alz4*aly8;
big_jac[1 ][ 17] = 0;
big_jac[2 ][ 17] = alx1*alz3 - alx3*alz1 - alx1*alz5 + alz1*alx5 + alx3*alz5 - alx5*alz3 + alz1*myco - alz3*myco;
big_jac[3 ][ 17] = 0;
big_jac[4 ][ 17] = alx1*alz2 - alx2*alz1 - alx1*alz3 + alx3*alz1 + alx2*alz3 - alx3*alz2 + alx1*myco - alx3*myco;
big_jac[5 ][ 17] = alx3*alz1 - alx1*alz3 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx5*alz7 + alx7*alz5 - alz1*myco + alz3*myco + alz5*myco - alz7*myco;
big_jac[6 ][ 17] = alx3*alz1 - alx1*alz3 + 2*alx1*alz5 - 2*alz1*alx5 - alx1*alz7 - 2*alx3*alz5 + alz1*alx7 + 2*alx5*alz3 + alx3*alz7 - alz3*alx7 - alz1*myco + alz3*myco;
big_jac[7 ][ 17] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 + alx1*alz5 - alz1*alx5 - alx2*alz5 + alz2*alx5 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 + alx3*myco + alz3*myco - alz4*myco - alx7*myco - myco*myco;
big_jac[8 ][ 17] = 0;
big_jac[9 ][ 17] = 2*alx2*alz1 - 2*alx1*alz2 + alx1*alz3 - alx3*alz1 + alx1*alz4 - 2*alx2*alz3 + 2*alx3*alz2 - alz1*alx4 - alx3*alz4 + alx4*alz3 - alx1*myco + alx3*myco;
big_jac[10 ][ 17] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - alx1*myco + alx2*myco + alx3*myco - alx4*myco;
big_jac[11 ][ 17] = alx1*alz3 - alx3*alz1 - alx1*alz7 + alx3*alz5 + alz1*alx7 - alx5*alz3 + alx5*alz7 - alx7*alz5 + alz1*myco - alz3*myco - alz5*myco + alz7*myco;
big_jac[12 ][ 17] = alx1*alz3 - alx3*alz1 - alx1*alz4 + alz1*alx4 - alx1*alz7 + alx3*alz5 + alz1*alx7 - alx5*alz3 + alx1*alz8 - alz1*alx8 - alx4*alz5 + alx5*alz4 + alx5*alz7 - alx7*alz5 - alx5*alz8 + alz5*alx8 - alz3*myco + alz4*myco + alz7*myco - alz8*myco;
big_jac[13 ][ 17] = alz1*alx5 - alx1*alz5 + alx1*alz7 + alx3*alz5 - alz1*alx7 - alx5*alz3 - alx3*alz7 + alz3*alx7;
big_jac[14 ][ 17] = 2*alx1*alz2 - 2*alx2*alz1 - 2*alx1*alz4 + 2*alx2*alz3 - 2*alx3*alz2 + 2*alz1*alx4 - 2*alx1*alz5 + 2*alz1*alx5 + 2*alx2*alz5 + 2*alx3*alz4 - 2*alx4*alz3 - 2*alz2*alx5 + 2*alx1*alz7 + 2*alx3*alz5 - 2*alz1*alx7 - 2*alx5*alz3 - 2*alx2*alz7 - 2*alx4*alz5 + 2*alz2*alx7 + 2*alx5*alz4 - 2*alx3*alz7 + 2*alz3*alx7 + 2*alx4*alz7 - 2*alz4*alx7;
big_jac[15 ][ 17] = alx3*alz1 - alx1*alz3 + alx1*alz4 + alx2*alz3 - alx3*alz2 - alz1*alx4 - alx2*alz4 + alx4*alz2 + alx1*alz7 - alz1*alx7 - alx1*alz8 - alx2*alz7 + alz1*alx8 + alz2*alx7 + alx2*alz8 - alz2*alx8 - alx3*myco + alx4*myco + alx7*myco - alx8*myco;
big_jac[16 ][ 17] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alx2*alz3 - alx3*alz2 + alz1*alx4 + alx3*alz4 - alx4*alz3;
big_jac[17 ][ 17] = alx3*alz1 - alx1*alz3 + alx1*alz4 + alx2*alz3 - alx3*alz2 - alz1*alx4 - alx2*alz4 + alx4*alz2 + alx1*myco - alx2*myco - alx3*myco + alx4*myco;
big_jac[18 ][ 17] = alx3*alz1 - alx1*alz3 + alx1*alz4 - alz1*alx4 - alx3*alz4 + alx4*alz3 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx1*alz8 + alz1*alx8 + alx4*alz5 - alx5*alz4 + alx3*alz8 - alx4*alz7 - alz3*alx8 + alz4*alx7 - alx5*alz7 + alx7*alz5 + alx5*alz8 - alz5*alx8 - alx7*alz8 + alx8*alz7;
big_jac[19 ][ 17] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 + alx1*alz5 - alz1*alx5 - alx2*alz5 - alx3*alz4 + alx4*alz3 + alz2*alx5 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx2*alz7 + alx4*alz5 - alz2*alx7 - alx5*alz4 + alx3*alz7 - alz3*alx7 - alx4*alz7 + alz4*alx7;
big_jac[20 ][ 17] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - alx1*alz7 + alz1*alx7 + alx1*alz8 + alx2*alz7 - alz1*alx8 - alz2*alx7 - alx2*alz8 + alx3*alz7 + alz2*alx8 - alz3*alx7 - alx3*alz8 - alx4*alz7 + alz3*alx8 + alz4*alx7 + alx4*alz8 - alz4*alx8;
big_jac[1 ][ 18] = 0;
big_jac[2 ][ 18] = aly1*alx3 - alx1*aly3 + alx1*aly5 - aly1*alx5 - alx3*aly5 + aly3*alx5 - alx1*myco - aly1*myco + aly3*myco + alx5*myco + myco*myco;
big_jac[3 ][ 18] = 0;
big_jac[4 ][ 18] = alx2*aly1 - alx1*aly2 + alx1*aly3 - aly1*alx3 - alx2*aly3 + alx3*aly2 + alx1*myco - alx2*myco;
big_jac[5 ][ 18] = alx1*aly3 - aly1*alx3 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx5*aly7 - aly5*alx7 + aly1*myco - aly3*myco - aly5*myco + aly7*myco;
big_jac[6 ][ 18] = alx1*aly3 - aly1*alx3 - 2*alx1*aly5 + 2*aly1*alx5 + alx1*aly7 - aly1*alx7 + 2*alx3*aly5 - 2*aly3*alx5 - alx3*aly7 + aly3*alx7 + 2*alx1*myco + aly1*myco - alx3*myco - aly3*myco - 2*alx5*myco + alx7*myco - myco*myco;
big_jac[7 ][ 18] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 - alx1*aly5 + aly1*alx5 + alx2*aly5 - aly2*alx5 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 - aly3*myco + aly4*myco;
big_jac[8 ][ 18] = 0;
big_jac[9 ][ 18] = 2*alx1*aly2 - 2*alx2*aly1 - alx1*aly3 + aly1*alx3 - alx1*aly4 + 2*alx2*aly3 + aly1*alx4 - 2*alx3*aly2 + alx3*aly4 - alx4*aly3 - 2*alx1*myco + 2*alx2*myco + alx3*myco - alx4*myco;
big_jac[10 ][ 18] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 - alx2*aly4 + aly2*alx4;
big_jac[11 ][ 18] = aly1*alx3 - alx1*aly3 + alx1*aly7 - aly1*alx7 - alx3*aly5 + aly3*alx5 - alx5*aly7 + aly5*alx7 - aly1*myco + aly3*myco + aly5*myco - aly7*myco;
big_jac[12 ][ 18] = aly1*alx3 - alx1*aly3 + alx1*aly4 - aly1*alx4 + alx1*aly7 - aly1*alx7 - alx3*aly5 + aly3*alx5 - alx1*aly8 + aly1*alx8 + alx4*aly5 - alx5*aly4 - alx5*aly7 + aly5*alx7 + alx5*aly8 - aly5*alx8 + aly3*myco - aly4*myco - aly7*myco + aly8*myco;
big_jac[13 ][ 18] = alx1*aly5 - aly1*alx5 - alx1*aly7 + aly1*alx7 - alx3*aly5 + aly3*alx5 + alx3*aly7 - aly3*alx7 - alx1*myco + alx3*myco + alx5*myco - alx7*myco;
big_jac[14 ][ 18] = 2*alx2*aly1 - 2*alx1*aly2 + 2*alx1*aly4 - 2*alx2*aly3 - 2*aly1*alx4 + 2*alx3*aly2 + 2*alx1*aly5 - 2*aly1*alx5 - 2*alx2*aly5 - 2*alx3*aly4 + 2*aly2*alx5 + 2*alx4*aly3 - 2*alx1*aly7 + 2*aly1*alx7 - 2*alx3*aly5 + 2*aly3*alx5 + 2*alx2*aly7 - 2*aly2*alx7 + 2*alx4*aly5 - 2*alx5*aly4 + 2*alx3*aly7 - 2*aly3*alx7 - 2*alx4*aly7 + 2*aly4*alx7;
big_jac[15 ][ 18] = alx1*aly3 - aly1*alx3 - alx1*aly4 - alx2*aly3 + aly1*alx4 + alx3*aly2 + alx2*aly4 - aly2*alx4 - alx1*aly7 + aly1*alx7 + alx1*aly8 + alx2*aly7 - aly1*alx8 - aly2*alx7 - alx2*aly8 + aly2*alx8;
big_jac[16 ][ 18] = alx2*aly1 - alx1*aly2 + alx1*aly4 - alx2*aly3 - aly1*alx4 + alx3*aly2 - alx3*aly4 + alx4*aly3 + alx1*myco - alx2*myco - alx3*myco + alx4*myco;
big_jac[17 ][ 18] = alx1*aly3 - aly1*alx3 - alx1*aly4 - alx2*aly3 + aly1*alx4 + alx3*aly2 + alx2*aly4 - aly2*alx4;
big_jac[18 ][ 18] = alx1*aly3 - aly1*alx3 - alx1*aly4 + aly1*alx4 + alx3*aly4 - alx4*aly3 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx1*aly8 - aly1*alx8 - alx4*aly5 + alx5*aly4 - alx3*aly8 + alx4*aly7 + aly3*alx8 - aly4*alx7 + alx5*aly7 - aly5*alx7 - alx5*aly8 + aly5*alx8 + alx7*aly8 - alx8*aly7;
big_jac[19 ][ 18] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 - alx1*aly5 + aly1*alx5 + alx2*aly5 + alx3*aly4 - aly2*alx5 - alx4*aly3 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx2*aly7 + aly2*alx7 - alx4*aly5 + alx5*aly4 - alx3*aly7 + aly3*alx7 + alx4*aly7 - aly4*alx7;
big_jac[20 ][ 18] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 - alx2*aly4 + aly2*alx4 + alx1*aly7 - aly1*alx7 - alx1*aly8 - alx2*aly7 + aly1*alx8 + aly2*alx7 + alx2*aly8 - alx3*aly7 - aly2*alx8 + aly3*alx7 + alx3*aly8 + alx4*aly7 - aly3*alx8 - aly4*alx7 - alx4*aly8 + aly4*alx8;

break;
case 7:

big_jac[1 ][ 19] = 0;
big_jac[2 ][ 19] = aly1*alz2 - aly2*alz1 - aly1*alz5 + alz1*aly5 + aly2*alz5 - alz2*aly5 + aly1*myco - aly5*myco;
big_jac[3 ][ 19] = aly2*alz1 - aly1*alz2 + aly1*alz3 - alz1*aly3 - aly2*alz3 + aly3*alz2 - aly1*myco - alz1*myco + aly3*myco + alz2*myco + myco*myco;
big_jac[4 ][ 19] = 0;
big_jac[5 ][ 19] = aly2*alz1 - aly1*alz2 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly5*alz6 + aly6*alz5;
big_jac[6 ][ 19] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 + aly1*alz5 - alz1*aly5 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 - aly3*alz5 + alz3*aly5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 - alz2*myco + alz6*myco;
big_jac[7 ][ 19] = aly2*alz1 - aly1*alz2 + 2*aly1*alz5 - 2*alz1*aly5 - aly1*alz6 - 2*aly2*alz5 + alz1*aly6 + 2*alz2*aly5 + aly2*alz6 - alz2*aly6 - 2*aly1*myco + aly2*myco + 2*aly5*myco - aly6*myco;
big_jac[8 ][ 19] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 + aly3*alz4 - aly4*alz3 + alz1*myco - alz2*myco - alz3*myco + alz4*myco;
big_jac[9 ][ 19] = aly1*alz2 - aly2*alz1 - 2*aly1*alz3 + 2*alz1*aly3 + aly1*alz4 + 2*aly2*alz3 - alz1*aly4 - 2*aly3*alz2 - aly2*alz4 + alz2*aly4 + 2*aly1*myco - aly2*myco + alz1*myco - 2*aly3*myco - alz2*myco + aly4*myco - myco*myco;
big_jac[10 ][ 19] = 0;
big_jac[11 ][ 19] = aly1*alz2 - aly2*alz1 - aly1*alz4 + alz1*aly4 - aly1*alz6 + aly2*alz5 + alz1*aly6 - alz2*aly5 + aly1*alz8 - alz1*aly8 - aly4*alz5 + aly5*alz4 + aly5*alz6 - aly6*alz5 - aly5*alz8 + alz5*aly8;
big_jac[12 ][ 19] = aly1*alz2 - aly2*alz1 - aly1*alz6 + aly2*alz5 + alz1*aly6 - alz2*aly5 + aly5*alz6 - aly6*alz5;
big_jac[13 ][ 19] = aly2*alz1 - aly1*alz2 + aly1*alz4 - aly2*alz3 - alz1*aly4 + aly3*alz2 + aly1*alz6 - alz1*aly6 - aly3*alz4 + aly4*alz3 - aly1*alz8 + alz1*aly8 - aly3*alz6 + alz3*aly6 + aly3*alz8 - alz3*aly8 + alz2*myco - alz4*myco - alz6*myco + alz8*myco;
big_jac[14 ][ 19] = 2*aly1*alz3 - 2*alz1*aly3 - 2*aly1*alz4 - 2*aly2*alz3 + 2*alz1*aly4 + 2*aly3*alz2 - 2*aly1*alz5 + 2*aly2*alz4 + 2*alz1*aly5 - 2*alz2*aly4 + 2*aly1*alz6 + 2*aly2*alz5 - 2*alz1*aly6 - 2*alz2*aly5 - 2*aly2*alz6 + 2*aly3*alz5 + 2*alz2*aly6 - 2*alz3*aly5 - 2*aly3*alz6 - 2*aly4*alz5 + 2*alz3*aly6 + 2*aly5*alz4 + 2*aly4*alz6 - 2*alz4*aly6;
big_jac[15 ][ 19] = alz1*aly5 - aly1*alz5 + aly1*alz6 + aly2*alz5 - alz1*aly6 - alz2*aly5 - aly2*alz6 + alz2*aly6 + aly1*myco - aly2*myco - aly5*myco + aly6*myco;
big_jac[16 ][ 19] = aly2*alz1 - aly1*alz2 + aly1*alz4 - aly2*alz3 - alz1*aly4 + aly3*alz2 - aly3*alz4 + aly4*alz3 - alz1*myco + alz2*myco + alz3*myco - alz4*myco;
big_jac[17 ][ 19] = aly1*alz3 - alz1*aly3 - aly1*alz4 - aly2*alz3 + alz1*aly4 + aly3*alz2 + aly2*alz4 - alz2*aly4 - aly1*myco + aly2*myco + aly3*myco - aly4*myco;
big_jac[18 ][ 19] = aly2*alz1 - aly1*alz2 + aly1*alz4 - alz1*aly4 - aly2*alz4 + alz2*aly4 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly1*alz8 + alz1*aly8 + aly4*alz5 - aly5*alz4 + aly2*alz8 - alz2*aly8 - aly4*alz6 + alz4*aly6 - aly5*alz6 + aly6*alz5 + aly5*alz8 - alz5*aly8 - aly6*alz8 + alz6*aly8;
big_jac[19 ][ 19] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 - aly1*alz6 + alz1*aly6 + aly3*alz4 - aly4*alz3 + aly2*alz6 - alz2*aly6 + aly1*alz8 - alz1*aly8 + aly3*alz6 - alz3*aly6 - aly2*alz8 + alz2*aly8 - aly4*alz6 + alz4*aly6 - aly3*alz8 + alz3*aly8 + aly4*alz8 - alz4*aly8;
big_jac[20 ][ 19] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 + aly1*alz5 - aly2*alz4 - alz1*aly5 + alz2*aly4 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 + aly2*alz6 - aly3*alz5 - alz2*aly6 + alz3*aly5 + aly3*alz6 + aly4*alz5 - alz3*aly6 - aly5*alz4 - aly4*alz6 + alz4*aly6;
big_jac[1 ][ 20] = 0;
big_jac[2 ][ 20] = alx2*alz1 - alx1*alz2 + alx1*alz5 - alz1*alx5 - alx2*alz5 + alz2*alx5 - alx1*myco - alz1*myco + alz2*myco + alx5*myco + myco*myco;
big_jac[3 ][ 20] = alx1*alz2 - alx2*alz1 - alx1*alz3 + alx3*alz1 + alx2*alz3 - alx3*alz2 + alx1*myco - alx3*myco;
big_jac[4 ][ 20] = 0;
big_jac[5 ][ 20] = alx1*alz2 - alx2*alz1 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx5*alz6 - alx6*alz5 + alz1*myco - alz2*myco - alz5*myco + alz6*myco;
big_jac[6 ][ 20] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 - alx1*alz5 + alz1*alx5 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 + alx3*alz5 - alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 - alz2*myco + alz4*myco;
big_jac[7 ][ 20] = alx1*alz2 - alx2*alz1 - 2*alx1*alz5 + 2*alz1*alx5 + alx1*alz6 + 2*alx2*alz5 - alz1*alx6 - 2*alz2*alx5 - alx2*alz6 + alz2*alx6 + 2*alx1*myco - alx2*myco + alz1*myco - alz2*myco - 2*alx5*myco + alx6*myco - myco*myco;
big_jac[8 ][ 20] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 - alx3*alz4 + alx4*alz3;
big_jac[9 ][ 20] = alx2*alz1 - alx1*alz2 + 2*alx1*alz3 - 2*alx3*alz1 - alx1*alz4 - 2*alx2*alz3 + 2*alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - 2*alx1*myco + alx2*myco + 2*alx3*myco - alx4*myco;
big_jac[10 ][ 20] = 0;
big_jac[11 ][ 20] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alz1*alx4 + alx1*alz6 - alx2*alz5 - alz1*alx6 + alz2*alx5 - alx1*alz8 + alz1*alx8 + alx4*alz5 - alx5*alz4 - alx5*alz6 + alx6*alz5 + alx5*alz8 - alz5*alx8 + alz2*myco - alz4*myco - alz6*myco + alz8*myco;
big_jac[12 ][ 20] = alx2*alz1 - alx1*alz2 + alx1*alz6 - alx2*alz5 - alz1*alx6 + alz2*alx5 - alx5*alz6 + alx6*alz5 - alz1*myco + alz2*myco + alz5*myco - alz6*myco;
big_jac[13 ][ 20] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alx2*alz3 - alx3*alz2 + alz1*alx4 - alx1*alz6 + alx3*alz4 + alz1*alx6 - alx4*alz3 + alx1*alz8 + alx3*alz6 - alz1*alx8 - alz3*alx6 - alx3*alz8 + alz3*alx8;
big_jac[14 ][ 20] = 2*alx3*alz1 - 2*alx1*alz3 + 2*alx1*alz4 + 2*alx2*alz3 - 2*alx3*alz2 - 2*alz1*alx4 + 2*alx1*alz5 - 2*alx2*alz4 - 2*alz1*alx5 + 2*alx4*alz2 - 2*alx1*alz6 - 2*alx2*alz5 + 2*alz1*alx6 + 2*alz2*alx5 + 2*alx2*alz6 - 2*alx3*alz5 - 2*alz2*alx6 + 2*alx5*alz3 + 2*alx3*alz6 + 2*alx4*alz5 - 2*alx5*alz4 - 2*alz3*alx6 - 2*alx4*alz6 + 2*alx6*alz4;
big_jac[15 ][ 20] = alx1*alz5 - alz1*alx5 - alx1*alz6 - alx2*alz5 + alz1*alx6 + alz2*alx5 + alx2*alz6 - alz2*alx6 - alx1*myco + alx2*myco + alx5*myco - alx6*myco;
big_jac[16 ][ 20] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alx2*alz3 - alx3*alz2 + alz1*alx4 + alx3*alz4 - alx4*alz3;
big_jac[17 ][ 20] = alx3*alz1 - alx1*alz3 + alx1*alz4 + alx2*alz3 - alx3*alz2 - alz1*alx4 - alx2*alz4 + alx4*alz2 + alx1*myco - alx2*myco - alx3*myco + alx4*myco;
big_jac[18 ][ 20] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alz1*alx4 + alx2*alz4 - alx4*alz2 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx1*alz8 - alz1*alx8 - alx4*alz5 + alx5*alz4 - alx2*alz8 + alx4*alz6 + alz2*alx8 - alx6*alz4 + alx5*alz6 - alx6*alz5 - alx5*alz8 + alz5*alx8 + alx6*alz8 - alx8*alz6;
big_jac[19 ][ 20] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 + alx1*alz6 - alx3*alz4 - alz1*alx6 + alx4*alz3 - alx2*alz6 + alz2*alx6 - alx1*alz8 - alx3*alz6 + alz1*alx8 + alz3*alx6 + alx2*alz8 + alx4*alz6 - alz2*alx8 - alx6*alz4 + alx3*alz8 - alz3*alx8 - alx4*alz8 + alz4*alx8;
big_jac[20 ][ 20] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 - alx1*alz5 + alx2*alz4 + alz1*alx5 - alx4*alz2 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 - alx2*alz6 + alx3*alz5 + alz2*alx6 - alx5*alz3 - alx3*alz6 - alx4*alz5 + alx5*alz4 + alz3*alx6 + alx4*alz6 - alx6*alz4;
big_jac[1 ][ 21] = 0;
big_jac[2 ][ 21] = alx1*aly2 - alx2*aly1 - alx1*aly5 + aly1*alx5 + alx2*aly5 - aly2*alx5 + aly1*myco - aly2*myco;
big_jac[3 ][ 21] = alx2*aly1 - alx1*aly2 + alx1*aly3 - aly1*alx3 - alx2*aly3 + alx3*aly2 + alx1*myco - alx2*myco;
big_jac[4 ][ 21] = 0;
big_jac[5 ][ 21] = alx2*aly1 - alx1*aly2 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx5*aly6 + alx6*aly5 - aly1*myco + aly2*myco + aly5*myco - aly6*myco;
big_jac[6 ][ 21] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 + alx1*aly5 - aly1*alx5 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 - alx3*aly5 + aly3*alx5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 + alx2*myco + aly2*myco - aly4*myco - alx6*myco - myco*myco;
big_jac[7 ][ 21] = alx2*aly1 - alx1*aly2 + 2*alx1*aly5 - 2*aly1*alx5 - alx1*aly6 - 2*alx2*aly5 + aly1*alx6 + 2*aly2*alx5 + alx2*aly6 - aly2*alx6 - aly1*myco + aly2*myco;
big_jac[8 ][ 21] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 + alx3*aly4 - alx4*aly3 - alx1*myco + alx2*myco + alx3*myco - alx4*myco;
big_jac[9 ][ 21] = alx1*aly2 - alx2*aly1 - 2*alx1*aly3 + 2*aly1*alx3 + alx1*aly4 + 2*alx2*aly3 - aly1*alx4 - 2*alx3*aly2 - alx2*aly4 + aly2*alx4 - alx1*myco + alx2*myco;
big_jac[10 ][ 21] = 0;
big_jac[11 ][ 21] = alx1*aly2 - alx2*aly1 - alx1*aly4 + aly1*alx4 - alx1*aly6 + alx2*aly5 + aly1*alx6 - aly2*alx5 + alx1*aly8 - aly1*alx8 - alx4*aly5 + alx5*aly4 + alx5*aly6 - alx6*aly5 - alx5*aly8 + aly5*alx8 - aly2*myco + aly4*myco + aly6*myco - aly8*myco;
big_jac[12 ][ 21] = alx1*aly2 - alx2*aly1 - alx1*aly6 + alx2*aly5 + aly1*alx6 - aly2*alx5 + alx5*aly6 - alx6*aly5 + aly1*myco - aly2*myco - aly5*myco + aly6*myco;
big_jac[13 ][ 21] = alx2*aly1 - alx1*aly2 + alx1*aly4 - alx2*aly3 - aly1*alx4 + alx3*aly2 + alx1*aly6 - aly1*alx6 - alx3*aly4 + alx4*aly3 - alx1*aly8 + aly1*alx8 - alx3*aly6 + aly3*alx6 + alx3*aly8 - aly3*alx8 - alx2*myco + alx4*myco + alx6*myco - alx8*myco;
big_jac[14 ][ 21] = 2*alx1*aly3 - 2*aly1*alx3 - 2*alx1*aly4 - 2*alx2*aly3 + 2*aly1*alx4 + 2*alx3*aly2 - 2*alx1*aly5 + 2*alx2*aly4 + 2*aly1*alx5 - 2*aly2*alx4 + 2*alx1*aly6 + 2*alx2*aly5 - 2*aly1*alx6 - 2*aly2*alx5 - 2*alx2*aly6 + 2*alx3*aly5 + 2*aly2*alx6 - 2*aly3*alx5 - 2*alx3*aly6 - 2*alx4*aly5 + 2*aly3*alx6 + 2*alx5*aly4 + 2*alx4*aly6 - 2*aly4*alx6;
big_jac[15 ][ 21] = aly1*alx5 - alx1*aly5 + alx1*aly6 + alx2*aly5 - aly1*alx6 - aly2*alx5 - alx2*aly6 + aly2*alx6;
big_jac[16 ][ 21] = alx2*aly1 - alx1*aly2 + alx1*aly4 - alx2*aly3 - aly1*alx4 + alx3*aly2 - alx3*aly4 + alx4*aly3 + alx1*myco - alx2*myco - alx3*myco + alx4*myco;
big_jac[17 ][ 21] = alx1*aly3 - aly1*alx3 - alx1*aly4 - alx2*aly3 + aly1*alx4 + alx3*aly2 + alx2*aly4 - aly2*alx4;
big_jac[18 ][ 21] = alx2*aly1 - alx1*aly2 + alx1*aly4 - aly1*alx4 - alx2*aly4 + aly2*alx4 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx1*aly8 + aly1*alx8 + alx4*aly5 - alx5*aly4 + alx2*aly8 - aly2*alx8 - alx4*aly6 + aly4*alx6 - alx5*aly6 + alx6*aly5 + alx5*aly8 - aly5*alx8 - alx6*aly8 + aly6*alx8;
big_jac[19 ][ 21] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 - alx1*aly6 + aly1*alx6 + alx3*aly4 - alx4*aly3 + alx2*aly6 - aly2*alx6 + alx1*aly8 - aly1*alx8 + alx3*aly6 - aly3*alx6 - alx2*aly8 + aly2*alx8 - alx4*aly6 + aly4*alx6 - alx3*aly8 + aly3*alx8 + alx4*aly8 - aly4*alx8;
big_jac[20 ][ 21] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 + alx1*aly5 - alx2*aly4 - aly1*alx5 + aly2*alx4 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 + alx2*aly6 - alx3*aly5 - aly2*alx6 + aly3*alx5 + alx3*aly6 + alx4*aly5 - aly3*alx6 - alx5*aly4 - alx4*aly6 + aly4*alx6;

break;
case 8:

big_jac[1 ][ 22] = 0;
big_jac[2 ][ 22] = 0;
big_jac[3 ][ 22] = 0;
big_jac[4 ][ 22] = 0;
big_jac[5 ][ 22] = 0;
big_jac[6 ][ 22] = alz1*aly3 - aly1*alz3 + aly1*alz5 - alz1*aly5 - aly3*alz5 + alz3*aly5 + alz1*myco - alz5*myco;
big_jac[7 ][ 22] = aly1*alz2 - aly2*alz1 - aly1*alz5 + alz1*aly5 + aly2*alz5 - alz2*aly5 + aly1*myco - aly5*myco;
big_jac[8 ][ 22] = 0;
big_jac[9 ][ 22] = aly2*alz1 - aly1*alz2 + aly1*alz3 - alz1*aly3 - aly2*alz3 + aly3*alz2 - aly1*myco - alz1*myco + aly3*myco + alz2*myco + myco*myco;
big_jac[10 ][ 22] = 0;
big_jac[11 ][ 22] = aly1*alz3 - alz1*aly3 - aly1*alz7 + alz1*aly7 + aly3*alz5 - alz3*aly5 + aly5*alz7 - alz5*aly7;
big_jac[12 ][ 22] = aly2*alz1 - aly1*alz2 + aly1*alz6 - aly2*alz5 - alz1*aly6 + alz2*aly5 - aly5*alz6 + aly6*alz5;
big_jac[13 ][ 22] = alz1*aly5 - aly1*alz5 + aly1*alz7 - alz1*aly7 + aly3*alz5 - alz3*aly5 - aly3*alz7 + alz3*aly7 - alz1*myco + alz3*myco + alz5*myco - alz7*myco;
big_jac[14 ][ 22] = 0;
big_jac[15 ][ 22] = aly1*alz5 - alz1*aly5 - aly1*alz6 - aly2*alz5 + alz1*aly6 + alz2*aly5 + aly2*alz6 - alz2*aly6 - aly1*myco + aly2*myco + aly5*myco - aly6*myco;
big_jac[16 ][ 22] = aly1*alz2 - aly2*alz1 - aly1*alz4 + aly2*alz3 + alz1*aly4 - aly3*alz2 + aly3*alz4 - aly4*alz3 + alz1*myco - alz2*myco - alz3*myco + alz4*myco;
big_jac[17 ][ 22] = alz1*aly3 - aly1*alz3 + aly1*alz4 + aly2*alz3 - alz1*aly4 - aly3*alz2 - aly2*alz4 + alz2*aly4 + aly1*myco - aly2*myco - aly3*myco + aly4*myco;
big_jac[18 ][ 22] = aly1*alz2 - aly2*alz1 - aly1*alz3 + alz1*aly3 + aly2*alz3 - aly3*alz2 - aly1*alz6 + aly2*alz5 + alz1*aly6 - alz2*aly5 + aly1*alz7 - alz1*aly7 - aly3*alz5 + alz3*aly5 - aly2*alz7 + aly3*alz6 + alz2*aly7 - alz3*aly6 + aly5*alz6 - aly6*alz5 - aly5*alz7 + alz5*aly7 + aly6*alz7 - aly7*alz6;
big_jac[19 ][ 22] = aly2*alz1 - aly1*alz2 + aly1*alz4 - aly2*alz3 - alz1*aly4 + aly3*alz2 + aly1*alz5 - alz1*aly5 - aly2*alz5 - aly3*alz4 + alz2*aly5 + aly4*alz3 - aly1*alz7 + alz1*aly7 - aly3*alz5 + alz3*aly5 + aly2*alz7 - alz2*aly7 + aly4*alz5 - aly5*alz4 + aly3*alz7 - alz3*aly7 - aly4*alz7 + alz4*aly7;
big_jac[20 ][ 22] = aly1*alz3 - alz1*aly3 - aly1*alz4 - aly2*alz3 + alz1*aly4 + aly3*alz2 - aly1*alz5 + aly2*alz4 + alz1*aly5 - alz2*aly4 + aly1*alz6 + aly2*alz5 - alz1*aly6 - alz2*aly5 - aly2*alz6 + aly3*alz5 + alz2*aly6 - alz3*aly5 - aly3*alz6 - aly4*alz5 + alz3*aly6 + aly5*alz4 + aly4*alz6 - alz4*aly6;
big_jac[1 ][ 23] = 0;
big_jac[2 ][ 23] = 0;
big_jac[3 ][ 23] = 0;
big_jac[4 ][ 23] = 0;
big_jac[5 ][ 23] = 0;
big_jac[6 ][ 23] = alx1*alz3 - alx3*alz1 - alx1*alz5 + alz1*alx5 + alx3*alz5 - alx5*alz3 + alz1*myco - alz3*myco;
big_jac[7 ][ 23] = alx2*alz1 - alx1*alz2 + alx1*alz5 - alz1*alx5 - alx2*alz5 + alz2*alx5 - alx1*myco - alz1*myco + alz2*myco + alx5*myco + myco*myco;
big_jac[8 ][ 23] = 0;
big_jac[9 ][ 23] = alx1*alz2 - alx2*alz1 - alx1*alz3 + alx3*alz1 + alx2*alz3 - alx3*alz2 + alx1*myco - alx3*myco;
big_jac[10 ][ 23] = 0;
big_jac[11 ][ 23] = alx3*alz1 - alx1*alz3 + alx1*alz7 - alx3*alz5 - alz1*alx7 + alx5*alz3 - alx5*alz7 + alx7*alz5 - alz1*myco + alz3*myco + alz5*myco - alz7*myco;
big_jac[12 ][ 23] = alx1*alz2 - alx2*alz1 - alx1*alz6 + alx2*alz5 + alz1*alx6 - alz2*alx5 + alx5*alz6 - alx6*alz5 + alz1*myco - alz2*myco - alz5*myco + alz6*myco;
big_jac[13 ][ 23] = alx1*alz5 - alz1*alx5 - alx1*alz7 - alx3*alz5 + alz1*alx7 + alx5*alz3 + alx3*alz7 - alz3*alx7;
big_jac[14 ][ 23] = 0;
big_jac[15 ][ 23] = alz1*alx5 - alx1*alz5 + alx1*alz6 + alx2*alz5 - alz1*alx6 - alz2*alx5 - alx2*alz6 + alz2*alx6 + alx1*myco - alx2*myco - alx5*myco + alx6*myco;
big_jac[16 ][ 23] = alx2*alz1 - alx1*alz2 + alx1*alz4 - alx2*alz3 + alx3*alz2 - alz1*alx4 - alx3*alz4 + alx4*alz3;
big_jac[17 ][ 23] = alx1*alz3 - alx3*alz1 - alx1*alz4 - alx2*alz3 + alx3*alz2 + alz1*alx4 + alx2*alz4 - alx4*alz2 - alx1*myco + alx2*myco + alx3*myco - alx4*myco;
big_jac[18 ][ 23] = alx2*alz1 - alx1*alz2 + alx1*alz3 - alx3*alz1 - alx2*alz3 + alx3*alz2 + alx1*alz6 - alx2*alz5 - alz1*alx6 + alz2*alx5 - alx1*alz7 + alx3*alz5 + alz1*alx7 - alx5*alz3 + alx2*alz7 - alx3*alz6 - alz2*alx7 + alz3*alx6 - alx5*alz6 + alx6*alz5 + alx5*alz7 - alx7*alz5 - alx6*alz7 + alx7*alz6;
big_jac[19 ][ 23] = alx1*alz2 - alx2*alz1 - alx1*alz4 + alx2*alz3 - alx3*alz2 + alz1*alx4 - alx1*alz5 + alz1*alx5 + alx2*alz5 + alx3*alz4 - alx4*alz3 - alz2*alx5 + alx1*alz7 + alx3*alz5 - alz1*alx7 - alx5*alz3 - alx2*alz7 - alx4*alz5 + alz2*alx7 + alx5*alz4 - alx3*alz7 + alz3*alx7 + alx4*alz7 - alz4*alx7;
big_jac[20 ][ 23] = alx3*alz1 - alx1*alz3 + alx1*alz4 + alx2*alz3 - alx3*alz2 - alz1*alx4 + alx1*alz5 - alx2*alz4 - alz1*alx5 + alx4*alz2 - alx1*alz6 - alx2*alz5 + alz1*alx6 + alz2*alx5 + alx2*alz6 - alx3*alz5 - alz2*alx6 + alx5*alz3 + alx3*alz6 + alx4*alz5 - alx5*alz4 - alz3*alx6 - alx4*alz6 + alx6*alz4;
big_jac[1 ][ 24] = 0;
big_jac[2 ][ 24] = 0;
big_jac[3 ][ 24] = 0;
big_jac[4 ][ 24] = 0;
big_jac[5 ][ 24] = 0;
big_jac[6 ][ 24] = aly1*alx3 - alx1*aly3 + alx1*aly5 - aly1*alx5 - alx3*aly5 + aly3*alx5 - alx1*myco - aly1*myco + aly3*myco + alx5*myco + myco*myco;
big_jac[7 ][ 24] = alx1*aly2 - alx2*aly1 - alx1*aly5 + aly1*alx5 + alx2*aly5 - aly2*alx5 + aly1*myco - aly2*myco;
big_jac[8 ][ 24] = 0;
big_jac[9 ][ 24] = alx2*aly1 - alx1*aly2 + alx1*aly3 - aly1*alx3 - alx2*aly3 + alx3*aly2 + alx1*myco - alx2*myco;
big_jac[10 ][ 24] = 0;
big_jac[11 ][ 24] = alx1*aly3 - aly1*alx3 - alx1*aly7 + aly1*alx7 + alx3*aly5 - aly3*alx5 + alx5*aly7 - aly5*alx7 + aly1*myco - aly3*myco - aly5*myco + aly7*myco;
big_jac[12 ][ 24] = alx2*aly1 - alx1*aly2 + alx1*aly6 - alx2*aly5 - aly1*alx6 + aly2*alx5 - alx5*aly6 + alx6*aly5 - aly1*myco + aly2*myco + aly5*myco - aly6*myco;
big_jac[13 ][ 24] = aly1*alx5 - alx1*aly5 + alx1*aly7 - aly1*alx7 + alx3*aly5 - aly3*alx5 - alx3*aly7 + aly3*alx7 + alx1*myco - alx3*myco - alx5*myco + alx7*myco;
big_jac[14 ][ 24] = 0;
big_jac[15 ][ 24] = alx1*aly5 - aly1*alx5 - alx1*aly6 - alx2*aly5 + aly1*alx6 + aly2*alx5 + alx2*aly6 - aly2*alx6;
big_jac[16 ][ 24] = alx1*aly2 - alx2*aly1 - alx1*aly4 + alx2*aly3 + aly1*alx4 - alx3*aly2 + alx3*aly4 - alx4*aly3 - alx1*myco + alx2*myco + alx3*myco - alx4*myco;
big_jac[17 ][ 24] = aly1*alx3 - alx1*aly3 + alx1*aly4 + alx2*aly3 - aly1*alx4 - alx3*aly2 - alx2*aly4 + aly2*alx4;
big_jac[18 ][ 24] = alx1*aly2 - alx2*aly1 - alx1*aly3 + aly1*alx3 + alx2*aly3 - alx3*aly2 - alx1*aly6 + alx2*aly5 + aly1*alx6 - aly2*alx5 + alx1*aly7 - aly1*alx7 - alx3*aly5 + aly3*alx5 - alx2*aly7 + alx3*aly6 + aly2*alx7 - aly3*alx6 + alx5*aly6 - alx6*aly5 - alx5*aly7 + aly5*alx7 + alx6*aly7 - alx7*aly6;
big_jac[19 ][ 24] = alx2*aly1 - alx1*aly2 + alx1*aly4 - alx2*aly3 - aly1*alx4 + alx3*aly2 + alx1*aly5 - aly1*alx5 - alx2*aly5 - alx3*aly4 + aly2*alx5 + alx4*aly3 - alx1*aly7 + aly1*alx7 - alx3*aly5 + aly3*alx5 + alx2*aly7 - aly2*alx7 + alx4*aly5 - alx5*aly4 + alx3*aly7 - aly3*alx7 - alx4*aly7 + aly4*alx7;
big_jac[20 ][ 24] = alx1*aly3 - aly1*alx3 - alx1*aly4 - alx2*aly3 + aly1*alx4 + alx3*aly2 - alx1*aly5 + alx2*aly4 + aly1*alx5 - aly2*alx4 + alx1*aly6 + alx2*aly5 - aly1*alx6 - aly2*alx5 - alx2*aly6 + alx3*aly5 + aly2*alx6 - aly3*alx5 - alx3*aly6 - alx4*aly5 + aly3*alx6 + alx5*aly4 + alx4*aly6 - aly4*alx6;

}												// fin switch



/*----------------------------------------------
  code matlab :

big_jac = big_jac * myco^(-3);

myP = big_jac(:,1+3*(indice-1):3+3*(indice-1));
pente = myP * dir_des;
------------------------------------------------*/

for (i=0; i<20; i++)
{ pente[i] = 0;
  for (j=0; j<3; j++)
    pente[i] += big_jac[i+1][1+3*(indice-1)+j] * pow(myco,-3) * dir_des[j];
}

}





/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*%-----------------------------------------------------------
%----- fichier : poly_init.m
%----- auteur  : C. Heinrich
%----- date    : juillet 2002
%----- objet   : initialisation du polynome
%-----------------------------------------------------------


function [alpha] = poly_init(alx,aly,alz,l);*/

void TOP_poly_init(double *alx, double *aly, double *alz, int l, double *alpha)
{
double alx1,alx2,alx3,alx4,alx5,alx6,alx7,alx8;
double aly1,aly2,aly3,aly4,aly5,aly6,aly7,aly8;
double alz1,alz2,alz3,alz4,alz5,alz6,alz7,alz8;




   //   myco = 2^(-5*l/2);
   // double myco = pow(2,-5*l/2);
/* Convention O.Musse : facteur de la fonction d'echelle = 1 */
//l=0;

  // double myco = 1.0;
   double myco=pow(2.0,-1.0*l);
   int i;


alx1 = alx[0];
alx2 = alx[1];
alx3 = alx[2];
alx4 = alx[3];
alx5 = alx[4];
alx6 = alx[5];
alx7 = alx[6];
alx8 = alx[7];

aly1 = aly[0];
aly2 = aly[1];
aly3 = aly[2];
aly4 = aly[3];
aly5 = aly[4];
aly6 = aly[5];
aly7 = aly[6];
aly8 = aly[7];

alz1 = alz[0];
alz2 = alz[1];
alz3 = alz[2];
alz4 = alz[3];
alz5 = alz[4];
alz6 = alz[5];
alz7 = alz[6];
alz8 = alz[7];


alpha[0] = alx1*aly2*alz3 - alx1*aly3*alz2 - alx2*aly1*alz3 + alx2*alz1*aly3 + aly1*alx3*alz2 - alx3*aly2*alz1 \
- alx1*aly2*alz5 + alx1*alz2*aly5 + alx2*aly1*alz5 - alx2*alz1*aly5 - aly1*alz2*alx5 + aly2*alz1*alx5 + alx1*aly3*alz5 \
- alx1*alz3*aly5 - aly1*alx3*alz5 + aly1*alx5*alz3 + alx3*alz1*aly5 - alz1*aly3*alx5 - alx2*aly3*alz5 + alx2*alz3*aly5 \
+ alx3*aly2*alz5 - alx3*alz2*aly5 - aly2*alx5*alz3 + aly3*alz2*alx5 - alx1*aly3*myco - alx1*alz2*myco + alx2*alz1*myco \
+ aly1*alx3*myco - aly1*alz2*myco + aly2*alz1*myco + alx1*aly5*myco - aly1*alx5*myco + aly1*alz3*myco - alz1*aly3*myco \
+ alx1*alz5*myco - aly2*alz3*myco - alz1*alx5*myco + aly3*alz2*myco - alx2*alz5*myco - alx3*aly5*myco + aly3*alx5*myco \
+ alz2*alx5*myco + myco*myco*myco - alx1*myco*myco - aly1*myco*myco - alz1*myco*myco + aly3*myco*myco + alz2*myco*myco + alx5*myco*myco;

alpha[1] = 2*alx1*aly3*alz2 - 2*alx1*aly2*alz3 + 2*alx2*aly1*alz3 - 2*alx2*alz1*aly3 - 2*aly1*alx3*alz2 + 2*alx3*aly2*alz1 \
+ alx1*aly2*alz5 - alx1*alz2*aly5 - alx2*aly1*alz5 + alx2*alz1*aly5 + aly1*alz2*alx5 - aly2*alz1*alx5 - alx1*aly3*alz5 \
+ alx1*alz3*aly5 + aly1*alx3*alz5 - aly1*alx5*alz3 - alx3*alz1*aly5 + alz1*aly3*alx5 + alx1*aly2*alz7 - alx1*aly3*alz6 \
- alx1*alz2*aly7 + alx1*alz3*aly6 - alx2*aly1*alz7 + alx2*alz1*aly7 + 2*alx2*aly3*alz5 - 2*alx2*alz3*aly5 + aly1*alx3*alz6 \
 + aly1*alz2*alx7 - aly1*alz3*alx6 - 2*alx3*aly2*alz5 - alx3*alz1*aly6 + 2*alx3*alz2*aly5 - aly2*alz1*alx7 + \
 2*aly2*alx5*alz3 + alz1*aly3*alx6 - 2*aly3*alz2*alx5 + alx1*aly5*alz6 - alx1*aly6*alz5 - aly1*alx5*alz6 + aly1*alx6*alz5 \
  + alz1*alx5*aly6 - alz1*alx6*aly5 - alx1*aly5*alz7 + alx1*alz5*aly7 + aly1*alx5*alz7 - aly1*alx7*alz5 - alz1*alx5*aly7 \
   + alz1*aly5*alx7 + alx2*aly5*alz7 - alx2*alz5*aly7 - alx3*aly5*alz6 + alx3*aly6*alz5 - aly2*alx5*alz7 + aly2*alx7*alz5 \
    + aly3*alx5*alz6 - aly3*alx6*alz5 + alz2*alx5*aly7 - alz2*aly5*alx7 - alx5*alz3*aly6 + alz3*alx6*aly5 + alx1*aly3*myco \
     + alx1*alz2*myco - alx2*alz1*myco - aly1*alx3*myco + 2*aly1*alz2*myco - 2*aly2*alz1*myco - 2*aly1*alz3*myco \
      + 2*alz1*aly3*myco + 2*aly2*alz3*myco - 2*aly3*alz2*myco - alx1*aly7*myco - alx1*alz6*myco + alx2*alz5*myco \
      + aly1*alx7*myco + alx3*aly5*myco + alz1*alx6*myco - aly3*alx5*myco - alz2*alx5*myco - aly1*alz6*myco \
      + aly2*alz5*myco + alz1*aly6*myco - alz2*aly5*myco + aly1*alz7*myco - alz1*aly7*myco - aly3*alz5*myco \
      + alz3*aly5*myco - aly2*alz7*myco + aly3*alz6*myco + alz2*aly7*myco - alz3*aly6*myco + alx5*aly7*myco \
      + alx5*alz6*myco - alx6*alz5*myco - aly5*alx7*myco + aly1*myco*myco + alz1*myco*myco - aly3*myco*myco - alz2*myco*myco \
      - aly5*myco*myco - alz5*myco*myco + aly7*myco*myco + alz6*myco*myco;

alpha[2] = alx1*aly3*alz2 - alx1*aly2*alz3 + alx2*aly1*alz3 - alx2*alz1*aly3 - aly1*alx3*alz2 + alx3*aly2*alz1 \
  + 2*alx1*aly2*alz5 - alx1*aly3*alz4 - 2*alx1*alz2*aly5 + alx1*aly4*alz3 - 2*alx2*aly1*alz5 + 2*alx2*alz1*aly5 \
  + aly1*alx3*alz4 - aly1*alx4*alz3 + 2*aly1*alz2*alx5 - alx3*alz1*aly4 - 2*aly2*alz1*alx5 + alz1*alx4*aly3 \
  - alx1*aly3*alz5 + alx1*alz3*aly5 + aly1*alx3*alz5 - aly1*alx5*alz3 - alx3*alz1*aly5 + alz1*aly3*alx5 \
  - alx1*aly2*alz7 + alx1*alz2*aly7 - alx1*aly4*alz5 + alx1*aly5*alz4 + alx2*aly1*alz7 - alx2*alz1*aly7 \
  + 2*alx2*aly3*alz5 - 2*alx2*alz3*aly5 + aly1*alx4*alz5 - aly1*alz2*alx7 - aly1*alx5*alz4 - 2*alx3*aly2*alz5 \
  + 2*alx3*alz2*aly5 + aly2*alz1*alx7 + 2*aly2*alx5*alz3 - alz1*alx4*aly5 + alz1*alx5*aly4 - 2*aly3*alz2*alx5 \
  + alx1*aly3*alz7 - alx1*alz3*aly7 - aly1*alx3*alz7 + aly1*alz3*alx7 + alx3*alz1*aly7 - alz1*aly3*alx7 \
  - alx2*aly3*alz7 + alx2*alz3*aly7 + alx3*aly2*alz7 - alx3*alz2*aly7 + alx3*aly4*alz5 - alx3*aly5*alz4 - aly2*alz3*alx7 \
  - alx4*aly3*alz5 + alx4*alz3*aly5 + aly3*alz2*alx7 + aly3*alx5*alz4 - alx5*aly4*alz3 + 2*alx1*alz2*myco \
  - 2*alx2*alz1*myco + aly1*alz2*myco - aly2*alz1*myco - alx1*aly5*myco - alx1*alz4*myco + alx2*alz3*myco \
  + aly1*alx5*myco - alx3*alz2*myco + alz1*alx4*myco - 2*alx1*alz5*myco - aly1*alz4*myco + aly2*alz3*myco \
  + 2*alz1*alx5*myco + alz1*aly4*myco - aly3*alz2*myco + alx1*aly7*myco + 2*alx2*alz5*myco - aly1*alx7*myco \
  + alx3*aly5*myco - aly3*alx5*myco - 2*alz2*alx5*myco + alx1*alz7*myco + alx3*alz5*myco - alz1*alx7*myco \
  + aly3*alz4*myco - alx5*alz3*myco - aly4*alz3*myco - alx2*alz7*myco - alx3*aly7*myco - alx4*alz5*myco + aly3*alx7*myco \
   + alz2*alx7*myco + alx5*alz4*myco + alx1*myco*myco - alx3*myco*myco + alz1*myco*myco - alz2*myco*myco - alx5*myco*myco \
   - alz3*myco*myco + alz4*myco*myco + alx7*myco*myco;

alpha[3] = alx1*aly3*alz2 - alx1*aly2*alz3 + alx2*aly1*alz3 - alx2*alz1*aly3 - aly1*alx3*alz2 + alx3*aly2*alz1 \
+ alx1*aly2*alz4 - alx1*alz2*aly4 - alx2*aly1*alz4 + alx2*alz1*aly4 + aly1*alx4*alz2 - aly2*alz1*alx4 + alx1*aly2*alz5 \
- alx1*alz2*aly5 - alx2*aly1*alz5 + alx2*alz1*aly5 + aly1*alz2*alx5 - aly2*alz1*alx5 - alx1*aly2*alz6 - 2*alx1*aly3*alz5 \
+ alx1*alz2*aly6 + 2*alx1*alz3*aly5 + alx2*aly1*alz6 - alx2*alz1*aly6 + 2*aly1*alx3*alz5 - aly1*alz2*alx6 \
- 2*aly1*alx5*alz3 - 2*alx3*alz1*aly5 + aly2*alz1*alx6 + 2*alz1*aly3*alx5 + alx1*aly3*alz6 + alx1*aly4*alz5 \
- alx1*alz3*aly6 - alx1*aly5*alz4 + 2*alx2*aly3*alz5 - 2*alx2*alz3*aly5 - aly1*alx3*alz6 - aly1*alx4*alz5 + aly1*alx5*alz4 \
 + aly1*alz3*alx6 - 2*alx3*aly2*alz5 + alx3*alz1*aly6 + 2*alx3*alz2*aly5 + 2*aly2*alx5*alz3 + alz1*alx4*aly5 \
 - alz1*aly3*alx6 - alz1*alx5*aly4 - 2*aly3*alz2*alx5 - alx2*aly3*alz6 - alx2*aly4*alz5 + alx2*alz3*aly6 + alx2*aly5*alz4 \
 + alx3*aly2*alz6 - alx3*alz2*aly6 + aly2*alx4*alz5 - aly2*alx5*alz4 - aly2*alz3*alx6 - alx4*alz2*aly5 + aly3*alz2*alx6 \
 + alz2*alx5*aly4 + 2*alx1*aly3*myco - 2*aly1*alx3*myco - alx1*aly4*myco - alx2*aly3*myco + aly1*alx4*myco \
 + alx3*aly2*myco - 2*alx1*aly5*myco + 2*aly1*alx5*myco - aly1*alz3*myco + alz1*aly3*myco + alx1*aly6*myco \
 - alx1*alz5*myco + alx2*aly5*myco - aly1*alx6*myco + aly1*alz4*myco - aly2*alx5*myco + aly2*alz3*myco + alz1*alx5*myco \
  - alz1*aly4*myco - aly3*alz2*myco + alx1*alz6*myco + alx2*alz5*myco + 2*alx3*aly5*myco - aly2*alz4*myco \
  - alz1*alx6*myco - 2*aly3*alx5*myco - alz2*alx5*myco + alz2*aly4*myco - alx2*alz6*myco - alx3*aly6*myco \
   - alx4*aly5*myco + aly3*alx6*myco + alz2*alx6*myco + alx5*aly4*myco + alx1*myco*myco - alx2*myco*myco + aly1*myco*myco \
   - aly2*myco*myco - aly3*myco*myco - alx5*myco*myco + aly4*myco*myco + alx6*myco*myco;

alpha[4] = alx1*aly2*alz3 - alx1*aly3*alz2 - alx2*aly1*alz3 + alx2*alz1*aly3 + aly1*alx3*alz2 - alx3*aly2*alz1 \
 - alx1*aly2*alz7 + alx1*aly3*alz6 + alx1*alz2*aly7 - alx1*alz3*aly6 + alx2*aly1*alz7 - alx2*alz1*aly7 - alx2*aly3*alz5 \
  + alx2*alz3*aly5 - aly1*alx3*alz6 - aly1*alz2*alx7 + aly1*alz3*alx6 + alx3*aly2*alz5 + alx3*alz1*aly6 - alx3*alz2*aly5 \
  + aly2*alz1*alx7 - aly2*alx5*alz3 - alz1*aly3*alx6 + aly3*alz2*alx5 + alx1*aly6*alz7 - alx1*aly7*alz6 - alx2*aly5*alz7 \
   + alx2*alz5*aly7 - aly1*alx6*alz7 + aly1*alx7*alz6 + alx3*aly5*alz6 - alx3*aly6*alz5 + aly2*alx5*alz7 - aly2*alx7*alz5 \
    + alz1*alx6*aly7 - alz1*alx7*aly6 - aly3*alx5*alz6 + aly3*alx6*alz5 - alz2*alx5*aly7 + alz2*aly5*alx7 + alx5*alz3*aly6 \
     - alz3*alx6*aly5 - alx5*aly6*alz7 + alx5*aly7*alz6 + alx6*aly5*alz7 - alx6*alz5*aly7 - aly5*alx7*alz6 + alx7*aly6*alz5 \
      - aly1*alz2*myco + aly2*alz1*myco + aly1*alz3*myco - alz1*aly3*myco - aly2*alz3*myco + aly3*alz2*myco \
      + aly1*alz6*myco - aly2*alz5*myco - alz1*aly6*myco + alz2*aly5*myco - aly1*alz7*myco + alz1*aly7*myco \
      + aly3*alz5*myco - alz3*aly5*myco + aly2*alz7*myco - aly3*alz6*myco - alz2*aly7*myco + alz3*aly6*myco \
      - aly5*alz6*myco + aly6*alz5*myco + aly5*alz7*myco - alz5*aly7*myco - aly6*alz7*myco + aly7*alz6*myco;

alpha[5] = 2*alx1*aly2*alz3 - 2*alx1*aly3*alz2 - 2*alx2*aly1*alz3 + 2*alx2*alz1*aly3 + 2*aly1*alx3*alz2 - 2*alx3*aly2*alz1 \
 - 2*alx1*aly2*alz5 + 2*alx1*aly3*alz4 + 2*alx1*alz2*aly5 - 2*alx1*aly4*alz3 + 2*alx2*aly1*alz5 - 2*alx2*alz1*aly5 \
 - 2*aly1*alx3*alz4 + 2*aly1*alx4*alz3 - 2*aly1*alz2*alx5 + 2*alx3*alz1*aly4 + 2*aly2*alz1*alx5 - 2*alz1*alx4*aly3 \
 + alx1*aly3*alz5 - alx1*alz3*aly5 - aly1*alx3*alz5 + aly1*alx5*alz3 + alx3*alz1*aly5 - alz1*aly3*alx5 + alx1*aly3*alz6 \
 + alx1*aly4*alz5 - alx1*alz3*aly6 - alx1*aly5*alz4 - 3*alx2*aly3*alz5 + 3*alx2*alz3*aly5 - aly1*alx3*alz6 - aly1*alx4*alz5 \
  + aly1*alx5*alz4 + aly1*alz3*alx6 + 3*alx3*aly2*alz5 + alx3*alz1*aly6 - 3*alx3*alz2*aly5 - 3*aly2*alx5*alz3 \
  + alz1*alx4*aly5 - alz1*aly3*alx6 - alz1*alx5*aly4 + 3*aly3*alz2*alx5 - alx1*aly3*alz7 + alx1*alz3*aly7 \
  + aly1*alx3*alz7 - aly1*alz3*alx7 - alx3*alz1*aly7 + alz1*aly3*alx7 - alx1*aly3*alz8 + alx1*aly4*alz7 \
  + alx1*alz3*aly8 - 2*alx1*aly5*alz6 - alx1*alz4*aly7 + 2*alx1*aly6*alz5 + alx2*aly3*alz7 - alx2*alz3*aly7 \
  + aly1*alx3*alz8 - aly1*alx4*alz7 + 2*aly1*alx5*alz6 - aly1*alz3*alx8 - 2*aly1*alx6*alz5 + aly1*alz4*alx7 \
  - alx3*aly2*alz7 - alx3*alz1*aly8 + alx3*alz2*aly7 - 2*alx3*aly4*alz5 + 2*alx3*aly5*alz4 + aly2*alz3*alx7 \
  + alz1*alx4*aly7 + alz1*aly3*alx8 - 2*alz1*alx5*aly6 - alz1*aly4*alx7 + 2*alz1*alx6*aly5 + 2*alx4*aly3*alz5 \
  - 2*alx4*alz3*aly5 - aly3*alz2*alx7 - 2*aly3*alx5*alz4 + 2*alx5*aly4*alz3 + alx1*aly5*alz7 - alx1*alz5*aly7 \
  - aly1*alx5*alz7 + aly1*alx7*alz5 + alz1*alx5*aly7 - alz1*aly5*alx7 + alx1*aly5*alz8 - alx1*aly6*alz7 \
  - alx1*alz5*aly8 + alx1*aly7*alz6 - alx2*aly5*alz7 + alx2*alz5*aly7 - aly1*alx5*alz8 + aly1*alx6*alz7 - aly1*alx7*alz6 \
   + aly1*alz5*alx8 + 2*alx3*aly5*alz6 - 2*alx3*aly6*alz5 + aly2*alx5*alz7 - aly2*alx7*alz5 + alz1*alx5*aly8 \
   - alz1*alx6*aly7 - alz1*aly5*alx8 + alz1*alx7*aly6 - 2*aly3*alx5*alz6 + 2*aly3*alx6*alz5 - alz2*alx5*aly7 \
   + alz2*aly5*alx7 + 2*alx5*alz3*aly6 - 2*alz3*alx6*aly5 - alx3*aly5*alz7 + alx3*alz5*aly7 + aly3*alx5*alz7 \
   - aly3*alx7*alz5 - alx5*alz3*aly7 + alz3*aly5*alx7 - alx3*aly5*alz8 + alx3*aly6*alz7 \
+ alx3*alz5*aly8 - alx3*aly7*alz6 + alx4*aly5*alz7 - alx4*alz5*aly7 + aly3*alx5*alz8 - aly3*alx6*alz7 + aly3*alx7*alz6 \
- aly3*alz5*alx8 - alx5*aly4*alz7 - alx5*alz3*aly8 + alx5*alz4*aly7 + aly4*alx7*alz5 + alz3*alx6*aly7 + alz3*aly5*alx8 \
- alz3*alx7*aly6 - aly5*alz4*alx7 - 2*alx1*alz2*myco + 2*alx2*alz1*myco - 2*aly1*alz2*myco + 2*aly2*alz1*myco \
+ alx1*alz4*myco - alx2*alz3*myco + alx3*alz2*myco - alz1*alx4*myco + 2*aly1*alz4*myco - 2*aly2*alz3*myco \
- 2*alz1*aly4*myco + 2*aly3*alz2*myco + 2*alx1*alz6*myco - 2*alx2*alz5*myco - 2*alz1*alx6*myco + 2*alz2*alx5*myco \
+ aly1*alz6*myco - aly2*alz5*myco - alz1*aly6*myco - 2*aly3*alz4*myco + alz2*aly5*myco + 2*aly4*alz3*myco \
- alx1*alz8*myco + alx2*alz7*myco - alx3*alz6*myco + alz1*alx8*myco + alx4*alz5*myco - alz2*alx7*myco - alx5*alz4*myco \
+ alz3*alx6*myco - aly1*alz8*myco + aly2*alz7*myco + alz1*aly8*myco - aly3*alz6*myco - alz2*aly7*myco + aly4*alz5*myco \
+ alz3*aly6*myco - aly5*alz4*myco - 2*alx5*alz6*myco + 2*alx6*alz5*myco + aly3*alz8*myco - aly4*alz7*myco \
- alz3*aly8*myco + alz4*aly7*myco + alx5*alz8*myco - alx6*alz7*myco + alx7*alz6*myco - alz5*alx8*myco - alz1*myco*myco \
+ alz2*myco*myco + alz3*myco*myco - alz4*myco*myco + alz5*myco*myco - alz6*myco*myco - alz7*myco*myco + alz8*myco*myco;

alpha[6] = 2*alx1*aly2*alz3 - 2*alx1*aly3*alz2 - 2*alx2*aly1*alz3 + 2*alx2*alz1*aly3 + 2*aly1*alx3*alz2 - 2*alx3*aly2*alz1\
 - 2*alx1*aly2*alz4 + 2*alx1*alz2*aly4 + 2*alx2*aly1*alz4 - 2*alx2*alz1*aly4 - 2*aly1*alx4*alz2 + 2*aly2*alz1*alx4 \
 - alx1*aly2*alz5 + alx1*alz2*aly5 + alx2*aly1*alz5 - alx2*alz1*aly5 - aly1*alz2*alx5 + aly2*alz1*alx5 + alx1*aly2*alz6 \
 + 2*alx1*aly3*alz5 - alx1*alz2*aly6 - 2*alx1*alz3*aly5 - alx2*aly1*alz6 + alx2*alz1*aly6 - 2*aly1*alx3*alz5 \
 + aly1*alz2*alx6 + 2*aly1*alx5*alz3 + 2*alx3*alz1*aly5 - aly2*alz1*alx6 - 2*alz1*aly3*alx5 - alx1*aly2*alz7 \
 + alx1*alz2*aly7 - alx1*aly4*alz5 + alx1*aly5*alz4 + alx2*aly1*alz7 - alx2*alz1*aly7 - 3*alx2*aly3*alz5 \
 + 3*alx2*alz3*aly5 + aly1*alx4*alz5 - aly1*alz2*alx7 - aly1*alx5*alz4 + 3*alx3*aly2*alz5 - 3*alx3*alz2*aly5 \
 + aly2*alz1*alx7 - 3*aly2*alx5*alz3 - alz1*alx4*aly5 + alz1*alx5*aly4 + 3*aly3*alz2*alx5 + alx1*aly2*alz8 \
 - alx1*alz2*aly8 - alx1*aly4*alz6 + alx1*alz4*aly6 - alx2*aly1*alz8 + alx2*alz1*aly8 + alx2*aly3*alz6 + 2*alx2*aly4*alz5 \
 - alx2*alz3*aly6 - 2*alx2*aly5*alz4 + aly1*alx4*alz6 + aly1*alz2*alx8 - aly1*alx6*alz4 - alx3*aly2*alz6 \
 + alx3*alz2*aly6 - aly2*alz1*alx8 - 2*aly2*alx4*alz5 + 2*aly2*alx5*alz4 + aly2*alz3*alx6 - alz1*alx4*aly6 \
 + alz1*aly4*alx6 + 2*alx4*alz2*aly5 - aly3*alz2*alx6 - 2*alz2*alx5*aly4 - alx1*aly5*alz6 + alx1*aly6*alz5 \
 + aly1*alx5*alz6 - aly1*alx6*alz5 - alz1*alx5*aly6 + alz1*alx6*aly5 + 2*alx1*aly5*alz7 - 2*alx1*alz5*aly7 \
 + alx2*aly5*alz6 - alx2*aly6*alz5 - 2*aly1*alx5*alz7 + 2*aly1*alx7*alz5 - aly2*alx5*alz6 + aly2*alx6*alz5 \
 + 2*alz1*alx5*aly7 - 2*alz1*aly5*alx7 + alz2*alx5*aly6 - alz2*alx6*aly5 - alx1*aly5*alz8 - alx1*aly6*alz7 \
 + alx1*alz5*aly8 + alx1*aly7*alz6 - 2*alx2*aly5*alz7 + 2*alx2*alz5*aly7 + aly1*alx5*alz8 + aly1*alx6*alz7 \
 - aly1*alx7*alz6 - aly1*alz5*alx8 + alx3*aly5*alz6 - alx3*aly6*alz5 + 2*aly2*alx5*alz7 - 2*aly2*alx7*alz5 \
 - alz1*alx5*aly8 - alz1*alx6*aly7 + alz1*aly5*alx8 + alz1*alx7*aly6 - aly3*alx5*alz6 + aly3*alx6*alz5 \
 - 2*alz2*alx5*aly7 + 2*alz2*aly5*alx7 + alx5*alz3*aly6 - alz3*alx6*aly5 + alx2*aly5*alz8 + alx2*aly6*alz7 \
- alx2*alz5*aly8 - alx2*aly7*alz6 - aly2*alx5*alz8 - aly2*alx6*alz7 + aly2*alx7*alz6 + aly2*alz5*alx8 \
- alx4*aly5*alz6 + alx4*aly6*alz5 + alz2*alx5*aly8 + alz2*alx6*aly7 - alz2*aly5*alx8 - alz2*alx7*aly6 \
+ alx5*aly4*alz6 - alx5*alz4*aly6 - aly4*alx6*alz5 + alx6*aly5*alz4 - 2*alx1*aly3*myco + 2*aly1*alx3*myco \
+ alx1*aly4*myco + alx2*aly3*myco - aly1*alx4*myco - alx3*aly2*myco + 2*aly1*alz3*myco - 2*alz1*aly3*myco \
- 2*aly1*alz4*myco - 2*aly2*alz3*myco + 2*alz1*aly4*myco + 2*aly3*alz2*myco + 2*alx1*aly7*myco - 2*aly1*alx7*myco \
- 2*alx3*aly5*myco + 2*aly2*alz4*myco + 2*aly3*alx5*myco - 2*alz2*aly4*myco - alx1*aly8*myco - alx2*aly7*myco \
+ aly1*alx8*myco + alx3*aly6*myco + aly2*alx7*myco + alx4*aly5*myco - aly3*alx6*myco - alx5*aly4*myco - aly1*alz7*myco \
 + alz1*aly7*myco + aly3*alz5*myco - alz3*aly5*myco + aly1*alz8*myco + aly2*alz7*myco - alz1*aly8*myco - aly3*alz6*myco \
 - alz2*aly7*myco - aly4*alz5*myco + alz3*aly6*myco + aly5*alz4*myco - aly2*alz8*myco + alz2*aly8*myco \
 - 2*alx5*aly7*myco + aly4*alz6*myco + 2*aly5*alx7*myco - alz4*aly6*myco + alx5*aly8*myco + alx6*aly7*myco \
 - aly5*alx8*myco - alx7*aly6*myco - aly1*myco*myco + aly2*myco*myco + aly3*myco*myco - aly4*myco*myco + aly5*myco*myco \
 - aly6*myco*myco - aly7*myco*myco + aly8*myco*myco;

alpha[7] = alx1*alz2*aly5 - alx1*aly2*alz5 + alx2*aly1*alz5 - alx2*alz1*aly5 - aly1*alz2*alx5 + aly2*alz1*alx5 \
+ alx1*aly2*alz7 - alx1*alz2*aly7 + alx1*aly4*alz5 - alx1*aly5*alz4 - alx2*aly1*alz7 + alx2*alz1*aly7 - alx2*aly3*alz5 \
+ alx2*alz3*aly5 - aly1*alx4*alz5 + aly1*alz2*alx7 + aly1*alx5*alz4 + alx3*aly2*alz5 - alx3*alz2*aly5 - aly2*alz1*alx7 \
- aly2*alx5*alz3 + alz1*alx4*aly5 - alz1*alx5*aly4 + aly3*alz2*alx5 - alx1*aly4*alz7 + alx1*alz4*aly7 + alx2*aly3*alz7 \
- alx2*alz3*aly7 + aly1*alx4*alz7 - aly1*alz4*alx7 - alx3*aly2*alz7 + alx3*alz2*aly7 - alx3*aly4*alz5 + alx3*aly5*alz4 \
+ aly2*alz3*alx7 - alz1*alx4*aly7 + alz1*aly4*alx7 + alx4*aly3*alz5 - alx4*alz3*aly5 - aly3*alz2*alx7 - aly3*alx5*alz4 \
+ alx5*aly4*alz3 + alx3*aly4*alz7 - alx3*alz4*aly7 - alx4*aly3*alz7 + alx4*alz3*aly7 + aly3*alz4*alx7 - aly4*alz3*alx7 \
- alx1*alz2*myco + alx2*alz1*myco + alx1*alz4*myco - alx2*alz3*myco + alx3*alz2*myco - alz1*alx4*myco + alx1*alz5*myco \
- alz1*alx5*myco - alx2*alz5*myco - alx3*alz4*myco + alx4*alz3*myco + alz2*alx5*myco - alx1*alz7*myco - alx3*alz5*myco \
+ alz1*alx7*myco + alx5*alz3*myco + alx2*alz7*myco + alx4*alz5*myco - alz2*alx7*myco - alx5*alz4*myco + alx3*alz7*myco \
- alz3*alx7*myco - alx4*alz7*myco + alz4*alx7*myco;

alpha[8] = alx1*aly2*alz3 - alx1*aly3*alz2 - alx2*aly1*alz3 + alx2*alz1*aly3 + aly1*alx3*alz2 - alx3*aly2*alz1 \
- alx1*aly2*alz4 + alx1*alz2*aly4 + alx2*aly1*alz4 - alx2*alz1*aly4 - aly1*alx4*alz2 + aly2*alz1*alx4 - 2*alx1*aly2*alz5\
 + alx1*aly3*alz4 + 2*alx1*alz2*aly5 - alx1*aly4*alz3 + 2*alx2*aly1*alz5 - 2*alx2*alz1*aly5 - aly1*alx3*alz4 \
 + aly1*alx4*alz3 - 2*aly1*alz2*alx5 + alx3*alz1*aly4 + 2*aly2*alz1*alx5 - alz1*alx4*aly3 + 2*alx1*aly2*alz6 \
 + 2*alx1*aly3*alz5 - 2*alx1*alz2*aly6 - 2*alx1*alz3*aly5 - 2*alx2*aly1*alz6 + 2*alx2*alz1*aly6 - alx2*aly3*alz4 \
 + alx2*aly4*alz3 - 2*aly1*alx3*alz5 + 2*aly1*alz2*alx6 + 2*aly1*alx5*alz3 + alx3*aly2*alz4 + 2*alx3*alz1*aly5 \
 - alx3*alz2*aly4 - 2*aly2*alz1*alx6 - aly2*alx4*alz3 - 2*alz1*aly3*alx5 + alx4*aly3*alz2 + alx1*aly2*alz7 \
 - alx1*aly3*alz6 - alx1*alz2*aly7 + alx1*alz3*aly6 - alx2*aly1*alz7 + alx2*alz1*aly7 - 3*alx2*aly3*alz5 \
 + 3*alx2*alz3*aly5 + aly1*alx3*alz6 + aly1*alz2*alx7 - aly1*alz3*alx6 + 3*alx3*aly2*alz5 - alx3*alz1*aly6 \
 - 3*alx3*alz2*aly5 - aly2*alz1*alx7 - 3*aly2*alx5*alz3 + alz1*aly3*alx6 + 3*aly3*alz2*alx5 - alx1*aly2*alz8 \
 - 2*alx1*aly3*alz7 + alx1*alz2*aly8 - alx1*aly4*alz6 + 2*alx1*alz3*aly7 + alx1*alz4*aly6 + alx2*aly1*alz8 \
 - alx2*alz1*aly8 + 2*alx2*aly3*alz6 + alx2*aly4*alz5 - 2*alx2*alz3*aly6 - alx2*aly5*alz4 + 2*aly1*alx3*alz7 \
 + aly1*alx4*alz6 - aly1*alz2*alx8 - 2*aly1*alz3*alx7 - aly1*alx6*alz4 - 2*alx3*aly2*alz6 - 2*alx3*alz1*aly7 \
 + 2*alx3*alz2*aly6 + aly2*alz1*alx8 - aly2*alx4*alz5 + aly2*alx5*alz4 + 2*aly2*alz3*alx6 - alz1*alx4*aly6 \
 + 2*alz1*aly3*alx7 + alz1*aly4*alx6 + alx4*alz2*aly5 - 2*aly3*alz2*alx6 - alz2*alx5*aly4 + alx1*aly3*alz8 \
 + alx1*aly4*alz7 - alx1*alz3*aly8 - alx1*alz4*aly7 + 2*alx2*aly3*alz7 - 2*alx2*alz3*aly7 - aly1*alx3*alz8 \
 - aly1*alx4*alz7 + aly1*alz3*alx8 + aly1*alz4*alx7 - 2*alx3*aly2*alz7 + alx3*alz1*aly8 + 2*alx3*alz2*aly7 \
 - alx3*aly4*alz5 + alx3*aly5*alz4 + 2*aly2*alz3*alx7 + alz1*alx4*aly7 - alz1*aly3*alx8 - alz1*aly4*alx7 \
 + alx4*aly3*alz5 - alx4*alz3*aly5 - 2*aly3*alz2*alx7 - aly3*alx5*alz4 + alx5*aly4*alz3 - alx2*aly3*alz8 \
 - alx2*aly4*alz7 + alx2*alz3*aly8 + alx2*alz4*aly7 + alx3*aly2*alz8 - alx3*alz2*aly8 + alx3*aly4*alz6 - alx3*alz4*aly6 \
  + aly2*alx4*alz7 - aly2*alz3*alx8 - aly2*alz4*alx7 - alx4*aly3*alz6 - alx4*alz2*aly7 + alx4*alz3*aly6 + aly3*alz2*alx8 \
  + aly3*alx6*alz4 + alz2*aly4*alx7 - aly4*alz3*alx6 + 2*alx1*aly5*myco - 2*aly1*alx5*myco - alx1*aly6*myco \
  + 2*alx1*alz5*myco - alx2*aly5*myco + aly1*alx6*myco + aly2*alx5*myco - 2*alz1*alx5*myco - 2*alx1*aly7*myco \
  - 2*alx1*alz6*myco - 2*alx2*alz5*myco + 2*aly1*alx7*myco - 2*alx3*aly5*myco + 2*alz1*alx6*myco + 2*aly3*alx5*myco \
  + 2*alz2*alx5*myco + alx1*aly8*myco - alx1*alz7*myco + alx2*aly7*myco + 2*alx2*alz6*myco - aly1*alx8*myco \
  + alx3*aly6*myco - alx3*alz5*myco - aly2*alx7*myco + alz1*alx7*myco + alx4*aly5*myco - aly3*alx6*myco - 2*alz2*alx6*myco \
  - alx5*aly4*myco + alx5*alz3*myco + alx1*alz8*myco + alx2*alz7*myco + 2*alx3*aly7*myco + alx3*alz6*myco - alz1*alx8*myco \
   + alx4*alz5*myco - 2*aly3*alx7*myco - alz2*alx7*myco - alx5*alz4*myco - alz3*alx6*myco - alx2*alz8*myco \
   - alx3*aly8*myco - alx4*aly7*myco - alx4*alz6*myco + aly3*alx8*myco + alz2*alx8*myco + aly4*alx7*myco \
   + alx6*alz4*myco - alx1*myco*myco + alx2*myco*myco + alx3*myco*myco - alx4*myco*myco + alx5*myco*myco - alx6*myco*myco - alx7*myco*myco + alx8*myco*myco;

alpha[9] = alx1*aly3*alz5 - alx1*alz3*aly5 - aly1*alx3*alz5 + aly1*alx5*alz3 + alx3*alz1*aly5 - alz1*aly3*alx5 \
- alx1*aly3*alz6 - alx1*aly4*alz5 + alx1*alz3*aly6 + alx1*aly5*alz4 - alx2*aly3*alz5 + alx2*alz3*aly5 + aly1*alx3*alz6 \
 + aly1*alx4*alz5 - aly1*alx5*alz4 - aly1*alz3*alx6 + alx3*aly2*alz5 - alx3*alz1*aly6 - alx3*alz2*aly5 - aly2*alx5*alz3 \
  - alz1*alx4*aly5 + alz1*aly3*alx6 + alz1*alx5*aly4 + aly3*alz2*alx5 + alx1*aly4*alz6 - alx1*alz4*aly6 + alx2*aly3*alz6 \
   + alx2*aly4*alz5 - alx2*alz3*aly6 - alx2*aly5*alz4 - aly1*alx4*alz6 + aly1*alx6*alz4 - alx3*aly2*alz6 + alx3*alz2*aly6 \
    - aly2*alx4*alz5 + aly2*alx5*alz4 + aly2*alz3*alx6 + alz1*alx4*aly6 - alz1*aly4*alx6 + alx4*alz2*aly5 - aly3*alz2*alx6 \
    - alz2*alx5*aly4 - alx2*aly4*alz6 + alx2*alz4*aly6 + aly2*alx4*alz6 - aly2*alx6*alz4 - alx4*alz2*aly6 + alz2*aly4*alx6 \
    - alx1*aly3*myco + aly1*alx3*myco + alx1*aly4*myco + alx2*aly3*myco - aly1*alx4*myco - alx3*aly2*myco + alx1*aly5*myco \
    - alx2*aly4*myco - aly1*alx5*myco + aly2*alx4*myco - alx1*aly6*myco - alx2*aly5*myco + aly1*alx6*myco + aly2*alx5*myco \
    + alx2*aly6*myco - alx3*aly5*myco - aly2*alx6*myco + aly3*alx5*myco + alx3*aly6*myco + alx4*aly5*myco - aly3*alx6*myco \
    - alx5*aly4*myco - alx4*aly6*myco + aly4*alx6*myco;

alpha[10] = alx1*aly3*alz2 - alx1*aly2*alz3 + alx2*aly1*alz3 - alx2*alz1*aly3 - aly1*alx3*alz2 + alx3*aly2*alz1 \
- alx1*aly3*alz4 + alx1*aly4*alz3 + aly1*alx3*alz4 - aly1*alx4*alz3 - alx3*alz1*aly4 + alz1*alx4*aly3 + alx1*aly2*alz7 \
 - alx1*aly3*alz6 - alx1*alz2*aly7 + alx1*alz3*aly6 - alx2*aly1*alz7 + alx2*alz1*aly7 + alx2*aly3*alz5 \
 - alx2*alz3*aly5 + aly1*alx3*alz6 + aly1*alz2*alx7 - aly1*alz3*alx6 - alx3*aly2*alz5 - alx3*alz1*aly6 + alx3*alz2*aly5 \
  - aly2*alz1*alx7 + aly2*alx5*alz3 + alz1*aly3*alx6 - aly3*alz2*alx5 + alx1*aly3*alz8 - alx1*aly4*alz7 - alx1*alz3*aly8 \
  + alx1*alz4*aly7 - aly1*alx3*alz8 + aly1*alx4*alz7 + aly1*alz3*alx8 - aly1*alz4*alx7 + alx3*alz1*aly8 + alx3*aly4*alz5 \
   - alx3*aly5*alz4 - alz1*alx4*aly7 - alz1*aly3*alx8 + alz1*aly4*alx7 - alx4*aly3*alz5 + alx4*alz3*aly5 + aly3*alx5*alz4 \
    - alx5*aly4*alz3 - alx1*aly6*alz7 + alx1*aly7*alz6 + alx2*aly5*alz7 - alx2*alz5*aly7 + aly1*alx6*alz7 - aly1*alx7*alz6 \
     - alx3*aly5*alz6 + alx3*aly6*alz5 - aly2*alx5*alz7 + aly2*alx7*alz5 - alz1*alx6*aly7 + alz1*alx7*aly6 + aly3*alx5*alz6 \
     - aly3*alx6*alz5 + alz2*alx5*aly7 - alz2*aly5*alx7 - alx5*alz3*aly6 + alz3*alx6*aly5 - alx1*aly7*alz8 \
     + alx1*aly8*alz7 + aly1*alx7*alz8 - aly1*alx8*alz7 + alx3*aly5*alz8 - alx3*alz5*aly8 - alz1*alx7*aly8 \
     + alz1*alx8*aly7 - alx4*aly5*alz7 + alx4*alz5*aly7 - aly3*alx5*alz8 + aly3*alz5*alx8 + alx5*aly4*alz7 \
     + alx5*alz3*aly8 - alx5*alz4*aly7 - aly4*alx7*alz5 - alz3*aly5*alx8 + aly5*alz4*alx7 + alx5*aly6*alz7 \
     - alx5*aly7*alz6 - alx6*aly5*alz7 + alx6*alz5*aly7 + aly5*alx7*alz6 - alx7*aly6*alz5 + alx5*aly7*alz8 \
     - alx5*aly8*alz7 - aly5*alx7*alz8 + aly5*alx8*alz7 + alx7*alz5*aly8 - alz5*alx8*aly7 + aly1*alz2*myco \
     - aly2*alz1*myco - aly1*alz4*myco + aly2*alz3*myco + alz1*aly4*myco - aly3*alz2*myco - aly1*alz6*myco \
     + aly2*alz5*myco + alz1*aly6*myco + aly3*alz4*myco - alz2*aly5*myco - aly4*alz3*myco + aly1*alz8*myco \
     - aly2*alz7*myco - alz1*aly8*myco + aly3*alz6*myco + alz2*aly7*myco - aly4*alz5*myco - alz3*aly6*myco \
     + aly5*alz4*myco - aly3*alz8*myco + aly4*alz7*myco + alz3*aly8*myco + aly5*alz6*myco - alz4*aly7*myco \
 - aly6*alz5*myco - aly5*alz8*myco + aly6*alz7*myco + alz5*aly8*myco - aly7*alz6*myco + aly7*alz8*myco - aly8*alz7*myco;

alpha[11] = alx1*aly3*alz2 - alx1*aly2*alz3 + alx2*aly1*alz3 - alx2*alz1*aly3 - aly1*alx3*alz2 + alx3*aly2*alz1 \
+ alx1*aly2*alz4 - alx1*alz2*aly4 - alx2*aly1*alz4 + alx2*alz1*aly4 + aly1*alx4*alz2 - aly2*alz1*alx4 + alx1*aly2*alz7 \
 - alx1*aly3*alz6 - alx1*alz2*aly7 + alx1*alz3*aly6 - alx2*aly1*alz7 + alx2*alz1*aly7 + alx2*aly3*alz5 - alx2*alz3*aly5 \
  + aly1*alx3*alz6 + aly1*alz2*alx7 - aly1*alz3*alx6 - alx3*aly2*alz5 - alx3*alz1*aly6 + alx3*alz2*aly5 - aly2*alz1*alx7 \
  + aly2*alx5*alz3 + alz1*aly3*alx6 - aly3*alz2*alx5 - alx1*aly2*alz8 + alx1*alz2*aly8 + alx1*aly4*alz6 \
  - alx1*alz4*aly6 + alx2*aly1*alz8 - alx2*alz1*aly8 - alx2*aly4*alz5 + alx2*aly5*alz4 - aly1*alx4*alz6 \
  - aly1*alz2*alx8 + aly1*alx6*alz4 + aly2*alz1*alx8 + aly2*alx4*alz5 - aly2*alx5*alz4 + alz1*alx4*aly6 \
  - alz1*aly4*alx6 - alx4*alz2*aly5 + alz2*alx5*aly4 - alx1*aly6*alz7 + alx1*aly7*alz6 + alx2*aly5*alz7 \
  - alx2*alz5*aly7 + aly1*alx6*alz7 - aly1*alx7*alz6 - alx3*aly5*alz6 + alx3*aly6*alz5 - aly2*alx5*alz7 \
  + aly2*alx7*alz5 - alz1*alx6*aly7 + alz1*alx7*aly6 + aly3*alx5*alz6 - aly3*alx6*alz5 + alz2*alx5*aly7 \
  - alz2*aly5*alx7 - alx5*alz3*aly6 + alz3*alx6*aly5 + alx1*aly6*alz8 - alx1*alz6*aly8 - alx2*aly5*alz8 \
  + alx2*alz5*aly8 - aly1*alx6*alz8 + aly1*alx8*alz6 + aly2*alx5*alz8 - aly2*alz5*alx8 + alz1*alx6*aly8 \
  - alz1*aly6*alx8 + alx4*aly5*alz6 - alx4*aly6*alz5 - alz2*alx5*aly8 + alz2*aly5*alx8 - alx5*aly4*alz6 \
  + alx5*alz4*aly6 + aly4*alx6*alz5 - alx6*aly5*alz4 + alx5*aly6*alz7 - alx5*aly7*alz6 - alx6*aly5*alz7 \
  + alx6*alz5*aly7 + aly5*alx7*alz6 - alx7*aly6*alz5 - alx5*aly6*alz8 + alx5*alz6*aly8 + alx6*aly5*alz8 \
  - alx6*alz5*aly8 - aly5*alx8*alz6 + aly6*alz5*alx8 - aly1*alz3*myco + alz1*aly3*myco + aly1*alz4*myco \
  + aly2*alz3*myco - alz1*aly4*myco - aly3*alz2*myco - aly2*alz4*myco + alz2*aly4*myco + aly1*alz7*myco \
  - alz1*aly7*myco - aly3*alz5*myco + alz3*aly5*myco - aly1*alz8*myco - aly2*alz7*myco + alz1*aly8*myco \
  + aly3*alz6*myco + alz2*aly7*myco + aly4*alz5*myco - alz3*aly6*myco - aly5*alz4*myco + aly2*alz8*myco \
  - alz2*aly8*myco - aly4*alz6*myco + alz4*aly6*myco - aly5*alz7*myco \
 + alz5*aly7*myco + aly5*alz8*myco + aly6*alz7*myco - alz5*aly8*myco - aly7*alz6*myco - aly6*alz8*myco + alz6*aly8*myco;

alpha[12] = alx1*aly2*alz5 - alx1*alz2*aly5 - alx2*aly1*alz5 + alx2*alz1*aly5 + aly1*alz2*alx5 - aly2*alz1*alx5 \
- alx1*aly2*alz7 + alx1*alz2*aly7 - alx1*aly4*alz5 + alx1*aly5*alz4 + alx2*aly1*alz7 - alx2*alz1*aly7 + alx2*aly3*alz5 \
 - alx2*alz3*aly5 + aly1*alx4*alz5 - aly1*alz2*alx7 - aly1*alx5*alz4 - alx3*aly2*alz5 + alx3*alz2*aly5 + aly2*alz1*alx7 \
 + aly2*alx5*alz3 - alz1*alx4*aly5 + alz1*alx5*aly4 - aly3*alz2*alx5 + alx1*aly4*alz7 + alx1*aly5*alz6 - alx1*alz4*aly7 \
 - alx1*aly6*alz5 - alx2*aly3*alz7 + alx2*alz3*aly7 - aly1*alx4*alz7 - aly1*alx5*alz6 + aly1*alx6*alz5 + aly1*alz4*alx7 \
 + alx3*aly2*alz7 - alx3*alz2*aly7 + alx3*aly4*alz5 - alx3*aly5*alz4 - aly2*alz3*alx7 + alz1*alx4*aly7 + alz1*alx5*aly6 \
 - alz1*aly4*alx7 - alz1*alx6*aly5 - alx4*aly3*alz5 + alx4*alz3*aly5 + aly3*alz2*alx7 + aly3*alx5*alz4 - alx5*aly4*alz3 \
 - alx1*aly5*alz8 + alx1*aly6*alz7 + alx1*alz5*aly8 - alx1*aly7*alz6 + aly1*alx5*alz8 - aly1*alx6*alz7 + aly1*alx7*alz6 \
 - aly1*alz5*alx8 - alx3*aly4*alz7 - alx3*aly5*alz6 + alx3*alz4*aly7 + alx3*aly6*alz5 - alz1*alx5*aly8 + alz1*alx6*aly7 \
 + alz1*aly5*alx8 - alz1*alx7*aly6 + alx4*aly3*alz7 - alx4*alz3*aly7 + aly3*alx5*alz6 - aly3*alx6*alz5 - aly3*alz4*alx7 \
 - alx5*alz3*aly6 + aly4*alz3*alx7 + alz3*alx6*aly5 + alx1*aly7*alz8 - alx1*aly8*alz7 - aly1*alx7*alz8 + aly1*alx8*alz7 \
 + alx3*aly5*alz8 - alx3*aly6*alz7 - alx3*alz5*aly8 + alx3*aly7*alz6 + alz1*alx7*aly8 - alz1*alx8*aly7 - aly3*alx5*alz8 \
 + aly3*alx6*alz7 - aly3*alx7*alz6 + aly3*alz5*alx8 + alx5*alz3*aly8 - alz3*alx6*aly7 - alz3*aly5*alx8 + alz3*alx7*aly6 \
 - alx3*aly7*alz8 + alx3*aly8*alz7 + aly3*alx7*alz8 - aly3*alx8*alz7 - alz3*alx7*aly8 + alz3*alx8*aly7 + alx1*alz2*myco \
 - alx2*alz1*myco - alx1*alz4*myco + alx2*alz3*myco - alx3*alz2*myco + alz1*alx4*myco - alx1*alz6*myco + alx2*alz5*myco \
 + alx3*alz4*myco + alz1*alx6*myco - alx4*alz3*myco - alz2*alx5*myco + alx1*alz8*myco - alx2*alz7*myco + alx3*alz6*myco \
 - alz1*alx8*myco - alx4*alz5*myco + alz2*alx7*myco + alx5*alz4*myco - alz3*alx6*myco - alx3*alz8*myco + alx4*alz7*myco \
 + alx5*alz6*myco + alz3*alx8*myco - alx6*alz5*myco - alz4*alx7*myco - alx5*alz8*myco + alx6*alz7*myco - alx7*alz6*myco \
 + alz5*alx8*myco + alx7*alz8*myco - alx8*alz7*myco;

alpha[13] = 2*alx1*aly3*alz2 - 2*alx1*aly2*alz3 + 2*alx2*aly1*alz3 - 2*alx2*alz1*aly3 - 2*aly1*alx3*alz2 + 2*alx3*aly2*alz1 \
 + 2*alx1*aly2*alz4 - 2*alx1*alz2*aly4 - 2*alx2*aly1*alz4 + 2*alx2*alz1*aly4 + 2*aly1*alx4*alz2 - 2*aly2*alz1*alx4 + \
 2*alx1*aly2*alz5 - 2*alx1*aly3*alz4 - 2*alx1*alz2*aly5 + 2*alx1*aly4*alz3 - 2*alx2*aly1*alz5 + 2*alx2*alz1*aly5 \
 + 2*aly1*alx3*alz4 - 2*aly1*alx4*alz3 + 2*aly1*alz2*alx5 - 2*alx3*alz1*aly4 - 2*aly2*alz1*alx5 + 2*alz1*alx4*aly3 \
 - 2*alx1*aly2*alz6 - 2*alx1*aly3*alz5 + 2*alx1*alz2*aly6 + 2*alx1*alz3*aly5 + 2*alx2*aly1*alz6 - 2*alx2*alz1*aly6 \
 + 2*alx2*aly3*alz4 - 2*alx2*aly4*alz3 + 2*aly1*alx3*alz5 - 2*aly1*alz2*alx6 - 2*aly1*alx5*alz3 - 2*alx3*aly2*alz4 \
 - 2*alx3*alz1*aly5 + 2*alx3*alz2*aly4 + 2*aly2*alz1*alx6 + 2*aly2*alx4*alz3 + 2*alz1*aly3*alx5 - 2*alx4*aly3*alz2 \
 + 4*alx2*aly3*alz5 - 4*alx2*alz3*aly5 - 4*alx3*aly2*alz5 + 4*alx3*alz2*aly5 + 4*aly2*alx5*alz3 - 4*aly3*alz2*alx5 \
 + 2*alx1*aly3*alz7 + 2*alx1*aly4*alz6 - 2*alx1*alz3*aly7 - 2*alx1*alz4*aly6 - 2*alx2*aly3*alz6 - 2*alx2*aly4*alz5 \
 + 2*alx2*alz3*aly6 + 2*alx2*aly5*alz4 - 2*aly1*alx3*alz7 - 2*aly1*alx4*alz6 + 2*aly1*alz3*alx7 + 2*aly1*alx6*alz4 \
 + 2*alx3*aly2*alz6 + 2*alx3*alz1*aly7 - 2*alx3*alz2*aly6 + 2*aly2*alx4*alz5 - 2*aly2*alx5*alz4 - 2*aly2*alz3*alx6 \
 + 2*alz1*alx4*aly6 - 2*alz1*aly3*alx7 - 2*alz1*aly4*alx6 - 2*alx4*alz2*aly5 + 2*aly3*alz2*alx6 + 2*alz2*alx5*aly4 \
 - 2*alx1*aly4*alz7 + 2*alx1*aly5*alz6 + 2*alx1*alz4*aly7 - 2*alx1*aly6*alz5 - 2*alx2*aly3*alz7 + 2*alx2*alz3*aly7 \
 + 2*aly1*alx4*alz7 - 2*aly1*alx5*alz6 + 2*aly1*alx6*alz5 - 2*aly1*alz4*alx7 + 2*alx3*aly2*alz7 - 2*alx3*alz2*aly7 \
 + 2*alx3*aly4*alz5 - 2*alx3*aly5*alz4 - 2*aly2*alz3*alx7 - 2*alz1*alx4*aly7\
 + 2*alz1*alx5*aly6 + 2*alz1*aly4*alx7 - 2*alz1*alx6*aly5 - 2*alx4*aly3*alz5 + 2*alx4*alz3*aly5 + 2*aly3*alz2*alx7 \
 + 2*aly3*alx5*alz4 - 2*alx5*aly4*alz3 - 2*alx1*aly5*alz7 + 2*alx1*alz5*aly7 + 2*alx2*aly4*alz7 - 2*alx2*aly5*alz6 \
 - 2*alx2*alz4*aly7 + 2*alx2*aly6*alz5 + 2*aly1*alx5*alz7 - 2*aly1*alx7*alz5 - 2*alx3*aly4*alz6 + 2*alx3*alz4*aly6 \
 - 2*aly2*alx4*alz7 + 2*aly2*alx5*alz6 - 2*aly2*alx6*alz5 + 2*aly2*alz4*alx7 - 2*alz1*alx5*aly7 + 2*alz1*aly5*alx7 \
 + 2*alx4*aly3*alz6 + 2*alx4*alz2*aly7 - 2*alx4*alz3*aly6 - 2*aly3*alx6*alz4 - 2*alz2*alx5*aly6 - 2*alz2*aly4*alx7 \
 + 2*alz2*alx6*aly5 + 2*aly4*alz3*alx6 + 2*alx1*aly6*alz7 - 2*alx1*aly7*alz6 + 2*alx2*aly5*alz7 - 2*alx2*alz5*aly7 \
 - 2*aly1*alx6*alz7 + 2*aly1*alx7*alz6 - 2*alx3*aly5*alz6 + 2*alx3*aly6*alz5 - 2*aly2*alx5*alz7 + 2*aly2*alx7*alz5 \
 + 2*alz1*alx6*aly7 - 2*alz1*alx7*aly6 + 2*aly3*alx5*alz6 - 2*aly3*alx6*alz5 + 2*alz2*alx5*aly7 - 2*alz2*aly5*alx7 \
 - 2*alx5*alz3*aly6 + 2*alz3*alx6*aly5 - 2*alx2*aly6*alz7 + 2*alx2*aly7*alz6 + 2*alx3*aly5*alz7 - 2*alx3*alz5*aly7 \
 + 2*aly2*alx6*alz7 - 2*aly2*alx7*alz6 + 2*alx4*aly5*alz6 - 2*alx4*aly6*alz5 - 2*aly3*alx5*alz7 + 2*aly3*alx7*alz5 \
 - 2*alz2*alx6*aly7 + 2*alz2*alx7*aly6 - 2*alx5*aly4*alz6 + 2*alx5*alz3*aly7 + 2*alx5*alz4*aly6 + 2*aly4*alx6*alz5 \
 - 2*alz3*aly5*alx7 - 2*alx6*aly5*alz4 - 2*alx3*aly6*alz7 + 2*alx3*aly7*alz6 - 2*alx4*aly5*alz7 + 2*alx4*alz5*aly7 \
 + 2*aly3*alx6*alz7 - 2*aly3*alx7*alz6 + 2*alx5*aly4*alz7 - 2*alx5*alz4*aly7 - 2*aly4*alx7*alz5 - 2*alz3*alx6*aly7 \
 + 2*alz3*alx7*aly6 + 2*aly5*alz4*alx7 + 2*alx4*aly6*alz7 - 2*alx4*aly7*alz6 - 2*aly4*alx6*alz7 + 2*aly4*alx7*alz6 \
 + 2*alx6*alz4*aly7 - 2*alz4*alx7*aly6;

alpha[14] = alx1*alz3*aly5 - alx1*aly3*alz5 + aly1*alx3*alz5 - aly1*alx5*alz3 - alx3*alz1*aly5 + alz1*aly3*alx5 \
+ alx1*aly3*alz6 + alx1*aly4*alz5 - alx1*alz3*aly6 - alx1*aly5*alz4 + alx2*aly3*alz5 - alx2*alz3*aly5 - aly1*alx3*alz6 \
- aly1*alx4*alz5 + aly1*alx5*alz4 + aly1*alz3*alx6 - alx3*aly2*alz5 + alx3*alz1*aly6 + alx3*alz2*aly5 + aly2*alx5*alz3 \
+ alz1*alx4*aly5 - alz1*aly3*alx6 - alz1*alx5*aly4 - aly3*alz2*alx5 - alx1*aly4*alz6 + alx1*alz4*aly6 - alx2*aly3*alz6 \
- alx2*aly4*alz5 + alx2*alz3*aly6 + alx2*aly5*alz4 + aly1*alx4*alz6 - aly1*alx6*alz4 + alx3*aly2*alz6 - alx3*alz2*aly6 \
+ aly2*alx4*alz5 - aly2*alx5*alz4 - aly2*alz3*alx6 - alz1*alx4*aly6 + alz1*aly4*alx6 - alx4*alz2*aly5 + aly3*alz2*alx6 \
+ alz2*alx5*aly4 + alx2*aly4*alz6 - alx2*alz4*aly6 - aly2*alx4*alz6 + aly2*alx6*alz4 + alx4*alz2*aly6 - alz2*aly4*alx6 \
- alx1*aly5*alz7 + alx1*alz5*aly7 + aly1*alx5*alz7 - aly1*alx7*alz5 - alz1*alx5*aly7 + alz1*aly5*alx7 + alx1*aly5*alz8 \
+ alx1*aly6*alz7 - alx1*alz5*aly8 - alx1*aly7*alz6 + alx2*aly5*alz7 - alx2*alz5*aly7 - aly1*alx5*alz8 - aly1*alx6*alz7 \
+ aly1*alx7*alz6 + aly1*alz5*alx8 - aly2*alx5*alz7 + aly2*alx7*alz5 + alz1*alx5*aly8 + alz1*alx6*aly7 - alz1*aly5*alx8 \
- alz1*alx7*aly6 + alz2*alx5*aly7 - alz2*aly5*alx7 - alx1*aly6*alz8 + alx1*alz6*aly8 - alx2*aly5*alz8 - alx2*aly6*alz7 \
+ alx2*alz5*aly8 + alx2*aly7*alz6 + aly1*alx6*alz8 - aly1*alx8*alz6 + aly2*alx5*alz8 + aly2*alx6*alz7 - aly2*alx7*alz6 \
- aly2*alz5*alx8 - alz1*alx6*aly8 + alz1*aly6*alx8 - alz2*alx5*aly8 - alz2*alx6*aly7 + alz2*aly5*alx8 + alz2*alx7*aly6 \
+ alx2*aly6*alz8 - alx2*alz6*aly8 - aly2*alx6*alz8 + aly2*alx8*alz6 + alz2*alx6*aly8 - alz2*aly6*alx8 + alx1*aly3*myco \
- aly1*alx3*myco - alx1*aly4*myco - alx2*aly3*myco + aly1*alx4*myco + alx3*aly2*myco + alx2*aly4*myco - aly2*alx4*myco \
- alx1*aly7*myco + aly1*alx7*myco + alx3*aly5*myco - aly3*alx5*myco + alx1*aly8*myco + alx2*aly7*myco - aly1*alx8*myco \
- alx3*aly6*myco - aly2*alx7*myco - alx4*aly5*myco + aly3*alx6*myco + alx5*aly4*myco - alx2*aly8*myco + aly2*alx8*myco \
+ alx4*aly6*myco - aly4*alx6*myco + alx5*aly7*myco - aly5*alx7*myco - alx5*aly8*myco - alx6*aly7*myco + aly5*alx8*myco \
+ alx7*aly6*myco + alx6*aly8*myco - aly6*alx8*myco;

alpha[15] = alx1*aly2*alz5 - alx1*alz2*aly5 - alx2*aly1*alz5 + alx2*alz1*aly5 + aly1*alz2*alx5 - aly2*alz1*alx5 \
- alx1*aly2*alz6 + alx1*alz2*aly6 + alx2*aly1*alz6 - alx2*alz1*aly6 - aly1*alz2*alx6 + aly2*alz1*alx6 - alx1*aly2*alz7 \
 + alx1*alz2*aly7 - alx1*aly4*alz5 + alx1*aly5*alz4 + alx2*aly1*alz7 - alx2*alz1*aly7 + alx2*aly3*alz5 - alx2*alz3*aly5 \
 + aly1*alx4*alz5 - aly1*alz2*alx7 - aly1*alx5*alz4 - alx3*aly2*alz5 + alx3*alz2*aly5 + aly2*alz1*alx7 + aly2*alx5*alz3 \
 - alz1*alx4*aly5 + alz1*alx5*aly4 - aly3*alz2*alx5 + alx1*aly2*alz8 - alx1*alz2*aly8 + alx1*aly4*alz6 - alx1*alz4*aly6 \
 - alx2*aly1*alz8 + alx2*alz1*aly8 - alx2*aly3*alz6 + alx2*alz3*aly6 - aly1*alx4*alz6 + aly1*alz2*alx8 + aly1*alx6*alz4 \
 + alx3*aly2*alz6 - alx3*alz2*aly6 - aly2*alz1*alx8 - aly2*alz3*alx6 + alz1*alx4*aly6 - alz1*aly4*alx6 + aly3*alz2*alx6 \
 + alx1*aly4*alz7 - alx1*alz4*aly7 - alx2*aly3*alz7 + alx2*alz3*aly7 - aly1*alx4*alz7 + aly1*alz4*alx7 + alx3*aly2*alz7 \
 - alx3*alz2*aly7 + alx3*aly4*alz5 - alx3*aly5*alz4 - aly2*alz3*alx7 + alz1*alx4*aly7 - alz1*aly4*alx7 - alx4*aly3*alz5 \
 + alx4*alz3*aly5 + aly3*alz2*alx7 + aly3*alx5*alz4 - alx5*aly4*alz3 - alx1*aly4*alz8 + alx1*alz4*aly8 + alx2*aly3*alz8 \
 - alx2*alz3*aly8 + aly1*alx4*alz8 - aly1*alz4*alx8 - alx3*aly2*alz8 + alx3*alz2*aly8 - alx3*aly4*alz6 + alx3*alz4*aly6 \
 + aly2*alz3*alx8 - alz1*alx4*aly8 + alz1*aly4*alx8 + alx4*aly3*alz6 - alx4*alz3*aly6 - aly3*alz2*alx8 - aly3*alx6*alz4 \
 + aly4*alz3*alx6 - alx3*aly4*alz7 + alx3*alz4*aly7 + alx4*aly3*alz7 - alx4*alz3*aly7 - aly3*alz4*alx7 + aly4*alz3*alx7 \
 + alx3*aly4*alz8 - alx3*alz4*aly8 - alx4*aly3*alz8 + alx4*alz3*aly8 + aly3*alz4*alx8 - aly4*alz3*alx8 - alx1*alz5*myco \
 + alz1*alx5*myco + alx1*alz6*myco + alx2*alz5*myco - alz1*alx6*myco - alz2*alx5*myco + alx1*alz7*myco - alx2*alz6*myco \
 + alx3*alz5*myco - alz1*alx7*myco + alz2*alx6*myco - alx5*alz3*myco - alx1*alz8*myco - alx2*alz7*myco - alx3*alz6*myco \
 + alz1*alx8*myco - alx4*alz5*myco + alz2*alx7*myco + alx5*alz4*myco + alz3*alx6*myco + alx2*alz8*myco - alx3*alz7*myco \
 + alx4*alz6*myco - alz2*alx8*myco + alz3*alx7*myco \
 - alx6*alz4*myco + alx3*alz8*myco + alx4*alz7*myco - alz3*alx8*myco - alz4*alx7*myco - alx4*alz8*myco + alz4*alx8*myco;

alpha[16] = alx1*alz3*aly5 - alx1*aly3*alz5 + aly1*alx3*alz5 - aly1*alx5*alz3 - alx3*alz1*aly5 + alz1*aly3*alx5 \
+ alx1*aly3*alz6 + alx1*aly4*alz5 - alx1*alz3*aly6 - alx1*aly5*alz4 + alx2*aly3*alz5 - alx2*alz3*aly5 - aly1*alx3*alz6 \
 - aly1*alx4*alz5 + aly1*alx5*alz4 + aly1*alz3*alx6 - alx3*aly2*alz5 + alx3*alz1*aly6 + alx3*alz2*aly5 + aly2*alx5*alz3 \
 + alz1*alx4*aly5 - alz1*aly3*alx6 - alz1*alx5*aly4 - aly3*alz2*alx5 + alx1*aly3*alz7 - alx1*aly4*alz6 - alx1*alz3*aly7 \
 + alx1*alz4*aly6 - alx2*aly3*alz6 - alx2*aly4*alz5 + alx2*alz3*aly6 + alx2*aly5*alz4 - aly1*alx3*alz7 + aly1*alx4*alz6 \
 + aly1*alz3*alx7 - aly1*alx6*alz4 + alx3*aly2*alz6 + alx3*alz1*aly7 - alx3*alz2*aly6 + aly2*alx4*alz5 - aly2*alx5*alz4 \
 - aly2*alz3*alx6 - alz1*alx4*aly6 - alz1*aly3*alx7 + alz1*aly4*alx6 - alx4*alz2*aly5 + aly3*alz2*alx6 + alz2*alx5*aly4 \
 - alx1*aly3*alz8 - alx1*aly4*alz7 + alx1*alz3*aly8 + alx1*alz4*aly7 - alx2*aly3*alz7 + alx2*aly4*alz6 + alx2*alz3*aly7 \
 - alx2*alz4*aly6 + aly1*alx3*alz8 + aly1*alx4*alz7 - aly1*alz3*alx8 - aly1*alz4*alx7 + alx3*aly2*alz7 - alx3*alz1*aly8 \
 - alx3*alz2*aly7 - aly2*alx4*alz6 - aly2*alz3*alx7 + aly2*alx6*alz4 - alz1*alx4*aly7 + alz1*aly3*alx8 + alz1*aly4*alx7 \
 + alx4*alz2*aly6 + aly3*alz2*alx7 - alz2*aly4*alx6 + alx1*aly4*alz8 - alx1*alz4*aly8 + alx2*aly3*alz8 + alx2*aly4*alz7 \
 - alx2*alz3*aly8 - alx2*alz4*aly7 - aly1*alx4*alz8 + aly1*alz4*alx8 - alx3*aly2*alz8 + alx3*alz2*aly8 - aly2*alx4*alz7 \
 + aly2*alz3*alx8 + aly2*alz4*alx7 + alz1*alx4*aly8 - alz1*aly4*alx8 + alx4*alz2*aly7 - aly3*alz2*alx8 - alz2*aly4*alx7 \
 - alx2*aly4*alz8 + alx2*alz4*aly8 + aly2*alx4*alz8 - aly2*alz4*alx8 - alx4*alz2*aly8 + alz2*aly4*alx8 - alx1*aly5*myco \
 + aly1*alx5*myco + alx1*aly6*myco + alx2*aly5*myco - aly1*alx6*myco - aly2*alx5*myco + alx1*aly7*myco - alx2*aly6*myco \
 - aly1*alx7*myco + alx3*aly5*myco + aly2*alx6*myco - aly3*alx5*myco - alx1*aly8*myco - alx2*aly7*myco + aly1*alx8*myco \
 - alx3*aly6*myco + aly2*alx7*myco - alx4*aly5*myco + aly3*alx6*myco + alx5*aly4*myco + alx2*aly8*myco - alx3*aly7*myco \
 - aly2*alx8*myco + alx4*aly6*myco + aly3*alx7*myco \
 - aly4*alx6*myco + alx3*aly8*myco + alx4*aly7*myco - aly3*alx8*myco - aly4*alx7*myco - alx4*aly8*myco + aly4*alx8*myco;

alpha[17] = alx1*aly2*alz3 - alx1*aly3*alz2 - alx2*aly1*alz3 + alx2*alz1*aly3 + aly1*alx3*alz2 - alx3*aly2*alz1 \
- alx1*aly2*alz4 + alx1*alz2*aly4 + alx2*aly1*alz4 - alx2*alz1*aly4 - aly1*alx4*alz2 + aly2*alz1*alx4 + alx1*aly3*alz4 \
- alx1*aly4*alz3 - aly1*alx3*alz4 + aly1*alx4*alz3 + alx3*alz1*aly4 - alz1*alx4*aly3 - alx2*aly3*alz4 + alx2*aly4*alz3 \
+ alx3*aly2*alz4 - alx3*alz2*aly4 - aly2*alx4*alz3 + alx4*aly3*alz2 - alx1*aly2*alz7 + alx1*aly3*alz6 + alx1*alz2*aly7 \
- alx1*alz3*aly6 + alx2*aly1*alz7 - alx2*alz1*aly7 - alx2*aly3*alz5 + alx2*alz3*aly5 - aly1*alx3*alz6 - aly1*alz2*alx7 \
+ aly1*alz3*alx6 + alx3*aly2*alz5 + alx3*alz1*aly6 - alx3*alz2*aly5 + aly2*alz1*alx7 - aly2*alx5*alz3 - alz1*aly3*alx6 \
+ aly3*alz2*alx5 + alx1*aly2*alz8 - alx1*alz2*aly8 - alx1*aly4*alz6 + alx1*alz4*aly6 - alx2*aly1*alz8 + alx2*alz1*aly8 \
+ alx2*aly4*alz5 - alx2*aly5*alz4 + aly1*alx4*alz6 + aly1*alz2*alx8 - aly1*alx6*alz4 - aly2*alz1*alx8 - aly2*alx4*alz5 \
+ aly2*alx5*alz4 - alz1*alx4*aly6 + alz1*aly4*alx6 + alx4*alz2*aly5 - alz2*alx5*aly4 - alx1*aly3*alz8 + alx1*aly4*alz7 \
+ alx1*alz3*aly8 - alx1*alz4*aly7 + aly1*alx3*alz8 - aly1*alx4*alz7 - aly1*alz3*alx8 + aly1*alz4*alx7 - alx3*alz1*aly8 \
- alx3*aly4*alz5 + alx3*aly5*alz4 + alz1*alx4*aly7 + alz1*aly3*alx8 - alz1*aly4*alx7 + alx4*aly3*alz5 - alx4*alz3*aly5 \
- aly3*alx5*alz4 + alx5*aly4*alz3 + alx2*aly3*alz8 - alx2*aly4*alz7 - alx2*alz3*aly8 + alx2*alz4*aly7 - alx3*aly2*alz8 \
+ alx3*alz2*aly8 + alx3*aly4*alz6 - alx3*alz4*aly6 + aly2*alx4*alz7 + aly2*alz3*alx8 - aly2*alz4*alx7 - alx4*aly3*alz6 \
- alx4*alz2*aly7 + alx4*alz3*aly6 - aly3*alz2*alx8 + aly3*alx6*alz4 + alz2*aly4*alx7 - aly4*alz3*alx6 + alx1*aly6*alz7 \
- alx1*aly7*alz6 - alx2*aly5*alz7 + alx2*alz5*aly7 - aly1*alx6*alz7 + aly1*alx7*alz6 + alx3*aly5*alz6 - alx3*aly6*alz5 \
+ aly2*alx5*alz7 - aly2*alx7*alz5 + alz1*alx6*aly7 - alz1*alx7*aly6 - aly3*alx5*alz6 + aly3*alx6*alz5 - alz2*alx5*aly7 \
+ alz2*aly5*alx7 + alx5*alz3*aly6 - alz3*alx6*aly5 - alx1*aly6*alz8 + alx1*alz6*aly8 + alx2*aly5*alz8 - alx2*alz5*aly8 \
+ aly1*alx6*alz8 - aly1*alx8*alz6 - aly2*alx5*alz8 \
 + aly2*alz5*alx8 - alz1*alx6*aly8 + alz1*aly6*alx8 - alx4*aly5*alz6 + alx4*aly6*alz5 + alz2*alx5*aly8 - alz2*aly5*alx8 \
 + alx5*aly4*alz6 - alx5*alz4*aly6 - aly4*alx6*alz5 + alx6*aly5*alz4 + alx1*aly7*alz8 - alx1*aly8*alz7 - aly1*alx7*alz8 \
 + aly1*alx8*alz7 - alx3*aly5*alz8 + alx3*alz5*aly8 + alz1*alx7*aly8 - alz1*alx8*aly7 + alx4*aly5*alz7 - alx4*alz5*aly7 \
 + aly3*alx5*alz8 - aly3*alz5*alx8 - alx5*aly4*alz7 - alx5*alz3*aly8 + alx5*alz4*aly7 + aly4*alx7*alz5 + alz3*aly5*alx8 \
 - aly5*alz4*alx7 - alx2*aly7*alz8 + alx2*aly8*alz7 + alx3*aly6*alz8 - alx3*alz6*aly8 + aly2*alx7*alz8 - aly2*alx8*alz7 \
 - alx4*aly6*alz7 + alx4*aly7*alz6 - aly3*alx6*alz8 + aly3*alx8*alz6 - alz2*alx7*aly8 + alz2*alx8*aly7 + aly4*alx6*alz7 \
 - aly4*alx7*alz6 + alz3*alx6*aly8 - alz3*aly6*alx8 - alx6*alz4*aly7 + alz4*alx7*aly6 - alx5*aly6*alz7 + alx5*aly7*alz6 \
 + alx6*aly5*alz7 - alx6*alz5*aly7 - aly5*alx7*alz6 + alx7*aly6*alz5 + alx5*aly6*alz8 - alx5*alz6*aly8 - alx6*aly5*alz8 \
 + alx6*alz5*aly8 + aly5*alx8*alz6 - aly6*alz5*alx8 - alx5*aly7*alz8 + alx5*aly8*alz7 + aly5*alx7*alz8 - aly5*alx8*alz7 \
 - alx7*alz5*aly8 + alz5*alx8*aly7 + alx6*aly7*alz8 - alx6*aly8*alz7 - alx7*aly6*alz8 + alx7*alz6*aly8 + aly6*alx8*alz7 \
 - alx8*aly7*alz6;

alpha[18] = alx1*alz2*aly5 - alx1*aly2*alz5 + alx2*aly1*alz5 - alx2*alz1*aly5 - aly1*alz2*alx5 + aly2*alz1*alx5 \
+ alx1*aly2*alz6 - alx1*alz2*aly6 - alx2*aly1*alz6 + alx2*alz1*aly6 + aly1*alz2*alx6 - aly2*alz1*alx6 + alx1*aly2*alz7 \
- alx1*alz2*aly7 + alx1*aly4*alz5 - alx1*aly5*alz4 - alx2*aly1*alz7 + alx2*alz1*aly7 - alx2*aly3*alz5 + alx2*alz3*aly5 \
- aly1*alx4*alz5 + aly1*alz2*alx7 + aly1*alx5*alz4 + alx3*aly2*alz5 - alx3*alz2*aly5 - aly2*alz1*alx7 - aly2*alx5*alz3 \
+ alz1*alx4*aly5 - alz1*alx5*aly4 + aly3*alz2*alx5 - alx1*aly2*alz8 + alx1*alz2*aly8 - alx1*aly4*alz6 + alx1*alz4*aly6 \
+ alx2*aly1*alz8 - alx2*alz1*aly8 + alx2*aly3*alz6 - alx2*alz3*aly6 + aly1*alx4*alz6 - aly1*alz2*alx8 - aly1*alx6*alz4 \
- alx3*aly2*alz6 + alx3*alz2*aly6 + aly2*alz1*alx8 + aly2*alz3*alx6 - alz1*alx4*aly6 + alz1*aly4*alx6 - aly3*alz2*alx6 \
- alx1*aly4*alz7 - alx1*aly5*alz6 + alx1*alz4*aly7 + alx1*aly6*alz5 + alx2*aly3*alz7 - alx2*alz3*aly7 + aly1*alx4*alz7 \
+ aly1*alx5*alz6 - aly1*alx6*alz5 - aly1*alz4*alx7 - alx3*aly2*alz7 + alx3*alz2*aly7 - alx3*aly4*alz5 + alx3*aly5*alz4 \
+ aly2*alz3*alx7 - alz1*alx4*aly7 - alz1*alx5*aly6 + alz1*aly4*alx7 + alz1*alx6*aly5 + alx4*aly3*alz5 - alx4*alz3*aly5 \
- aly3*alz2*alx7 - aly3*alx5*alz4 + alx5*aly4*alz3 + alx1*aly4*alz8 - alx1*alz4*aly8 - alx2*aly3*alz8 + alx2*alz3*aly8 \
+ alx2*aly5*alz6 - alx2*aly6*alz5 - aly1*alx4*alz8 + aly1*alz4*alx8 + alx3*aly2*alz8 - alx3*alz2*aly8 + alx3*aly4*alz6 \
- alx3*alz4*aly6 - aly2*alx5*alz6 - aly2*alz3*alx8 + aly2*alx6*alz5 + alz1*alx4*aly8 - alz1*aly4*alx8 - alx4*aly3*alz6 \
+ alx4*alz3*aly6 + aly3*alz2*alx8 + aly3*alx6*alz4 + alz2*alx5*aly6 - alz2*alx6*aly5 - aly4*alz3*alx6 + alx1*aly5*alz8 \
- alx1*aly6*alz7 - alx1*alz5*aly8 + alx1*aly7*alz6 - aly1*alx5*alz8 + aly1*alx6*alz7 - aly1*alx7*alz6 + aly1*alz5*alx8 \
+ alx3*aly4*alz7 + alx3*aly5*alz6 - alx3*alz4*aly7 - alx3*aly6*alz5 + alz1*alx5*aly8 - alz1*alx6*aly7 - alz1*aly5*alx8 \
+ alz1*alx7*aly6 - alx4*aly3*alz7 + alx4*alz3*aly7 - aly3*alx5*alz6 + aly3*alx6*alz5 + aly3*alz4*alx7 + alx5*alz3*aly6 \
- aly4*alz3*alx7 - alz3*alx6*aly5 - alx2*aly5*alz8 \
 + alx2*aly6*alz7 + alx2*alz5*aly8 - alx2*aly7*alz6 - alx3*aly4*alz8 + alx3*alz4*aly8 + aly2*alx5*alz8 - aly2*alx6*alz7 \
 + aly2*alx7*alz6 - aly2*alz5*alx8 + alx4*aly3*alz8 - alx4*alz3*aly8 - alx4*aly5*alz6 + alx4*aly6*alz5 - aly3*alz4*alx8 \
 - alz2*alx5*aly8 + alz2*alx6*aly7 + alz2*aly5*alx8 - alz2*alx7*aly6 + alx5*aly4*alz6 - alx5*alz4*aly6 + aly4*alz3*alx8 \
 - aly4*alx6*alz5 + alx6*aly5*alz4 - alx1*aly7*alz8 + alx1*aly8*alz7 + aly1*alx7*alz8 - aly1*alx8*alz7 - alx3*aly5*alz8 \
 + alx3*aly6*alz7 + alx3*alz5*aly8 - alx3*aly7*alz6 - alz1*alx7*aly8 + alz1*alx8*aly7 + aly3*alx5*alz8 - aly3*alx6*alz7 \
 + aly3*alx7*alz6 - aly3*alz5*alx8 - alx5*alz3*aly8 + alz3*alx6*aly7 + alz3*aly5*alx8 - alz3*alx7*aly6 + alx2*aly7*alz8 \
 - alx2*aly8*alz7 - aly2*alx7*alz8 + aly2*alx8*alz7 + alx4*aly5*alz8 - alx4*aly6*alz7 - alx4*alz5*aly8 + alx4*aly7*alz6 \
 + alz2*alx7*aly8 - alz2*alx8*aly7 - alx5*aly4*alz8 + alx5*alz4*aly8 + aly4*alx6*alz7 - aly4*alx7*alz6 + aly4*alz5*alx8 \
 - alx6*alz4*aly7 - aly5*alz4*alx8 + alz4*alx7*aly6 + alx3*aly7*alz8 - alx3*aly8*alz7 - aly3*alx7*alz8 + aly3*alx8*alz7 \
 + alz3*alx7*aly8 - alz3*alx8*aly7 - alx4*aly7*alz8 + alx4*aly8*alz7 + aly4*alx7*alz8 - aly4*alx8*alz7 - alz4*alx7*aly8 \
 + alz4*alx8*aly7;

alpha[19] = alx1*aly3*alz5 - alx1*alz3*aly5 - aly1*alx3*alz5 + aly1*alx5*alz3 + alx3*alz1*aly5 - alz1*aly3*alx5 \
- alx1*aly3*alz6 - alx1*aly4*alz5 + alx1*alz3*aly6 + alx1*aly5*alz4 - alx2*aly3*alz5 + alx2*alz3*aly5 + aly1*alx3*alz6 \
+ aly1*alx4*alz5 - aly1*alx5*alz4 - aly1*alz3*alx6 + alx3*aly2*alz5 - alx3*alz1*aly6 - alx3*alz2*aly5 - aly2*alx5*alz3 \
- alz1*alx4*aly5 + alz1*aly3*alx6 + alz1*alx5*aly4 + aly3*alz2*alx5 - alx1*aly3*alz7 + alx1*aly4*alz6 + alx1*alz3*aly7 \
- alx1*alz4*aly6 + alx2*aly3*alz6 + alx2*aly4*alz5 - alx2*alz3*aly6 - alx2*aly5*alz4 + aly1*alx3*alz7 - aly1*alx4*alz6 \
 - aly1*alz3*alx7 + aly1*alx6*alz4 - alx3*aly2*alz6 - alx3*alz1*aly7 + alx3*alz2*aly6 - aly2*alx4*alz5 + aly2*alx5*alz4 \
 + aly2*alz3*alx6 + alz1*alx4*aly6 + alz1*aly3*alx7 - alz1*aly4*alx6 + alx4*alz2*aly5 - aly3*alz2*alx6 - alz2*alx5*aly4 \
 + alx1*aly3*alz8 + alx1*aly4*alz7 - alx1*alz3*aly8 - alx1*alz4*aly7 + alx2*aly3*alz7 - alx2*aly4*alz6 - alx2*alz3*aly7 \
 + alx2*alz4*aly6 - aly1*alx3*alz8 - aly1*alx4*alz7 + aly1*alz3*alx8 + aly1*alz4*alx7 - alx3*aly2*alz7 + alx3*alz1*aly8 \
 + alx3*alz2*aly7 + aly2*alx4*alz6 + aly2*alz3*alx7 - aly2*alx6*alz4 + alz1*alx4*aly7 - alz1*aly3*alx8 - alz1*aly4*alx7 \
 - alx4*alz2*aly6 - aly3*alz2*alx7 + alz2*aly4*alx6 - alx1*aly4*alz8 + alx1*aly5*alz7 + alx1*alz4*aly8 - alx1*alz5*aly7 \
 - alx2*aly3*alz8 - alx2*aly4*alz7 + alx2*alz3*aly8 + alx2*alz4*aly7 + aly1*alx4*alz8 - aly1*alx5*alz7 - aly1*alz4*alx8 \
 + aly1*alx7*alz5 + alx3*aly2*alz8 - alx3*alz2*aly8 + aly2*alx4*alz7 - aly2*alz3*alx8 - aly2*alz4*alx7 - alz1*alx4*aly8 \
 + alz1*alx5*aly7 + alz1*aly4*alx8 - alz1*aly5*alx7 \
 - alx4*alz2*aly7 + aly3*alz2*alx8 + alz2*aly4*alx7 - alx1*aly5*alz8 - alx1*aly6*alz7 + alx1*alz5*aly8 + alx1*aly7*alz6 \
 + alx2*aly4*alz8 - alx2*aly5*alz7 - alx2*alz4*aly8 + alx2*alz5*aly7 + aly1*alx5*alz8 + aly1*alx6*alz7 - aly1*alx7*alz6 \
 - aly1*alz5*alx8 - aly2*alx4*alz8 + aly2*alx5*alz7 + aly2*alz4*alx8 - aly2*alx7*alz5 - alz1*alx5*aly8 - alz1*alx6*aly7 \
 + alz1*aly5*alx8 + alz1*alx7*aly6 + alx4*alz2*aly8 - alz2*alx5*aly7 - alz2*aly4*alx8 + alz2*aly5*alx7 + alx1*aly6*alz8 \
 - alx1*alz6*aly8 + alx2*aly5*alz8 + alx2*aly6*alz7 - alx2*alz5*aly8 - alx2*aly7*alz6 - aly1*alx6*alz8 + aly1*alx8*alz6 \
 - alx3*aly5*alz7 + alx3*alz5*aly7 - aly2*alx5*alz8 - aly2*alx6*alz7 + aly2*alx7*alz6 + aly2*alz5*alx8 + alz1*alx6*aly8 \
 - alz1*aly6*alx8 + aly3*alx5*alz7 - aly3*alx7*alz5 + alz2*alx5*aly8 + alz2*alx6*aly7 - alz2*aly5*alx8 - alz2*alx7*aly6 \
 - alx5*alz3*aly7 + alz3*aly5*alx7 - alx2*aly6*alz8 + alx2*alz6*aly8 + alx3*aly5*alz8 + alx3*aly6*alz7 - alx3*alz5*aly8 \
 - alx3*aly7*alz6 + aly2*alx6*alz8 - aly2*alx8*alz6 + alx4*aly5*alz7 - alx4*alz5*aly7 - aly3*alx5*alz8 - aly3*alx6*alz7 \
 + aly3*alx7*alz6 + aly3*alz5*alx8 - alz2*alx6*aly8 + alz2*aly6*alx8 - alx5*aly4*alz7 + alx5*alz3*aly8 + alx5*alz4*aly7 \
 + aly4*alx7*alz5 + alz3*alx6*aly7 - alz3*aly5*alx8 - alz3*alx7*aly6 - aly5*alz4*alx7 - alx3*aly6*alz8 + alx3*alz6*aly8 \
 - alx4*aly5*alz8 - alx4*aly6*alz7 + alx4*alz5*aly8 + alx4*aly7*alz6 + aly3*alx6*alz8 - aly3*alx8*alz6 + alx5*aly4*alz8 \
 - alx5*alz4*aly8 + aly4*alx6*alz7 - aly4*alx7*alz6 - aly4*alz5*alx8 - alz3*alx6*aly8 + alz3*aly6*alx8 - alx6*alz4*aly7 \
 + aly5*alz4*alx8 + alz4*alx7*aly6 + alx4*aly6*alz8 - alx4*alz6*aly8 - aly4*alx6*alz8 + aly4*alx8*alz6 + alx6*alz4*aly8 \
 - alz4*aly6*alx8;

for (i=0; i<20; i++) alpha[i] = alpha[i] *  pow(myco,-3.0);

}
