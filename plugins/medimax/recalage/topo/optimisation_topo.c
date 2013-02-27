/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>

#include "math/imx_matrix.h"
#include "noyau/io/imx_log.h"

#include "recalage/chps_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/optimisation_topo.h"
#include "recalage/topo/distance_topo.h"
#include "recalage/topo/gradient_distance_topo.h"
#include "recalage/topo/hessien_distance_topo.h"

/*******************************************************************************
********************************************************************************
*************************** MINIMISATIONS **************************************
********************************************************************************
*******************************************************************************/


 
/*******************************************************************************
**     TOP_linemin_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter
**                   ,dist,nb_param,min_param,max_param,prec_param,param,direc)
**
**     minimisation unidimensionnelle suivant la direction direc
**  par decoupage de l'intervalle de recherche
**
**     utilisation:
**                  imref,imreca,imres: images
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**          direc: vecteur donnant la direction de minimisation
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retournes dans param
*******************************************************************************/

double TOP_linemin_locale_3d    (grphic3d *imref, grphic3d *imreca, grphic3d *imres,  
                                                            field3d *champ_fin, transf_func_t the_transf, InterpolationFct inter,
                                                            dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, 
                                                            double *param_norm,double *moins_direc,
                                                            int topi, int topj, int topk, double Jm, double JM,TSlpqr *Slpqr)
{
  double cf,cfmax,cfopt,Emin,E,*p,*p_norm;
  double max_grad_cfmax;
  int i,j,l,nb_pas,topD,continuer;
  int width,height,depth,resol; 
  int niter;
  double x0,x1,x2,x3=0,f1,f2,aux,dGOLD=0.61803;
  double precis;
  int i_min;
 //---Variables dediees au recalage symetrique -------
  double opp_cfmax,*opp_param_norm,*opp_moins_direc;
  double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;

  
  width=BASE3D.width;height=BASE3D.height;depth=BASE3D.depth;
  resol=BASE3D.resol; 
 
  p=CALLOC(nb_param,double);
  p_norm=CALLOC(nb_param,double);

  topD=(int)pow(2.0,1.0*BASE3D.resol)-1;
  l = TOP_conv_ind(topi,topj,topk,nb_param);
    
 
    
  //---------------------------------------------------------------
  //----- determination de cfmin et cfmax -------------------------
  //---------------------------------------------------------------

  moins_direc[l] = 1.0*moins_direc[l] / width;
  moins_direc[l+1] = 1.0*moins_direc[l+1] / height;
  moins_direc[l+2] = 1.0*moins_direc[l+2] / depth;

  
    cfmax = TOP_linemin_maxpas_3d(nb_param,param_norm,moins_direc,topi,topj,topk,Jm,JM,Slpqr);



// ------------------------------------------------------------
//----- Dans le cas du recalage symetrique, on calcul aussi ---
//----- le pas max pour la transfo symetrique et on prend le --
//----- le min des 2 ------------------------------------------
//-------------------------------------------------------------


    if((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
        || ((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
    {
    opp_param_norm=CALLOC(nb_param,double);
    opp_moins_direc=CALLOC(nb_param,double);
    
    for(i=0;i<nb_param;i++)
        {
        opp_param_norm[i]=-1.0*lambda_ref*param_norm[i];
        opp_moins_direc[i]=-1.0*moins_direc[i];
        }
    opp_cfmax = TOP_linemin_maxpas_3d(nb_param,opp_param_norm,opp_moins_direc,topi,topj,topk,Jm,JM,Slpqr);

    if (lambda_ref>0)
    opp_cfmax=opp_cfmax/lambda_ref;
    else
    opp_cfmax=HUGE_VAL;
    
    
    cfmax=MINI(cfmax,opp_cfmax);
    }



  //----- Recherche du deplacement relatif max ----
  max_grad_cfmax=0;
  if (fabs(moins_direc[l])>max_grad_cfmax) max_grad_cfmax=fabs(moins_direc[l]);
  if (fabs(moins_direc[l+1])>max_grad_cfmax) max_grad_cfmax=fabs(moins_direc[l+1]);
  if (fabs(moins_direc[l+2])>max_grad_cfmax) max_grad_cfmax=fabs(moins_direc[l+2]);
 
 
  moins_direc[l] = 1.0*moins_direc[l] * width;
  moins_direc[l+1] = 1.0*moins_direc[l+1] * height;
  moins_direc[l+2] = 1.0*moins_direc[l+2] * depth;
  
  
  nb_pas = NB_PAS_LINEMIN_TOPO;

cfopt=0.0;

if ((max_grad_cfmax*cfmax)>(nb_pas*TOPO_QUANTIFICATION_PARAM)) // on ne fait une recherche lineaire que si les parametres sont susceptibles de bouger de plus de TOPO_QUANTIFICATION_PARAM
 {
 if (max_grad_cfmax*cfmax<(1./topD))
    {  
        //---------------------------------------------------------------
        //---------- Recherche lineaire ---------------------------------
        //---------------------------------------------------------------
    
     
    Emin=dist(imref,imreca,nb_param,param,param_norm,topi,topj,topk,Slpqr,regularisation);
    
    cfopt=0.0;
    continuer = 1;
   
    for (i=0; i<nb_param; i++) {p[i] = param[i];p_norm[i]=param_norm[i];}
        i_min=0;
    for (i=1; (i<=nb_pas)/*&&(continuer==1)*/; i++)
        { 
        cf=1.0*i*cfmax/(1.0*nb_pas);
        for (j=l; j<l+3; j++) p[j] = param[j] - cf*moins_direc[j];
            p_norm[l]    =1.0*p[l]  /width;
            p_norm[l+1]=1.0*p[l+1]/height;
            p_norm[l+2]=1.0*p[l+2]/depth;
            
        E=dist(imref, imreca,nb_param,p,p_norm,topi,topj,topk,Slpqr,regularisation);
        
        if (E<Emin) { Emin=E; cfopt=cf;i_min=i;}
        if (E>Emin) {continuer=0;}
        
        } 
    i--;
        for (j=l; j<l+3; j++) p[j] = param[j] - cfopt*moins_direc[j];
        
        for (j=l; j<l+3; j++) 
                if (fabs(p[j])<TOPO_QUANTIFICATION_PARAM)  /* C'est pour en finir avec les bug num�riques ... */
                    p[j]=0.0;
        
        
            p_norm[l]    =1.0*p[l]  /width;
            p_norm[l+1]=1.0*p[l+1]/height;
            p_norm[l+2]=1.0*p[l+2]/depth;

    i=i_min;
    while ((TOP_verification_locale(nb_param,p,p_norm,Jm,JM,topi,topj,topk,Slpqr)==0)&&(i>0))
        { 
                i--;
            cf=1.0*i*cfmax/(1.0*nb_pas);
                cfopt = cf;
            
                for (j=l; j<l+3; j++) p[j] = param[j] - cf*moins_direc[j];
                
                p_norm[l]    =1.0*p[l]  /width;
                p_norm[l+1]=1.0*p[l+1]/height;
                p_norm[l+2]=1.0*p[l+2]/depth;
        #ifndef TOPO_VERBOSE
        printf("Avertissements ci-dessus pris en compte....\n\n");
        #endif  
                
        }
  
    }
else
    {
        //---------------------------------------------------------------
        //----- Recherche par nombre d'or -------------------------------
        //---------------------------------------------------------------
 
    
    for (j=0; j<nb_param; j++) {p[j] = param[j];p_norm[j]=param_norm[j];}
    Emin=dist(imref,imreca,nb_param,param,param_norm,topi,topj,topk,Slpqr,regularisation);
    i_min=0;
    
    // Initialisation
    for (i=1;i<nb_pas;i++)
    {
    cf=0.1*i*cfmax;
    for (j=l; j<l+3; j++) p[j] = param[j] - cf*moins_direc[j];
        p_norm[l]    =1.0*p[l]  /width;
        p_norm[l+1]=1.0*p[l+1]/height;
        p_norm[l+2]=1.0*p[l+2]/depth;
    E=dist(imref, imreca,nb_param,p,p_norm,topi,topj,topk,Slpqr,regularisation);
    if (E<Emin) {Emin = E; i_min = i;}
    }
    
    if (i_min == 0)
        {
        x0=0;
        x3=0.1*cfmax;
        x1=dGOLD*x3;    
        x2=x1+(1-dGOLD)*(x3-x1);
        }
    else
        {
        x0=0.1*(i_min-1)*cfmax;
        x1=0.1*i_min*cfmax;
        x2=x1+(1-dGOLD)*(x3-x1);
        x3=0.1*(i_min+1)*cfmax;
        }
            
    
    for (j=l; j<l+3; j++) p[j] = param[j] - x1*moins_direc[j];
    p_norm[l]    =1.0*p[l]  /width;
    p_norm[l+1]=1.0*p[l+1]/height;
    p_norm[l+2]=1.0*p[l+2]/depth;
    f1=dist(imref,imreca,nb_param,p,p_norm,topi,topj,topk,Slpqr,regularisation);
    
    for (j=l; j<l+3; j++) p[j] = param[j] - x2*moins_direc[j];
    p_norm[l]    =1.0*p[l]  /width;
    p_norm[l+1]=1.0*p[l+1]/height;
    p_norm[l+2]=1.0*p[l+2]/depth;
    f2=dist(imref,imreca,nb_param,p,p_norm,topi,topj,topk,Slpqr,regularisation);

    niter=0;
    aux=fabs(1.0*max_grad_cfmax*(x0-x3));
    precis= 0.1/topD;

    while (aux>precis)
        {
        if (f2<f1)
            {
            aux=dGOLD*x2+(1-dGOLD)*x3;
            x0=x1;x1=x2;x2=aux;
            f1=f2;
            for (j=l; j<l+3; j++) p[j] = param[j] - x2*moins_direc[j];
            p_norm[l]    =1.0*p[l]  /width;
            p_norm[l+1]=1.0*p[l+1]/height;
            p_norm[l+2]=1.0*p[l+2]/depth;
            f2=dist(imref,imreca,nb_param,p,p_norm,topi,topj,topk,Slpqr,regularisation);
            }
        else
            {
            aux=dGOLD*x1+(1-dGOLD)*x0;
            x3=x2;x2=x1;x1=aux;
            f2=f1;
            for (j=l; j<l+3; j++) p[j] = param[j] - x1*moins_direc[j];
            p_norm[l]    =1.0*p[l]  /width;
            p_norm[l+1]=1.0*p[l+1]/height;
            p_norm[l+2]=1.0*p[l+2]/depth;
            f1=dist(imref,imreca,nb_param,p,p_norm,topi,topj,topk,Slpqr,regularisation);
            }
        niter++;
        aux=fabs(1.0*max_grad_cfmax*(x0-x3));
        } 


    if ((Emin<f1)&&(Emin<f2)){cfopt=0;}
    else if (f1<f2) { cfopt = x1; Emin = f1;}
        else {cfopt = x2; Emin = f2;}
    
        
        
    i=nb_pas;   
    while ((TOP_verification_locale(nb_param,p,p_norm,Jm,JM,topi,topj,topk,Slpqr)==0)&&(i>0))
        { 
                i--;
            cf=1.0*i*cfopt/(1.0*nb_pas);
                cfopt = cf;
            
                for (j=l; j<l+3; j++) p[j] = param[j] - cf*moins_direc[j];
                
                p_norm[l]    =1.0*p[l]  /width;
                p_norm[l+1]=1.0*p[l+1]/height;
                p_norm[l+2]=1.0*p[l+2]/depth;
        #ifndef TOPO_VERBOSE
        printf("Avertissements ci-dessus pris en compte....\n");
        #endif
        }
        
    }
}
  
    // ------------------------------------------------------------
    //----- On verifie que le pas est aussi acceptable pour la ----
    //----- transfo symetrique ------------------------------------
    //-------------------------------------------------------------
    if(((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d))&&(cfopt != 0))
    {
    
        opp_param_norm[l] = -1.0*lambda_ref*(param[l] + cf*moins_direc[l])/width;
        opp_param_norm[l+1] = -1.0* lambda_ref*(param[l+1] + cf*moins_direc[l+1])/height;
        opp_param_norm[l+2] = -1.0* lambda_ref*(param[l+2] + cf*moins_direc[l+2])/depth;
        
        
    i=nb_pas;   
    while ((TOP_verification_locale(nb_param,p,opp_param_norm,Jm,JM,topi,topj,topk,Slpqr)==0)&&(i>0))
        { 
                i--;
            cf=1.0*i*cfopt/(1.0*nb_pas);
                cfopt = cf;
        
                opp_param_norm[l] = -1.0* lambda_ref*(param[l] + cf*moins_direc[l])/width;
                opp_param_norm[l+1] = -1.0* lambda_ref*(param[l+1] + cf*moins_direc[l+1])/height;
                opp_param_norm[l+2] = -1.0* lambda_ref*(param[l+2] + cf*moins_direc[l+2])/depth;
            
            }
        }

    
    //---------------------------------------------------------------
  //----- mise a jour des parametres ------------------------------
  //---------------------------------------------------------------
  if (cfopt != 0)
  { 
    for (j=l; j<l+3; j++) param[j] = param[j] - cfopt*moins_direc[j];    
    }
       
    param_norm[l]  =1.0*param[l]  /width;
    param_norm[l+1]=1.0*param[l+1]/height;
    param_norm[l+2]=1.0*param[l+2]/depth;

    Emin=dist(imref,imreca,nb_param,param,param_norm,topi,topj,topk,Slpqr,regularisation);
 
 
// ------------------------------------------------------------
//----- Liberation de variables uniquement allouees lors du----
//----- recalage symetrique------------------------------------
//-------------------------------------------------------------
    if((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d))
    {
    free(opp_param_norm);
    free(opp_moins_direc);
    }
 
 
 
 free(p);free(p_norm);
 return(Emin);
}


/*******************************************************************************
**     TOP_min_desc_grad_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     recalage par descente de gradient                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_desc_grad_locale_3d  (grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                            field3d *champ_fin,transf_func_t the_transf, InterpolationFct inter,
                                            dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i; 
 field3d  *gradIref,*maskgradIref=NULL;
 double *grad,*p,*p_i,gradl[3],*param_norm;
 double Edebut,Efin,E1,E2,E,Eprec;
 int    nb,nbreduc,stop;
 double prreduc,maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;   
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt,topD1;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11,l;
 double xm,xM,ym,yM,zm,zM;
 int (*gradient)();   
 reg_grad_locale_t gradient_reg;   
 
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif

    
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
  
 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_Lp_locale_3d) gradient=gradient_base_Lp_locale_3d;
 if (dist==Energie_Lp_sym_locale_3d) gradient=gradient_base_Lp_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_L1norm_locale_3d) gradient=gradient_base_L1norm_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) gradient=gradient_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) gradient=gradient_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) gradient=gradient_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) gradient=gradient_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) gradient=gradient_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) gradient=gradient_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) gradient=gradient_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) gradient=gradient_quad_locale_symetrique_coupe_3d;
 if (dist==Energie_quad_sous_champ_locale_3d) gradient=gradient_quad_sous_champ_locale_3d;
 if (dist==Energie_groupwise_variance_locale_3d) gradient=gradient_groupwise_variance_locale_3d;
 if (dist==Energie_groupwise_variance_locale_nonsym_3d) gradient=gradient_groupwise_variance_locale_nonsym_3d;

 gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) gradient_reg=gradient_regularisation_log_jacobien_local;
 if (regularisation==regularisation_energie_membrane_jacobien_local) gradient_reg=gradient_regularisation_energie_membrane_jacobien_local;
 if (regularisation==regularisation_log_jacobien_centre_local) gradient_reg=gradient_regularisation_log_jacobien_centre_local;
 if (regularisation==regularisation_dist_identite_local) gradient_reg=gradient_regularisation_dist_identite_local;
 if (regularisation==regularisation_patch_imagebased_local) gradient_reg=gradient_regularisation_patch_imagebased_local;
 
 
 /*allocation memoire des variables*/
 //p=alloc_dvector(nb_param);
 p_i=alloc_dvector(nb_param);
 param_norm=alloc_dvector(nb_param);
 
 grad=alloc_dvector(nb_param);
 
 
 if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
 gradIref=cr_field3d(wdth,hght,dpth);
 else
 gradIref=NULL;
 
 
 if (dist==Energie_groupwise_variance_locale_3d)
    {
     _ParamRecalageBspline.grad_reca_groupwise=cr_field3d(wdth,hght,dpth);
     _ParamRecalageBspline.grad_ref_groupwise=cr_field3d(wdth,hght,dpth);
     _ParamRecalageBspline.grad_ref2_groupwise=cr_field3d(wdth,hght,dpth);
     }
 

 masque_bloc=alloc_imatrix_3d(topD,topD,topD);

if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
 maskgradIref=cr_field3d(wdth,hght,dpth);


 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


/*Initialisation de masque_bloc*/
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;
    
for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

/* initialisation des parametres normalises sur [0,1] */
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }
    
//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 

 /*calcul de l'energie de depart*/
    Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  E2=Edebut;

 /* calcul du gradient de l'image a recaler */
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d))
    imx_gradient_3d_p(imreca,gradIref,0,4);
else
  if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
    imx_gradient_3d_p(imreca,gradIref,2,4);
 
 if (dist==Energie_ICP_sym_locale_3d)
    imx_gradient_3d_p(imreca->mask,maskgradIref,0,4);

//-------------------------------------------------------
//----- Pour le recalage symetrique, on calcule aussi----
//----- le gradient de imref ----------------------------
//-------------------------------------------------------

 if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
    imx_gradient_3d_p(imref,maskgradIref,2,4);





 if (dist==Energie_groupwise_variance_locale_3d)
    {
    
     imx_gradient_3d_p(_ParamRecalageBspline.reca_groupwise,_ParamRecalageBspline.grad_reca_groupwise,2,4);
     mul_field_3d(_ParamRecalageBspline.grad_reca_groupwise,_ParamRecalageBspline.grad_reca_groupwise, _ParamRecalageBspline.reca_groupwise->rcoeff);

     imx_gradient_3d_p(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.grad_ref_groupwise,2,4);
     mul_field_3d(_ParamRecalageBspline.grad_ref_groupwise, _ParamRecalageBspline.grad_ref_groupwise, _ParamRecalageBspline.ref_groupwise->rcoeff);
     
     imx_gradient_3d_p(_ParamRecalageBspline.ref2_groupwise,_ParamRecalageBspline.grad_ref2_groupwise,2,4);
     mul_field_3d(_ParamRecalageBspline.grad_ref2_groupwise, _ParamRecalageBspline.grad_ref2_groupwise, _ParamRecalageBspline.ref2_groupwise->rcoeff);
     
      /*// calcul plus chiade du gradient de ref2
     field3d  *grad_tmp;
     grad_tmp=cr_field3d(wdth,hght,dpth);
     
     for(i=0; i<wdth; i++)
     for(j=0; j<hght; j++)
     for(k=0; k<dpth; k++)
        {
        _ParamRecalageBspline.grad_ref2_groupwise->raw[i][j][k].x=0.0;
        _ParamRecalageBspline.grad_ref2_groupwise->raw[i][j][k].y=0.0;
        _ParamRecalageBspline.grad_ref2_groupwise->raw[i][j][k].z=0.0;
        }
     
     for (l=0; l< _ParamRecalageBspline.nb_tot_groupwise; l++)
        if (l!=_ParamRecalageBspline.nb_reca_groupwise)
        {
        rcoef2=_ParamRecalageBspline.serie_groupwise[l]->rcoeff*_ParamRecalageBspline.serie_groupwise[l]->rcoeff;
        imx_gradient_3d_p(_ParamRecalageBspline.serie_groupwise[l],grad_tmp,2,4);
      
        for(i=0; i<wdth; i++)
        for(j=0; j<hght; j++)
        for(k=0; k<dpth; k++)
            {
            tmp=_ParamRecalageBspline.serie_groupwise[l]->mri[i][j][k]*rcoef2;
            _ParamRecalageBspline.grad_ref2_groupwise->raw[i][j][k].x+=grad_tmp->raw[i][j][k].x*tmp;
            _ParamRecalageBspline.grad_ref2_groupwise->raw[i][j][k].y+=grad_tmp->raw[i][j][k].y*tmp;
            _ParamRecalageBspline.grad_ref2_groupwise->raw[i][j][k].z+=grad_tmp->raw[i][j][k].z*tmp;
            }
     
        
        }
     
     tmp=2.0;
     mul_field_3d(_ParamRecalageBspline.grad_ref2_groupwise, _ParamRecalageBspline.grad_ref2_groupwise, tmp);
    
     
      free_field3d(grad_tmp);*/
      }

 
 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;


 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 if ((dist==Energie_atrophie_log_jacobien_locale_3d)|| (dist==Energie_atrophie_jacobien_locale_3d)) precis=pow(0.1,PRECISION_SIMULATION)*precis; // on augmente un peu la precision pour la simulation d'atrophie
 
  
 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
    #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. debut resolution..... "); 
    #endif
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
    #endif

E=0;
 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;

for(pari=MINI(2,resol)-1;pari>=0;pari--)
for(parj=MINI(2,resol)-1;parj>=0;parj--)
for(park=MINI(2,resol)-1;park>=0;park--)
    {
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
            
        //---------------------------------------------------------------
        //----- initialisation Slpqr            -------------------------
        //---------------------------------------------------------------
        l=TOP_conv_ind(topi,topj,topk,nb_param);
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

        /* Calcul du gradient*/
        
        gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
        
                
        
            
        topD1 = l;
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
        grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];   
            
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
            {
            if (dist==Energie_IM_locale_3d)
                init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
                    
                
            E2=TOP_linemin_locale_3d(imref,imreca,imres,champ_fin,the_transf,inter,dist,regularisation,nb_param,p,param_norm,grad,topi,topj,topk,Jmin,Jmax,Slpqr);
            E=E+E2; 
        
            
            if (dist==Energie_IM_locale_3d)
             init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            }
                //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
            
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
             stop++;
        }
      
    }
    
    for (i=0;i<nb_param;i++) p_i[i]=p[i];
    
    printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
  nb++;
 } while (nb<50 && Eprec!=E);

 /*calcul de l'energie finale*/
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];


  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         /*si l'extension .trf existe on la supprime*/
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 /*liberation memoire des variables*/
 //free_dvector(p,nb_param);
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);
if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
 free_field3d(gradIref);

 if (dist==Energie_groupwise_variance_locale_3d)
    {
     free_field3d(_ParamRecalageBspline.grad_reca_groupwise);
     free_field3d(_ParamRecalageBspline.grad_ref_groupwise);
     free_field3d(_ParamRecalageBspline.grad_ref2_groupwise);
     }


 free_transf3d(transfo);


 if  ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        free_field3d(maskgradIref);
 
 free_imatrix_3d(masque_bloc);
 return(Efin);
}

/*******************************************************************************
**     TOP_min_icm_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     recalage par ICM                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_icm_locale_3d    (grphic3d *imref, grphic3d *imreca, grphic3d *imres,
                                    field3d *champ_fin, transf_func_t the_transf, InterpolationFct inter,
                                dist_func_locale_t dist, reg_func_locale_t regularisation, int nb_param, double *param, int ***masque_param,
                                                            double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i,j; 
 double *grad,*p,*p_i,*param_norm;
 double Edebut,Efin,E1,E2,E,Eprec;
 int    nb,nbreduc,stop;
 double prreduc,maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;   
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt,topD1;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11,l;
 double xm,xM,ym,yM,zm,zM;   
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif

    
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
  
  
 /*allocation memoire des variables*/
 //p=alloc_dvector(nb_param);
 p_i=alloc_dvector(nb_param);
 param_norm=alloc_dvector(nb_param);
 
 grad=alloc_dvector(nb_param);
 masque_bloc=alloc_imatrix_3d(topD,topD,topD);

 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


/*Initialisation de masque_bloc*/
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;
    
  for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];
 
 /* initialisation des parametres normalises sur [0,1] */
    for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }
 
//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    }   
 /*calcul de l'energie de depart*/
  Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  E2=Edebut;

 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;

 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
    #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. debut resolution..... "); 
    #endif
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
    #endif

E=0;
 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
/*for (pari=0;pari<MINI(2,resol);pari++)
 for (parj=0;parj<MINI(2,resol);parj++)
  for (park=0;park<MINI(2,resol);park++)*/
for(pari=MINI(2,resol)-1;pari>=0;pari--)
    for(parj=MINI(2,resol)-1;parj>=0;parj--)
  for(park=MINI(2,resol)-1;park>=0;park--)
    
    {
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
        
        //---------------------------------------------------------------
        //----- initialisation Slpqr            -------------------------
        //---------------------------------------------------------------
            l=TOP_conv_ind(topi,topj,topk,nb_param);
            
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

        /* Calcul du gradient*/
        topD1 = l;;
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
        
        if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);


                
        for (i=0;i<3;i++)
            for (j=-1;j<2;j=j+2)
            {
            grad[topD1]=0;grad[topD1+1]=0;grad[topD1+2]=0;    
            grad[topD1+i]=j;
            E2=TOP_linemin_locale_3d(imref,imreca,imres,champ_fin,the_transf,inter,dist,regularisation,nb_param,p,param_norm,grad,topi,topj,topk,Jmin,Jmax,Slpqr);
            E=E+E2; 
            }

        if (dist==Energie_IM_locale_3d)
         init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        
                //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
            /*for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
            for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
            for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
                masque_bloc[i][j][k]=1;*/
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
            stop++;
        }
      
    }
    
    for (i=0;i<nb_param;i++) p_i[i]=p[i];
    
    printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
        
  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE 
    aff_log("Verif. fin iteration ..... "); 
    #endif
    TOP_verification_all(nb_param,p,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE 
    aff_log("finie \n");
    #endif
    }
  #endif

    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
  nb++;
 } while (nb<50 && Eprec!=E);

 /*calcul de l'energie finale*/
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         /*si l'extension .trf existe on la supprime*/
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 /*liberation memoire des variables*/
 //free_dvector(p,nb_param);
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);
 free_imatrix_3d(masque_bloc);
 free_transf3d(transfo);

 return(Efin);
}


double TOP_min_marquardt_locale_3d_bloc(grphic3d *imref, grphic3d *imreca, grphic3d *imres, double *p, double *p_i, double *param_norm, double *opp_param_norm, int nb_param, double *grad, dist_func_locale_t dist, reg_func_locale_t regularisation, reg_grad_locale_t gradient_reg, field3d *gradIref, field3d *maskgradIref, reg_hessien_locale_t hessien_reg,
 int (*gradient)(),  
 int (*func_hessien)(),  
 hessien3d *hessienIref, hessien3d *maskhessienIref, double Jmin, double Jmax, int ***masque_param,  int ***masque_bloc, int topi, int topj, int topk)
{
int topD1,i,j,itmp,imin;
double xm,xM,ym,yM,zm,zM,p_ini[3], Etmp,Emin,gradl[3],p_tmp[3], E;
int wdth,hght,dpth; 
TSlpqr Slpqr[8];
int *x00,*x11,*y00,*y11,*z00,*z11;
int res_inv=0,continu,nb_iter,tutu;
double nuk,ck=10,ck1,seuil_nuk;
mat3d pseudo_hessien,hessien,inv_hessien;
double energie_precision = 0.01;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
double maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;    
int nb_pas = NB_PAS_LINEMIN_TOPO; 
int resol=BASE3D.resol;
int topD=(int)pow(2.0,1.0*resol)-1;
double max_delta_p=0, delta_p;  
            

x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
wdth=imref->width;hght=imref->height;dpth=imref->depth;
 
 // --- On fixe la précision à 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la précision en recalage symétrique soit équivalent que celle obtenue en non symetrique

 if ((dist==Energie_atrophie_log_jacobien_locale_3d)||(dist==Energie_atrophie_jacobien_locale_3d)) {precis=pow(0.1,PRECISION_SIMULATION)*precis;}
 
 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 

    topD1=TOP_conv_ind(topi,topj,topk,nb_param);
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
        
            
                
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
        
        for (i=0;i<3;i++)
            p_ini[i]=p[topD1+i];
    
        if (dist==Energie_IM_locale_3d)
                init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);

        
        if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
    
    
    
        continu=0;
        nuk=0;
        ck=10;
        ck1=10;
        seuil_nuk=0;
        nb_iter=0;
        tutu=1;
        if (Etmp==0) continu=1;
        do{
        Emin=Etmp;
        if (tutu>0)
            {
            // Calcul du gradient

    
            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            func_hessien(imref, imreca, nb_param,p,param_norm,gradIref,maskgradIref,hessienIref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
    
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=hessien.H[i][j];
        
        if ((nuk==0)&&(seuil_nuk==0))
            {
            nuk=0.01*sqrt(hessien.H[0][0]*hessien.H[0][0]+hessien.H[1][1]*hessien.H[1][1]+hessien.H[2][2]*hessien.H[2][2]);
            seuil_nuk=100.0*sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);
            
            
                if ((seuil_nuk<100*nuk)||(nuk<0.0001))
                {
                nuk=0.01;seuil_nuk=1000000; //on fixe des valeurs arbitraires
                }
            }
            }   
            
        

            
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
        {       
        // Calcul du pseudo hessien Hk=H+nuk*I jusqu'a ce qu'il soit defini positif
        res_inv=0;
        while (res_inv==0)
            {
                for (i=0;i<3;i++)
                    pseudo_hessien.H[i][i]=hessien.H[i][i]+nuk; 
            
                res_inv=inversion_mat3d(&pseudo_hessien,&inv_hessien);
                
                if(res_inv==0)
                    nuk=nuk*ck;
            }
        
        if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        
        for (i=0;i<3;i++)
            p_tmp[i]=p[topD1+i];
                
            max_delta_p=0;  
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    {
                    delta_p=inv_hessien.H[i][j]*gradl[j];
                    max_delta_p=MAXI(max_delta_p,fabs(delta_p));
                    p[topD1+i]=p[topD1+i]-delta_p;
                    }
        
        
        if (max_delta_p>TOPO_QUANTIFICATION_PARAM)  // on ne calcul le critere que si la modif des params est superieur a TOPO_QUANTIFICATION_PARAM
            {
            for (i=0;i<3;i++)  
                    if (fabs(p[topD1+i])<TOPO_QUANTIFICATION_PARAM)  /* C'est pour en finir avec les bug num�riques ... */
                        p[topD1+i]=0.0;
        
        
        
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;    

                    
                                        
                    Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
            }
        else
            Etmp=HUGE_VAL;
    
    
        if (Etmp>Emin)
                {
                for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
                    
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                nuk=nuk*ck; 
                Etmp=Emin;
                if (nuk>seuil_nuk) {continu=1;}
                tutu=0;
                nb_iter=0;  
                }
        else{
                {
                nuk=1.0*nuk/ck1; tutu=2;
                nb_iter++;
                if (((1.0*fabs((Emin-Etmp)/Emin))<energie_precision)||(nb_iter>10)||(Emin==0))
                    {continu=1;Emin=Etmp;}
                }
                                
                }
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
                    
        }
        else
        {
        continu=1;
        }           
        }while (continu==0);
        
        
        // V�rification de la conservation de la topologie pour les nouveaux param
        

        if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
         || ((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
        {
        for(i=0;i<nb_param;i++)
            opp_param_norm[i]=-1.0*lambda_ref*param_norm[i];
        
        }


        if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)||(res_inv==-1)
                        ||(TOP_verification_locale(nb_param,p,opp_param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0))
        {
            for (i=0;i<3;i++)
                p_tmp[i]=p[topD1+i];


        if (dist==Energie_IM_locale_3d)
            init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            Emin=HUGE_VAL;imin=0;
            continu=1;
            for (itmp=0;((itmp<nb_pas)&&(continu==1));itmp++)
                {
                
                for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i]+1.0*itmp/nb_pas*(p_tmp[i]-p_ini[i]);
            
            
                for (i=0;i<3;i++)
                    if (fabs(p[topD1+i])<TOPO_QUANTIFICATION_PARAM) /* C'est pour en finir avec les bug num�riques ... */
                        p[topD1+i]=0.0;
        
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    

                Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
    
                if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
                ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
                    {
                    for(i=0;i<3;i++)
                            opp_param_norm[topD1+i]=-1.0*lambda_ref*param_norm[topD1+i];
        
                    }

        
                if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)
                    ||(TOP_verification_locale(nb_param,p,opp_param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0))
                 continu=0;
                else if (Etmp<Emin) {Emin=Etmp;imin=itmp;}
                }
                for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i]+1.0*imin/nb_pas*(p_tmp[i]-p_ini[i]);
            
                
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                
                if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
                ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
                    {
                    for(i=0;i<3;i++)
                            opp_param_norm[topD1+i]=-1.0*lambda_ref*param_norm[topD1+i];
        
                    }

                
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);


        }
        
        if (Emin==HUGE_VAL)
            {//printf("c'est un bug issu des pb numerique de topologie ...\n");
            
            if (dist==Energie_IM_locale_3d)
                init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    
                
                Emin=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
        
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            }
        E=E+Emin;
        
    
        //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
        

return(Emin);
}   
/*******************************************************************************
**     TOP_min_marquardt_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     Minimisation par la m�thode de Levenberg-Marquardt                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_marquardt_locale_3d_parallel (grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,
                                            transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,
                                                                        int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i; 
 field3d *gradIref,*maskgradIref=NULL;
 double *grad,*p,*p_i,*param_norm;
 double Edebut,Efin,E1,E2,E,Eprec,Etmp;
 int    nb,nbreduc,stop,nb_tmp;
 double prreduc,precis,precisx,precisy,precisz; 
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11;
 int (*gradient)();  
 int (*func_hessien)();  
 reg_grad_locale_t gradient_reg;  
 reg_hessien_locale_t hessien_reg;
 hessien3d *hessienIref,*maskhessienIref=NULL; 
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif
double energie_precision = 0.01;
double *opp_param_norm=NULL;
int nb_iter_max=50;

 if (dist==Energie_IM_locale_3d) energie_precision=0.0001;  /* les variations relatives de l'IM sont en g�n�ral bcp plus faibles */


 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
 

 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) gradient=gradient_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) gradient=gradient_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) gradient=gradient_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) gradient=gradient_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) gradient=gradient_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) gradient=gradient_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) gradient=gradient_quad_locale_symetrique_coupe_3d;
 if (dist==Energie_quad_sous_champ_locale_3d) gradient=gradient_quad_sous_champ_locale_3d;
 if (dist==Energie_groupwise_variance_locale_3d) gradient=gradient_groupwise_variance_locale_3d;
 if (dist==Energie_groupwise_variance_locale_nonsym_3d) gradient=gradient_groupwise_variance_locale_nonsym_3d;
 
 func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) func_hessien=hessien_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) func_hessien=hessien_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) func_hessien=hessien_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) func_hessien=hessien_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) func_hessien=hessien_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) func_hessien=hessien_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) func_hessien=hessien_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) func_hessien=hessien_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) func_hessien=hessien_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) func_hessien=hessien_quad_locale_symetrique_coupe_3d;
 if (dist==Energie_quad_sous_champ_locale_3d) func_hessien=hessien_quad_sous_champ_locale_3d;
 if (dist==Energie_groupwise_variance_locale_3d) func_hessien=hessien_groupwise_variance_locale_3d;
 if (dist==Energie_groupwise_variance_locale_nonsym_3d) func_hessien=hessien_groupwise_variance_locale_nonsym_3d;

 
 gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) gradient_reg=gradient_regularisation_log_jacobien_local;
 if (regularisation==regularisation_dist_identite_local) gradient_reg=gradient_regularisation_dist_identite_local;
 if (regularisation==regularisation_patch_imagebased_local) gradient_reg=gradient_regularisation_patch_imagebased_local;

 hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) hessien_reg=hessien_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) hessien_reg=hessien_regularisation_log_jacobien_local;
 if (regularisation==regularisation_dist_identite_local) hessien_reg=hessien_regularisation_dist_identite_local;
 if (regularisation==regularisation_patch_imagebased_local) hessien_reg=hessien_regularisation_patch_imagebased_local;

 
    //allocation memoire des variables
    //p=alloc_dvector(nb_param);
    p_i=alloc_dvector(nb_param);
    param_norm=alloc_dvector(nb_param);
 
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
    ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
    opp_param_norm=alloc_dvector(nb_param);
    
 
 
    grad=alloc_dvector(nb_param);
    //gradIref=cr_field3d(wdth,hght,dpth);   // pour economiser de la place memoire
    gradIref=champ_fin;
    masque_bloc=alloc_imatrix_3d(topD,topD,topD);

    
    if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
    hessienIref=cr_hessien3d(wdth,hght,dpth);
    else
    hessienIref=NULL;
    
    if (dist==Energie_groupwise_variance_locale_3d)
        {
        _ParamRecalageBspline.grad_reca_groupwise=cr_field3d(wdth,hght,dpth);
        _ParamRecalageBspline.grad_ref_groupwise=cr_field3d(wdth,hght,dpth);
        _ParamRecalageBspline.grad_ref2_groupwise=cr_field3d(wdth,hght,dpth);
        
        _ParamRecalageBspline.hessien_reca_groupwise=cr_hessien3d(wdth,hght,dpth);
        _ParamRecalageBspline.hessien_ref_groupwise=cr_hessien3d(wdth,hght,dpth);
        _ParamRecalageBspline.hessien_ref2_groupwise=cr_hessien3d(wdth,hght,dpth);
        }
        
        
    
    
     if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        {   
        maskhessienIref=cr_hessien3d(wdth,hght,dpth);
        maskgradIref=cr_field3d(wdth,hght,dpth);
        }
        
 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


//Initialisation de masque_bloc
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;

for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

// initialisation des parametres normalises sur [0,1] 
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 
        
 //calcul de l'energie de depart
    Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  E2=Edebut;

 // calcul du gradient de l'image a recaler 
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d))
        imx_gradient_3d_p(imreca,gradIref,0,4);
 else
if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
    imx_gradient_3d_p(imreca,gradIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_gradient_3d_p(imreca->mask,maskgradIref,0,4);
    
 
 //--- Pour le recalage symetrique, on calcule le gradient de imref -----
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        imx_gradient_3d_p(imref,maskgradIref,2,4);
 
 
 // calcul du Hessien de l'image 
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        imx_hessien_3d_p(imreca,hessienIref,0,4);
 else
  if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
        imx_hessien_3d_p(imreca,hessienIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_hessien_3d_p(imreca->mask,maskhessienIref,0,4);


 //--- Pour le recalage symetrique, on calcule le hessien de imref -----
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
    imx_hessien_3d_p(imref,maskhessienIref,2,4);


//---- Pour le recalage groupwise, calcul des gradients et hessiens
if (dist==Energie_groupwise_variance_locale_3d)
        {
        imx_gradient_3d_p(_ParamRecalageBspline.reca_groupwise,_ParamRecalageBspline.grad_reca_groupwise,2,4);
        mul_field_3d(_ParamRecalageBspline.grad_reca_groupwise,_ParamRecalageBspline.grad_reca_groupwise, _ParamRecalageBspline.reca_groupwise->rcoeff);

        imx_gradient_3d_p(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.grad_ref_groupwise,2,4);
        mul_field_3d(_ParamRecalageBspline.grad_ref_groupwise,_ParamRecalageBspline.grad_ref_groupwise, _ParamRecalageBspline.ref_groupwise->rcoeff);
        
        imx_gradient_3d_p(_ParamRecalageBspline.ref2_groupwise,_ParamRecalageBspline.grad_ref2_groupwise,2,4);
        mul_field_3d(_ParamRecalageBspline.grad_ref2_groupwise,_ParamRecalageBspline.grad_ref2_groupwise, _ParamRecalageBspline.ref2_groupwise->rcoeff);

        imx_hessien_3d_p(_ParamRecalageBspline.reca_groupwise,_ParamRecalageBspline.hessien_reca_groupwise,2,4);
        mul_hessien_3d(_ParamRecalageBspline.hessien_reca_groupwise, _ParamRecalageBspline.hessien_reca_groupwise, _ParamRecalageBspline.reca_groupwise->rcoeff);

        imx_hessien_3d_p(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.hessien_ref_groupwise,2,4);
        mul_hessien_3d(_ParamRecalageBspline.hessien_ref_groupwise, _ParamRecalageBspline.hessien_ref_groupwise, _ParamRecalageBspline.ref_groupwise->rcoeff);

        imx_hessien_3d_p(_ParamRecalageBspline.ref2_groupwise,_ParamRecalageBspline.hessien_ref2_groupwise,2,4);
        mul_hessien_3d(_ParamRecalageBspline.hessien_ref2_groupwise, _ParamRecalageBspline.hessien_ref2_groupwise, _ParamRecalageBspline.ref2_groupwise->rcoeff);
        
        }



 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);


 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;
 
 
 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 if ((dist==Energie_atrophie_log_jacobien_locale_3d)||(dist==Energie_atrophie_jacobien_locale_3d)) {precis=pow(0.1,PRECISION_SIMULATION)*precis; nb_iter_max=200;}
 
 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
 

E=0;

 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
    nb_tmp=0;


for(pari=MINI(2,resol)-1;pari>=0;pari--)
    for(parj=MINI(2,resol)-1;parj>=0;parj--)
  for(park=MINI(2,resol)-1;park>=0;park--)
    {
 #pragma omp parallel for private(topi,topj,topk,Etmp) shared(E) schedule(dynamic)
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
                
        Etmp=TOP_min_marquardt_locale_3d_bloc(imref, imreca,imres, p, p_i, param_norm, opp_param_norm,  nb_param, grad,  dist,  regularisation,  gradient_reg, gradIref, maskgradIref,  hessien_reg, gradient,func_hessien, hessienIref, maskhessienIref, Jmin, Jmax, masque_param,  masque_bloc, topi, topj, topk);
                    
        E=E+Etmp;
        }
#pragma end parallel for

    stop++;
    }
        for (i=0;i<nb_param;i++) p_i[i]=p[i];
        printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
        
  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE 
    aff_log("Verif. fin iteration ..... "); 
    #endif
    TOP_verification_all(nb_param,p,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE 
    aff_log("finie \n");
    #endif
    }
  #endif
    
    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
    
            
  nb++;
 } while (nb<nb_iter_max && Eprec!=E);

 //calcul de l'energie finale
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 //liberation memoire des variables
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);

 if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
 ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
    free_dvector(opp_param_norm,nb_param);

 free_imatrix_3d(masque_bloc);
 
 if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d)) free_hessien3d(hessienIref);

 free_transf3d(transfo);

 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
     {free_hessien3d(maskhessienIref);
        free_field3d(maskgradIref);
        }


if (dist==Energie_groupwise_variance_locale_3d)
        {
        free_field3d(_ParamRecalageBspline.grad_reca_groupwise);
        free_field3d(_ParamRecalageBspline.grad_ref_groupwise);
        free_field3d(_ParamRecalageBspline.grad_ref2_groupwise);
        
        free_hessien3d(_ParamRecalageBspline.hessien_reca_groupwise);
        free_hessien3d(_ParamRecalageBspline.hessien_ref_groupwise);
        free_hessien3d(_ParamRecalageBspline.hessien_ref2_groupwise);
        }
        


 return(Efin);
}


/*******************************************************************************
**     TOP_min_marquardt_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     Minimisation par la méthode de Levenberg-Marquardt                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_marquardt_locale_3d  (grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,
                                            transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,
                                                                        int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i,j/*,k*/; 
 field3d *gradIref,*maskgradIref=NULL;
 double *grad,*p,*p_i,gradl[3],*param_norm,p_ini[3],p_tmp[3];
 double Edebut,Efin,E1,E2,E,Eprec,Emin,Etmp;
 int    nb,nbreduc,stop,nb_iter,nb_tmp;
 double prreduc,maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;   
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt,topD1;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11;
 int (*gradient)();  
 int (*func_hessien)();  
 reg_grad_locale_t gradient_reg;  
 reg_hessien_locale_t hessien_reg;
 hessien3d *hessienIref,*maskhessienIref=NULL; 
 mat3d pseudo_hessien,hessien,inv_hessien;
 //double cfmax;
 double xm,xM,ym,yM,zm,zM;
 int res_inv=0,continu;
 double nuk,ck=10,ck1,seuil_nuk;
 int tutu;
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif
int nb_pas = NB_PAS_LINEMIN_TOPO,itmp,imin; 
double energie_precision = 0.01;
double *opp_param_norm=NULL;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;
int nb_iter_max=50;
double max_delta_p=0, delta_p;  
            
 if (dist==Energie_IM_locale_3d) energie_precision=0.0001;  /* les variations relatives de l'IM sont en général bcp plus faibles */

 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
 

 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) gradient=gradient_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) gradient=gradient_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) gradient=gradient_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) gradient=gradient_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) gradient=gradient_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) gradient=gradient_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) gradient=gradient_quad_locale_symetrique_coupe_3d;
 if (dist==Energie_quad_sous_champ_locale_3d) gradient=gradient_quad_sous_champ_locale_3d;
 if (dist==Energie_groupwise_variance_locale_3d) gradient=gradient_groupwise_variance_locale_3d;
 if (dist==Energie_groupwise_variance_locale_nonsym_3d) gradient=gradient_groupwise_variance_locale_nonsym_3d;
 
 func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) func_hessien=hessien_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) func_hessien=hessien_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) func_hessien=hessien_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) func_hessien=hessien_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) func_hessien=hessien_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) func_hessien=hessien_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) func_hessien=hessien_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) func_hessien=hessien_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) func_hessien=hessien_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) func_hessien=hessien_quad_locale_symetrique_coupe_3d;
 if (dist==Energie_quad_sous_champ_locale_3d) func_hessien=hessien_quad_sous_champ_locale_3d;
 if (dist==Energie_groupwise_variance_locale_3d) func_hessien=hessien_groupwise_variance_locale_3d;
 if (dist==Energie_groupwise_variance_locale_nonsym_3d) func_hessien=hessien_groupwise_variance_locale_nonsym_3d;

 
 gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) gradient_reg=gradient_regularisation_log_jacobien_local;
 if (regularisation==regularisation_dist_identite_local) gradient_reg=gradient_regularisation_dist_identite_local;
 if (regularisation==regularisation_patch_imagebased_local) gradient_reg=gradient_regularisation_patch_imagebased_local;

 hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) hessien_reg=hessien_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) hessien_reg=hessien_regularisation_log_jacobien_local;
 if (regularisation==regularisation_dist_identite_local) hessien_reg=hessien_regularisation_dist_identite_local;
 if (regularisation==regularisation_patch_imagebased_local) hessien_reg=hessien_regularisation_patch_imagebased_local;

 
    //allocation memoire des variables
    //p=alloc_dvector(nb_param);
    p_i=alloc_dvector(nb_param);
    param_norm=alloc_dvector(nb_param);
 
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
    ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
    opp_param_norm=alloc_dvector(nb_param);
    
 
 
    grad=alloc_dvector(nb_param);
    //gradIref=cr_field3d(wdth,hght,dpth);   // pour economiser de la place memoire
    gradIref=champ_fin;
    masque_bloc=alloc_imatrix_3d(topD,topD,topD);

    
    if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
    hessienIref=cr_hessien3d(wdth,hght,dpth);
    else
    hessienIref=NULL;
    
    if (dist==Energie_groupwise_variance_locale_3d)
        {
        _ParamRecalageBspline.grad_reca_groupwise=cr_field3d(wdth,hght,dpth);
        _ParamRecalageBspline.grad_ref_groupwise=cr_field3d(wdth,hght,dpth);
        _ParamRecalageBspline.grad_ref2_groupwise=cr_field3d(wdth,hght,dpth);
        
        _ParamRecalageBspline.hessien_reca_groupwise=cr_hessien3d(wdth,hght,dpth);
        _ParamRecalageBspline.hessien_ref_groupwise=cr_hessien3d(wdth,hght,dpth);
        _ParamRecalageBspline.hessien_ref2_groupwise=cr_hessien3d(wdth,hght,dpth);
        }
        
        
    
    
     if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        {   
        maskhessienIref=cr_hessien3d(wdth,hght,dpth);
        maskgradIref=cr_field3d(wdth,hght,dpth);
        }
        
 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


//Initialisation de masque_bloc
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;

for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

// initialisation des parametres normalises sur [0,1] 
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 
        
 //calcul de l'energie de depart
    Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  E2=Edebut;

 // calcul du gradient de l'image a recaler 
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d))
        imx_gradient_3d_p(imreca,gradIref,0,4);
 else
if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
    imx_gradient_3d_p(imreca,gradIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_gradient_3d_p(imreca->mask,maskgradIref,0,4);
    
 
 //--- Pour le recalage symetrique, on calcule le gradient de imref -----
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        imx_gradient_3d_p(imref,maskgradIref,2,4);
 
 
 // calcul du Hessien de l'image 
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        imx_hessien_3d_p(imreca,hessienIref,0,4);
 else
  if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d))
        imx_hessien_3d_p(imreca,hessienIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_hessien_3d_p(imreca->mask,maskhessienIref,0,4);


 //--- Pour le recalage symetrique, on calcule le hessien de imref -----
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
    imx_hessien_3d_p(imref,maskhessienIref,2,4);


//---- Pour le recalage groupwise, calcul des gradients et hessiens
if (dist==Energie_groupwise_variance_locale_3d)
        {
        imx_gradient_3d_p(_ParamRecalageBspline.reca_groupwise,_ParamRecalageBspline.grad_reca_groupwise,2,4);
        mul_field_3d(_ParamRecalageBspline.grad_reca_groupwise,_ParamRecalageBspline.grad_reca_groupwise, _ParamRecalageBspline.reca_groupwise->rcoeff);

        imx_gradient_3d_p(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.grad_ref_groupwise,2,4);
        mul_field_3d(_ParamRecalageBspline.grad_ref_groupwise,_ParamRecalageBspline.grad_ref_groupwise, _ParamRecalageBspline.ref_groupwise->rcoeff);
        
        imx_gradient_3d_p(_ParamRecalageBspline.ref2_groupwise,_ParamRecalageBspline.grad_ref2_groupwise,2,4);
        mul_field_3d(_ParamRecalageBspline.grad_ref2_groupwise,_ParamRecalageBspline.grad_ref2_groupwise, _ParamRecalageBspline.ref2_groupwise->rcoeff);

        imx_hessien_3d_p(_ParamRecalageBspline.reca_groupwise,_ParamRecalageBspline.hessien_reca_groupwise,2,4);
        mul_hessien_3d(_ParamRecalageBspline.hessien_reca_groupwise, _ParamRecalageBspline.hessien_reca_groupwise, _ParamRecalageBspline.reca_groupwise->rcoeff);

        imx_hessien_3d_p(_ParamRecalageBspline.ref_groupwise,_ParamRecalageBspline.hessien_ref_groupwise,2,4);
        mul_hessien_3d(_ParamRecalageBspline.hessien_ref_groupwise, _ParamRecalageBspline.hessien_ref_groupwise, _ParamRecalageBspline.ref_groupwise->rcoeff);

        imx_hessien_3d_p(_ParamRecalageBspline.ref2_groupwise,_ParamRecalageBspline.hessien_ref2_groupwise,2,4);
        mul_hessien_3d(_ParamRecalageBspline.hessien_ref2_groupwise, _ParamRecalageBspline.hessien_ref2_groupwise, _ParamRecalageBspline.ref2_groupwise->rcoeff);
        
        }



 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);


 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;
 
 
 // --- On fixe la précision à 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la précision en recalage symétrique soit équivalent que celle obtenue en non symetrique

 if ((dist==Energie_atrophie_log_jacobien_locale_3d)||(dist==Energie_atrophie_jacobien_locale_3d)) {precis=pow(0.1,PRECISION_SIMULATION)*precis; nb_iter_max=200;}
 
 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
 

E=0;

 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
    nb_tmp=0;


for(pari=MINI(2,resol)-1;pari>=0;pari--)
    for(parj=MINI(2,resol)-1;parj>=0;parj--)
  for(park=MINI(2,resol)-1;park>=0;park--)
    
    {
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
                
        topD1=TOP_conv_ind(topi,topj,topk,nb_param);
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
        
            
                
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
        
        for (i=0;i<3;i++)
            p_ini[i]=p[topD1+i];
    
        if (dist==Energie_IM_locale_3d)
                init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);

        
        if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
    
    
    
        continu=0;
        nuk=0;
        ck=10;
        ck1=10;
        seuil_nuk=0;
        nb_iter=0;
        tutu=1;
        if (Etmp==0) continu=1;
        do{
        Emin=Etmp;
        if (tutu>0)
            {
            // Calcul du gradient

    
            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            func_hessien(imref, imreca, nb_param,p,param_norm,gradIref,maskgradIref,hessienIref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
    
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=hessien.H[i][j];
        
        if ((nuk==0)&&(seuil_nuk==0))
            {
            nuk=0.01*sqrt(hessien.H[0][0]*hessien.H[0][0]+hessien.H[1][1]*hessien.H[1][1]+hessien.H[2][2]*hessien.H[2][2]);
            seuil_nuk=100.0*sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);
            
            
                if ((seuil_nuk<100*nuk)||(nuk<0.0001))
                {
                nuk=0.01;seuil_nuk=1000000; //on fixe des valeurs arbitraires
                }
            }
            }   
            
        

            
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
        {       
        // Calcul du pseudo hessien Hk=H+nuk*I jusqu'a ce qu'il soit defini positif
        res_inv=0;
        while (res_inv==0)
            {
                for (i=0;i<3;i++)
                    pseudo_hessien.H[i][i]=hessien.H[i][i]+nuk; 
            
                res_inv=inversion_mat3d(&pseudo_hessien,&inv_hessien);
                
                if(res_inv==0)
                    nuk=nuk*ck;
            }
        
        if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        
        for (i=0;i<3;i++)
            p_tmp[i]=p[topD1+i];
                
            max_delta_p=0;  
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    {
                    delta_p=inv_hessien.H[i][j]*gradl[j];
                    max_delta_p=MAXI(max_delta_p,fabs(delta_p));
                    p[topD1+i]=p[topD1+i]-delta_p;
                    }
        
        
        if (max_delta_p>TOPO_QUANTIFICATION_PARAM)  // on ne calcul le critere que si la modif des params est superieur a TOPO_QUANTIFICATION_PARAM
            {
            for (i=0;i<3;i++)  
                    if (fabs(p[topD1+i])<TOPO_QUANTIFICATION_PARAM)  /* C'est pour en finir avec les bug numériques ... */
                        p[topD1+i]=0.0;
        
        
        
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;    

                    
                                        
                    Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
            }
        else
            Etmp=HUGE_VAL;
    
    
        if (Etmp>Emin)
                {
                for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
                    
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                nuk=nuk*ck; 
                Etmp=Emin;
                if (nuk>seuil_nuk) {continu=1;}
                tutu=0;
                nb_iter=0;  
                }
        else{
                {
                nuk=1.0*nuk/ck1; tutu=2;
                nb_iter++;
                if (((1.0*fabs((Emin-Etmp)/Emin))<energie_precision)||(nb_iter>10)||(Emin==0))
                    {continu=1;Emin=Etmp;}
                }
                                
                }
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
                    
        }
        else
        {
        continu=1;
        }           
        }while (continu==0);
        
        
        // Vérification de la conservation de la topologie pour les nouveaux param
        

        if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
         || ((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
        {
        for(i=0;i<nb_param;i++)
            opp_param_norm[i]=-1.0*lambda_ref*param_norm[i];
        
        }


        if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)||(res_inv==-1)
                        ||(TOP_verification_locale(nb_param,p,opp_param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0))
        {
            for (i=0;i<3;i++)
                p_tmp[i]=p[topD1+i];


        if (dist==Energie_IM_locale_3d)
            init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            Emin=HUGE_VAL;imin=0;
            continu=1;
            for (itmp=0;((itmp<nb_pas)&&(continu==1));itmp++)
                {
                
                for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i]+1.0*itmp/nb_pas*(p_tmp[i]-p_ini[i]);
            
            
                for (i=0;i<3;i++)
                    if (fabs(p[topD1+i])<TOPO_QUANTIFICATION_PARAM) /* C'est pour en finir avec les bug numériques ... */
                        p[topD1+i]=0.0;
        
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    

                Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
    
                if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
                ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
                    {
                    for(i=0;i<3;i++)
                            opp_param_norm[topD1+i]=-1.0*lambda_ref*param_norm[topD1+i];
        
                    }

        
                if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)
                    ||(TOP_verification_locale(nb_param,p,opp_param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0))
                 continu=0;
                else if (Etmp<Emin) {Emin=Etmp;imin=itmp;}
                }
                for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i]+1.0*imin/nb_pas*(p_tmp[i]-p_ini[i]);
            
                
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                
                if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
                ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
                    {
                    for(i=0;i<3;i++)
                            opp_param_norm[topD1+i]=-1.0*lambda_ref*param_norm[topD1+i];
        
                    }

                
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);


        }
        
        if (Emin==HUGE_VAL)
            {//printf("c'est un bug issu des pb numerique de topologie ...\n");
            
            if (dist==Energie_IM_locale_3d)
                init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    
                
                Emin=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
        
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            }
        E=E+Emin;
        
    
        //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
        stop++;
        
        }
        
                          
    }
    
        for (i=0;i<nb_param;i++) p_i[i]=p[i];
        printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
        
  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE 
    aff_log("Verif. fin iteration ..... "); 
    #endif
    TOP_verification_all(nb_param,p,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE 
    aff_log("finie \n");
    #endif
    }
  #endif
    
    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
    
            
  nb++;
 } while (nb<nb_iter_max && Eprec!=E);

 //calcul de l'energie finale
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 //liberation memoire des variables
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);

 if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d)||(dist==Energie_groupwise_variance_locale_3d)
 ||((dist==Energie_quad_sous_champ_locale_3d)&&(TOPO_PONDERATION_RECALAGE_SYMETRIQUE > 0.0)))
    free_dvector(opp_param_norm,nb_param);

 free_imatrix_3d(masque_bloc);
 
 if ((dist!=Energie_quad_sous_champ_locale_3d)&&(dist!=Energie_groupwise_variance_locale_3d)) free_hessien3d(hessienIref);

 free_transf3d(transfo);

 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
     {free_hessien3d(maskhessienIref);
        free_field3d(maskgradIref);
        }


if (dist==Energie_groupwise_variance_locale_3d)
        {
        free_field3d(_ParamRecalageBspline.grad_reca_groupwise);
        free_field3d(_ParamRecalageBspline.grad_ref_groupwise);
        free_field3d(_ParamRecalageBspline.grad_ref2_groupwise);
        
        free_hessien3d(_ParamRecalageBspline.hessien_reca_groupwise);
        free_hessien3d(_ParamRecalageBspline.hessien_ref_groupwise);
        free_hessien3d(_ParamRecalageBspline.hessien_ref2_groupwise);
        }
        


 return(Efin);
}








/*******************************************************************************
**     TOP_min_marquardt_locale_lent_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     Minimisation par la m�thode de Levenberg-Marquardt                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_marquardt_locale_lent_3d (grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,
                                            transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,
                                                                        int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i,j/*,k*/; 
 field3d *gradIref,*maskgradIref=NULL;
 double *grad,*p,*p_i,gradl[3],*param_norm,p_ini[3],p_tmp[3];
 double Edebut,Efin,E1,E2,E,Eprec,Emin,Etmp;
 int    nb,nbreduc,stop,nb_iter,nb_tmp;
 double prreduc,maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;   
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt,topD1;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11;
 int (*gradient)();  
 int (*func_hessien)();  
 reg_grad_locale_t gradient_reg;  
 reg_hessien_locale_t hessien_reg;
 hessien3d *hessienIref,*maskhessienIref=NULL; 
 mat3d pseudo_hessien,hessien,inv_hessien;
 //double cfmax;
 double xm,xM,ym,yM,zm,zM;
 int res_inv=0,continu,pas;
 double nuk,ck=10,ck1,seuil_nuk;
 int tutu;
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif
int nb_pas = NB_PAS_LINEMIN_TOPO,itmp,imin; 
double energie_precision = 0.01;
double *opp_param_norm=NULL;
double lambda_ref = (1.0-TOPO_PONDERATION_RECALAGE_SYMETRIQUE)/TOPO_PONDERATION_RECALAGE_SYMETRIQUE;

 if (dist==Energie_IM_locale_3d) energie_precision=0.0001;  /* les variations relatives de l'IM sont en g�n�ral bcp plus faibles */

 
 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
 

 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) gradient=gradient_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) gradient=gradient_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) gradient=gradient_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) gradient=gradient_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) gradient=gradient_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) gradient=gradient_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) gradient=gradient_quad_locale_symetrique_coupe_3d;
 
 func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) func_hessien=hessien_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) func_hessien=hessien_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) func_hessien=hessien_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) func_hessien=hessien_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) func_hessien=hessien_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) func_hessien=hessien_base_ICP_sym_locale_3d;
 if (dist==Energie_atrophie_jacobien_locale_3d) func_hessien=hessien_atrophie_jacobien_locale_3d;
 if (dist==Energie_atrophie_log_jacobien_locale_3d) func_hessien=hessien_atrophie_log_jacobien_locale_3d;
 if (dist==Energie_quad_locale_symetrique_3d) func_hessien=hessien_quad_locale_symetrique_3d;
 if (dist==Energie_quad_locale_symetrique_coupe_3d) func_hessien=hessien_quad_locale_symetrique_coupe_3d;
 
 
  gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) gradient_reg=gradient_regularisation_log_jacobien_local;

 hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) hessien_reg=hessien_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 
    //allocation memoire des variables
    //p=alloc_dvector(nb_param);
    p_i=alloc_dvector(nb_param);
    param_norm=alloc_dvector(nb_param);
 
 if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
    opp_param_norm=alloc_dvector(nb_param);
    
 
 
    grad=alloc_dvector(nb_param);
  //gradIref=cr_field3d(wdth,hght,dpth);   // pour economiser de la place memoire
    gradIref=champ_fin;
    masque_bloc=alloc_imatrix_3d(topD,topD,topD);

    
    hessienIref=cr_hessien3d(wdth,hght,dpth);

     if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        {   
        maskhessienIref=cr_hessien3d(wdth,hght,dpth);
        maskgradIref=cr_field3d(wdth,hght,dpth);
        }
        
 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


//Initialisation de masque_bloc
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;

for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

// initialisation des parametres normalises sur [0,1] 
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 
        
 //calcul de l'energie de depart
    Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  E2=Edebut;

 // calcul du gradient de l'image a recaler 
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d))
        imx_gradient_3d_p(imreca,gradIref,0,4);
 else
        imx_gradient_3d_p(imreca,gradIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_gradient_3d_p(imreca->mask,maskgradIref,0,4);
    
 
 //--- Pour le recalage symetrique, on calcule le gradient de imref -----
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        imx_gradient_3d_p(imref,maskgradIref,2,4);
 
 
 // calcul du Hessien de l'image 
 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_ICP_locale_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        imx_hessien_3d_p(imreca,hessienIref,0,4);
 else
        imx_hessien_3d_p(imreca,hessienIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_hessien_3d_p(imreca->mask,maskhessienIref,0,4);


 //--- Pour le recalage symetrique, on calcule le hessien de imref -----
    if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
    imx_hessien_3d_p(imref,maskhessienIref,2,4);


 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);


 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;

E=0;

//E=Edebut;
 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
    nb_tmp=0;

/*for (pari=0;pari<MINI(2,resol);pari++)
 for (parj=0;parj<MINI(2,resol);parj++)
  for (park=0;park<MINI(2,resol);park++)*/
/*if (dist==Energie_geman_locale_3d)
    update_TOPO_ALPHA_ROBUST_global(imref,imreca, nb_param, p);*/
for(pas=0; pas<11; pas++)
    {
    
 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 

for(pari=MINI(2,resol)-1;pari>=0;pari--)
    for(parj=MINI(2,resol)-1;parj>=0;parj--)
  for(park=MINI(2,resol)-1;park>=0;park--)
    
    {
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
            
            /*borne_jacobien_bspline1_local_3d( nb_param,  p, param_norm, topi,  topj,  topk, &min1, &max1);
            regularisation_energie_membrane_local(nb_param,  p, topi,  topj,  topk, &min2, &max2);

            if (((min1-min2)/min1>0.0001)||((max1-max2)/max1>0.0001))
                    printf ("Aie un bug : min1 : %f   min2 : %f   max1 : %f   max2 : %f \n",min1,min2,max1,max2);
            else
                printf ("OK : min1 : %f   min2 : %f   max1 : %f   max2 : %f \n",min1,min2,max1,max2);*/
                
            topD1=TOP_conv_ind(topi,topj,topk,nb_param);
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
        
            
                
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
        
            for (i=0;i<3;i++)
                p_ini[i]=p[topD1+i];
    
    if (dist==Energie_IM_locale_3d)
                init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);

        
    if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
    
    
    
        continu=0;
        nuk=0;
        ck=10;
        ck1=10;
        seuil_nuk=0;
        nb_iter=0;
        tutu=1;
        if (Etmp==0) continu=1;
        do{
        Emin=Etmp;
        if (tutu>0)
            {
            // Calcul du gradient

    
            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            func_hessien(imref, imreca, nb_param,p,param_norm,gradIref,maskgradIref,hessienIref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
    
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=hessien.H[i][j];
        
        if ((nuk==0)&&(seuil_nuk==0))
            {
            nuk=0.01*sqrt(hessien.H[0][0]*hessien.H[0][0]+hessien.H[1][1]*hessien.H[1][1]+hessien.H[2][2]*hessien.H[2][2]);
            seuil_nuk=100.0*sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);
            /*normG=sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);
            normH=sqrt(hessien.H[0][0]*hessien.H[0][0]+hessien.H[1][1]*hessien.H[1][1]+hessien.H[2][2]*hessien.H[2][2]);
            nuk=MINI(pow(0.001*normG,2)/normH,0.01*normH);
            seuil_nuk=MAXI(pow(1000*normG,2)/normH,100*normG);*/
            
            
                if ((seuil_nuk<100*nuk)||(nuk<0.0001))
                {
                nuk=0.01;seuil_nuk=1000000; //on fixe des valeurs arbitraires
                }
            }
            }   
            
        

            
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
        {       
        // Calcul du pseudo hessien Hk=H+nuk*I jusqu'a ce qu'il soit defini positif
        res_inv=0;
        while (res_inv==0)
            {
                for (i=0;i<3;i++)
                    pseudo_hessien.H[i][i]=hessien.H[i][i]+nuk; 
            
                res_inv=inversion_mat3d(&pseudo_hessien,&inv_hessien);
                
                if(res_inv==0)
                    nuk=nuk*ck;
            }
        
        if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        
        for (i=0;i<3;i++)
            p_tmp[i]=p[topD1+i];
                
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    p[topD1+i]=p[topD1+i]-0.1*pas*inv_hessien.H[i][j]*gradl[j];
        
            for (i=0;i<3;i++)  
                    if (fabs(p[topD1+i])<0.000001)  /* C'est pour en finir avec les bug num�riques ... */
                        p[topD1+i]=0.0;
        
        
        
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;    

                    
                                        
        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
    
    
    /*if (isinf(p[topD1])||isinf(p[topD1+1])||isinf(p[topD1+2]))
        {printf("Y a un bug \n");
        for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    printf("hessien.H[%d][%d]=%f\n",i,j,hessien.H[i][j]);
                    
                        for (i=0;i<3;i++)
                        printf("grad[%d]=%f\n",i,gradl[i]);
                    
                    }*/
    
        if (Etmp>Emin)
                {
                for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
                    
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                nuk=nuk*ck; 
                Etmp=Emin;
                if (nuk>seuil_nuk) {continu=1;}
                tutu=0;
                nb_iter=0;  
                }
        else{
                    /*if (TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)
                        {
                        for (i=0;i<3;i++)
                            p[topD1+i]=p_tmp[i];
                    
                        param_norm[topD1]  =p[topD1]  /wdth;    
                        param_norm[topD1+1]=p[topD1+1]/hght;    
                        param_norm[topD1+2]=p[topD1+2]/dpth;
                        Etmp=Emin;
                        continu=1;
                        res_inv=-1;
                        }
                    else*/
                        {
                        nuk=1.0*nuk/ck1; tutu=2;
                        nb_iter++;
                            if (((1.0*fabs((Emin-Etmp)/Emin))<energie_precision)||(nb_iter>0)||(Emin==0))
                                {continu=1;Emin=Etmp;}
                                
                            }
                                
                }
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
                    
        }
        else
        {
        continu=1;
        }           
        }while (continu==0);
        
        
        // V�rification de la conservation de la topologie pour les nouveaux param
        
/*      if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)||(res_inv==-1))
            {
            // on fait une recherche lin�aire suivant la descente de gradient
            for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
    
            param_norm[topD1]  =p[topD1]  /wdth;    
            param_norm[topD1+1]=p[topD1+1]/hght;    
            param_norm[topD1+2]=p[topD1+2]/dpth;    
        
            gradient(imref,imreca,nb_param,p,param_norm,gradIref,gradl,topi,topj,topk,Slpqr);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            //printf("\n On fait de la descente de gradient !!!! %f \n",nuk);
            
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
                Emin=TOP_linemin_locale_3d(imref,imreca,imres,champ_fin,the_transf,inter,dist,nb_param,p,param_norm,grad,topi,topj,topk,Jmin,Jmax,Slpqr);
            nb_tmp++;   
    
            }
*/

 if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
        {
        for(i=0;i<nb_param;i++)
            opp_param_norm[i]=-1.0*lambda_ref*param_norm[i];
        
        }


if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)||(res_inv==-1)
                        ||(TOP_verification_locale(nb_param,p,opp_param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0))
        {
            for (i=0;i<3;i++)
                p_tmp[i]=p[topD1+i];


    if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            Emin=HUGE_VAL;imin=0;
            continu=1;
            for (itmp=0;((itmp<nb_pas)&&(continu==1));itmp++)
                {
                
                for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i]+1.0*itmp/nb_pas*(p_tmp[i]-p_ini[i]);
            
            
            for (i=0;i<3;i++)
                    if (fabs(p[topD1+i])<0.000001) /* C'est pour en finir avec les bug num�riques ... */
                        p[topD1+i]=0.0;
        
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    

                Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
    
                if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
                    {
                    for(i=0;i<3;i++)
                            opp_param_norm[topD1+i]=-1.0*lambda_ref*param_norm[topD1+i];
        
                    }

        
                if ((TOP_verification_locale(nb_param,p,param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0)
                    ||(TOP_verification_locale(nb_param,p,opp_param_norm,Jmin,Jmax,topi,topj,topk,Slpqr)==0))
                 continu=0;
                else if (Etmp<Emin) {Emin=Etmp;imin=itmp;}
                }
                for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i]+1.0*imin/nb_pas*(p_tmp[i]-p_ini[i]);
            
                
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                
                if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
                    {
                    for(i=0;i<3;i++)
                            opp_param_norm[topD1+i]=-1.0*lambda_ref*param_norm[topD1+i];
        
                    }

                
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);


        }
        
        if (Emin==HUGE_VAL)
            {//printf("c'est un bug issu des pb numerique de topologie ...\n");
            
            if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    
                
                Emin=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
        
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            }
        E=E+Emin;
        
            //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
            /*for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
            for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
            for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
                masque_bloc[i][j][k]=1;*/
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
        stop++;
        
        }
    }
        
                          
    }
    
        for (i=0;i<nb_param;i++) p_i[i]=p[i];
        printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
        
  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE 
    aff_log("Verif. fin iteration ..... "); 
    #endif
    TOP_verification_all(nb_param,p,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE 
    aff_log("finie \n");
    #endif
    }
  #endif
    
    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
    
/*          base_to_field_3d(nb_param,p,champ_fin,imref,imreca);
            inter_labeled_3d(imreca->mask,champ_fin,imres); 
            save_mri_ipb_3d_p("test_anim",imres);*/
            
  nb++;
 } while (nb<50 && Eprec!=E);

 //calcul de l'energie finale
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 //liberation memoire des variables
 //free_dvector(p,nb_param);
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);

 if ((dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
    free_dvector(opp_param_norm,nb_param);

 //free_field3d(gradIref);
 free_imatrix_3d(masque_bloc);
 free_hessien3d(hessienIref);

 free_transf3d(transfo);

 if ((dist==Energie_ICP_sym_locale_3d)||(dist==Energie_quad_locale_symetrique_3d)||(dist==Energie_quad_locale_symetrique_coupe_3d))
     {free_hessien3d(maskhessienIref);
        free_field3d(maskgradIref);
        }

 return(Efin);
}

/*******************************************************************************
**     TOP_min_marquardt_penal_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     Minimisation par la m�thode de Levenberg-Marquardt                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_marquardt_penal_locale_3d    (grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,
                                                        transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist, reg_func_locale_t regularisation,
                                                                                    int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i,j; 
 field3d *gradIref,*maskgradIref=NULL;
 double *grad,*p,*p_i,gradl[3],*param_norm,p_ini[3],p_tmp[3];
 double Edebut,Efin,E1,E2,E,Eprec,Emin,Etmp;
 int    nb,nbreduc,stop,nb_iter,nb_tmp;
 double prreduc,maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;   
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt,topD1;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11;
 int (*gradient)();  
 reg_grad_locale_t gradient_reg;  
 int (*func_hessien)();  
 reg_hessien_locale_t hessien_reg;
 hessien3d *hessienIref,*maskhessienIref=NULL; 
 mat3d pseudo_hessien,hessien,inv_hessien;
 //double cfmax;
 double xm,xM,ym,yM,zm,zM;
 int res_inv,continu;
 double nuk,ck=100,ck1=10,seuil_nuk;
 int tutu;
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif
double energie_precision = 0.01;

 if (dist==Energie_IM_locale_3d) energie_precision=0.0001;  /* les variations relatives de l'IM sont en g�n�ral bcp plus faibles */

 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
  
 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) gradient=gradient_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) gradient=gradient_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) gradient=gradient_base_ICP_sym_locale_3d;


 func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) func_hessien=hessien_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) func_hessien=hessien_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) func_hessien=hessien_base_geman_locale_3d;
 if (dist==Energie_IM_locale_3d) func_hessien=hessien_base_IM_locale_3d;
 if (dist==Energie_ICP_locale_3d) func_hessien=hessien_base_ICP_locale_3d;
 if (dist==Energie_ICP_sym_locale_3d) func_hessien=hessien_base_ICP_sym_locale_3d;

  gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) gradient_reg=gradient_regularisation_log_jacobien_local;

 hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) hessien_reg=hessien_regularisation_energie_membrane_Lp_local;
 if (regularisation==regularisation_log_jacobien_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 

    //allocation memoire des variables
    p=alloc_dvector(nb_param);
    p_i=alloc_dvector(nb_param);
    param_norm=alloc_dvector(nb_param);
 
    grad=alloc_dvector(nb_param);
  //gradIref=cr_field3d(wdth,hght,dpth);   // pour economiser de la place memoire
    gradIref=champ_fin;
    masque_bloc=alloc_imatrix_3d(topD,topD,topD);
    
    
    hessienIref=cr_hessien3d(wdth,hght,dpth);

 if (dist==Energie_ICP_sym_locale_3d)
    {maskhessienIref=cr_hessien3d(wdth,hght,dpth);
    maskgradIref=cr_field3d(wdth,hght,dpth);  
    }

 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


//Initialisation de masque_bloc
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;

            
for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

// initialisation des parametres normalises sur [0,1] 
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 
        
 //calcul de l'energie de depart
    Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  //Etmp=Energie_globale_sym_3d(imref,imreca,nb_param,p,param_norm,dist,masque_param);
  E2=Edebut;

 // calcul du gradient de l'image a recaler 
 imx_gradient_3d_p(imreca,gradIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_gradient_3d_p(imreca->mask,maskgradIref,2,4);
    
 
 // calucl du Hessiend e l'image 
    imx_hessien_3d_p(imreca,hessienIref,2,4);

 if (dist==Energie_ICP_sym_locale_3d)
    imx_hessien_3d_p(imreca->mask,maskhessienIref,2,4);

 aff_log("Edebut: %.2f    nb_param: %d \n",Edebut,nb_param);

 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;

 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
    
E=0;
 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
    nb_tmp=0;
/*for (pari=0;pari<MINI(2,resol);pari++)
 for (parj=0;parj<MINI(2,resol);parj++)
  for (park=0;park<MINI(2,resol);park++)*/
for(pari=MINI(2,resol)-1;pari>=0;pari--)
    for(parj=MINI(2,resol)-1;parj>=0;parj--)
  for(park=MINI(2,resol)-1;park>=0;park--)
    
    {
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
            
            topD1=TOP_conv_ind(topi,topj,topk,nb_param);
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
                
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
        
            for (i=0;i<3;i++)
                p_ini[i]=p[topD1+i];
    
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;
    
        if (dist==Energie_IM_locale_3d)
            init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);
        
        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
        
        if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

    
        continu=0;
        nuk=0;
        ck=10;
        ck1=10;
        seuil_nuk=0;
        nb_iter=0;
        tutu=1;
        if (Etmp==0) continu=1;
        do{
        Emin=Etmp;
        if (tutu>0)
            {
            // Calcul du gradient

            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            func_hessien(imref, imreca, nb_param,p,param_norm,gradIref,maskgradIref,hessienIref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=hessien.H[i][j];
            
            // Estimation de l'espace de recherche de nuk
         if ((nuk==0)&&(seuil_nuk==0))
            {
            nuk=0.01*sqrt(hessien.H[0][0]*hessien.H[0][0]+hessien.H[1][1]*hessien.H[1][1]+hessien.H[2][2]*hessien.H[2][2]);
            seuil_nuk=100.0*sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);
                if ((seuil_nuk<100*nuk)||(nuk<0.0001))
                {
                nuk=0.01;seuil_nuk=1000000; //on fixe des valeurs arbitraires
                }
            }
            
            }   
            
        
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
        {
        // Calcul du pseudo hessien Hk=H+nuk*I jusqu'a ce qu'il soit defini positif
        res_inv=0; 
         
        while (res_inv==0)
            {
                for (i=0;i<3;i++)
                    pseudo_hessien.H[i][i]=hessien.H[i][i]+nuk; 
            
                res_inv=inversion_mat3d(&pseudo_hessien,&inv_hessien);
                
                if(res_inv==0)
                    nuk=nuk*ck;
            }
        
        if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

        
        for (i=0;i<3;i++)
            p_tmp[i]=p[topD1+i];
                
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    p[topD1+i]=p[topD1+i]-inv_hessien.H[i][j]*gradl[j];
        
        
        for (i=0;i<3;i++)
                    if (fabs(p[topD1+i])<0.000001) /* C'est pour en finir avec les bug num�riques ... */
                        p[topD1+i]=0.0;
        
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;    
        
        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
      
    /*  borne_jacobien_bspline1_local_3d(nb_param,p,param_norm, topi,  topj,  topk, &minJ, &maxJ);
        if ((minJ<0)&&(Etmp<HUGE_VAL))
            printf("violation de la topologie discrete\n");*/
        /*if (Etmp==HUGE_VAL)
            {
            px=p[topD1];
            py=p[topD1+1];
            pz=p[topD1+2];
            
            for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
            
            param_norm[topD1]  =p[topD1]  /wdth;    
            param_norm[topD1+1]=p[topD1+1]/hght;    
            param_norm[topD1+2]=p[topD1+2]/dpth;
            
            p[topD1]=px;
            param_norm[topD1]  =p[topD1]  /wdth;
            Etmpx=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr);
            
            p[topD1]=p_tmp[0];
            param_norm[topD1]  =p[topD1]  /wdth;
            
            p[topD1+1]=py;
            param_norm[topD1+1]  =p[topD1+1]  /hght;
            Etmpy=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr);
            
            p[topD1+1]=p_tmp[1];
            param_norm[topD1+1]  =p[topD1+1]  /hght;
            
            p[topD1+2]=pz;
            param_norm[topD1+2]  =p[topD1+2]  /dpth;
            Etmpz=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr);
            
            p[topD1+2]=p_tmp[2];
            param_norm[topD1+2]  =p[topD1+2]  /dpth;
            
            if ((Etmpx<Etmpy)&&(Etmpx<Etmpz))
                {Etmp=Etmpx;
                p[topD1]=px;
                param_norm[topD1]  =p[topD1]  /wdth;
                }
            
            if ((Etmpy<Etmpx)&&(Etmpy<Etmpz))
                {Etmp=Etmpy;
                p[topD1+1]=py;
                param_norm[topD1+1]  =p[topD1+1]  /hght;
                }
                
            if ((Etmpz<Etmpx)&&(Etmpz<Etmpy))
                {Etmp=Etmpz;
                p[topD1+2]=pz;
                param_norm[topD1+2]  =p[topD1+2]  /dpth;
                }
            }*/
                
        if (Etmp>Emin)
                {
                for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
                    
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                nuk=nuk*ck; 
                Etmp=Emin;
                if (nuk>seuil_nuk) {continu=1;}
                tutu=0;
                nb_iter=0;  
                }
        else{
                        {
                        nuk=1.0*nuk/ck1; tutu=2;
                        nb_iter++;
                            if (((1.0*(Emin-Etmp)/Emin)<energie_precision)||(nb_iter>10)||(Emin==0))
                                {continu=1;Emin=Etmp;}
                        }
                                
                }
            if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

            }
        else
            {continu=1;}
                    
        }while (continu==0);
        
        
        
        if (Emin==HUGE_VAL)
            {for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    
        
        if (dist==Energie_IM_locale_3d)
        init_IM_locale_before_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);

                Emin=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);

        if (dist==Energie_IM_locale_3d)
            init_IM_locale_after_optimization_3d(imref, imreca, nb_param, p, topi, topj, topk);


            }
        E=E+Emin;
        
        
                //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
            /*for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
            for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
            for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
                masque_bloc[i][j][k]=1;*/
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
        stop++;
        }
        
            
              
    }
    
    for (i=0;i<nb_param;i++) p_i[i]=p[i];
    
    printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
//  printf("Descente de gradient sur %f %% des blocs \n",100.0*nb_tmp/stop);
        
  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE 
    aff_log("Verif. fin iteration ..... "); 
    #endif
    TOP_verification_all(nb_param,p,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE 
    aff_log("finie \n");
    #endif
    }
  #endif
    
    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
    
/*          base_to_field_3d(nb_param,p,champ_fin,imref,imreca);
            inter_labeled_3d(imreca->mask,champ_fin,imres); 
            save_mri_ipb_3d_p("test_anim",imres);*/
            
  nb++;
 } while (nb<50 && Eprec!=E);

 //calcul de l'energie finale
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 //liberation memoire des variables
 free_dvector(p,nb_param);
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);
 //free_field3d(gradIref);
 free_imatrix_3d(masque_bloc);
 free_hessien3d(hessienIref);

 if (dist==Energie_ICP_sym_locale_3d)
    { free_hessien3d(maskhessienIref);
        free_field3d(maskgradIref);
    }
 return(Efin);
}

/*******************************************************************************
**     TOP_min_marquardt_penal_locale_gradient_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     Minimisation par la m�thode de Levenberg-Marquardt                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_marquardt_penal_locale_gradient_3d   (grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,
                                                                        transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist,reg_func_locale_t regularisation,
                                                                                                    int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i,j,k; 
 field3d *gradIref,*gradIreca,*gradGxref,*gradGyref,*gradGzref,*maskgradIref=NULL;
 double *grad,*p,*p_i,gradl[3],*param_norm,p_ini[3],p_tmp[3];
 double Edebut,Efin,E1,E2,E,Eprec,Emin,Etmp;
 int    nb,nbreduc,stop,nb_iter,nb_tmp;
 double prreduc,maxpx,maxpy,maxpz,precis,precisx,precisy,precisz;   
 int topi,topj,topk,topD,resol;
 int pari,parj,park;
 int compt,topD1;
 int ***masque_bloc;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11;
 int (*gradient)();  
 reg_grad_locale_t gradient_reg;  
 int (*func_hessien)();  
 reg_hessien_locale_t hessien_reg;
 hessien3d *hessienIref,*hessienGxref,*hessienGyref,*hessienGzref,*maskhessienIref=NULL; 
 mat3d pseudo_hessien,hessien,inv_hessien;
 //double cfmax;
 double xm,xM,ym,yM,zm,zM;
 int res_inv,continu;
 double nuk,ck=100,ck1=10,seuil_nuk;
 int tutu;
 grphic3d *gxref,*gyref,*gzref,*gxreca,*gyreca,*gzreca;
 double alpha=0.5;
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif

 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
  
 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 
 func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) func_hessien=hessien_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) func_hessien=hessien_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
if (dist==Energie_geman_locale_3d) func_hessien=hessien_base_geman_locale_3d;
 
 
  gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;

 hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) hessien_reg=hessien_regularisation_energie_membrane_Lp_local;
 
    //allocation memoire des variables
    p=alloc_dvector(nb_param);
    p_i=alloc_dvector(nb_param);
    param_norm=alloc_dvector(nb_param);
 
    grad=alloc_dvector(nb_param);
  //gradIref=cr_field3d(wdth,hght,dpth);   // pour economiser de la place memoire
    gradIref=champ_fin;
    gradIreca=cr_field3d(wdth,hght,dpth);
    gradGxref=cr_field3d(wdth,hght,dpth);
    gradGyref=cr_field3d(wdth,hght,dpth);
    gradGzref=cr_field3d(wdth,hght,dpth);
    
    masque_bloc=alloc_imatrix_3d(topD,topD,topD);

    hessienIref=cr_hessien3d(wdth,hght,dpth);
    hessienGxref=cr_hessien3d(wdth,hght,dpth);
    hessienGyref=cr_hessien3d(wdth,hght,dpth);
    hessienGzref=cr_hessien3d(wdth,hght,dpth);


gxref=cr_grphic3d(imref);
gyref=cr_grphic3d(imref);
gzref=cr_grphic3d(imref);
gxreca=cr_grphic3d(imreca);
gyreca=cr_grphic3d(imreca);
gzreca=cr_grphic3d(imreca);

 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


//Initialisation de masque_bloc
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;

for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

// initialisation des parametres normalises sur [0,1] 
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 
        

 // calcul du gradient de l'image a recaler 
 imx_gradient_3d_p(imreca,gradIref,2,4);
 imx_gradient_3d_p(imref,gradIreca,2,4);
 
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
    {
    gxref->mri[i][j][k]=(int)gradIreca->raw[i][j][k].x;
    gyref->mri[i][j][k]=(int)gradIreca->raw[i][j][k].y;
    gzref->mri[i][j][k]=(int)gradIreca->raw[i][j][k].z;
    gxreca->mri[i][j][k]=(int)gradIref->raw[i][j][k].x;
    gyreca->mri[i][j][k]=(int)gradIref->raw[i][j][k].y;
    gzreca->mri[i][j][k]=(int)gradIref->raw[i][j][k].z;
    }     
 
 imx_gradient_3d_p(gxreca,gradGxref,2,4);
 imx_gradient_3d_p(gyreca,gradGyref,2,4);
 imx_gradient_3d_p(gzreca,gradGzref,2,4);
 
    imx_hessien_3d_p(gxreca,hessienGxref,2,4);
    imx_hessien_3d_p(gyreca,hessienGyref,2,4);
    imx_hessien_3d_p(gzreca,hessienGzref,2,4);

 // calucl du Hessiend e l'image 
 imx_hessien_3d_p(imreca,hessienIref,2,4);

 //calcul de l'energie de depart
    Edebut=(1.0-alpha)*Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  Edebut=Edebut+alpha*Energie_globale_3d(gxref,gxreca,nb_param,p,param_norm,dist,regularisation,masque_param)/3.0;
  Edebut=Edebut+alpha*Energie_globale_3d(gyref,gyreca,nb_param,p,param_norm,dist,regularisation,masque_param)/3.0;
  Edebut=Edebut+alpha*Energie_globale_3d(gzref,gzreca,nb_param,p,param_norm,dist,regularisation,masque_param)/3.0;
  
    E2=Edebut;


 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;

 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
    
E=0;
 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
    nb_tmp=0;
/*for (pari=0;pari<MINI(2,resol);pari++)
 for (parj=0;parj<MINI(2,resol);parj++)
  for (park=0;park<MINI(2,resol);park++)*/
for(pari=MINI(2,resol)-1;pari>=0;pari--)
    for(parj=MINI(2,resol)-1;parj>=0;parj--)
  for(park=MINI(2,resol)-1;park>=0;park--)
    
    {
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
       if ((topi%2==pari)&&(topj%2==parj)&&(topk%2==park))
        if ((masque_param[topi][topj][topk]!=0)&&(masque_bloc[topi][topj][topk]!=0))
        {
            
            topD1=TOP_conv_ind(topi,topj,topk,nb_param);
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
            
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
        
            for (i=0;i<3;i++)
                p_ini[i]=p[topD1+i];
    
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;
                    
        Etmp=(1.0-alpha)*dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
        Etmp=Etmp+alpha*dist(gxref,gxreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
        Etmp=Etmp+alpha*dist(gyref,gyreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
        Etmp=Etmp+alpha*dist(gzref,gzreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
        
            
        continu=0;
        nuk=0;
        ck=10;
        ck1=10;
        seuil_nuk=0;
        nb_iter=0;
        tutu=1;
        if (Etmp==0) continu=1;
        do{
        Emin=Etmp;
        if (tutu>0)
            {
            // Calcul du gradient

            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=(1.0-alpha)*gradl[0];grad[topD1+1]=(1.0-alpha)*gradl[1];grad[topD1+2]=(1.0-alpha)*gradl[2];    
            
            gradient(gxref,gxreca,nb_param,p,param_norm,gradGxref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=grad[topD1]+alpha*gradl[0]/3.0;grad[topD1+1]=grad[topD1+1]+alpha*gradl[1]/3.0;grad[topD1+2]=grad[topD1+2]+alpha*gradl[2]/3.0;
            
            gradient(gyref,gyreca,nb_param,p,param_norm,gradGyref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=grad[topD1]+alpha*gradl[0]/3.0;grad[topD1+1]=grad[topD1+1]+alpha*gradl[1]/3.0;grad[topD1+2]=grad[topD1+2]+alpha*gradl[2]/3.0;
            
            gradient(gzref,gzreca,nb_param,p,param_norm,gradGzref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=grad[topD1]+alpha*gradl[0]/3.0;grad[topD1+1]=grad[topD1+1]+alpha*gradl[1]/3.0;grad[topD1+2]=grad[topD1+2]+alpha*gradl[2]/3.0;
            
            gradl[0]=grad[topD1];gradl[1]=grad[topD1+1];gradl[2]=grad[topD1+2];
            
            
            
            func_hessien(imref, imreca, nb_param,p,param_norm,gradIref,maskgradIref,hessienIref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=(1.0-alpha)*hessien.H[i][j];
            
            
            func_hessien(gxref, gxreca, nb_param,p,param_norm,gradGxref,maskgradIref,hessienGxref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=pseudo_hessien.H[i][j]+alpha*hessien.H[i][j]/3.0;
            
            func_hessien(gyref, gyreca, nb_param,p,param_norm,gradGyref,maskgradIref,hessienGyref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=pseudo_hessien.H[i][j]+alpha*hessien.H[i][j]/3.0;
            
            func_hessien(gzref, gzreca, nb_param,p,param_norm,gradGzref,maskgradIref,hessienGzref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=pseudo_hessien.H[i][j]+alpha*hessien.H[i][j]/3.0;
            
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    hessien.H[i][j]=pseudo_hessien.H[i][j];
            
            // Estimation de l'espace de recherche de nuk
         if ((nuk==0)&&(seuil_nuk==0))
            {
            nuk=0.01*sqrt(hessien.H[0][0]*hessien.H[0][0]+hessien.H[1][1]*hessien.H[1][1]+hessien.H[2][2]*hessien.H[2][2]);
            seuil_nuk=100.0*sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);
                if ((seuil_nuk<100*nuk)||(nuk<0.0001))
                {
                nuk=0.01;seuil_nuk=1000000; //on fixe des valeurs arbitraires
                }
            }
            
            }   
            
        
        if((gradl[0]!=0)||(gradl[1]!=0)||(gradl[2]!=0))
        {
        // Calcul du pseudo hessien Hk=H+nuk*I jusqu'a ce qu'il soit defini positif
        res_inv=0; 
         
        while (res_inv==0)
            {
                for (i=0;i<3;i++)
                    pseudo_hessien.H[i][i]=hessien.H[i][i]+nuk; 
            
                res_inv=inversion_mat3d(&pseudo_hessien,&inv_hessien);
                
                if(res_inv==0)
                    nuk=nuk*ck;
            }
        
        
        for (i=0;i<3;i++)
            p_tmp[i]=p[topD1+i];
                
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    p[topD1+i]=p[topD1+i]-inv_hessien.H[i][j]*gradl[j];
                
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;    
        
        Etmp=(1.0-alpha)*dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);  
        Etmp=Etmp+alpha*dist(gxref,gxreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
        Etmp=Etmp+alpha*dist(gyref,gyreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
        Etmp=Etmp+alpha*dist(gzref,gzreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;

        if (isnan(Etmp))
            Etmp=HUGE_VAL;
                            
    /*  borne_jacobien_bspline1_local_3d(nb_param,p,param_norm, topi,  topj,  topk, &minJ, &maxJ);
        if ((minJ<0)&&(Etmp<HUGE_VAL))
            printf("violation de la topologie discrete\n");*/
        /*if (Etmp==HUGE_VAL)
            {
            px=p[topD1];
            py=p[topD1+1];
            pz=p[topD1+2];
            
            for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
            
            param_norm[topD1]  =p[topD1]  /wdth;    
            param_norm[topD1+1]=p[topD1+1]/hght;    
            param_norm[topD1+2]=p[topD1+2]/dpth;
            
            p[topD1]=px;
            param_norm[topD1]  =p[topD1]  /wdth;
            Etmpx=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr);
            
            p[topD1]=p_tmp[0];
            param_norm[topD1]  =p[topD1]  /wdth;
            
            p[topD1+1]=py;
            param_norm[topD1+1]  =p[topD1+1]  /hght;
            Etmpy=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr);
            
            p[topD1+1]=p_tmp[1];
            param_norm[topD1+1]  =p[topD1+1]  /hght;
            
            p[topD1+2]=pz;
            param_norm[topD1+2]  =p[topD1+2]  /dpth;
            Etmpz=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr);
            
            p[topD1+2]=p_tmp[2];
            param_norm[topD1+2]  =p[topD1+2]  /dpth;
            
            if ((Etmpx<Etmpy)&&(Etmpx<Etmpz))
                {Etmp=Etmpx;
                p[topD1]=px;
                param_norm[topD1]  =p[topD1]  /wdth;
                }
            
            if ((Etmpy<Etmpx)&&(Etmpy<Etmpz))
                {Etmp=Etmpy;
                p[topD1+1]=py;
                param_norm[topD1+1]  =p[topD1+1]  /hght;
                }
                
            if ((Etmpz<Etmpx)&&(Etmpz<Etmpy))
                {Etmp=Etmpz;
                p[topD1+2]=pz;
                param_norm[topD1+2]  =p[topD1+2]  /dpth;
                }
            }*/
                
        if (Etmp>Emin)
                {
                for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
                    
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                nuk=nuk*ck; 
                Etmp=Emin;
                if (nuk>seuil_nuk) {continu=1;}
                tutu=0;
                nb_iter=0;  
                }
        else{
                        {
                        nuk=1.0*nuk/ck1; tutu=2;
                        nb_iter++;
                            if (((1.0*(Emin-Etmp)/Emin)<0.01)||(nb_iter>10))
                                {continu=1;Emin=Etmp;}
                        }
                                
                }
            
            }
        else
            {continu=1;}
                    
        }while (continu==0);
        
        
        
        if (Emin==HUGE_VAL)
            {for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    
                
                Emin=(1.0-alpha)*dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
                Emin=Emin+alpha*dist(gxref,gxreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
                Emin=Emin+alpha*dist(gyref,gyreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
                Emin=Emin+alpha*dist(gzref,gzreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation)/3.0;
                
            }
        E=E+Emin;
        
                //-------Evaluation du deplacement max------------
        maxpx=1.0*fabs(p_i[topD1]-p[topD1]);
        maxpy=1.0*fabs(p_i[topD1+1]-p[topD1+1]);
        maxpz=1.0*fabs(p_i[topD1+2]-p[topD1+2]);
        
        if ((maxpx<precisx)&&(maxpy<precisy)&&(maxpz<precisz)) masque_bloc[topi][topj][topk]=0;
        else {/* on debloque les blocs voisins */
            /*for (i=MAXI(topi-1,0);i<=MINI(topi+1,topD-1);i++)
            for (j=MAXI(topj-1,0);j<=MINI(topj+1,topD-1);j++)
            for (k=MAXI(topk-1,0);k<=MINI(topk+1,topD-1);k++)
                masque_bloc[i][j][k]=1;*/
                if (maxpx>precisx) 
                    {   if (topi>0) masque_bloc[topi-1][topj][topk]=1;
                        if (topi<topD-1) masque_bloc[topi+1][topj][topk]=1;
                    }
                if (maxpy>precisy) 
                    {   if (topj>0) masque_bloc[topi][topj-1][topk]=1;
                        if (topj<topD-1) masque_bloc[topi][topj+1][topk]=1;
                    }
                if (maxpz>precisz) 
                    {   if (topk>0) masque_bloc[topi][topj][topk-1]=1;
                        if (topk<topD-1) masque_bloc[topi][topj][topk+1]=1;
                    }
             }
        stop++;
        }
        
            
              
    }
    
    for (i=0;i<nb_param;i++) p_i[i]=p[i];
    
    printf("Minimisation sur %f %% des blocs Eprec : %f    E : %f\n",300.0*stop/nb_param,Eprec,E);
//  printf("Descente de gradient sur %f %% des blocs \n",100.0*nb_tmp/stop);
        
  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE 
    aff_log("Verif. fin iteration ..... "); 
    #endif
    TOP_verification_all(nb_param,p,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE 
    aff_log("finie \n");
    #endif
    }
  #endif
    
    if (nomfichiertrf!=NULL)
    {
    #ifndef TOPO_COMPARE_CRITERE 
    compare_critere(imref,imreca,imres,champ_fin,p,param_norm,nb_param,masque_param,nb);
    #endif
    }
    
/*          base_to_field_3d(nb_param,p,champ_fin,imref,imreca);
            inter_labeled_3d(imreca->mask,champ_fin,imres); 
            save_mri_ipb_3d_p("test_anim",imres);*/
            
  nb++;
 } while (nb<50 && Eprec!=E);

 //calcul de l'energie finale
  Efin=(1.0-alpha)*Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  Efin=Efin+alpha*Energie_globale_3d(gxref,gxreca,nb_param,p,param_norm,dist,regularisation,masque_param)/3.0;
  Efin=Efin+alpha*Energie_globale_3d(gyref,gyreca,nb_param,p,param_norm,dist,regularisation,masque_param)/3.0;
  Efin=Efin+alpha*Energie_globale_3d(gzref,gzreca,nb_param,p,param_norm,dist,regularisation,masque_param)/3.0;
aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 //liberation memoire des variables
 free_dvector(p,nb_param);
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);
 //free_field3d(gradIref);
 free_field3d(gradIreca);
 free_field3d(gradGxref);free_field3d(gradGyref);free_field3d(gradGzref);

 free_imatrix_3d(masque_bloc);
 free_hessien3d(hessienIref);
 free_hessien3d(hessienGxref);free_hessien3d(hessienGyref);free_hessien3d(hessienGzref);

free_grphic3d(gxref);free_grphic3d(gyref);free_grphic3d(gzref);
free_grphic3d(gxreca);free_grphic3d(gyreca);free_grphic3d(gzreca);
    
 return(Efin);
}


/*******************************************************************************
**     TOP_min_marquardt_adapt_penal_locale_3d(imref,imreca,imres,champ_ini,champ_fin,the_transf,inter  
**                   ,dist,nb_param,min_param,max_param,prec_param,param)
**                                                                       
**     Minimisation par la m�thode de Levenberg-Marquardt                      
**                                                                       
**     utilisation:                                                      
**                  imref,imreca,imres: images       
**                  champ_ini : champ initial
**                              NULL = pas de champ initial
**                  champ_fin : champ final resultat de la transfo
**                  the_transf : pointeur sur la fonction de transformation
**                           int the_transf(nb_param,param,champ)
**                  inter : pointeur sur la fonction d'interpolation
**                           int inter(imdeb,champ,imres)
**                  dist : pointeur sur la fonction de distance
**                         double dist(im1,im2)
**                  nb_param: nombre de parametre a minimiser
**                  min_param: bornes inferieures des parametres
**                  max_param: bornes superieures des parametres
**                  prec_param: precision demandee pour chaque parametre
**                  param: resultats donnant la valeur des parametres au
**                         point definissant le minimum
**      la fonction retourne un double contenant la valeur du minimum
**      Les parametres correspondant a ce minimum sont retourne dans param
*******************************************************************************/
double TOP_min_marquardt_adapt_penal_locale_3d  (grphic3d *imref, grphic3d *imreca, grphic3d *imres, field3d *champ_fin,
                                                                    transf_func_t the_transf, InterpolationFct inter, dist_func_locale_t dist,reg_func_locale_t regularisation,
                                                                                                int nb_param,  double *param, int ***masque_param,double Jmin, double Jmax,char *nomfichiertrf)
{
 int wdth,hght,dpth,i,j; 
 field3d *gradIref,*maskgradIref=NULL;
 double *grad,*p,*p_i,gradl[3],*param_norm,p_ini[3],p_tmp[3];
 double Edebut,Efin,E1,E2,E,Eprec,Emin,Etmp;
 int    nb,nbreduc,stop,nb_iter,nb_tmp;
 double prreduc,precis,precisx,precisy,precisz; 
 int topi,topj,topk,topD,resol;
 int compt,topD1;
 int ***masque_bloc;
 double ***norme_gradient;
 transf3d *transfo;
 TSlpqr Slpqr[8];
 int *x00,*x11,*y00,*y11,*z00,*z11;
 int (*gradient)();  
 reg_grad_locale_t gradient_reg;  
 int (*func_hessien)();  
 reg_hessien_locale_t hessien_reg;
 hessien3d *hessienIref, *maskhessienIref=NULL; 
 mat3d pseudo_hessien,hessien,inv_hessien;
 //double cfmax;
 double xm,xM,ym,yM,zm,zM;
 int res_inv,continu;
 double nuk,ck=10;
 int tutu/*,nbg=0*/;
#ifndef SAVE_INTERMEDIAIRE
 char nomfichier[255];
 char temp[2];
 char *ext;
#endif
double max_gradient,seuil_gradient,min_gradient;
int imax=0,jmax=0,kmax=0;


 wdth=imref->width;hght=imref->height;dpth=imref->depth;
 resol=BASE3D.resol;
 topD=(int)pow(2.0,1.0*resol)-1;
 x00=BASE3D.x0;x11=BASE3D.x1;y00=BASE3D.y0;y11=BASE3D.y1;z00=BASE3D.z0;z11=BASE3D.z1;
  
 gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) gradient=gradient_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) gradient=gradient_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) gradient=gradient_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) gradient=gradient_base_quad_sym_locale_3d;
if (dist==Energie_geman_locale_3d) gradient=gradient_base_geman_locale_3d;
 
 func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_L1L2_locale_3d) func_hessien=hessien_base_L1L2_locale_3d;
 if (dist==Energie_L1L2_sym_locale_3d) func_hessien=hessien_base_L1L2_sym_locale_3d;
 if (dist==Energie_quad_topo_locale_3d) func_hessien=hessien_base_quad_locale_3d;
 if (dist==Energie_quad_sym_topo_locale_3d) func_hessien=hessien_base_quad_sym_locale_3d;
 if (dist==Energie_geman_locale_3d) func_hessien=hessien_base_geman_locale_3d;


  gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) gradient_reg=gradient_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) gradient_reg=gradient_regularisation_energie_membrane_Lp_local;

 hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_local) hessien_reg=hessien_regularisation_energie_membrane_local;
 if (regularisation==regularisation_energie_membrane_Lp_local) hessien_reg=hessien_regularisation_energie_membrane_Lp_local;
 
    //allocation memoire des variables
    p=alloc_dvector(nb_param);
    p_i=alloc_dvector(nb_param);
    param_norm=alloc_dvector(nb_param);
 
    grad=alloc_dvector(nb_param);
  //gradIref=cr_field3d(wdth,hght,dpth);   // pour economiser de la place memoire
    gradIref=champ_fin;
    masque_bloc=alloc_imatrix_3d(topD,topD,topD);
 norme_gradient=alloc_dmatrix_3d(topD,topD,topD);
    hessienIref=cr_hessien3d(wdth,hght,dpth);


 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->nb_param=nb_param;transfo->param=CALLOC(nb_param,double);
 transfo->typetrans=BSPLINE3D;
 transfo->resol=BASE3D.resol;transfo->degre=1;   // Valeur codee en dur car ca ne marche que pour les splines de degre 1
 p = transfo->param; 


//Initialisation de masque_bloc
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        masque_bloc[topi][topj][topk]=1;

for (i=0;i<nb_param;i++) p_i[i]=p[i]=param[i];

// initialisation des parametres normalises sur [0,1] 
 for (i=0;i<nb_param/3;i++)
    {
    param_norm[3*i]  =1.0*p[3*i]  /wdth;
    param_norm[3*i+1]=1.0*p[3*i+1]/hght;
    param_norm[3*i+2]=1.0*p[3*i+2]/dpth;
    }

//Mise a jour de Jm et JM dans les Slpqr
for (i=0;i<8;i++)
    {
    Slpqr[i].Jm=Jmin;
    Slpqr[i].JM=Jmax;
    } 
        
 //calcul de l'energie de depart
    Edebut=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  E2=Edebut;

 // calcul du gradient de l'image a recaler 
 imx_gradient_3d_p(imreca,gradIref,2,4);
 
 // calucl du Hessiend e l'image 
 imx_hessien_3d_p(imreca,hessienIref,2,4);

 aff_log("Edebut: %.2f   nb_param: %d \n",Edebut,nb_param);

 nb=nbreduc=0;prreduc=0.1/100.0;
 stop=1;

 // --- On fixe la pr�cision � 10% de la taille du support de spline
 precis=MAXI(TOPO_PONDERATION_RECALAGE_SYMETRIQUE,1.0);  // pour que la pr�cision en recalage sym�trique soit �quivalent que celle obtenue en non symetrique

 precisx=CRITERE_ARRET_PARAM*precis*(x11[0]-x00[0]);
 precisy=CRITERE_ARRET_PARAM*precis*(y11[0]-y00[0]);
 precisz=CRITERE_ARRET_PARAM*precis*(z11[0]-z00[0]);
 
/* calcul des normes de gradient */
  for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        {
            topD1=TOP_conv_ind(topi,topj,topk,nb_param);
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
        
            gradient(imref,imreca,nb_param,p,param_norm,gradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            norme_gradient[topi][topj][topk]=sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);           
        }

max_gradient=0;
/* recherche du gradient max */ 
 for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
            {
            if(norme_gradient[topi][topj][topk]>max_gradient)
                {max_gradient=norme_gradient[topi][topj][topk];
                imax=topi;jmax=topj;kmax=topk;
                }
            }

min_gradient=HUGE_VAL;
/* recherche du gradient max */ 
 for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
            {
            if(norme_gradient[topi][topj][topk]<min_gradient)
                {min_gradient=norme_gradient[topi][topj][topk];
                }
            }

seuil_gradient=min_gradient+0.1*(max_gradient-min_gradient);

/* seuillage du gradient */
 for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
        {
        if(norme_gradient[topi][topj][topk]<seuil_gradient)
            norme_gradient[topi][topj][topk]=0;
        }
        
        
E=0;
 do
 {
  Eprec=E;
  E=0;
  E1=E2;
  compt=1;
  stop=0;
    nb_tmp=0;

    topi=imax;topj=jmax;topk=kmax;
            
            topD1=TOP_conv_ind(topi,topj,topk,nb_param);
        // printf("%d ieme bloc sur %d \r",topD1/3,nb_param/3);
            
    
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
        
            for (i=0;i<3;i++)
                p_ini[i]=p[topD1+i];
    
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;
                    
        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
        
            
        continu=0;
        nuk=10000;
        nb_iter=0;
        tutu=1;
        if (Etmp==0) continu=1;
        do{
        Emin=Etmp;
        if (tutu>0)
            {
            // Calcul du gradient

            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            func_hessien(imref, imreca, nb_param,p,param_norm,gradIref,maskgradIref,hessienIref,maskhessienIref,& hessien,topi,topj,topk,Slpqr,hessien_reg);
            
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    pseudo_hessien.H[i][j]=hessien.H[i][j];
            }   
        // Calcul du pseudo hessien Hk=H+nuk*I jusqu'a ce qu'il soit defini positif
        res_inv=0;
        
        while (res_inv==0)
            {
                for (i=0;i<3;i++)
                    pseudo_hessien.H[i][i]=hessien.H[i][i]+nuk; 
            
                res_inv=inversion_mat3d(&pseudo_hessien,&inv_hessien);
                
                if(res_inv==0)
                    nuk=nuk*ck;
            }
        
        
        for (i=0;i<3;i++)
            p_tmp[i]=p[topD1+i];
                
            for (i=0;i<3;i++)
                for (j=0;j<3;j++)
                    p[topD1+i]=p[topD1+i]-inv_hessien.H[i][j]*gradl[j];
        
                    param_norm[topD1]  =p[topD1]  /wdth;    
                    param_norm[topD1+1]=p[topD1+1]/hght;    
                    param_norm[topD1+2]=p[topD1+2]/dpth;    
        
        Etmp=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
      
                
        if (Etmp>Emin)
                {
                for (i=0;i<3;i++)
                    p[topD1+i]=p_tmp[i];
                    
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;
                
                nuk=nuk*ck; 
                Etmp=Emin;
                if (nuk>100000000) {continu=1;}
                tutu=0;
                nb_iter=0;  
                }
        else{
                        {
                        nuk=1.0*nuk/ck; tutu=2;
                        nb_iter++;
                            if (((1.0*(Emin-Etmp)/Emin)<0.01)||(nb_iter>10))
                                {continu=1;Emin=Etmp;}
                        }
                                
                }
            
                    
        }while (continu==0);
        
        
        
        if (Emin==HUGE_VAL)
            {for (i=0;i<3;i++)
                p[topD1+i]=p_ini[i];
            
                param_norm[topD1]  =p[topD1]  /wdth;    
                param_norm[topD1+1]=p[topD1+1]/hght;    
                param_norm[topD1+2]=p[topD1+2]/dpth;    
                Emin=dist(imref,imreca,nb_param,p,param_norm,topi,topj,topk,Slpqr,regularisation);
            }
        E=E+Emin;
        
/* mise a jour des normes de gradient */
  for (topi=MAXI(imax-1,0); topi<MINI(imax+1,topD); topi++)
    for (topj=MAXI(jmax-1,0); topj<MINI(jmax+1,topD); topj++)
      for (topk=MAXI(kmax-1,0); topk<MINI(kmax+1,topD); topk++)
        {
            topD1=TOP_conv_ind(topi,topj,topk,nb_param);
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
        
            gradient(imref,imreca,nb_param,p,param_norm,gradIref,maskgradIref,gradl,topi,topj,topk,Slpqr,gradient_reg);
            grad[topD1]=gradl[0];grad[topD1+1]=gradl[1];grad[topD1+2]=gradl[2];    
            
            norme_gradient[topi][topj][topk]=sqrt(gradl[0]*gradl[0]+gradl[1]*gradl[1]+gradl[2]*gradl[2]);   
            if (norme_gradient[topi][topj][topk]<seuil_gradient) norme_gradient[topi][topj][topk]=0;
            if ((topi==imax)&&(topj==jmax)&&(topk==kmax)&&(norme_gradient[topi][topj][topk]==max_gradient)) {norme_gradient[topi][topj][topk]=0; /*nbg++;max_gradient=0;*/}     
        }

/*if (max_gradient>0)
    nbg=0;

if (nbg<2)
{*/
max_gradient=0;
/* recherche du gradient max */ 
 for (topi=0; topi<topD; topi++)
    for (topj=0; topj<topD; topj++)
      for (topk=0; topk<topD; topk++)
            {
            if(norme_gradient[topi][topj][topk]>max_gradient)
                {max_gradient=norme_gradient[topi][topj][topk];
                imax=topi;jmax=topj;kmax=topk;
                }
            }
            
    //}
} while ((max_gradient>0)/*&&(Eprec!=E)*/);

 //calcul de l'energie finale
  Efin=Energie_globale_3d(imref,imreca,nb_param,p,param_norm,dist,regularisation,masque_param);
  aff_log("nbiter: %d    Efin: %.2f  reduction: %.2f %%  \n",nb,Efin,(Edebut-Efin)/Edebut*100.0);

  for (i=0;i<nb_param;i++) param[i]=p[i];

  #ifndef TOPO_DEBUG
    {
    #ifndef TOPO_VERBOSE
    aff_log("Verif. fin resolution..... ");
    #endif 
    TOP_verification_all(nb_param,param,masque_param,Jmin,Jmax); 
    #ifndef TOPO_VERBOSE
    aff_log("finie \n");
    #endif
    }
  #endif

    #ifndef SAVE_INTERMEDIAIRE
    if (nomfichiertrf!=NULL)
        {
         //si l'extension .trf existe on la supprime
        strcpy(nomfichier,nomfichiertrf);
        ext=strstr(nomfichier,".trf");
        if (ext!=NULL) *ext='\0';

        strcat(nomfichier,"_");
        temp[0]=(char)(resol+48);
        temp[1]=0;
        strcat(nomfichier,temp);
    
        save_transf_3d(transfo,nomfichier);
        }
    #endif


 //liberation memoire des variables
 free_dvector(p,nb_param);
 free_dvector(p_i,nb_param);
 free_dvector(param_norm,nb_param);
 free_dvector(grad,nb_param);
 //free_field3d(gradIref);
 free_imatrix_3d(masque_bloc);
 free_dmatrix_3d(norme_gradient);
 free_hessien3d(hessienIref);

 return(Efin);
}






