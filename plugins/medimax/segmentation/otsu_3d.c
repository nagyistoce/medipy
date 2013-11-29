/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!************************************************************************
***	
***	\file:		otsu_3d.c
***
***	project:	Imagix 2.01
***			
***
***	\brief description:    Fichier source pour segmentation 3D
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	Copyright (c) 1997, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***
*************************************************************************/
#include <config.h>
#include <string.h> 	/* Needed for strcpy	*/

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "noyau/imx_lang.h"
#include "segmentation/otsu_3d.h"
#include "noyau/io/imx_export_file.h"

/*  Defini dans segmentation.c */
double calc_omega   (double *prob, int k1, int k2);
double calc_mu      (double *prob, int k1, int k2);


/********************************************************************
**   int otsu_thresholding_3d()
*/
/*! seuillage 3d
**
*****************************************************/
void otsu_thresholding_3d(void)
{
    int im_deb,im_res,thr_nb;
    int i;
    char *quest[16];

    im_deb=GET_PLACE3D("Original Image ? ");
    im_res=GET_PLACE3D("Result Image ? ");

    for(i=0;i<6;i++)
        quest[i]=CALLOC(80,char);

    strcpy(quest[1]," 1 threshold (2 regions) ");
    strcpy(quest[2]," 2 threshold (3 regions) ");
    strcpy(quest[3]," 3 threshold (4 regions) ");
    strcpy(quest[0]," Cancel                  ");
    strcpy(quest[4],"\0");

    thr_nb=GETV_QCM("How many thresholds ?",(char **)quest);

    if(thr_nb>MAX_THRESHOLD_NB || thr_nb==0)
        return;

    imx_otsu_thresholding_3d(im_deb,im_res,thr_nb);

    show_picture_3d(im_res);

}

/********************************************************************
**   int imx_otsu_thresholding_3d()
*/
/*!  Seuillage
**	\param im_deb : numero de l'image source
**	\param im_res : numero de l'image destination
**  \param thr_nb : nbre de seuils
**
*****************************************************/


void imx_otsu_thresholding_3d(int im_deb, int im_res, int thr_nb)
{
    grphic3d *imdeb,*imres;

    imdeb=ptr_img_3d(im_deb);
    imres=ptr_img_3d(im_res);

    imx_otsu_thresholding_3d_p(imdeb,imres,thr_nb);

}


/************************************************************
**  imx_otsu_thresholding_3d_p    
*/
/*!  \brief Seuillage 3d 
**	\param imdeb : image source
**	\param imres : image destination (E/S)
**  \param thr_nb ; nbre de seuils
**
*****************************************************************/

void imx_otsu_thresholding_3d_p(grphic3d *imdeb, grphic3d *imres, int thr_nb)
{
    int threshold[MAX_THRESHOLD_NB];
    int i,j,k,m;
    int width,height,depth;
    int flag=0;

    width=imdeb->width;
    height=imdeb->height;
    depth=imdeb->depth;

    imx_copie_param_3d_p(imdeb,imres);

    find_threshold_3d(imdeb,threshold,thr_nb,TRUE);



    for(i=0;i<width;i++)
        for(j=0;j<height;j++)
            for(k=0;k<depth;k++)
            {
                m=0;
                flag=0;

                while(m<thr_nb && flag==0)
                {
                    if(imdeb->mri[i][j][k]<=threshold[m])
                    {
                        imres->mri[i][j][k]=m;
                        flag=1;
                    }


                    else if(m==thr_nb-1)
                        imres->mri[i][j][k]=m+1;

                    m++;
                }
            }


    /**** Pour l'affichage ***/

    for(i=0;i<thr_nb;i++)
    {
        threshold[i]=(int)(0.5+(double)threshold[i]*(double)MAX_GREY_LEVEL/(double)imdeb->max_pixel);
        PUTI("Threshold  =",threshold[i]," ");
    }


    imx_inimaxpixel_3d_p(imres);
    imres->icomp=0;
    imres->rcoeff=1;

}

/********************************************************************/

void otsu_brain_3d(void)
{
    int im_deb,im_res;

    im_deb=GET_PLACE3D("Original Image ? ");
    im_res=GET_PLACE3D("Result Image ? ");

    imx_otsu_brain_3d(im_deb,im_res);

    imx_copie_param_3d(im_deb,im_res);
    imx_iniparaimg_3d(im_res);
    show_picture_3d(im_res);


}

/********************************************************************/

void imx_otsu_brain_3d(int im_deb, int im_res)
{
    grphic3d *imdeb,*imres;

    imdeb=ptr_img_3d(im_deb);
    imres=ptr_img_3d(im_res);

    imx_otsu_brain_3d_p(imdeb,imres);


}

/********************************************************************/

void imx_otsu_brain_3d_p(grphic3d *imdeb, grphic3d *imres)
{
    int threshold[MAX_THRESHOLD_NB];
    int i,j,k;
    int width,height,depth;
    int thr_nb;

    width=imdeb->width;
    height=imdeb->height;
    depth=imdeb->depth;

    thr_nb=2;


    find_threshold_3d(imdeb,threshold,thr_nb,TRUE);


    for(i=0;i<width;i++)
        for(j=0;j<height;j++)
            for(k=0;k<depth;k++)
            {
                if(imdeb->mri[i][j][k]>threshold[1])
                    imres->mri[i][j][k]=1;

                else 
                    imres->mri[i][j][k]=0;

            }

}

/**********************************************************************************/

void otsu_ventricules_3d(void)
{
    int im_deb,im_res;

    im_deb=GET_PLACE3D("Original Image ? ");
    im_res=GET_PLACE3D("Result Image ? ");

    imx_otsu_ventricules_3d(im_deb,im_res);

    imx_copie_param_3d(im_deb,im_res);
    imx_iniparaimg_3d(im_res);
    show_picture_3d(im_res);

}

/********************************************************************/

void imx_otsu_ventricules_3d(int im_deb, int im_res)
{
    grphic3d *imdeb,*imres;

    imdeb=ptr_img_3d(im_deb);
    imres=ptr_img_3d(im_res);

    imx_otsu_ventricules_3d_p(imdeb,imres);


}

/********************************************************************/

void imx_otsu_ventricules_3d_p(grphic3d *imdeb, grphic3d *imres)
{
    int threshold[MAX_THRESHOLD_NB];
    int i,j,k;
    int width,height,depth;
    int thr_nb;

    width=imdeb->width;
    height=imdeb->height;
    depth=imdeb->depth;

    thr_nb=2;


    find_threshold_3d(imdeb,threshold,thr_nb,TRUE);


    for(i=0;i<width;i++)
        for(j=0;j<height;j++)
            for(k=0;k<depth;k++)
            {
                if(imdeb->mri[i][j][k]>threshold[0] && imdeb->mri[i][j][k]<threshold[1])
                    imres->mri[i][j][k]=1;

                else 
                    imres->mri[i][j][k]=0;

            }


}

/**********************************************************************/
/*!
**  calcul d'un distribution de proba. de imdeb 
**
**  \param imdeb : image source
**	\param prob : tableau  
**	\pre prob doit etre alloue avant l'appel de la fonction
**
********************************************************************/
void calc_prob_table_3d(grphic3d *imdeb, double *prob, BOOL bWithZeroValue)
{
	int i,j,k,val = 0;
	int width,height,depth;
	double sum=0.0,max,min;
	
	width=imdeb->width;
	height=imdeb->height;
	depth=imdeb->depth;
	max=(double)imdeb->max_pixel;
	min=(double)imdeb->min_pixel;
	
	for(i=1;i<GREY_LEVEL_NB+1;i++)
		prob[i]=0.0;
	
	
	/* debut modif emily */
	if (max == min) {
		PUT_ERROR("ERREUR DANS LA FONCTION calc_prob_table_3d DE imx_3d.c, max=min");
		return;
	} else {
		/* fon modif emily */
		for(i=0;i<width;i++)
			for(j=0;j<height;j++)
				for(k=0;k<depth;k++)
				{
					if (bWithZeroValue || imdeb->mri[i][j][k])
					{
						val=(int)(0.5+((double)imdeb->mri[i][j][k]-min)*(double)MAX_GREY_LEVEL/(max-min)); 
					
						if(val>GREY_LEVEL_NB-1 || val<0) {
							printf("ERREUR calc_prob_table_3d =%d\n",val);
						}
						else  
							prob[val+1]+=1.0;
					}
				}
	}
	for(i=1;i<GREY_LEVEL_NB+1;i++)
		sum+=prob[i];
	
	for(i=1;i<GREY_LEVEL_NB+1;i++)
		prob[i]/=sum;
}

/*********************************
prob[val+1]+=1.0;
}

  for(i=1;i<GREY_LEVEL_NB+1;i++)
  sum+=prob[i];
  
	for(i=1;i<GREY_LEVEL_NB+1;i++)
    prob[i]/=sum;
	}

**********************************/

/**********************************************/
/*!
**	trouve les niveau de seuillages
**
**  \param imdeb : image source
**	\param threshold : le tableau qui contiendra les seuils a la fin de la fct
**	\param thr_nb : nbre de seuil
**  \pre : threshold doit etre alloue avant l'appel
** 
**********************************************/
void find_threshold_3d(grphic3d *imdeb, int *threshold, int thr_nb, BOOL bWithZeroValue)
{

    int i,k1,k2,k3,k4,k5;
//	int k6,j;				unused variable
    int width,height,depth;
//    double sum=0.0;
    double prob[GREY_LEVEL_NB+1];
    double omega1,omega2,omega3,omega4,omega5,omega6;
    double mu1,mu2,mu3,mu4,mu5,mu6;
    double dmu1,dmu2,dmu3,dmu4,dmu5,dmu6;
    double sigma,sigma_max;
    double mt;
//    int best;				unused variable
//   double max;			unused variable
//    int val;				unused variable

    width=imdeb->width;
    height=imdeb->height;
    depth=imdeb->depth;

    sigma_max=-1.0;

    calc_prob_table_3d(imdeb,prob,bWithZeroValue);

    mt=calc_mu(prob,1,GREY_LEVEL_NB);

    if(thr_nb==1)
    {
        threshold[0]=0;

        for(k1=1;k1<GREY_LEVEL_NB+1;k1++)
        {
            omega1=calc_omega(prob,1,k1);
            omega2=calc_omega(prob,k1+1,GREY_LEVEL_NB);
            mu1=calc_mu(prob,1,k1);
            mu2=calc_mu(prob,k1+1,GREY_LEVEL_NB);
            dmu1=(mu1-mt)*(mu1-mt);
            dmu2=(mu2-mt)*(mu2-mt);
            sigma=omega1*dmu1+omega2*dmu2;

            if(sigma>sigma_max)
            {
                sigma_max=sigma;
                threshold[0]=k1-1;
            }
        }
    }
    else if(thr_nb==2)
    {
        threshold[0]=0;
        threshold[1]=0;

        for(k2=1;k2<GREY_LEVEL_NB+1;k2++)
            for(k1=1;k1<k2;k1++)
            {
                omega1=calc_omega(prob,1,k1);
                omega2=calc_omega(prob,k1+1,k2);
                omega3=calc_omega(prob,k2+1,GREY_LEVEL_NB);
                mu1=calc_mu(prob,1,k1);
                mu2=calc_mu(prob,k1+1,k2);
                mu3=calc_mu(prob,k2+1,GREY_LEVEL_NB);
                dmu1=(mu1-mt)*(mu1-mt);
                dmu2=(mu2-mt)*(mu2-mt);
                dmu3=(mu3-mt)*(mu3-mt);
                sigma=omega1*dmu1+omega2*dmu2+omega3*dmu3;

                if(sigma>sigma_max)
                {
                    sigma_max=sigma;
                    threshold[0]=k1-1;
                    threshold[1]=k2-1;
                }
            }

    }
    else if(thr_nb==3)
    {
        threshold[0]=0;
        threshold[1]=0;
        threshold[2]=0;


        for(k3=1;k3<GREY_LEVEL_NB+1;k3++)
            for(k2=1;k2<k3;k2++)
                for(k1=1;k1<k2;k1++)
                {
                    omega1=calc_omega(prob,1,k1);
                    omega2=calc_omega(prob,k1+1,k2);
                    omega3=calc_omega(prob,k2+1,k3);
                    omega4=calc_omega(prob,k3+1,GREY_LEVEL_NB);
                    mu1=calc_mu(prob,1,k1);
                    mu2=calc_mu(prob,k1+1,k2);
                    mu3=calc_mu(prob,k2+1,k3);
                    mu4=calc_mu(prob,k3+1,GREY_LEVEL_NB);
                    dmu1=(mu1-mt)*(mu1-mt);
                    dmu2=(mu2-mt)*(mu2-mt);
                    dmu3=(mu3-mt)*(mu3-mt);
                    dmu4=(mu4-mt)*(mu4-mt);
                    sigma=omega1*dmu1+omega2*dmu2+omega3*dmu3+omega4*dmu4;

                    if(sigma>sigma_max)
                    {
                        sigma_max=sigma;
                        threshold[0]=k1-1;
                        threshold[1]=k2-1;
                        threshold[2]=k3-1;
                    }
                }

    }
    else if(thr_nb==4)
    {
        threshold[0]=0;
        threshold[1]=0;
        threshold[2]=0;
        threshold[3]=0;


        for(k4=1;k4<GREY_LEVEL_NB+1;k4++)
            for(k3=1;k3<k4;k3++)
                for(k2=1;k2<k3;k2++)
                    for(k1=1;k1<k2;k1++)
                    {
                        omega1=calc_omega(prob,1,k1);
                        omega2=calc_omega(prob,k1+1,k2);
                        omega3=calc_omega(prob,k2+1,k3);
                        omega4=calc_omega(prob,k3+1,k4);
                        omega5=calc_omega(prob,k4+1,GREY_LEVEL_NB);
                        mu1=calc_mu(prob,1,k1);
                        mu2=calc_mu(prob,k1+1,k2);
                        mu3=calc_mu(prob,k2+1,k3);
                        mu4=calc_mu(prob,k3+1,k4);
                        mu5=calc_mu(prob,k4+1,GREY_LEVEL_NB);
                        dmu1=(mu1-mt)*(mu1-mt);
                        dmu2=(mu2-mt)*(mu2-mt);
                        dmu3=(mu3-mt)*(mu3-mt);
                        dmu4=(mu4-mt)*(mu4-mt);
                        dmu5=(mu5-mt)*(mu5-mt);
                        sigma=omega1*dmu1+omega2*dmu2+omega3*dmu3+omega4*dmu4+omega5*dmu5;

                        if(sigma>sigma_max)
                        {
                            sigma_max=sigma;
                            threshold[0]=k1-1;
                            threshold[1]=k2-1;
                            threshold[2]=k3-1;
                            threshold[3]=k4-1;
                        }
                    }

    }
    else if(thr_nb==5)
    {
        threshold[0]=0;
        threshold[1]=0;
        threshold[2]=0;
        threshold[3]=0;
        threshold[4]=0;


        for(k5=1;k5<GREY_LEVEL_NB+1;k5++)
            for(k4=1;k4<k5;k4++)
                for(k3=1;k3<k4;k3++)
                    for(k2=1;k2<k3;k2++)
                        for(k1=1;k1<k2;k1++)
                        {
                            omega1=calc_omega(prob,1,k1);
                            omega2=calc_omega(prob,k1+1,k2);
                            omega3=calc_omega(prob,k2+1,k3);
                            omega4=calc_omega(prob,k3+1,k4);
                            omega5=calc_omega(prob,k4+1,k5);
                            omega6=calc_omega(prob,k5+1,GREY_LEVEL_NB);
                            mu1=calc_mu(prob,1,k1);
                            mu2=calc_mu(prob,k1+1,k2);
                            mu3=calc_mu(prob,k2+1,k3);
                            mu4=calc_mu(prob,k3+1,k4);
                            mu5=calc_mu(prob,k4+1,k5);
                            mu6=calc_mu(prob,k5+1,GREY_LEVEL_NB);
                            dmu1=(mu1-mt)*(mu1-mt);
                            dmu2=(mu2-mt)*(mu2-mt);
                            dmu3=(mu3-mt)*(mu3-mt);
                            dmu4=(mu4-mt)*(mu4-mt);
                            dmu5=(mu5-mt)*(mu5-mt);
                            dmu6=(mu6-mt)*(mu6-mt);
                            sigma=omega1*dmu1+omega2*dmu2+omega3*dmu3+omega4*dmu4+omega5*dmu5+omega6*dmu6;

                            if(sigma>sigma_max)
                            {
                                sigma_max=sigma;
                                threshold[0]=k1-1;
                                threshold[1]=k2-1;
                                threshold[2]=k3-1;
                                threshold[3]=k4-1;
                                threshold[4]=k5-1;
                            }
                        }

    }


    for(i=0;i<thr_nb;i++)
		threshold[i]=(int)floor((double)(0.5+threshold[i]*(double)imdeb->max_pixel/(double)MAX_GREY_LEVEL));



}

