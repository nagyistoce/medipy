/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!**********************************************************************
*** 
*** \file       correction.c
***
*** project:    Imagix 2.01
***         
***
*** \brief description: correction d'images
*** 
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.
*** All rights are reserved.
***
***     Last user action by: H. Boisgontier (mars 2007)
***
*************************************************************************/

#include "correction.h"


/************************ corr_sagittal_3d() ********************************/ 
/*  Correction coupe par coupe ( coupes sagittales ) d'une image 3D     */
/*                                      */
/*      Normalisation des moyennes et ecarts types                  */ 
/*                                      */
/*          im_ori  :image a normaliser                         */ 
/*          im_res  :image resultat                             */
/*                                                                          */
/****************************************************************************/
void corr_sagittal_3d(void) {
  int im_ori, im_res;
  
  im_ori=GET_PLACE3D(TEXT0231);
  im_res=GET_PLACE3D(TEXT0006);
  
  imx_corr_sagittal_3d(im_ori, im_res);
  
  show_picture_3d(im_res);
}

/****************************************************************************/

void imx_corr_sagittal_3d(int im_ori, int im_res) { 
  grphic3d *imori, *imres;
  
  imori=ptr_img_3d(im_ori); 
  imres=ptr_img_3d(im_res);
  
  imx_corr_sagittal_3d_p(imori, imres); 
}

/****************************************************************************/
/* La correction se fait uniquement suivant l'axe Y  */
/****************************************************************************/

int imx_corr_sagittal_3d_p(grphic3d *imori, grphic3d *imres)
{
  int x,y,z, cpt;
  long  moyP, moyI, moyM;
  float etyP, etyI, etyM, diff;

 moyP=0;
  cpt=0;
  // pour chaque coupe pair
  for(y=0; y<imori->height; y+=2)
    for(x=0; x<imori->width; x++)
      for(z=0; z<imori->depth; z++)
        if(imori->mri[x][y][z]!=0) {
          moyP += imori->mri[x][y][z];
          cpt++;
        }
  moyP/=cpt;

  etyP=0;

  // pour chaque coupe pair
  for(y=0; y<imori->height; y+=2)
    for(x=0; x<imori->width; x++)
      for(z=0; z<imori->depth; z++)
        if(imori->mri[x][y][z]!=0) {
          diff = imori->mri[x][y][z] - moyP;
          etyP += (diff*diff);
        }
  etyP=sqrt(etyP/(cpt-1));


 moyI=0;
  cpt=0;
  // pour chaque coupe impair
  for(y=1; y<imori->height; y+=2)
    for(x=0; x<imori->width; x++)
      for(z=0; z<imori->depth; z++)
        if(imori->mri[x][y][z]!=0) {
          moyI += imori->mri[x][y][z];
          cpt++;
        }
  moyI/=cpt;

  etyI=0;
  // pour chaque coupe pair
  for(y=1; y<imori->height; y+=2)
    for(x=0; x<imori->width; x++)
      for(z=0; z<imori->depth; z++)
        if(imori->mri[x][y][z]!=0) {
          diff = imori->mri[x][y][z] - moyI;
          etyI += (diff*diff);
        }
  etyI=sqrt(etyI/(cpt-1));

  moyM = (moyI + moyP)/2;
  etyM = (etyI + etyP)/2;
  // pour chaque coupe
  for(x=0; x<imori->width; x++)
    for(y=0; y<imori->height; y++)
      for(z=0; z<imori->depth; z++)
        if(imori->mri[x][y][z]==0) {
          imres->mri[x][y][z] = 0;
        } else {
          if(y%2==0) { // si coupe pair
            imres->mri[x][y][z] = (imori->mri[x][y][z] - moyP) / etyP * etyM + moyM;
          } else { // si coupe impair
            imres->mri[x][y][z] = (imori->mri[x][y][z] - moyI) / etyI * etyM + moyM;
          }
        }
  

  // recopie des attributs de l'image d'origine dans l'image finale
  imx_copie_param_3d_p(imori,imres);
 
 return(1);
}
