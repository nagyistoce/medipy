/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <sys/stat.h>   
#include <string.h>
#include "detection_fonctions_wrappees.h"

#include "noyau/imx_3d.h"

/*******************************************************************************
**	ComputeRocCurve
*******************************************************************************/
int ComputeRocCurve(grphic3d *imDetection,grphic3d *maskimDetection,grphic3d *maskBonnesDetec, int nbPointsH,char *fichier)
{

double s;
int p, fp, vp, pos, neg,  exists;
FILE *fd;
char  f_option[4];
struct  stat	buf;
int x,y,z;
float smin,smax;

/*//-// 07/07/2010 : ajout de #include "noyau/imx_3d.h" et  i
 * mx_inimaxminpixel_3d_p(imDetection); car problèmes de mise à jour min et
 * max des pixels (FELDER Christophe)*/
imx_inimaxminpixel_3d_p(imDetection);
/*//-//*/

if((fichier == NULL) || strlen(fichier) <=0) {
	PUT_WARN("Invalid Filename !");
    return(0);
	}
  
exists = (stat(fichier, &buf) == 0);
sprintf(f_option,"w");
/*if(exists) { 
    sprintf(f, "Erase old file ?");
    
    if(GET_EXIST("File exist, append ?", &err) != 1) {//le fichier existe, on ne veut pas ajouter
      if(XGET_EXIST(f, &err) == 0) {//on ne veut pas ecraser
        return(0);
      }
    }
    else sprintf(f_option,"a");
  }
 */ 
fd = fopen(fichier, f_option);
  
/* Calcul de la courbe ROC */
smin = imDetection->min_pixel * imDetection->rcoeff;
smax = imDetection->max_pixel * imDetection->rcoeff;

for(p=0; p<nbPointsH; p++) {
   	fp=0;vp=0;neg=0;pos=0;
    s = smin + (float)p/(float)nbPointsH * (smax - smin);
    for(x=0; x<imDetection->width; x++)
      for(y=0; y<imDetection->height; y++)
        for(z=0; z<imDetection->depth; z++)
          if(!imDetection->mask || imDetection->mask->mri[x][y][z]!=0) {
            if(maskBonnesDetec->mri[x][y][z]!=0) {// "positif"
              pos++;
              if(imDetection->mri[x][y][z] * imDetection->rcoeff >= s)
                vp++;
            } else { // "negatif"
              neg++;
              if(imDetection->mri[x][y][z] * imDetection->rcoeff >= s)
                fp++;
            }
          }
    fprintf(fd, "%f %f\n", (float)fp/(float)neg, (float)vp/(float)pos);
  }
fclose(fd);
return(1);
}
