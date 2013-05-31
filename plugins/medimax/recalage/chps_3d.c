/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/******************************************************************************
**/ 
/*! \file:      chps_3d.c
***
*** project:    Imagix 1.01
***         
***
*** \brief description:    Fichier de gestion et de traitement des champs 3D
*** 
*** 
*** Copyright (c) 1993, ULP-IPB Strasbourg.                                            
*** All rights are reserved.
***
***
********************************************************************************/
#include <config.h>
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>     
#include <time.h>


#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_lang.h"
#include "noyau/imx_3d.h"
#include "noyau/io/imx_head.h"

#include "noyau/imx_2d.h"
#include "noyau/extern.h"
#include "noyau/io/imx_log.h"

#include "recalage/chps_3d.h"
#include "recalage/mtch_3d.h"
#include "math/oper_3d.h"
#include "math/imx_matrix.h"
#include "math/imx_bspline.h"
#include "outils/imx_misc.h"
#include "noyau/mani_3d.h"
#include "noyau/gui/imx_picture3d.h"
#include "recalage/transformations_3d.h"

#define TOP_D(nb_param)  ( floor (pow(1.0*(nb_param)/3.0,1.0/3.0) + 0.1 ) )
#define TOP_conv_ind(i,j,k,nb_param) ( 3*( (i)*(TOP_D(nb_param))*(TOP_D(nb_param)) + (j)*(TOP_D(nb_param)) + (k)) )



/*******************************************************************************
********************************************************************************
****************** ALLOCATION ET LIBERATION MEMOIRE ****************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**     cr_field3d(wdth,hght,dpth)                                   
*/
/*!                                                                   
**     Creation et allocation memoire d'une structure champ pour les 
**  tailles specifiees en parametres  
**  \param  wdth, hght, dpth : taille du champ
**  \retval field3d, le champ alloue, exit(1) si echec
*******************************************************************************/
field3d *cr_field3d(int wdth, int hght, int dpth)
{ 
  vector3d ***data; 
  vector3d *mt;
  field3d  *ch;
  int i,j;
  
  /*allocation memoire de la structure champ*/
  ch=(field3d*)malloc(sizeof(field3d));
  if (ch==NULL) {printf("allocation failed in step 1 of cr_field3d \n");exit(1);}

  /*allocation memoire du tableau de donnee data*/
  data=(vector3d ***) malloc((size_t)((wdth+1)*sizeof(vector3d **)));
  if (data==NULL) {printf("allocation failed in step 2 of cr_field3d \n");exit(1);}

  data[0]=(vector3d **)malloc((size_t)((wdth*hght)*sizeof(vector3d *)));
  if (data[0]==NULL) {printf("allocation failed in step 3 of cr_field3d \n");exit(1);}
  
  for(i=1;i<wdth;i++)
     data[i]=data[i-1]+hght; 
 
  data[0][0]=(vector3d *)malloc((size_t)((wdth*hght*dpth)*sizeof(vector3d)));
  if (data[0][0]==NULL) {printf("\n allocation failed in step4  of cr_field3d \n"); exit(1);}
  mt=data[0][0];
  
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
   {data[i][j]=mt;mt=mt+dpth;}
  
  ch->raw=data;
  ch->width=wdth;ch->height=hght;ch->depth=dpth;
  ch->dx=ch->dy=ch->dz=0;
  
  return (ch);
 }


/*******************************************************************************
**    free_field3d(ch)                                   
*/
/*!                                                                   


**     liberation memoire d'une structure champ:
**   \param ch: champ a liberer
**   \retval 1
*******************************************************************************/
 int free_field3d(field3d *ch)
{
  vector3d ***data;
  
  data=ch->raw;

  /*liberation memoire des donnees data*/
  free((data[0][0]));
  free((data[0]));
  free(data);

  /*liberation memoire de la structure*/
  free(ch);
  
  return (1);
 }


/*******************************************************************************
**     cr_transf3d(wdth,hght,dpth,t_ini)                                   
*/                                                                   
/*!     Creation et allocation memoire d'une structure transf3d
**     si t_ini!=NULL allocation et remplissage de la chaine de caractere trans_ini
**
**  \param wdth, hght,dpth : taille
**  \param t_ini : transformation initial ?
**  \retval transf3d : la structure transf3d CALLOC'e, liberable avec free_transf3d()
*******************************************************************************/
transf3d *cr_transf3d(int wdth, int hght, int dpth, char *t_ini)
{ 
  transf3d *tr;
  
  
  tr=CALLOC(1,transf3d);
  
  tr->width=wdth;tr->height=hght;tr->depth=dpth;
  tr->dx=0.0;tr->dy=0.0;tr->dz=0.0;
  tr->param=NULL;
  
  if (t_ini!=NULL) 
    tr->trans_ini=strdup(t_ini);
  else
    tr->trans_ini=NULL;
  
  return (tr);
}

/*******************************************************************************
**     cr_transf3d(imref,t_ini)
*/
/*!     Creation et allocation memoire d'une structure transf3d
**     si t_ini!=NULL allocation et remplissage de la chaine de caractere trans_ini
**
**  \param imref : taille et dimension
**  \param t_ini : transformation initial ?
**  \retval transf3d : la structure transf3d CALLOC'e, liberable avec free_transf3d()
*******************************************************************************/
transf3d *cr_transf3d_p(grphic3d *imref, char *t_ini)
{
  transf3d *tr;


  tr=CALLOC(1,transf3d);

  tr->width=imref->width;tr->height=imref->height;tr->depth=imref->depth;
  tr->dx=imref->dx;tr->dy=imref->dy;tr->dz=imref->dz;
  tr->param=NULL;

  if (t_ini!=NULL)
    tr->trans_ini=strdup(t_ini);
  else
    tr->trans_ini=NULL;

  return (tr);
}

/*******************************************************************************
**     cr_copy_transf3d(transfo_src)
*/
/*!     creation d'une transformation par copie
**
**  \param transfo_src : transformation source
**  \retval transf3d : la transformation copie de transfo_src
*******************************************************************************/
transf3d *cr_copy_transf3d(transf3d *transfo_src)
{
  transf3d *transfo_dest;
  int res=0;

  //allocation de la transformation
  transfo_dest=CALLOC(1,transf3d);
  if (!transfo_dest)
  { fprintf(stderr, "l'allocation memoire a echoue dans cr_copy_transf3d\n"); return NULL; }
  //copie de la transformation source dans celle destination
  res=copie_transf_3d(transfo_src, transfo_dest);
  if (res)
  {
   fprintf (stderr, "la copie de la transformation a echoue dans cr_copy_transf3d\n");
   free_transf3d(transfo_dest);
   return NULL;
  }

  return transfo_dest;  
}


/*******************************************************************************
**     free_transf3d(tr)                                   
*/                                                                   
/*!
**     liberation memoire d'une structure transf3d
**  \param tr : la strucuture a liberer
**  \retval 1 si succes, 0 si echec.
*******************************************************************************/
int free_transf3d(transf3d *tr)
{
  if (!tr)
    return 0;
  switch (tr->typetrans)
  {
   case RIGID3D: free(tr->param);
   break;
   case RIGIDZOOM3D: free(tr->param);
   break;
   case AFFINE3D: free(tr->param);
   break;
   case RIGIDGLOBALZOOM3D: free(tr->param);
   break;
   case AFFINEDECOUPLE3D: free(tr->param);
   break;
   case BSPLINE3D: free(tr->param);
   break;
   case CHAMP3D: free(tr->ch);
   break;
   case CHAMPVOXEL3D: free(tr->ch);
   break;
  }
  if (tr->trans_ini!=NULL) free(tr->trans_ini);
  free(tr);
  return (1);
 }


/*******************************************************************************
********************************************************************************
*************** OPERATIONS MATHEMATIQUES SUR LES CHAMPS ************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**  oper_field_3d(champ1,champ2,champres)                                   
**                                                                   
**  Cette fonction permet de calculer le champ resultat d'un operation entre
**  deux champs contenus dans des fichiers (trf ou chp)
**  Les resultat est mis dans un fichier chp       
*******************************************************************************/
void    oper_field_3d(void)
{
char  nomfich1[500],nomfich2[500],nomfichres[500],str[500],*quest[10];
field3d *ch1=NULL,*ch2=NULL,*chres;
transf3d *transfo1,*transfo2,*transfores;
int e=0,oper_type,i;
TDimension wdth,hght,dpth;
float coeff=0;
double dx,dy,dz;
int res=0;

 /*choix de l'operation*/
  for(i=0;i<10;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],TEXT0075);
  strcpy(quest[1],TEXT0076);
  strcpy(quest[2],"Combinaison");
    strcpy(quest[3],"Division par un scalaire");
    strcpy(quest[4],"Carre");
    strcpy(quest[5],"Racine Carre");
    strcpy(quest[6],"\0");
  
    oper_type=GETV_QCM("Format d'enregistrement",(char **)quest);
  for(i=0;i<10;i++)
    free(quest[i]);

 /*lecture du premier champ*/
 strcpy(nomfich1,GET_FILE("Premier fichier trf",&e));
 if(e != 0)
     return ;
 put_file_extension(nomfich1,".trf",nomfich1);

  if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}

if (oper_type<3)
    {
 /*lecture du deuxieme champ*/
 strcpy(nomfich2,GET_FILE("Deuxieme fichier trf",&e));
 if(e != 0)
     return ;
 put_file_extension(nomfich2,".trf",nomfich2);

 if ((test_fichier_transf_3d(nomfich2))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}
    }


if (oper_type==3)
    {
    while (coeff==0)
        coeff=GET_FLOAT("coefficient", 1, &e);
    }
    
    
 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&e));
 put_file_extension(nomfichres,".trf",nomfichres);


  switch (oper_type)
  {
   case 0: 
   case 1:  /*  Adddition et soutraction de champs */
     /*lecture des champ*/
     transfo1=load_transf_3d(nomfich1);
     transfo2=load_transf_3d(nomfich2);
     ch1=transf_to_field_3d(transfo1,NULL,NULL);
     ch2=transf_to_field_3d(transfo2,NULL,NULL);

     dx=transfo1->dx;
         dy=transfo1->dy;
         dz=transfo1->dz;
    
         /*verification de la compatibilite de taille des champs*/
     wdth=ch1->width;hght=ch1->height;dpth=ch1->depth;
        if (transfo2->dx!=transfo1->dx || transfo2->dy!=transfo1->dy || transfo2->dz!=transfo1->dz)
     {
      PUT_ERROR("The two fields have not the same dx,dy,dz");
      free_field3d(ch1);free_field3d(ch2);
      return;
     }
  
         
         free_transf3d(transfo1);
     free_transf3d(transfo2);


         /*verification de la compatibilite des dx dy dz*/
    

     
        
        
     if (ch2->width!=(int)wdth || ch2->height!=(int)hght || ch2->depth!=(int)dpth)
     {
      PUT_ERROR("The two fields have not the same size");
      free_field3d(ch1);free_field3d(ch2);
      return;
     }
     else
     {
     /*allocation du champ resultat*/
      chres=cr_field3d(wdth,hght,dpth);
        
      switch (oper_type)
      {
       case 0: add_field_3d(ch1,ch2,chres);break;
       case 1: sub_field_3d(ch1,ch2,chres);break;
       default: add_field_3d(ch1,ch2,chres);break;
      } 
     }
     
     /*enregistrement du champ resultat*/
     transfo1=field_to_transf_3d(chres,NULL,NULL);
     transfo1->dx=dx;
         transfo1->dy=dy;
         transfo1->dz=dz;
         
         save_transf_3d(transfo1,nomfichres);
     free_transf3d(transfo1);
     
     /* Libere memoire allouee */
     free_field3d(ch1);free_field3d(ch2);free_field3d(chres);
     break;
   case 2 :   /*   Combinaison de champs */
   default:
     /*lecture des champ*/
     transfo1=load_transf_3d(nomfich1);
     transfo2=load_transf_3d(nomfich2);

     wdth=transfo1->width;hght=transfo1->height;dpth=transfo1->depth;
     // si transfo pas bspline, alors la taille de transfo1 n'a pas d'importance
     if (transfo1->typetrans==BSPLINE3D)
     {
       if (transfo2->width!=wdth || transfo2->height!=hght || transfo2->depth!=dpth)
       {
        PUT_ERROR("The two fields have not the same size");
        return;
       }
            
     }
     transfores=cr_transf3d(transfo2->width,transfo2->height,transfo2->depth,NULL);

     /*verification de la compatibilite de taille des champs*/
//     if (transfo2->width!=wdth || transfo2->height!=hght || transfo2->depth!=dpth)
//     {
//      PUT_ERROR("The two fields have not the same size");
//      return;
//     }
//     else
//     {
     /*allocation du champ resultat*/
//     }
       
                
                
            if((transfo1->typetrans == CHAMP3D)||(transfo1->typetrans == CHAMPVOXEL3D))
                {char *query[100];  
                
                query[0]="nearest";
                query[1]="linear";
                query[2]="bspline 2";
                query[3]="bspline 3";
                query[4]="bspline 4";
                query[5]="bspline 5";
                query[6]="\0";
                query[7]=NULL;
                res=GETV_QCM("Interpolation",(char **)query);
                }
                 
     comb_transf_3d(transfo1,transfo2,transfores,res);
     
     /*enregistrement du champ resultat*/
     save_transf_3d(transfores,nomfichres);
     
     /* Libere memoire allouee */
     free_transf3d(transfo1); free_transf3d(transfo2);free_transf3d(transfores);
     break;
    case 3 :   /*   division par un coeff */
     /*lecture des champ*/
     transfo1=load_transf_3d(nomfich1);
     ch1=transf_to_field_3d(transfo1,NULL,NULL);
            
          dx=transfo1->dx;
         dy=transfo1->dy;
         dz=transfo1->dz;  
     wdth=transfo1->width;hght=transfo1->height;dpth=transfo1->depth;
     chres=cr_field3d(wdth,hght,dpth);
    
         div_field_3d(ch1,chres,coeff);
    
         transfo1=field_to_transf_3d(chres,NULL,NULL);
     transfo1->dx=dx;
         transfo1->dy=dy;
         transfo1->dz=dz;
         
     /*enregistrement du champ resultat*/
     save_transf_3d(transfo1,nomfichres);
     
     /* Libere memoire allouee */
     free_transf3d(transfo1); free_field3d(ch1);free_field3d(chres);
     break;

        case 4 :   /*   mise au carr� */
     
         /*lecture des champ*/
     transfo1=load_transf_3d(nomfich1);
     ch1=transf_to_field_3d(transfo1,NULL,NULL);
        dx=transfo1->dx;
         dy=transfo1->dy;
         dz=transfo1->dz;
     wdth=transfo1->width;hght=transfo1->height;dpth=transfo1->depth;
     chres=cr_field3d(wdth,hght,dpth);
    
         square_field_3d(ch1,chres);
    
         transfo1=field_to_transf_3d(chres,NULL,NULL);
     
         transfo1->dx=dx;
         transfo1->dy=dy;
         transfo1->dz=dz;
         
     /*enregistrement du champ resultat*/
     save_transf_3d(transfo1,nomfichres);
     
     /* Libere memoire allouee */
     free_transf3d(transfo1); free_field3d(ch1);free_field3d(chres);
     break;
         
            case 5 :   /*   racine carr� */
     
         /*lecture des champ*/
     transfo1=load_transf_3d(nomfich1);
     ch1=transf_to_field_3d(transfo1,NULL,NULL);
        dx=transfo1->dx;
         dy=transfo1->dy;
         dz=transfo1->dz;
     wdth=transfo1->width;hght=transfo1->height;dpth=transfo1->depth;
     chres=cr_field3d(wdth,hght,dpth);
    
         rootsquare_field_3d(ch1,chres);
    
         transfo1=field_to_transf_3d(chres,NULL,NULL);
     
         transfo1->dx=dx;
         transfo1->dy=dy;
         transfo1->dz=dz;
         
     /*enregistrement du champ resultat*/
     save_transf_3d(transfo1,nomfichres);
     
     /* Libere memoire allouee */
     free_transf3d(transfo1); free_field3d(ch1);free_field3d(chres);
     break;
    
  } 
}   



/*******************************************************************************
**  add_field_3d(champ1,champ2,champres)                                   
*/                                                                   
/*! Addition de deux champs de deformations
**  Les champs sont supposes alloues avant l'appel de la fonction          
**  \param champ1 : champres = champs1+champ2
**  \param champ2 : champres = champs1+champ2
**  \param champres : champres = champs1+champ2
**  \retval 1 si reussite, 0 si echec
*******************************************************************************/
int add_field_3d(field3d *champ1, field3d *champ2, field3d *champres)
{
 int    wdth,hght,dpth,i,j,k;
 vector3d ***data1,***data2,***datares;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 
 if (champ2->width!=wdth || champ2->height!=hght || champ2->depth!=dpth || 
 champres->width!=wdth || champres->height!=hght || champres->depth!=dpth)
 {printf("les champs ne sont pas de la bonne taille dans add_field_3d\n");return(0);}
 
 data1=champ1->raw;data2=champ2->raw;datares=champres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    (datares[i][j][k]).x=(data1[i][j][k]).x+(data2[i][j][k]).x;
    (datares[i][j][k]).y=(data1[i][j][k]).y+(data2[i][j][k]).y;
    (datares[i][j][k]).z=(data1[i][j][k]).z+(data2[i][j][k]).z;
   }

 return(1);
}

/*******************************************************************************
**  sub_field_3d(champ1,champ2,champres)                                   
*/
/*!                                                                   
**  Soustraction de deux champs de deformations champres=champ1-champ2
**  Les champs sont supposes alloues avant l'appel de la fonction 
**  \param champ1, champ2 : champ1 - champ2
**  \param champres : le champ resultant (E/S)
**  \retval 1 si reussite , 0 si echec         
*******************************************************************************/
int sub_field_3d(field3d *champ1, field3d *champ2, field3d *champres)
{
 int    wdth,hght,dpth,i,j,k;
 vector3d ***data1,***data2,***datares;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 
 if (champ2->width!=wdth || champ2->height!=hght || champ2->depth!=dpth || 
 champres->width!=wdth || champres->height!=hght || champres->depth!=dpth)
 {printf("les champs ne sont pas de la bonne taille dans add_field_3d\n");return(0);}
 
 data1=champ1->raw;data2=champ2->raw;datares=champres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    (datares[i][j][k]).x=(data1[i][j][k]).x-(data2[i][j][k]).x;
    (datares[i][j][k]).y=(data1[i][j][k]).y-(data2[i][j][k]).y;
    (datares[i][j][k]).z=(data1[i][j][k]).z-(data2[i][j][k]).z;
   }


 return(1);
}

/*******************************************************************************
**  div_field_3d(champ1,champres,coeff)                                   
*/                                                                   
/*! division d'un champ de transfo par u coeff
**  Les champs sont supposes alloues avant l'appel de la fonction          
**  \retval 1 si reussite, 0 si echec
*******************************************************************************/
int div_field_3d(field3d *champ1, field3d *champres, float coeff)
{
 int    wdth,hght,dpth,i,j,k;
 vector3d ***data1,***datares;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 
 
 data1=champ1->raw;datares=champres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    (datares[i][j][k]).x=(data1[i][j][k]).x/coeff;
    (datares[i][j][k]).y=(data1[i][j][k]).y/coeff;
    (datares[i][j][k]).z=(data1[i][j][k]).z/coeff;
   }
  
 return(1);
}

/*******************************************************************************
**  mul_field_3d(champ1,champres,coeff)                                   
*/                                                                   
/*! Multiplication d'un champ de transfo par u coeff
**  Les champs sont supposes alloues avant l'appel de la fonction          
**  \retval 1 si reussite, 0 si echec
*******************************************************************************/
int mul_field_3d(field3d *champ1, field3d *champres, float coeff)
{
 int    wdth,hght,dpth,i,j,k;
 vector3d ***data1,***datares;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 
 
 data1=champ1->raw;datares=champres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    (datares[i][j][k]).x=(data1[i][j][k]).x*coeff;
    (datares[i][j][k]).y=(data1[i][j][k]).y*coeff;
    (datares[i][j][k]).z=(data1[i][j][k]).z*coeff;
   }
  
 return(1);
}


/*******************************************************************************
**  square_field_3d(champ1,champres,coeff)                                   
*/                                                                   
/*! division d'un champ de transfo par u coeff
**  Les champs sont supposes alloues avant l'appel de la fonction          
**  \retval 1 si reussite, 0 si echec
*******************************************************************************/
int square_field_3d(field3d *champ1, field3d *champres)
{
 int    wdth,hght,dpth,i,j,k;
 vector3d ***data1,***datares;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 
 
 data1=champ1->raw;datares=champres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    (datares[i][j][k]).x=(data1[i][j][k]).x*(data1[i][j][k]).x;
    (datares[i][j][k]).y=(data1[i][j][k]).y*(data1[i][j][k]).y;
    (datares[i][j][k]).z=(data1[i][j][k]).z*(data1[i][j][k]).z;
   }
  
 return(1);
}

/*******************************************************************************
**  rootsquare_field_3d(champ1,champres,coeff)                                   
*/                                                                   
/*! division d'un champ de transfo par u coeff
**  Les champs sont supposes alloues avant l'appel de la fonction          
**  \retval 1 si reussite, 0 si echec
*******************************************************************************/
int rootsquare_field_3d(field3d *champ1, field3d *champres)
{
 int    wdth,hght,dpth,i,j,k;
 vector3d ***data1,***datares;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 
 
 data1=champ1->raw;datares=champres->raw;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    (datares[i][j][k]).x=sqrt((data1[i][j][k]).x);
    (datares[i][j][k]).y=sqrt((data1[i][j][k]).y);
    (datares[i][j][k]).z=sqrt((data1[i][j][k]).z);
   }
  
 return(1);
}

/*******************************************************************************
**  decomposer_rotation(matrice_rot, angles_rotation)
*/
/*!
**  retrouve les angles d'une rotation selon Ox, Oy et Oz
**  a partir d'un matrice de rotation
**  /!\ il n'y a pas de controle sur le fait que ce soit bien une matrice de rotation
**  \param matrice_rot : la matrice de rotation
**  \param angles_rotation : le vecteur des angles alloue
**  \retval 0 si reussite
*******************************************************************************/
int decomposer_rotation(double ** matrice_rot, double *angles_rotation)
{
 double tetay=0.0, tetax=0.0, tetaz=0.0;
 double eps=1e-4;

 //si pas de gimbal lock, on peut d�terminer tous les angles
 if (sqrt(matrice_rot[0][1]*matrice_rot[0][1]+matrice_rot[0][0]*matrice_rot[0][0])>eps)
 {
  tetaz=atan2(-matrice_rot[0][1], matrice_rot[0][0]);
  tetax=atan2(-matrice_rot[1][2], matrice_rot[2][2]);
  if (sin(tetaz)>eps) tetay=atan2(-matrice_rot[0][2],-matrice_rot[0][1]/sin(tetaz));
  else tetay=atan2(-matrice_rot[0][2],matrice_rot[0][0]/cos(tetaz));
 }
 //gimbal lock => on pose tetax=0
 else
 {
  tetax=0.0;
  tetaz=atan2(matrice_rot[1][0], matrice_rot[1][1]);
  tetay=-asin(matrice_rot[0][2]);
 }

 angles_rotation[0]=tetax; angles_rotation[1]=tetay; angles_rotation[2]=tetaz;

 return 0;
}

/*******************************************************************************
**  composer_rotation (rotation1, rotation2, rotation_res)
*/
/*!
**  retrouve les angles d'une rotation selon Ox, Oy et Oz
**  issue de la composee de 2 rotations
**  \param rotation1, rotation2 : les rotations a composer
**  \param rotation_res : le vecteur des angles resultant
**  \retval 0 si reussite
*******************************************************************************/
int composer_rotation (double *rotation1, double *rotation2, double *rotation_res)
{
 double **mat_rot1, **mat_rot2, **mat_rot3;
 double rot_res[3]={0.0};
 double tetax,tetay,tetaz;
 double cxx,sxx,cyy,syy,czz,szz;
 int i;

 mat_rot1=alloc_dmatrix(3,3);
 mat_rot2=alloc_dmatrix(3,3);

 tetax=rotation1[0]; tetay=rotation1[1]; tetaz=rotation1[2];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 mat_rot1[0][0]=cyy*czz;              mat_rot1[0][1]=-cyy*szz;             mat_rot1[0][2]=-syy;
 mat_rot1[1][0]=cxx*szz-sxx*syy*czz;  mat_rot1[1][1]=cxx*czz+sxx*syy*szz;  mat_rot1[1][2]=-sxx*cyy;
 mat_rot1[2][0]=sxx*szz+cxx*syy*czz;  mat_rot1[2][1]=sxx*czz-cxx*syy*szz;  mat_rot1[2][2]=cxx*cyy;

 tetax=rotation2[0]; tetay=rotation2[1]; tetaz=rotation2[2];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 mat_rot2[0][0]=cyy*czz;              mat_rot2[0][1]=-cyy*szz;             mat_rot2[0][2]=-syy;
 mat_rot2[1][0]=cxx*szz-sxx*syy*czz;  mat_rot2[1][1]=cxx*czz+sxx*syy*szz;  mat_rot2[1][2]=-sxx*cyy;
 mat_rot2[2][0]=sxx*szz+cxx*syy*czz;  mat_rot2[2][1]=sxx*czz-cxx*syy*szz;  mat_rot2[2][2]=cxx*cyy;

 mat_rot3=multiply_matrix(mat_rot1, 3, 3, mat_rot2, 3, 3);
 decomposer_rotation(mat_rot3, rot_res);

 for (i=0;i<3;i++) rotation_res[i]=rot_res[i];

 free_dmatrix(mat_rot1,3,3); free_dmatrix(mat_rot2,3,3); free_dmatrix(mat_rot3,3,3);

 return 0;
}

/*******************************************************************************
**  comb_transf_3d(transfo2,transfo1,transfores)                                   

**                                                                   
**  \brief Combinaison de deux champs de deformations
**  Les champs sont supposes alloues avant l'appel de la fonction          
**
**  \param transfo2,transfo1 : transformation a combiner
**  \param transfores : la transformation resultante (doit etre allouee)
**  \retval  int : 1
*******************************************************************************/
int comb_transf_3d (transf3d *transfo2,transf3d *transfo1,transf3d *transfores, int interpolation)
{
double rotation1 [3]={0.0}, rotation2 [3] = {0.0}, rotation_res [3] = {0.0};
field3d *ch1,*ch2,*chres;
vector3d ***data1,***data2,***datares;
transf3d *transfotmp;
int wdth,hght,dpth,i,j,k,l,itmp,jtmp,ktmp;
int ix,iy,iz,nb_param;
TRANSFO3D typetrans1,typetrans2;
double tetax1,tetay1,tetaz1,sx1,sy1,sz1,tx1,ty1,tz1,xc1,yc1,zc1,tetax2,tetay2,tetaz2,tx2,ty2,tz2,xc2,yc2,zc2;
double cxx,sxx,cyy,syy,czz,szz;
double a11,a12,a13,a21,a22,a23,a31,a32,a33,b11,b12,b13,b21,b22,b23,b31,b32,b33; 
double tx,ty,tz;//,xc,yc,zc;
double tmpx,tmpy,tmpz,x,y,z;
double *param,*param1,*param2,*p; 
scal_func_t scal_func;
int *x0,*x1,*y0,*y1,*z0,*z1,D;
double *fx,*fy,*fz,px,py,pz,f;
grphic3d *imref,*imreca;
int wdth2,hght2,dpth2;
typetrans1=transfo1->typetrans;
typetrans2=transfo2->typetrans;

imref = cr_grphic3d_modif(1,1,1,0.0,1.0,0);
imref->dx=transfo1->dx ; imref->dy=transfo1->dy ; imref->dz=transfo1->dz;
imref->width=transfo1->width ; imref->height=transfo1->height ; imref->depth=transfo1->depth ;
imreca = cr_grphic3d_modif(1,1,1,0.0,1.0,0);
imreca->dx=transfo2->dx ; imreca->dy=transfo2->dy ; imreca->dz=transfo2->dz;
imreca->width=transfo2->width ; imreca->height=transfo2->height ; imreca->depth=transfo2->depth ;


if (typetrans1==RIGID3D && typetrans2==RIGID3D)
    {
    // ***************************************************
    // Le resultat est une transformation de type RIGID3D 
    // ***************************************************
    
    // recuperation des parametres de la transformation
    tetax1=transfo1->param[0];tetay1=transfo1->param[1];tetaz1=transfo1->param[2];
    tx1=transfo1->param[3];ty1=transfo1->param[4];tz1=transfo1->param[5];
    xc1=transfo1->param[6];yc1=transfo1->param[7];zc1=transfo1->param[8];
    
    tetax2=transfo2->param[0];tetay2=transfo2->param[1];tetaz2=transfo2->param[2];
    tx2=transfo2->param[3];ty2=transfo2->param[4];tz2=transfo2->param[5];
    xc2=transfo2->param[6];yc2=transfo2->param[7];zc2=transfo2->param[8];
        
    // Mise a jour du nombre de parametres et du type de transformation
    transfores->nb_param=9;
    transfores->typetrans=RIGID3D;
  transfores->dx = transfo1->dx;transfores->dy = transfo1->dy;transfores->dz = transfo1->dz;

  if (transfores->param) FREE(transfores->param);
    param=CALLOC(9,double);
    transfores->param=param;
    
    // Mise a jour des angles
 rotation1[0]=tetax1; rotation1[1]=tetay1; rotation1[2]=tetaz1;
 rotation2[0]=tetax2; rotation2[1]=tetay2; rotation2[2]=tetaz2;
 composer_rotation (rotation2, rotation1, rotation_res);
 param[0]=rotation_res[0]; param[1]=rotation_res[1]; param[2]=rotation_res[2];
    
    // Mise a jour du centre de rotation 
    param[6]=xc1;
    param[7]=yc1;
    param[8]=zc1;
    
    // Calcul de la translation resultante
    cxx=cos(tetax2);sxx=sin(tetax2);cyy=cos(tetay2);syy=sin(tetay2);czz=cos(tetaz2);szz=sin(tetaz2);
    a11=cyy*czz;a12=-cyy*szz;a13=-syy;
    a21=cxx*szz-sxx*syy*czz;
    a22=cxx*czz+sxx*syy*szz;
    a23=-sxx*cyy;
    a31=sxx*szz+cxx*syy*czz;
    a32=sxx*czz-cxx*syy*szz;
    a33=cxx*cyy;
    
    tx=tx1+xc1-xc2;
    ty=ty1+yc1-yc2;
    tz=tz1+zc1-zc2;

    param[3]=a11*tx+a12*ty+a13*tz+tx2+xc2-xc1;
    param[4]=a21*tx+a22*ty+a23*tz+ty2+yc2-yc1;
    param[5]=a31*tx+a32*ty+a33*tz+tz2+zc2-zc1;
    
    }
    
else if (typetrans1==RIGIDZOOM3D && typetrans2==RIGID3D)
    {
    // *******************************************************
    // Le resultat est une transformation de type RIGIDZOOM3D 
    // *******************************************************
    
    // recuperation des parametres de la transformation
    tetax1=transfo1->param[0];tetay1=transfo1->param[1];tetaz1=transfo1->param[2];
    sx1=transfo1->param[3];sy1=transfo1->param[4];sz1=transfo1->param[5];
    tx1=transfo1->param[6];ty1=transfo1->param[7];tz1=transfo1->param[8];
    xc1=transfo1->param[9];yc1=transfo1->param[10];zc1=transfo1->param[11];
    
    tetax2=transfo2->param[0];tetay2=transfo2->param[1];tetaz2=transfo2->param[2];
    tx2=transfo2->param[3];ty2=transfo2->param[4];tz2=transfo2->param[5];
    xc2=transfo2->param[6];yc2=transfo2->param[7];zc2=transfo2->param[8];
        
    // Mise a jour du nombre de parametres et du type de transformation
    transfores->nb_param=12;
    transfores->typetrans=RIGIDZOOM3D;
  transfores->dx = transfo1->dx;transfores->dy = transfo1->dy;transfores->dz = transfo1->dz;

  if (transfores->param) FREE(transfores->param);
  param=CALLOC(12,double);
    transfores->param=param;
    
    // Mise a jour des angles
 rotation1[0]=tetax1; rotation1[1]=tetay1; rotation1[2]=tetaz1;
 rotation2[0]=tetax2; rotation2[1]=tetay2; rotation2[2]=tetaz2;
 composer_rotation (rotation2, rotation1, rotation_res);
 param[0]=rotation_res[0]; param[1]=rotation_res[1]; param[2]=rotation_res[2];
    
    // Mise a jour du centre de rotation 
    param[9]=xc1;
    param[10]=yc1;
    param[11]=zc1;
    
    // Mise a jour des coeeficients d'echelle
    param[3]=sx1;
    param[4]=sy1;
    param[5]=sz1;
        
    // Calcul de la translation resultante
    cxx=cos(tetax2);sxx=sin(tetax2);cyy=cos(tetay2);syy=sin(tetay2);czz=cos(tetaz2);szz=sin(tetaz2);
    a11=cyy*czz;a12=-cyy*szz;a13=-syy;
    a21=cxx*szz-sxx*syy*czz;
    a22=cxx*czz+sxx*syy*szz;
    a23=-sxx*cyy;
    a31=sxx*szz+cxx*syy*czz;
    a32=sxx*czz-cxx*syy*szz;
    a33=cxx*cyy;
    
    tx=tx1+xc1-xc2;
    ty=ty1+yc1-yc2;
    tz=tz1+zc1-zc2;
    
    param[6]=a11*tx+a12*ty+a13*tz+tx2+xc2-xc1;
    param[7]=a21*tx+a22*ty+a23*tz+ty2+yc2-yc1;
    param[8]=a31*tx+a32*ty+a33*tz+tz2+zc2-zc1;
    
    }
    
else if ((typetrans1==RIGID3D && typetrans2==RIGIDZOOM3D)||(typetrans1==RIGID3D && typetrans2==AFFINE3D)||(typetrans1==RIGIDZOOM3D &&
typetrans2==RIGIDZOOM3D)||(typetrans1==RIGIDZOOM3D && typetrans2==AFFINE3D)||(typetrans1==AFFINE3D && typetrans2==RIGID3D)||(typetrans1==AFFINE3D &&
typetrans2==RIGIDZOOM3D)||(typetrans1==AFFINE3D && typetrans2==AFFINE3D))   
    {
    // ****************************************************
    // Le resultat est une transformation de type AFFINE3D 
    // ****************************************************
    
    // On se ramene a la composition de deux transformations affines par conversion des differents types de transformations 
    

        param1=CALLOC(15,double);
        for (i=0;i<transfo1->nb_param;i++)
        param1[i]=transfo1->param[i];
    
     switch (typetrans1)
        {
        case RIGID3D: rigid_to_rigidz_3d(param1);rigidz_to_affine_3d(param1);
        break;
        case RIGIDZOOM3D: rigidz_to_affine_3d(param1);
        break;
        //case AFFINE3D:
        default:
        break;
        }
    
        param2=CALLOC(15,double);
        for (i=0;i<transfo2->nb_param;i++)
        param2[i]=transfo2->param[i];
    
    switch (typetrans2)
        {
        case RIGID3D: rigid_to_rigidz_3d(param2);rigidz_to_affine_3d(param2);
        break;
        case RIGIDZOOM3D: rigidz_to_affine_3d(param2);
        break;
        //case AFFINE3D:
        default:
        break;
        }
        
    // recuperation des parametres de la transformation
    a11=param1[0];a12=param1[1];a13=param1[2];
    a21=param1[3];a22=param1[4];a23=param1[5];
    a31=param1[6];a32=param1[7];a33=param1[8];
    tx1=param1[9];ty1=param1[10];tz1=param1[11];
    xc1=param1[12];yc1=param1[13];zc1=param1[14];
    
    
    b11=param2[0];b12=param2[1];b13=param2[2];
    b21=param2[3];b22=param2[4];b23=param2[5];
    b31=param2[6];b32=param2[7];b33=param2[8];
    tx2=param2[9];ty2=param2[10];tz2=param2[11];
    xc2=param2[12];yc2=param2[13];zc2=param2[14];
    
    // Mise a jour du nombre de parametres et du type de transformation
    transfores->nb_param=15;
    transfores->typetrans=AFFINE3D;
    transfores->dx = transfo1->dx;transfores->dy = transfo1->dy;transfores->dz = transfo1->dz;

    if (transfores->param)
        FREE(transfores->param);
    transfores->param=CALLOC(15,double);
    param=transfores->param;
    
    // Mise a jour des coefficients aij
    param[0]=b11*a11+b12*a21+b13*a31;
    param[1]=b11*a12+b12*a22+b13*a32;

    param[2]=b11*a13+b12*a23+b13*a33;
    param[3]=b21*a11+b22*a21+b23*a31;
    param[4]=b21*a12+b22*a22+b23*a32;
    param[5]=b21*a13+b22*a23+b23*a33;
    param[6]=b31*a11+b32*a21+b33*a31;
    param[7]=b31*a12+b32*a22+b33*a32;
    param[8]=b31*a13+b32*a23+b33*a33;
    
    
    // Mise a jour du centre de rotation 
    param[12]=xc1;
    param[13]=yc1;
    param[14]=zc1;
    
    // Calcul de la translation resultante
    tx=tx1+xc1-xc2;
    ty=ty1+yc1-yc2;
    tz=tz1+zc1-zc2;
    
    param[9]=b11*tx+b12*ty+b13*tz+tx2+xc2-xc1;
    param[10]=b21*tx+b22*ty+b23*tz+ty2+yc2-yc1;
    param[11]=b31*tx+b32*ty+b33*tz+tz2+zc2-zc1;
    
    FREE(param1);FREE(param2);

  }
else
    {
    // ***************************************************
    // Le resultat est une transformation de type CHAMP3D 
    // ***************************************************
  //  !attention : gestion de l'anisotropie inexistante dx,dy,dz n'a pas de sens.
  
  // On convertit la premiere transformation en un champ de deplacement (en voxel)
    ch1=transf_to_field_3d(transfo1,imref,imreca);
    
    // Allocation du champ resultat
    wdth=transfo1->width;hght=transfo1->height;dpth=transfo1->depth;
    chres=cr_field3d(wdth,hght,dpth);
    data1=ch1->raw;
    datares=chres->raw;
    
/*  // on converti ch1 en mm
    for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                {
                    data1[i][j][k].x=data1[i][j][k].x*transfo1->dx;
                    data1[i][j][k].y=data1[i][j][k].y*transfo1->dy;
                    data1[i][j][k].z=data1[i][j][k].z*transfo1->dz;
                    
                    }
*/
    
    if ((typetrans2==RIGID3D)||(typetrans2==RIGIDZOOM3D)||(typetrans2==AFFINE3D)||(typetrans2==BSPLINE3D))
        {
        if ((typetrans2==RIGID3D)||(typetrans2==RIGIDZOOM3D)||(typetrans2==AFFINE3D))
        {
    // On traite le cas ou la seconde transformation est rigide, rigidezoom ou affine (on se ramene au cas affine)
      param2=CALLOC(15,double);
      for (i=0;i<transfo2->nb_param;i++)
        param2[i]=transfo2->param[i];

      switch (typetrans2)
      {
       case RIGID3D: rigid_to_rigid_global_zoom_3d(param2);
       case RIGIDGLOBALZOOM3D: rigid_global_zoom_to_rigidz_3d(param2);
       case RIGIDZOOM3D: rigidz_to_affine_decouple_3d(param2);
       case AFFINEDECOUPLE3D: affine_decouple_to_affine_3d(param2);
       case AFFINE3D : affine_to_affinesscg_3d(param2); break;
       default: fprintf (stderr, "SHOULD NOT BE THERE IN comb_transf_3d\n"); return 2;
     }
    
      a11=param2[0] ; a12=param2[1] ; a13=param2[2];
      a21=param2[3] ; a22=param2[4] ; a23=param2[5];
      a31=param2[6] ; a32=param2[7] ; a33=param2[8];
      tx =param2[9] ; ty =param2[10]; tz =param2[11];
//      xc =param2[12]; yc =param2[13]; zc =param2[14];      
    // Calcul du champ de deformation resultant   

      for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                {
            tmpx=(i+(data1[i][j][k]).x)*transfo2->dx;
                tmpy=(j+(data1[i][j][k]).y)*transfo2->dy;
                tmpz=(k+(data1[i][j][k]).z)*transfo2->dz;

                (datares[i][j][k]).x=(float)(a11*tmpx+a12*tmpy+a13*tmpz+tx-i*transfo1->dx);
                (datares[i][j][k]).y=(float)(a21*tmpx+a22*tmpy+a23*tmpz+ty-j*transfo1->dy);
                (datares[i][j][k]).z=(float)(a31*tmpx+a32*tmpy+a33*tmpz+tz-k*transfo1->dz);

                }


            // TODO ACHTUNG!!!!!
    //calcul du champ millimetrique a sauvegarder
        transfotmp=field_to_transf_3d(chres,NULL,NULL);

        transfores->ch=transfotmp->ch;
        transfores->typetrans = CHAMP3D;
        transfores->dx = transfo1->dx;
        transfores->dy = transfo1->dy;
        transfores->dz = transfo1->dz;
        free_field3d(chres);
        free_field3d(ch1);
        free(param2);


        }
        else
        {
        // On traite le cas ou la seconde transformation est du type Bspline
        switch (transfo2->degre)
        {
        case 0: scal_func=Bsplined0;break;
        case 1: scal_func=Bsplined1;break;
        case 2: scal_func=Bsplined2;break;
        default: scal_func=Bsplined0;break;
        }   
        p = CALLOC(transfo2->nb_param,double);  
        nb_param = init_base_3d(wdth, hght, dpth, scal_func);
        for (i = BASE3D.resol; i<transfo2->resol; i++)
          nb_param = base_resol_up_3d(p,nb_param);
        if (nb_param != transfo2->nb_param)
          printf("ERREUR dans transf_to_field_3d\n");
    
         x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
         fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;
         D=(int)floor(TOP_D(nb_param));
        
        for(i=0;i<nb_param;i++)
            p[i]=transfo2->param[i];
    
         
                for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
            for (k=0;k<dpth;k++)
                {
                tmpx=(i+(data1[i][j][k]).x)*transfo1->dx/transfo2->dx;
                    tmpy=(j+(data1[i][j][k]).y)*transfo1->dy/transfo2->dy;
                    tmpz=(k+(data1[i][j][k]).z)*transfo1->dz/transfo2->dz;
                
                    datares[i][j][k].x=(data1[i][j][k]).x*transfo1->dx;
                    datares[i][j][k].y=(data1[i][j][k]).y*transfo1->dy;
                    datares[i][j][k].z=(data1[i][j][k]).z*transfo1->dz;
//              printf("%d %d %d\r",i,j,k);                 
                    
                    ix=(int)floor(1.0*tmpx*(D+1)/wdth);
                    iy=(int)floor(1.0*tmpy*(D+1)/hght);
                    iz=(int)floor(1.0*tmpz*(D+1)/dpth);
                    
                    for (itmp=MAXI(ix-1,0);itmp<MINI(ix+1,D);itmp++)
                    for (jtmp=MAXI(iy-1,0);jtmp<MINI(iy+1,D);jtmp++)
                    for (ktmp=MAXI(iz-1,0);ktmp<MINI(iz+1,D);ktmp++)
                        {
                        l=(int)floor(TOP_conv_ind(itmp,jtmp,ktmp,nb_param)/3.);
                        px=p[3*l];py=p[3*l+1];pz=p[3*l+2];
                 
f=Bspline(transfo2->degre,2.0*((tmpx-x0[l])/(x1[l]-x0[l])))*Bspline(transfo2->degre,2.0*((tmpy-y0[l])/(y1[l]-y0[l])))*Bspline(transfo2->degre,2.0*((tmpz-z0[l])/(z1[l]-z0[l])));

                datares[i][j][k].x=(float)(datares[i][j][k].x+px*f);
                datares[i][j][k].y=(float)(datares[i][j][k].y+py*f);
                datares[i][j][k].z=(float)(datares[i][j][k].z+pz*f);       
                        }
                        
                    }
                
        end_base_3d();
    // TODO ACHTUNG!!!!
        transfotmp=field_to_transf_3d(chres,NULL,NULL);
        transfores->ch=transfotmp->ch;
    transfores->typetrans = CHAMP3D;
    transfores->dx = transfo1->dx;
    transfores->dy = transfo1->dy;
    transfores->dz = transfo1->dz;
        free_field3d(ch1);
        free_field3d(chres);
        free(p);
    }
        }
    else
        {
        // Cas Cham3D + Champ3D
        
            // On convertit la seconde transformation en un champ de deplacement
            ch2=transf_to_field_3d(transfo2,NULL,NULL);
            data2=ch2->raw;
            
            wdth2=ch2->width;
            hght2=ch2->height;
            dpth2=ch2->depth;
                 
        
            switch (interpolation)
            {
            case 0:
            // Methode du plus proche voisin
        for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                {
                //ix=(int)floor((i+data1[i][j][k].x)*transfo1->dx/transfo2->dx);
        //iy=(int)floor((j+data1[i][j][k].y)*transfo1->dy/transfo2->dy);
        //iz=(int)floor((k+data1[i][j][k].z)*transfo1->dz/transfo2->dz);
                
                ix=(int)floor(i+data1[i][j][k].x);
        iy=(int)floor(j+data1[i][j][k].y);
        iz=(int)floor(k+data1[i][j][k].z);
                
                
                if ((ix>=0)&&(ix<=wdth2-1)&&(iy>=0)&&(iy<=hght2-1)&&(iz>=0)&&(iz<=dpth2-1))
                    {
                    /*datares[i][j][k].x=data1[i][j][k].x*transfo1->dx+data2[ix][iy][iz].x;
                    datares[i][j][k].y=data1[i][j][k].y*transfo1->dy+data2[ix][iy][iz].y;
                    datares[i][j][k].z=data1[i][j][k].z*transfo1->dz+data2[ix][iy][iz].z;*/
                
                    //datares[i][j][k].x=(ch1->raw[i][j][k].x+i)*transfo2->dx-i*transfo1->dx+data2[ix][iy][iz].x;
                    //datares[i][j][k].y=(ch1->raw[i][j][k].y+j)*transfo2->dy-j*transfo1->dy+data2[ix][iy][iz].y;
                    //datares[i][j][k].z=(ch1->raw[i][j][k].z+k)*transfo2->dz-k*transfo1->dz+data2[ix][iy][iz].z;
                    
                    datares[i][j][k].x=(ch1->raw[i][j][k].x+i)*transfo2->dx-i*transfo1->dx+data2[ix][iy][iz].x;
                    datares[i][j][k].y=(ch1->raw[i][j][k].y+j)*transfo2->dy-j*transfo1->dy+data2[ix][iy][iz].y;
                    datares[i][j][k].z=(ch1->raw[i][j][k].z+k)*transfo2->dz-k*transfo1->dz+data2[ix][iy][iz].z;
                    
                    }
                else 
                    {
                    /*datares[i][j][k].x=data1[i][j][k].x*transfo1->dx;
                    datares[i][j][k].y=data1[i][j][k].y*transfo1->dy;
                    datares[i][j][k].z=data1[i][j][k].z*transfo1->dz;*/
                
                    //datares[i][j][k].x=(ch1->raw[i][j][k].x+i)*transfo2->dx-i*transfo1->dx;
                    //datares[i][j][k].y=(ch1->raw[i][j][k].y+j)*transfo2->dy-j*transfo1->dy;
                    //datares[i][j][k].z=(ch1->raw[i][j][k].z+k)*transfo2->dz-k*transfo1->dz;
                    
                    
                    // C'est un peu de la bidouille, mais au moins on est sur que la transformation sort du FOV
                    datares[i][j][k].x=10*wdth; 
                    datares[i][j][k].y=10*hght;
                    datares[i][j][k].z=10*dpth;
                    
                    
                    }
                }
                break;
                case 1:
                // Interpolation lineaire
        /*for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                {
                x=(i+data1[i][j][k].x)*transfo1->dx/transfo2->dx;y=(j+data1[i][j][k].y)*transfo1->dy/transfo2->dy;z=(k+data1[i][j][k].z)*transfo1->dz/transfo2->dz;
                ix=(int)floor(x);
        iy=(int)floor(y);
        iz=(int)floor(z);
                tmpx=x-ix;tmpy=y-iy;tmpz=z-iz;
                if ((ix>=0)&&(ix<=wdth-2)&&(iy>=0)&&(iy<=hght-2)&&(iz>=0)&&(iz<=dpth-2))
                    {
                    datares[i][j][k].x=(float)(data1[i][j][k].x*transfo1->dx+tmpx*tmpy*tmpz*data2[ix][iy][iz].x
                                                       +(1-tmpx)*tmpy*tmpz*data2[ix+1][iy][iz].x
                                                       +(1-tmpx)*(1-tmpy)*tmpz*data2[ix+1][iy+1][iz].x
                                                       +(1-tmpx)*tmpy*(1-tmpz)*data2[ix+1][iy][iz+1].x
                                                       +(1-tmpx)*(1-tmpy)*(1-tmpz)*data2[ix+1][iy+1][iz+1].x
                                                       +tmpx*(1-tmpy)*(1-tmpz)*data2[ix][iy+1][iz+1].x
                                                       +tmpx*tmpy*(1-tmpz)*data2[ix][iy][iz+1].x
                                                       +tmpx*(1-tmpy)*tmpz*data2[ix][iy+1][iz].x);
                    
                    datares[i][j][k].y=(float)(data1[i][j][k].y*transfo1->dy+tmpx*tmpy*tmpz*data2[ix][iy][iz].y
                                                       +(1-tmpx)*tmpy*tmpz*data2[ix+1][iy][iz].y
                                                       +(1-tmpx)*(1-tmpy)*tmpz*data2[ix+1][iy+1][iz].y
                                                       +(1-tmpx)*tmpy*(1-tmpz)*data2[ix+1][iy][iz+1].y
                                                       +(1-tmpx)*(1-tmpy)*(1-tmpz)*data2[ix+1][iy+1][iz+1].y
                                                       +tmpx*(1-tmpy)*(1-tmpz)*data2[ix][iy+1][iz+1].y
                                                       +tmpx*tmpy*(1-tmpz)*data2[ix][iy][iz+1].y
                                                       +tmpx*(1-tmpy)*tmpz*data2[ix][iy+1][iz].y);
                                                
                    datares[i][j][k].z=(float)(data1[i][j][k].z*transfo1->dz+tmpx*tmpy*tmpz*data2[ix][iy][iz].z
                                                       +(1-tmpx)*tmpy*tmpz*data2[ix+1][iy][iz].z
                                                       +(1-tmpx)*(1-tmpy)*tmpz*data2[ix+1][iy+1][iz].z
                                                       +(1-tmpx)*tmpy*(1-tmpz)*data2[ix+1][iy][iz+1].z
                                                       +(1-tmpx)*(1-tmpy)*(1-tmpz)*data2[ix+1][iy+1][iz+1].z
                                                       +tmpx*(1-tmpy)*(1-tmpz)*data2[ix][iy+1][iz+1].z
                                                       +tmpx*tmpy*(1-tmpz)*data2[ix][iy][iz+1].z
                                                       +tmpx*(1-tmpy)*tmpz*data2[ix][iy+1][iz].z);
                    }
                else 
                    {
                    datares[i][j][k].x=data1[i][j][k].x*transfo1->dx;
                    datares[i][j][k].y=data1[i][j][k].y*transfo1->dy;
                    datares[i][j][k].z=data1[i][j][k].z*transfo1->dz;
                    }*/
                    for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                {
                x=i+data1[i][j][k].x;
        y=j+data1[i][j][k].y;
        z=k+data1[i][j][k].z;
                
                ix=(int)floor(x);
        iy=(int)floor(y);
        iz=(int)floor(z);
                
                tmpx=x-ix;tmpy=y-iy;tmpz=z-iz;
                if ((ix>=0)&&(ix<=wdth2-2)&&(iy>=0)&&(iy<=hght2-2)&&(iz>=0)&&(iz<=dpth2-2))
                    {
                    datares[i][j][k].x=(float)((ch1->raw[i][j][k].x+i)*transfo2->dx-i*transfo1->dx+tmpx*tmpy*tmpz*data2[ix][iy][iz].x
                                                       +(1-tmpx)*tmpy*tmpz*data2[ix+1][iy][iz].x
                                                       +(1-tmpx)*(1-tmpy)*tmpz*data2[ix+1][iy+1][iz].x
                                                       +(1-tmpx)*tmpy*(1-tmpz)*data2[ix+1][iy][iz+1].x
                                                       +(1-tmpx)*(1-tmpy)*(1-tmpz)*data2[ix+1][iy+1][iz+1].x
                                                       +tmpx*(1-tmpy)*(1-tmpz)*data2[ix][iy+1][iz+1].x
                                                       +tmpx*tmpy*(1-tmpz)*data2[ix][iy][iz+1].x
                                                       +tmpx*(1-tmpy)*tmpz*data2[ix][iy+1][iz].x);
                    
                    datares[i][j][k].y=(float)((ch1->raw[i][j][k].y+j)*transfo2->dy-j*transfo1->dy+tmpx*tmpy*tmpz*data2[ix][iy][iz].y
                                                       +(1-tmpx)*tmpy*tmpz*data2[ix+1][iy][iz].y
                                                       +(1-tmpx)*(1-tmpy)*tmpz*data2[ix+1][iy+1][iz].y
                                                       +(1-tmpx)*tmpy*(1-tmpz)*data2[ix+1][iy][iz+1].y
                                                       +(1-tmpx)*(1-tmpy)*(1-tmpz)*data2[ix+1][iy+1][iz+1].y
                                                       +tmpx*(1-tmpy)*(1-tmpz)*data2[ix][iy+1][iz+1].y
                                                       +tmpx*tmpy*(1-tmpz)*data2[ix][iy][iz+1].y
                                                       +tmpx*(1-tmpy)*tmpz*data2[ix][iy+1][iz].y);
                                                
                    datares[i][j][k].z=(float)((ch1->raw[i][j][k].z+k)*transfo2->dz-k*transfo1->dz+tmpx*tmpy*tmpz*data2[ix][iy][iz].z
                                                       +(1-tmpx)*tmpy*tmpz*data2[ix+1][iy][iz].z
                                                       +(1-tmpx)*(1-tmpy)*tmpz*data2[ix+1][iy+1][iz].z
                                                       +(1-tmpx)*tmpy*(1-tmpz)*data2[ix+1][iy][iz+1].z
                                                       +(1-tmpx)*(1-tmpy)*(1-tmpz)*data2[ix+1][iy+1][iz+1].z
                                                       +tmpx*(1-tmpy)*(1-tmpz)*data2[ix][iy+1][iz+1].z
                                                       +tmpx*tmpy*(1-tmpz)*data2[ix][iy][iz+1].z
                                                       +tmpx*(1-tmpy)*tmpz*data2[ix][iy+1][iz].z);
                    }
                else 
                    {
                    //datares[i][j][k].x=(ch1->raw[i][j][k].x+i)*transfo2->dx-i*transfo1->dx;
                    //datares[i][j][k].y=(ch1->raw[i][j][k].y+j)*transfo2->dy-j*transfo1->dy;
                    //datares[i][j][k].z=(ch1->raw[i][j][k].z+k)*transfo2->dz-k*transfo1->dz;
                    
                    // C'est un peu de la bidouille, mais au moins on est sur que la transformation sort du FOV
                    datares[i][j][k].x=10*wdth; 
                    datares[i][j][k].y=10*hght;
                    datares[i][j][k].z=10*dpth;
    
                    }
                }
                break;
                case 2:
                Inter_Chps_Bspline_3d(ch2,ch1,chres,2,transfo2->dx,transfo2->dy,transfo2->dz,transfo1->dx,transfo1->dy,transfo1->dz);
                break;
                case 3:
                Inter_Chps_Bspline_3d(ch2,ch1,chres,3,transfo2->dx,transfo2->dy,transfo2->dz,transfo1->dx,transfo1->dy,transfo1->dz);
                break;
                case 4:
                Inter_Chps_Bspline_3d(ch2,ch1,chres,4,transfo2->dx,transfo2->dy,transfo2->dz,transfo1->dx,transfo1->dy,transfo1->dz);
                break;
                case 5:
                Inter_Chps_Bspline_3d(ch2,ch1,chres,5,transfo2->dx,transfo2->dy,transfo2->dz,transfo1->dx,transfo1->dy,transfo1->dz);
                break;
                default:break;
            }
                /*
                for (i=0;i<wdth;i++)
            for (j=0;j<hght;j++)
                for (k=0;k<dpth;k++)
                        {
                        if(abs(datares[i][j][k].x)>wdth*transfo1->dx)
                         datares[i][j][k].x=wdth*transfo1->dx;
                        if(abs(datares[i][j][k].y)>hght*transfo1->dy)
                         datares[i][j][k].y=hght*transfo1->dy;
                         if(abs(datares[i][j][k].z)>dpth*transfo1->dz)
                         datares[i][j][k].z=dpth*transfo1->dz;
                        }
                */      
                        
            transfotmp=field_to_transf_3d(chres,NULL,NULL);
            transfores->ch=transfotmp->ch;
      transfores->typetrans = CHAMP3D;
      transfores->dx = transfo1->dx;
      transfores->dy = transfo1->dy;
      transfores->dz = transfo1->dz;
            free_field3d(chres);
            free_field3d(ch1);
            free_field3d(ch2);
            

        }
        
        
    
    }
free_grphic3d(imref); 
free_grphic3d(imreca); 

    return(1);
}
/*******************************************************************************
**  Inter_Chps_Bspline_3d(ch1,ch2,chres,degre)                                   
**       Fais la composition de deux champs de transformation 
**      avec un interpolation Bspline de de degre n                            
*******************************************************************************/
void Inter_Chps_Bspline_3d(field3d *ch1, field3d *ch2, field3d *chres,int degre,double dx1, double dy1, double dz1, double dx2, double dy2, double dz2)
{
int wdth,hght,dpth,wdth2,hght2,dpth2,nbpoles;
int i,j,k,l,m,n;
double x,y,z,w,w2,w4,t,t0,t1;
vector3d vint,sumx,sumxy;
field3d *coeff;
double *row1Dx,*row1Dy,*row1Dz;

double epsilon=2.2204460492503131E-16;
double poles[2],xWeight[6],yWeight[6],zWeight[6]; 
int xIndex[6],yIndex[6],zIndex[6];

    
wdth=ch1->width;hght=ch1->height;dpth=ch1->depth;
wdth2=ch2->width;hght2=ch2->height;dpth2=ch2->depth;

// Calcul des coefficient d'interpolation a partir des echantillons de l'image de depart
coeff=cr_field3d(wdth, hght, dpth);
    
    // Recuperation de la valeur des poles en z
    switch (degre) 
    {
    case 2:
        poles[0]=(sqrt(8.0) - 3.0);
        nbpoles=1;
        break;
    case 3:
        poles[0]=(sqrt(3.0) - 2.0);
        nbpoles=1;
        break;
    case 4:
        poles[0]=(sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0);
        poles[1]=(sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0);
        nbpoles=2;
        break;
    case 5:
        poles[0]=(sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)- 13.0 / 2.0);
        poles[1]=(sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)- 13.0 / 2.0);
        nbpoles=2;
        break;
    default:
        {printf("Degre de spline invalide\n");return;}
    }

        
        // Calcul des coeeficients suivant x
        row1Dx=(double*)malloc(wdth*sizeof(double));
        row1Dy=(double*)malloc(wdth*sizeof(double));
        row1Dz=(double*)malloc(wdth*sizeof(double));
        
        
        for (k=0; k<dpth; k++)
        {
            for (j=0; j<hght; j++)
            {
                for (i=0; i<wdth; i++)
                    {
                    row1Dx[i]=ch1->raw[i][j][k].x;
                    row1Dy[i]=ch1->raw[i][j][k].y;
                    row1Dz[i]=ch1->raw[i][j][k].z;
                    }
                GetInterpolationCoefficients(row1Dx,wdth,poles,nbpoles,epsilon);
                GetInterpolationCoefficients(row1Dy,wdth,poles,nbpoles,epsilon);
                GetInterpolationCoefficients(row1Dz,wdth,poles,nbpoles,epsilon);                            
                for (i=0; i<wdth; i++)
                    {
                    coeff->raw[i][j][k].x=(float)row1Dx[i];
                    coeff->raw[i][j][k].y=(float)row1Dy[i];
                    coeff->raw[i][j][k].z=(float)row1Dz[i];
                    }
            }
        }
        free(row1Dx);free(row1Dy);free(row1Dz);
    
    
        // Calcul des coeeficients suivant y

        row1Dx=(double*)malloc(hght*sizeof(double));
        row1Dy=(double*)malloc(hght*sizeof(double));
        row1Dz=(double*)malloc(hght*sizeof(double));
        
                
        for (k=0; k<dpth; k++)
        {

            for (i=0; i<wdth; i++)
            {
                for (j=0; j<hght; j++)
                    {
                    row1Dx[j]=coeff->raw[i][j][k].x;
                    row1Dy[j]=coeff->raw[i][j][k].y;
                    row1Dz[j]=coeff->raw[i][j][k].z;
                    }
                GetInterpolationCoefficients(row1Dx,hght,poles,nbpoles,epsilon);
                GetInterpolationCoefficients(row1Dy,hght,poles,nbpoles,epsilon);
                GetInterpolationCoefficients(row1Dz,hght,poles,nbpoles,epsilon);
                for (j=0; j<hght; j++)
                {
                coeff->raw[i][j][k].x=(float)row1Dx[j];
                coeff->raw[i][j][k].y=(float)row1Dy[j];
                coeff->raw[i][j][k].z=(float)row1Dz[j];
                }
            }
        }
        free(row1Dx);free(row1Dy);free(row1Dz);
    
        // Calcul des coeeficients suivant z
        
        row1Dx=(double*)malloc(dpth*sizeof(double));
        row1Dy=(double*)malloc(dpth*sizeof(double));
        row1Dz=(double*)malloc(dpth*sizeof(double));
        
        
        for (i=0; i<wdth; i++)
        {
            for (j=0; j<hght; j++)
            {
                for (k=0; k<dpth; k++)
                {
                row1Dx[k]=coeff->raw[i][j][k].x;
                row1Dy[k]=coeff->raw[i][j][k].y;
                row1Dz[k]=coeff->raw[i][j][k].z;
                }
                GetInterpolationCoefficients(row1Dx,dpth,poles,nbpoles,epsilon);
                GetInterpolationCoefficients(row1Dy,dpth,poles,nbpoles,epsilon);
                GetInterpolationCoefficients(row1Dz,dpth,poles,nbpoles,epsilon);
                for (k=0; k<dpth; k++)
                {
                coeff->raw[i][j][k].x=(float)row1Dx[k];
                coeff->raw[i][j][k].y=(float)row1Dy[k];
                coeff->raw[i][j][k].z=(float)row1Dz[k];
                }
            }
        }
        free(row1Dx);free(row1Dy);free(row1Dz);
        
        for (i=0;i<wdth2;i++)
            for (j=0;j<hght2;j++)
                for (k=0;k<dpth2;k++)
                {
                /*x=(i+ch2->raw[i][j][k].x)*dx2/dx1;
                y=(j+ch2->raw[i][j][k].y)*dy2/dy1;
                z=(k+ch2->raw[i][j][k].z)*dz2/dz1;*/
                x=i+ch2->raw[i][j][k].x;
        y=j+ch2->raw[i][j][k].y;
        z=k+ch2->raw[i][j][k].z;
                
                if (x<0 || x>=wdth || y<0 || y>=hght || z<0 || z>=dpth)  
                    {
                    /*chres->raw[i][j][k].x=(float)(ch2->raw[i][j][k].x*dx2);
                    chres->raw[i][j][k].y=(float)(ch2->raw[i][j][k].y*dy2);
                    chres->raw[i][j][k].z=(float)(ch2->raw[i][j][k].z*dz2);*/
                    //chres->raw[i][j][k].x=(float)((ch2->raw[i][j][k].x+i)*dx1-i*dx2);
                    //chres->raw[i][j][k].y=(float)((ch2->raw[i][j][k].y+j)*dy1-j*dy2);
                    //chres->raw[i][j][k].z=(float)((ch2->raw[i][j][k].z+k)*dz1-k*dz2);
                    
                    // C'est un peu de la bidouille, mais au moins on est sur que la transformation sort du FOV
                    chres->raw[i][j][k].x=10*wdth; 
                    chres->raw[i][j][k].y=10*hght;
                    chres->raw[i][j][k].z=10*dpth;

                    }
                else
                    {
                    
                    //***********************************************
                    // Calcul de la valeur interpolee vint en (x,y,z)
                    //**********************************************/
                        
                    
                    xIndex[0]=(int)floor(x+0.5)-degre/2;
                    yIndex[0]=(int)floor(y+0.5)-degre/2;
                    zIndex[0]=(int)floor(z+0.5)-degre/2;
                    
                    for (l=0;l<degre;l++)
                        {
                        xIndex[l+1]=xIndex[l]+1;
                        yIndex[l+1]=yIndex[l]+1;
                        zIndex[l+1]=zIndex[l]+1;
                        }
                        
                    
                    // Calcul des poids de la fonction d'interpolation
                    
                    switch (degre) 
                    {
                    case 2:
                        /* x */
                        w = x - (double)xIndex[1];
                        xWeight[1] = 3.0 / 4.0 - w * w;
                        xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
                        xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
                        /* y */
                        w = y - (double)yIndex[1];
                        yWeight[1] = 3.0 / 4.0 - w * w;
                        yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
                        yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
                        /* z */
                        w = z - (double)zIndex[1];
                        zWeight[1] = 3.0 / 4.0 - w * w;
                        zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
                        zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
                        break;
                    case 3:
                        /* x */
                        w = x - (double)xIndex[1];
                        xWeight[3] = (1.0 / 6.0) * w * w * w;
                        xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
                        xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
                        xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
                        /* y */
                        w = y - (double)yIndex[1];
                        yWeight[3] = (1.0 / 6.0) * w * w * w;
                        yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
                        yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
                        yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
                        /* z */
                        w = z - (double)zIndex[1];
                        zWeight[3] = (1.0 / 6.0) * w * w * w;
                        zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
                        zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
                        zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
                        break;
                    case 4:
                        /* x */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        t = (1.0 / 6.0) * w2;
                        xWeight[0] = 1.0 / 2.0 - w;
                        xWeight[0] *= xWeight[0];
                        xWeight[0] *= (1.0 / 24.0) * xWeight[0];
                        t0 = w * (t - 11.0 / 24.0);
                        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
                        xWeight[1] = t1 + t0;
                        xWeight[3] = t1 - t0;
                        xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
                        xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
                        /* y */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        t = (1.0 / 6.0) * w2;
                        yWeight[0] = 1.0 / 2.0 - w;
                        yWeight[0] *= yWeight[0];
                        yWeight[0] *= (1.0 / 24.0) * yWeight[0];
                        t0 = w * (t - 11.0 / 24.0);
                        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
                        yWeight[1] = t1 + t0;
                        yWeight[3] = t1 - t0;
                        yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
                        yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
                        /* z */
                        w = z - (double)zIndex[2];
                        w2 = w * w;
                        t = (1.0 / 6.0) * w2;
                        zWeight[0] = 1.0 / 2.0 - w;
                        zWeight[0] *= zWeight[0];
                        zWeight[0] *= (1.0 / 24.0) * zWeight[0];
                        t0 = w * (t - 11.0 / 24.0);
                        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
                        zWeight[1] = t1 + t0;
                        zWeight[3] = t1 - t0;
                        zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
                        zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
                        break;
                    case 5:
                        /* x */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
                        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
                        t1 = (-1.0 / 12.0) * w * (t + 4.0);
                        xWeight[2] = t0 + t1;
                        xWeight[3] = t0 - t1;
                        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
                        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
                        xWeight[1] = t0 + t1;
                        xWeight[4] = t0 - t1;
                        /* y */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
                        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
                        t1 = (-1.0 / 12.0) * w * (t + 4.0);
                        yWeight[2] = t0 + t1;
                        yWeight[3] = t0 - t1;
                        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
                        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
                        yWeight[1] = t0 + t1;
                        yWeight[4] = t0 - t1;
                        /* z */
                        w = z - (double)zIndex[2];
                        w2 = w * w;
                        zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
                        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
                        t1 = (-1.0 / 12.0) * w * (t + 4.0);
                        zWeight[2] = t0 + t1;
                        zWeight[3] = t0 - t1;
                        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
                        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
                        zWeight[1] = t0 + t1;
                        zWeight[4] = t0 - t1;
                        break;
                    default:
                    printf("le degre de la spline est invalide\n");
                    }
                    
                    // Condition limite (mirror boundaries)
                    /*for (l=0;l<=degre;l++)
                        {
                        if (xIndex[l]<0)
                            xIndex[l]=0;
                        if (xIndex[l]>=wdth)
                            xIndex[l]=0;
                        if (yIndex[l]<0)
                            yIndex[l]=0;
                        if (yIndex[l]>=hght)
                            yIndex[l]=0;
                        if (zIndex[l]<0)
                            zIndex[l]=0;
                        if (zIndex[l]>=dpth)
                            zIndex[l]=0;
                        }*/
                    for (l=0;l<=degre;l++)
                        {
                        if (xIndex[l]<0)
                            xIndex[l]=-xIndex[l];
                        if (xIndex[l]>=wdth)
                            xIndex[l]=(int)floor(2.0*(wdth-1))-xIndex[l];
                        if (yIndex[l]<0)
                            yIndex[l]=-yIndex[l];
                        if (yIndex[l]>=hght)
                            yIndex[l]=(int)floor(2.0*(hght-1))-yIndex[l];
                        if (zIndex[l]<0)
                            zIndex[l]=-zIndex[l];
                        if (zIndex[l]>=dpth)
                            zIndex[l]=(int)floor(2.0*(dpth-1))-zIndex[l];
                        }
                    
                    
                    vint.x=0;
                    vint.y=0;
                    vint.z=0;
                    
                    for (n=0;n<=degre;n++) 
                        {
                        sumxy.x=0;
                        sumxy.y=0;
                        sumxy.z=0;
                        for (m=0;m<=degre;m++)
                            {
                            sumx.x=0;
                            sumx.y=0;
                            sumx.z=0;
                            for (l=0;l<=degre;l++)
                                {
                                sumx.x+=(float)(xWeight[l]*coeff->raw[xIndex[l]][yIndex[m]][zIndex[n]].x);
                                sumx.y+=(float)(xWeight[l]*coeff->raw[xIndex[l]][yIndex[m]][zIndex[n]].y);
                                sumx.z+=(float)(xWeight[l]*coeff->raw[xIndex[l]][yIndex[m]][zIndex[n]].z);
                                }
                            sumxy.x+=(float)(yWeight[m]*sumx.x);
                            sumxy.y+=(float)(yWeight[m]*sumx.y);
                            sumxy.z+=(float)(yWeight[m]*sumx.z);
                            }
                        vint.x+=(float)(zWeight[n]*sumxy.x);
                        vint.y+=(float)(zWeight[n]*sumxy.y);
                        vint.z+=(float)(zWeight[n]*sumxy.z);
                        
                        } 
                                        
                    
                    /*chres->raw[i][j][k].x=(float)(ch2->raw[i][j][k].x*dx2+vint.x);
                    chres->raw[i][j][k].y=(float)(ch2->raw[i][j][k].y*dy2+vint.y);
                    chres->raw[i][j][k].z=(float)(ch2->raw[i][j][k].z*dz2+vint.z);*/
                    chres->raw[i][j][k].x=(float)((ch2->raw[i][j][k].x+i)*dx1-i*dx2+vint.x);
                    chres->raw[i][j][k].y=(float)((ch2->raw[i][j][k].y+j)*dy1-j*dy2+vint.y);
                    chres->raw[i][j][k].z=(float)((ch2->raw[i][j][k].z+k)*dz1-k*dz2+vint.z);
                        }
                }
 
free_field3d(coeff);
}



/*******************************************************************************
**  Bspline(degre,valeur)                                   
**                                                                   
**  Renvoie la valeur en continu de la spline en une valeur donnee
**  (la spline de degre 3 correspond a Bspline0d1)
**   La valeuren laquelle la fonction est evaluee est une valeur normalisee dans
**  l'intervalle [0 degre+1]         
*******************************************************************************/

double Bspline(int degre,double valeur)
{
double res;
switch (degre)
    {
    case 0:
        if ((valeur>=0)&&(valeur<=1)) 
        res=1;
        else
        res=0;   
    break;
    case 1:
        if ((valeur>=0)&&(valeur<=1)) 
        res=valeur;
        else if ((valeur>1)&&(valeur<=2))
        res=2-valeur;
        else 
        res=0;
    break;
    case 2:
        if ((valeur>=0)&&(valeur<=1))
        res=valeur*valeur/2;
        else if ((valeur>1)&&(valeur<=2))
        res=3*valeur-valeur*valeur-3/2;
        else if ((valeur>2)&&(valeur<=3))
        res=valeur*valeur/2-3*valeur+9/2;
        else res=0;   
    break;
    case 3:
        if ((valeur>=0)&&(valeur<=0.5)) 
        res=-2*valeur;
        else if ((valeur>0.5)&&(valeur<=1))
        res=4*valeur-3;
        if ((valeur>1)&&(valeur<=1.5)) 
        res=5-4*valeur;
        else if ((valeur>1.5)&&(valeur<=2))
        res=2*valeur-4;
        else 
        res=0;
    break;
    default: // par defaut Bspline d'ordre 1
        if ((valeur>=0)&&(valeur<=1)) 
        res=valeur;
        else if ((valeur>1)&&(valeur<=2))
        res=2-valeur;
        else 
        res=0;
    break;
    }  
return(res);
}

/*******************************************************************************
********************************************************************************
************************** COMPARER CHAMP    ***********************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**     compare_field_3d()               
**                                                                  
**     Compare deux champ contenu dans des fichiers et affiche la difference
**  de module et la difference de phase dans deux image                          
*******************************************************************************/
void compare_field_3d(void)
{
 field3d *ch1,*ch2;
 transf3d *transfo;
 char   nomfich1[500],nomfich2[500];
 int   im_mod,im_phase,wdth,hght,dpth,e=0   ;


 /*lecture du premier champ*/
 strcpy(nomfich1,GET_FILE("*.trf",&e));
 put_file_extension(nomfich1,".trf",nomfich1);

 if ((test_fichier_transf_3d(nomfich1))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}

 /*lecture du deuxieme champ*/
 strcpy(nomfich2,GET_FILE("*.trf",&e));
 put_file_extension(nomfich2,".trf",nomfich2);

 if ((test_fichier_transf_3d(nomfich2))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}

 /*image de difference en module*/
 im_mod=GET_PLACE3D(TEXT0129);
 im_phase=GET_PLACE3D(TEXT0130);

 /*lecture des champ*/
 transfo=load_transf_3d(nomfich1);
 ch1=transf_to_field_3d(transfo,NULL,NULL);
 free_transf3d(transfo);

 transfo=load_transf_3d(nomfich2);
 ch2=transf_to_field_3d(transfo,NULL,NULL);
 free_transf3d(transfo);

 /*verification de la compatibilite de taille des champs*/
 wdth=ch1->width;hght=ch1->height;dpth=ch1->depth;
 if (ch2->width!=wdth || ch2->height!=hght || ch2->depth!=dpth)
 {
  PUT_ERROR("The two fields have not the same size");
  free_field3d(ch1);free_field3d(ch2);
  return;
 }
 else
 {
  imx_compare_field_3d_p(ch1,ch2,ptr_img_3d(im_mod),ptr_img_3d(im_phase));
  show_picture_3d(im_mod);show_picture_3d(im_phase);
 }
 free_field3d(ch1);free_field3d(ch2);
}

/*******************************************************************************
**     imx_compare_field_3d_p()             
**                                                                  
**     Compare deux champ contenu et donne le resultat de la difference en module
**  et en phase sous forme de deux images                          
*******************************************************************************/
int imx_compare_field_3d_p(field3d *ch1, field3d *ch2, grphic3d *immod, grphic3d *imphase)
{
 vector3d ***data1,***data2;
 vector3d v1,v2;
 double   norm1,norm2,maxmod,maxphase,rmod,rphase,c;
 double   ***module,***phase;
 int    wdth,hght,dpth,i,j,k;
 
 c=180.0/(double)PI;
 
 data1=ch1->raw;data2=ch2->raw;
 wdth=ch1->width;hght=ch1->height;dpth=ch1->depth;
 
 module=alloc_dmatrix_3d(wdth,hght,dpth);
 phase=alloc_dmatrix_3d(wdth,hght,dpth);
 
 maxmod=maxphase=0.0; 
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     v1=data1[i][j][k];v2=data2[i][j][k];
     norm1=sqrt(v1.x*v1.x+v1.y*v1.y+v1.z*v1.z);
     norm2=sqrt(v2.x*v2.x+v2.y*v2.y+v2.z*v2.z);
     module[i][j][k]=fabs(norm1-norm2);
     phase[i][j][k]=c*acos((v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)/norm1/norm2);
     if (module[i][j][k]>maxmod) maxmod=module[i][j][k];
     if (phase[i][j][k]>maxphase) maxphase=phase[i][j][k];
    }
 
 imx_brukermax_3d((float)maxmod,0.0,immod);
 imx_brukermax_3d((float)maxphase,0.0,imphase);
 immod->width=wdth;immod->height=hght;immod->depth=dpth;
 imphase->width=wdth;imphase->height=hght;imphase->depth=dpth;
 rmod=immod->rcoeff;rphase=imphase->rcoeff;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     immod->mri[i][j][k]=(int)(module[i][j][k]/rmod);
     imphase->mri[i][j][k]=(int)(phase[i][j][k]/rphase);
    }
 
 
 free_dmatrix_3d(module);free_dmatrix_3d(phase);
 return(1);
}


/*******************************************************************************
********************************************************************************
************************** CHAMPS  -> IMAGES ***********************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**     module_field_to_image_3d(ch,im)              
**                                                                  
**     module d'un champ dans une image                             
*******************************************************************************/
int module_field_to_image_3d(field3d *ch, grphic3d *im)
{
  int i,j,k,wdth,hght,dpth;
  double maxmod,minmod,m,cx,cy,cz;
  float rcoeff;
  vector3d ***data;
  
  wdth=ch->width;hght=ch->height;dpth=ch->depth;
  data=ch->raw;
  
  /*calcul de la norme max dans le champ*/
  cx=data[0][0][0].x;cy=data[0][0][0].y;cz=data[0][0][0].z;
  maxmod=sqrt(cx*cx+cy*cy+cz*cz);
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     {
      cx=data[i][j][k].x;cy=data[i][j][k].y;cz=data[i][j][k].z;
      m=sqrt(cx*cx+cy*cy+cz*cz);
      
            if (isinf(m))
                m=0;
        
         if (m>maxmod) maxmod=m;
     }
  /*calcul de maxpixel, rcoeff et icomp correspondant a maxmod*/  
  minmod=0.0;
  imx_brukermax_3d((float)maxmod,(float)minmod,im);
  rcoeff=im->rcoeff;
  
  /*mise a jour des caracteristique de l'image*/
  im->width=wdth;im->height=hght;im->depth=dpth;
  
  /*calcul de l'image*/
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
     cx=data[i][j][k].x;cy=data[i][j][k].y;cz=data[i][j][k].z;
     m=sqrt(cx*cx+cy*cy+cz*cz);
         
         if (isinf(m))
            m=-1.0*rcoeff;
         
     im->mri[i][j][k]=(int)(m/rcoeff);
    }
  
  return(1);
 }

/*******************************************************************************
**     composante_field_to_image_3d(ch,im,type)             
**                                                                  
**     composante d'un champ dans une image
**          si type =0  => composante suivant x
**          si type =1  => composante suivant y
**          si type =2  => composante suivant z
**                       
*******************************************************************************/
int composante_field_to_image_3d(field3d *ch, grphic3d *im, int type)
{
  int i,j,k,wdth,hght,dpth;
  double max,min;
  float rcoeff;
  vector3d ***data;
  
  wdth=ch->width;hght=ch->height;dpth=ch->depth;
  data=ch->raw;
  
  /*calcul la composante max et min du champ*/
  max  = - HUGE_VAL;
    min  =   HUGE_VAL;
    
    for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     {
      if (type==0)
                {
                max=MAXI(max,data[i][j][k].x);
                min=MINI(min,data[i][j][k].x);
                }
            
            if (type==1)
                {
                max=MAXI(max,data[i][j][k].y);
                min=MINI(min,data[i][j][k].y);
                }
            
            if (type==2)
                {
                max=MAXI(max,data[i][j][k].z);
                min=MINI(min,data[i][j][k].z);
                }
            }
  
    /*calcul de maxpixel, rcoeff et icomp correspondant a maxmod*/  
  imx_brukermax_3d((float)max,(float)min,im);
  rcoeff=im->rcoeff;
  
  /*mise a jour des caracteristique de l'image*/
  im->width=wdth;im->height=hght;im->depth=dpth;
  
  /*calcul de l'image*/
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
 
            if (type==0)
                im->mri[i][j][k]=(int)(data[i][j][k].x/rcoeff);

            if (type==1)
                im->mri[i][j][k]=(int)(data[i][j][k].y/rcoeff);
            
            if (type==2)
                im->mri[i][j][k]=(int)(data[i][j][k].z/rcoeff);
    }
  
  return(1);
 }

/*******************************************************************************
**     jacobien_field_to_image_3d(ch,im)                
**                                                                  
**     jacobien d'un champ dans une image                             
*******************************************************************************/
int jacobien_field_to_image_3d(field3d *ch, grphic3d *im)
{
  float rcoeff;
  int i,j,k,l,wdth,hght,dpth;
  double maxmod,m,J1,J2,J3,J4,J5,J6,J7,J8,J9,minmod,filt[5],singularite,Jtot=0;
  vector3d ***data;
  double   ***jacobien;
    
  wdth=ch->width;hght=ch->height;dpth=ch->depth;
  data=ch->raw;
  
  jacobien=alloc_dmatrix_3d(wdth,hght,dpth);
  
  //filt[1]=-0.833812;filt[2]=0.229945;filt[3]=0.0422064;
  filt[1]=0.3326;filt[2]=0.0612;filt[3]=0.015;
  //filt[1]=1;filt[2]=0;filt[3]=0;
  
  minmod=10000.0;
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     jacobien[i][j][k]=0.0;
  
  maxmod=0.0;
  for (i=3;i<wdth-4;i++)
   for (j=3;j<hght-4;j++)
    for (k=3;k<dpth-4;k++)
    {
     J1=J5=J9=1.0;J2=J3=J4=J6=J7=J8=0.0;
     for (l=1;l<=3;l++)
     {
         J1+=filt[l]*(data[i+l][j][k].x-data[i-l][j][k].x);
     J2+=filt[l]*(data[i][j+l][k].x-data[i][j-l][k].x);
     J3+=filt[l]*(data[i][j][k+l].x-data[i][j][k-l].x);
     J4+=filt[l]*(data[i+l][j][k].y-data[i-l][j][k].y);
     J5+=filt[l]*(data[i][j+l][k].y-data[i][j-l][k].y);
     J6+=filt[l]*(data[i][j][k+l].y-data[i][j][k-l].y);
     J7+=filt[l]*(data[i+l][j][k].z-data[i-l][j][k].z);
     J8+=filt[l]*(data[i][j+l][k].z-data[i][j-l][k].z);
     J9+=filt[l]*(data[i][j][k+l].z-data[i][j][k-l].z); 
          
         /*J1+=filt[l]*(data[i+l][j][k].x-data[i][j][k].x);
     J2+=filt[l]*(data[i][j+l][k].x-data[i][j][k].x);
     J3+=filt[l]*(data[i][j][k+l].x-data[i][j][k].x);
     J4+=filt[l]*(data[i+l][j][k].y-data[i][j][k].y);
     J5+=filt[l]*(data[i][j+l][k].y-data[i][j][k].y);
     J6+=filt[l]*(data[i][j][k+l].y-data[i][j][k].y);
     J7+=filt[l]*(data[i+l][j][k].z-data[i][j][k].z);
     J8+=filt[l]*(data[i][j+l][k].z-data[i][j][k].z);
     J9+=filt[l]*(data[i][j][k+l].z-data[i][j][k].z); */
           
     }
     m=J1*(J5*J9-J6*J8)-J4*(J2*J9-J3*J8)+J7*(J2*J6-J3*J5);
     jacobien[i][j][k]=m;
     Jtot=Jtot+m;
     if (m<minmod) minmod=m;
     if (m>maxmod) maxmod=m;
    }
    
  imx_brukermax_3d((float)maxmod,(float)minmod,im);
  rcoeff=im->rcoeff;
  
  im->width=wdth;im->height=hght;im->depth=dpth;

  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
     im->mri[i][j][k]=(int)(jacobien[i][j][k]/rcoeff);
    }

    singularite=0;  
    for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    if (jacobien[i][j][k]<0.0)
            singularite++;
            
   printf("Jacobien minimum: %.2f \n",minmod);
   printf("Nombre de singulairte :  %f \n",1.0*singularite/(wdth*hght*dpth));
   
   Jtot=Jtot/(wdth*hght*dpth);
   printf("Moyenne du Jacobien : %f \n",Jtot);
   
   
   
   free_dmatrix_3d(jacobien);
  
  imx_inimaxminpixel_3d_p(im);
 
  return(1);
 }


/*******************************************************************************
********************************************************************************
******************* CALCUL DU GRADIENT D'UNE IMAGE   ***************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**  int imx_gradient_3d_p(im,grad,method,t)
*/
/*! Calcul le gradient d'une image et met le resultat dans un champ de 
**  vecteurs
**      \param im:  image 
**      \param grad:    champ resultat
**      \param method:  methode de calcul 
**              \li  0: classic
**              \li  1: taylor
**              \li  2: leat square
**      \param t:   taille du filtre si necessaire 
**
**  \attention le champ est suppose alloue a la bonne taille avant l'appel
**  a cette procedure
**
**  \remark  pour plus de details sur les methodes voir le livre:
**  practical handbook on image processing for scientific applications
**  Bernd Jahne, crc press
**
*******************************************************************************/
int imx_gradient_3d_p(grphic3d *im, field3d *grad, int method, int t)
{
 int    i,j,k,l,wdth,hght,dpth;
 vector3d ***data,nulvec;
 double *filt,gx,gy,gz;
 
 wdth=im->width;hght=im->height;dpth=im->depth;
 data=grad->raw;
 nulvec.x=nulvec.y=nulvec.z=0.0;
 
 /*verification que le champ est de la bonne taille*/
 if (grad==NULL || grad->width!=wdth || grad->height!=hght || grad->depth!=dpth)
 {printf("Dans imx_gradient_3d_p, le champ n'est pas allouer ou n'a pas la bonne taille pour le calcul\n");return(1);}
 
 /*generation du filtre de convolution*/
  /*allocation du tableau filt*/
   if (t<2) t=2;
   filt=CALLOC(t+1,double); /*on considere que le filtre est antisymetrique*/
  /*remplissage du filtre en fonction du type de method et de la taille*/
   filt[0]=0.0;
   switch (method)
   {
    case 1: /*taylor expansion*/
     switch (t)
    {case 3: filt[1]=-45.0/60.0;filt[2]=9.0/60.0;filt[3]=-1.0/60.0; break;  
     case 4: filt[1]=-672.0/840.0;filt[2]=168.0/840.0;filt[3]=-32.0/840.0;filt[4]=3.0/840.0; break;
         default:filt[1]=-8.0/12.0;filt[2]=1/12.0;t=2; break;}
    break;
    case 2: /*least-square derivation (non recursive)*/
     switch (t)
    {case 3: filt[1]=-0.833812;filt[2]=0.229945;filt[3]=0.0422064; break;   
     case 4: filt[1]=-0.88464;filt[2]=0.298974;filt[3]=-0.0949175;filt[4]=0.0178608; break;
     case 5: filt[1]=-0.914685;filt[2]=0.346228;filt[3]=-0.138704;filt[4]=0.0453905;filt[5]=-0.0086445; break;
         default:filt[1]=-0.74038;filt[2]=0.12019;t=2; break;} 
    break;
    default:
      filt[1]=-1.0;t=1;
    break;

   }

 /*application du filtre de convolution a l'image: resultat dans grad*/
  /*les effets de bords sont traite en mettant le gradient a zero au bord*/
   for (i=0;i<wdth;i++) for(j=0;j<hght;j++) for (k=0;k<t;k++) {data[i][j][k]=data[i][j][dpth-k-1]=nulvec;}
   for (i=0;i<wdth;i++) for(j=0;j<t;j++) for (k=0;k<dpth;k++) {data[i][j][k]=data[i][hght-j-1][k]=nulvec;}
   for (i=0;i<t;i++) for(j=0;j<hght;j++) for (k=0;k<dpth;k++) {data[i][j][k]=data[wdth-i-1][j][k]=nulvec;}
  /*on applique le filtre de convolution*/
   for (i=t;i<wdth-t;i++)
    for (j=t;j<hght-t;j++)
     for (k=t;k<dpth-t;k++)
      {
       gx=gy=gz=0.0;
       for (l=1;l<=t;l++)
       {
        gx+=filt[l]*(double)(-im->mri[i+l][j][k]+im->mri[i-l][j][k]);
        gy+=filt[l]*(double)(-im->mri[i][j+l][k]+im->mri[i][j-l][k]);
        gz+=filt[l]*(double)(-im->mri[i][j][k+l]+im->mri[i][j][k-l]);
       }
       data[i][j][k].x=(float)gx;
       data[i][j][k].y=(float)gy;
       data[i][j][k].z=(float)gz;
      }
    
 free(filt);   
    
 return(1);
}


/*******************************************************************************
********************************************************************************
********************** CARACTERISTIQUES CHAMP **********************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**     erreur_field_3d(champ1,champ2)             
**                                                                   
**     calcule l'erreur quadratique entre deux champs                        
*******************************************************************************/
double erreur_field_3d(field3d *champ1, field3d *champ2)
{
 int i,j,k,wdth,hght,dpth;
 double som,diff;
 vector3d ***data1,***data2;
 
 wdth=champ1->width;hght=champ1->height;;dpth=champ1->depth;
 
 if (champ2->width!=wdth || champ2->height!=hght || champ2->depth!=dpth)
  {fprintf(stderr,"champ de taille differente dans erreur_field_3d \n");return(1);}
 
 data1=champ1->raw;data2=champ2->raw;
 
 som=0.0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;j++)
   {
    diff=data2[i][j][k].x-data1[i][j][k].x;som=som+diff*diff;
    diff=data2[i][j][k].y-data1[i][j][k].y;som=som+diff*diff;
    diff=data2[i][j][k].z-data1[i][j][k].z;som=som+diff*diff;
   } 
  
 som=sqrt(som)/(double)(2*wdth*hght*dpth);
 
 return(som);
}


/*******************************************************************************
**   carac_field_3d()                                  
*/
/*! fonction d'affichage dans le log du min et du max du champ suivant x, y, et z
**  \param champ : le champ
**  \retval 1                                                                       
**                
*******************************************************************************/
int carac_field_3d(field3d *champ)
{
 vector3d ***data;
 int    i,j,k,wdth,hght,dpth;
 double minx,maxx,miny,maxy,minz,maxz;
 
 
 wdth=champ->width;hght=champ->height;dpth=champ->depth;
 data=champ->raw;
 
 minx=maxx=data[0][0][0].x;miny=maxy=data[0][0][0].y;minz=maxz=data[0][0][0].z;
 
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    if (data[i][j][k].x<minx) minx=data[i][j][k].x;
    if (data[i][j][k].y<miny) miny=data[i][j][k].y;
    if (data[i][j][k].z<minz) minz=data[i][j][k].z;
    if (data[i][j][k].x>maxx) maxx=data[i][j][k].x;
    if (data[i][j][k].y>maxy) maxy=data[i][j][k].y;
    if (data[i][j][k].z>maxz) maxz=data[i][j][k].z;
   } 
 
 aff_log("minx: %.2f   maxx: %.2f   miny: %.2f   maxy: %.2f   minz: %.2f   maxz: %.2f\n",minx,maxx,miny,maxy,minz,maxz); 
 return(1); 
}


/*******************************************************************************
********************************************************************************
********************** METHODES D'INTERPOLATIONS *******************************
********************************************************************************
*******************************************************************************/
/*
ATTENTION DANS TOUTES LES METHODES D'INTERPOLATION ON SUPPOSE QUE IMRES EST
DIFFERENT DE IMDEB (GAIN DE TEMPS DE CALCUL)*/

/*! \addtogroup InterpolationFonction 
**
@{
*/
/*******************************************************************************
**    inter_nearest_3d(imdeb,champ,imres)                                   
*/                                                                       
/*!    fonction d'interpolation par la methode du plus proche voisin      
**    imres est l'image de im1 transformee par le champ                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
*******************************************************************************/
int inter_nearest_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 int wdthDeb,hghtDeb,dpthDeb;
 int wdthField,hghtField,dpthField;
 int i,j,k,x,y,z;
 vector3d   ***data;
 
 wdthDeb=imdeb->width;hghtDeb=imdeb->height;dpthDeb=imdeb->depth;
 wdthField=champ->width;hghtField=champ->height;dpthField=champ->depth;
 if (champ->width!=(int)imres->width || champ->height!=(int)imres->height || champ->depth!=(int)imres->depth)
  {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_nearest_3d\n");return(1);}
  
 data=champ->raw;
/* les parametres de imres ne sont pas forcement
ceux de imdeb (dx,dy,dz)
 imx_copie_param_3d_p(imdeb,imres);
 */
 imx_copie_image_params_3d_p(imdeb,imres);
 
 for (i=0;i<wdthField;i++)
  for (j=0;j<hghtField;j++)
   for (k=0;k<dpthField;k++)
    {
    //on prend les valeurs arrondies de x par un floor(x+0.5), pareil pour y et z
     x=(int)floor((double)i+data[i][j][k].x+0.5);
     y=(int)floor((double)j+data[i][j][k].y+0.5);
     z=(int)floor((double)k+data[i][j][k].z+0.5);
     //si le pixel le plus proche est dans l'image on prend sa valeur sinon on le met a 0 dans imres
     if (x>=0 && x<wdthDeb && y>=0 && y<hghtDeb && z>=0 && z<dpthDeb)
       {imres->mri[i][j][k]=imdeb->mri[x][y][z];}
     else
       {imres->mri[i][j][k]=0;}
    }
   
 return(0);  
}

/*******************************************************************************
**    inter_linear_3d(imdeb,champ,imres)                                  
*/
/*!                                                                       
**    fonction d'interpolation bilineaire                                
**    imres est l'image de im1 transformee par le champ                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
*******************************************************************************/
/*int inter_linear_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 int wdth,hght,dpth;
 int wdthDeb,hghtDeb,dpthDeb;
 int i,j,k,i2,j2,k2,i3,j3,k3;
 double x,y,z,xrel,yrel,zrel;
 double va,vb,vc,vd,ve,vf,vg,vh,val;
 vector3d ***data;
 TYPEMRI3D ***imresMRI=NULL, ***imdebMRI=NULL;
 
 
 wdthDeb=imdeb->width;hghtDeb=imdeb->height;dpthDeb=imdeb->depth;
 wdth=champ->width;hght=champ->height;dpth=champ->depth;

 if (champ->width!=(int)imres->width || champ->height!=(int)imres->height || champ->depth!=(int)imres->depth)
  {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_linear_3d\n");return(1);}  


 
 data=champ->raw;

 imx_copie_image_params_3d_p(imdeb,imres);
 imresMRI=imres->mri; imdebMRI=imdeb->mri;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    x=i+data[i][j][k].x;
    y=j+data[i][j][k].y;
    z=k+data[i][j][k].z;
    i2=(int)floor(x);j2=(int)floor(y);k2=(int)floor(z);
    i3=i2+1;j3=j2+1;k3=k2+1;

    if (i2>=0 && i3<wdthDeb && j2>=0 && j3<hghtDeb && k2>=0 && k3<dpthDeb)
      {// tout le voisinage est dans l'image
        xrel=x-i2;yrel=y-j2;zrel=z-k2;
    
        va=(double)imdebMRI[i2][j2][k2];
        vb=(double)imdebMRI[i2][j2][k3];
        vc=(double)imdebMRI[i2][j3][k2];
        vd=(double)imdebMRI[i2][j3][k3];
        ve=(double)imdebMRI[i3][j2][k2];
        vf=(double)imdebMRI[i3][j2][k3];
        vg=(double)imdebMRI[i3][j3][k2];
        vh=(double)imdebMRI[i3][j3][k3];

        val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        imresMRI[i][j][k]=(int)floor(val+0.5);
      }
      else
      {
       if (i3>=0 && i2<wdthDeb && j3>=0 && j2<hghtDeb && k3>=0 && k2<dpthDeb)
       {//une partie du voisinage est dans l'image
    
        xrel=x-i2;yrel=y-j2;zrel=z-k2;
    if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imdebMRI[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==dpthDeb) vb=0.0; else vb=(double)imdebMRI[i2][j2][k3];
        if (i2==-1 || j3==hghtDeb || k2==-1) vc=0.0; else vc=(double)imdebMRI[i2][j3][k2];
        if (i2==-1 || j3==hghtDeb || k3==dpthDeb) vd=0.0; else vd=(double)imdebMRI[i2][j3][k3];
        if (i3==wdthDeb || j2==-1 || k2==-1) ve=0.0; else ve=(double)imdebMRI[i3][j2][k2];
        if (i3==wdthDeb || j2==-1 || k3==dpthDeb) vf=0.0; else vf=(double)imdebMRI[i3][j2][k3];
        if (i3==wdthDeb || j3==hghtDeb || k2==-1) vg=0.0; else vg=(double)imdebMRI[i3][j3][k2];
        if (i3==wdthDeb || j3==hghtDeb || k3==dpthDeb) vh=0.0; else vh=(double)imdebMRI[i3][j3][k3];

        val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        imresMRI[i][j][k]=(int)floor(val+0.5);
       }
       else
       {//le voisinage est totalement en dehors de l'image
        imresMRI[i][j][k]=0;
       }
      }
   }
   
 return(0);  
} */

int inter_linear_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 int wdth,hght,dpth;
 int wdthDeb,hghtDeb,dpthDeb;
 int wdthDeb_m,hghtDeb_m,dpthDeb_m;
 int wdthDeb_M,hghtDeb_M,dpthDeb_M;
 int i,j,k,i2,j2,k2,i3,j3,k3;
 double x,y,z,xrel,yrel,zrel;
 double va,vb,vc,vd,ve,vf,vg,vh,val;
 vector3d ***data;
 TYPEMRI3D ***imresMRI=NULL, ***imdebMRI=NULL;


 wdthDeb=(int)imdeb->width;hghtDeb=(int)imdeb->height;dpthDeb=(int)imdeb->depth;
 wdth=(int)champ->width;hght=(int)champ->height;dpth=(int)champ->depth;

 if (champ->width!=(int)imres->width || champ->height!=(int)imres->height || champ->depth!=(int)imres->depth)
  {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_linear_3d\n");return(1);}

 data=champ->raw;

 imx_copie_image_params_3d_p(imdeb,imres);
 imresMRI=imres->mri; imdebMRI=imdeb->mri;
 wdthDeb_m=wdthDeb-1,hghtDeb_m=hghtDeb-1,dpthDeb_m=dpthDeb-1;
 wdthDeb_M=wdthDeb,hghtDeb_M=hghtDeb,dpthDeb_M=dpthDeb;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    x=i+data[i][j][k].x;
    y=j+data[i][j][k].y;
    z=k+data[i][j][k].z;
    
    i2=(int)floor(x);j2=(int)floor(y);k2=(int)floor(z);
    i3=i2+1;j3=j2+1;k3=k2+1;
    
         
      
      if ( i2>=0 && i2<wdthDeb_m && j2>=0 && j2<hghtDeb_m && k2>=0 && k2<dpthDeb_m )
      {
      // tout le voisinage est dans l'image
        i3=i2+1;j3=j2+1;k3=k2+1;
        xrel=x-i2;yrel=y-j2;zrel=z-k2;

        va=(double)imdebMRI[i2][j2][k2];
        vb=(double)imdebMRI[i2][j2][k3];
        vc=(double)imdebMRI[i2][j3][k2];
        vd=(double)imdebMRI[i2][j3][k3];
        ve=(double)imdebMRI[i3][j2][k2];
        vf=(double)imdebMRI[i3][j2][k3];
        vg=(double)imdebMRI[i3][j3][k2];
        vh=(double)imdebMRI[i3][j3][k3];

        val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        imresMRI[i][j][k]=(int)floor(val+0.5);
      }
      else if ( i2>=-1 && i3<=wdthDeb_M && j2>=-1 && j3<=hghtDeb_M && k2>=-1 && k3<=dpthDeb_M )
      {//une partie du voisinage est dans l'image
            
      
        xrel=x-i2;yrel=y-j2;zrel=z-k2;
        if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)imdebMRI[i2][j2][k2];
        if (i2==-1 || j2==-1 || k3==dpthDeb) vb=0.0; else vb=(double)imdebMRI[i2][j2][k3];
        if (i2==-1 || j3==hghtDeb || k2==-1) vc=0.0; else vc=(double)imdebMRI[i2][j3][k2];
        if (i2==-1 || j3==hghtDeb || k3==dpthDeb) vd=0.0; else vd=(double)imdebMRI[i2][j3][k3];
        if (i3==wdthDeb || j2==-1 || k2==-1) ve=0.0; else ve=(double)imdebMRI[i3][j2][k2];
        if (i3==wdthDeb || j2==-1 || k3==dpthDeb) vf=0.0; else vf=(double)imdebMRI[i3][j2][k3];
        if (i3==wdthDeb || j3==hghtDeb || k2==-1) vg=0.0; else vg=(double)imdebMRI[i3][j3][k2];
        if (i3==wdthDeb || j3==hghtDeb || k3==dpthDeb) vh=0.0; else vh=(double)imdebMRI[i3][j3][k3];

        val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
        +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
        +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
        imresMRI[i][j][k]=(int)floor(val+0.5);
       }
       else
       {//le voisinage est totalement en dehors de l'image
        imresMRI[i][j][k]=0;
       }
   }

 return(0);
}


/*******************************************************************************
**    inter_sinc_3d(imdeb,champ,imres)                                   
*/
/*!                                                                      
**    fonction d'interpolation par la methode du sinus cardinal      
**    imres est l'image de im1 transformee par le champ
**    apodisation de Hanning                
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
*******************************************************************************/
int inter_sinc_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 int wdthDeb,hghtDeb,dpthDeb;
 int wdthField,hghtField,dpthField;
 int   i,j,k,l,i1,j1,k1,i2,j2,k2,i1b,j1b,k1b,i2b,j2b,k2b,ic,jc,kc;
 double x,y,z,d,DX,DY,DZ,D,fx[10],fy[10],fz[10];
 double val,pix;
 vector3d   ***data=NULL;
 double ***tmp_res=NULL;
 TYPEMRI3D ***imdebMRI=NULL;
 const int tf=2;  //taille du filtre
 const double TX=1.0/(double)tf;
 const double TY=1.0/(double)tf;
 const double TZ=1.0/(double)tf;

 wdthDeb=imdeb->width;hghtDeb=imdeb->height;dpthDeb=imdeb->depth;
 wdthField=champ->width;hghtField=champ->height;dpthField=champ->depth;
 if (wdthField!=(int)imres->width || hghtField!=(int)imres->height || dpthField!=(int)imres->depth)
  {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_sinc_3d\n");return(1);}

 // image temporaire en double
 tmp_res=alloc_dmatrix_3d(wdthField, hghtField, dpthField);

 data=champ->raw;

 imdebMRI=imdeb->mri;

 for (i=0;i<wdthField;i++)
  for (j=0;j<hghtField;j++)
   for (k=0;k<dpthField;k++)
    {
     x=(double)i+data[i][j][k].x;
     y=(double)j+data[i][j][k].y;
     z=(double)k+data[i][j][k].z;

     //indices d'application du filtre
     i1=(int)floor(x)+1-tf;j1=(int)floor(y)+1-tf;k1=(int)floor(z)+1-tf;
     i2=(int)floor(x)+tf;j2=(int)floor(y)+tf;k2=(int)floor(z)+tf;

     if (i2>=0 && i1<wdthDeb && j2>=0 && j1<hghtDeb && k2>=0 && k1<dpthDeb) //une partie du voisinage est dans l'image: on calcule le filtre
     {
      DX=DY=DZ=0.0;
       //suivant x
       for (l=i1;l<=i2;l++)
        {d=x-(double)l;d=PI*sqrt(d*d);
         if (d!=0.0) fx[l-i1]=sin(d)/d*0.5*(1.0+cos(d*TX));
         else fx[l-i1]=1.0;
         DX=DX+fx[l-i1];
        }
       //suivant y
       for (l=j1;l<=j2;l++)
        {d=y-(double)l;d=PI*sqrt(d*d);
         if (d!=0.0) fy[l-j1]=sin(d)/d*0.5*(1.0+cos(d*TY));
         else fy[l-j1]=1.0;
         DY=DY+fy[l-j1];
        }
       //suivant z
       for (l=k1;l<=k2;l++)
        {d=z-(double)l;d=PI*sqrt(d*d);
         if (d!=0.0) fz[l-k1]=sin(d)/d*0.5*(1.0+cos(d*TZ));
         else fz[l-k1]=1.0;
         DZ=DZ+fz[l-k1];
        }
       D=DX*DY*DZ;

       //calcul du champ d'application du filtre
       i1b=MAXI(0, i1);
       j1b=MAXI(0, j1);
       k1b=MAXI(0, k1);
       i2b=MINI(wdthDeb-1, i2);
       j2b=MINI(hghtDeb-1, j2);
       k2b=MINI(dpthDeb-1, k2);

       //calcul de la valeur val du pixel
       val=0.0;
       for (ic=i1b;ic<=i2b;ic++)
        for (jc=j1b;jc<=j2b;jc++)
         for (kc=k1b;kc<=k2b;kc++)
         {
          pix=(double)(imdebMRI[ic][jc][kc]);
          val+=fx[ic-i1]*fy[jc-j1]*fz[kc-k1]/D*pix;
         }
       tmp_res[i][j][k]=val;

       
      }
      else //aucune partie du voisinage est dans l'image
      {
       tmp_res[i][j][k]=0.0;
      }
      
     }

 //normalisation
 normaliser_3d(imdeb, imres, tmp_res, wdthField, hghtField, dpthField);

 if (tmp_res) free_dmatrix_3d(tmp_res);

 return(0);
}


/*******************************************************************************
**    inter_qsinc2_3d(imdeb,champ,imres)                                   
*/                                                                       
/*!    fonction d'interpolation par la methode du sinus cardinal rapide     
**    imres est l'image de im1 transformee par le champ  
**  La taille du filtre est 2                
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
*******************************************************************************/
int inter_qsinc2_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 return inter_qsinc_3d(imdeb, champ, imres, 2);
}

/*******************************************************************************
**    inter_qsinc3_3d(imdeb,champ,imres)                                   
*/
/*!                                                                      
**    fonction d'interpolation par la methode du sinus cardinal rapide     
**    imres est l'image de im1 transformee par le champ  
**  La taille du filtre est 3                
**   \param imdeb : l'image de depart
**   \param champ : le champ
**   \param imres : l'image resultat (E/S)
**   \retval 1
*******************************************************************************/
int inter_qsinc3_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 return inter_qsinc_3d(imdeb, champ, imres, 3); 
}

/*******************************************************************************
**    inter_qsinc_3d(imdeb,champ,imres)
*/
/*!
**    fonction d'interpolation par la methode du sinus cardinal rapide
**    imres est l'image de im1 transformee par le champ
**    apodisation de Hanning 
**   \param imdeb : l'image de depart
**   \param champ : le champ
**   \param imres : l'image resultat (E/S)
**   \param tf : la taille du filtre
**   \retval 1
*******************************************************************************/
int inter_qsinc_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres, int tf)
{
 int wdth,hght,dpth;
 int wdthDeb,hghtDeb,dpthDeb;
 int   i,j,k,nbp,ic,jc,kc,ix0,iy0,iz0,ix1,iy1,iz1;
 int    ix,iy,iz,v1,v2,v3,v4,v5,v6,v7,v8;
 int   il,jl,kl,ir,jr,kr;
 double x,y,z,d,D;
 double fr[100][10],fl[100][10];
 double val,pix,pas;
 vector3d   ***data=NULL;
 double ***tmp_res=NULL;
 TYPEMRI3D ***imdebMRI=NULL;
 
 wdthDeb=imdeb->width;hghtDeb=imdeb->height;dpthDeb=imdeb->depth;
 wdth = champ->width;hght = champ->height;dpth = champ->depth;

 if (((int)imres->width != wdth) ||  ((int)imres->height != hght) || ((int)imres->depth != dpth))
  {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_qsinc_3d\n");return(1);}

 // image temporaire en double
 tmp_res=alloc_dmatrix_3d(wdth, hght, dpth);
  
 data=champ->raw;

 /*nombre de points de precision*/
 nbp=90;pas=1.0/(double)nbp;
 /*precalcul des filtres*/
 for (i=0;i<=nbp;i++)
 {
  D=0.0;
  for (j=0;j<tf;j++)
   {
    d=1.0-pas*(double)i+(double)j;
    d=PI*sqrt(d*d);
    if (d!=0.0) fr[i][j]=sin(d)/d*0.5*(1.0+cos(2.0*d/(double)(2*tf)));
      else  fr[i][j]=1.0;
    D=D+fr[i][j];
    d=pas*(double)i+(double)j;
    d=PI*sqrt(d*d);
    if (d!=0.0) fl[i][j]=sin(d)/d*0.5*(1.0+cos(2.0*d/(double)(2*tf)));
      else  fl[i][j]=1.0;
    D=D+fl[i][j];
   }
  for (j=0;j<tf;j++) {fr[i][j]=fr[i][j]/D;fl[i][j]=fl[i][j]/D;}
 }

 imdebMRI=imdeb->mri;

 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
    {
     x=(double)i+data[i][j][k].x;ix0=(int)floor(x);ix1=ix0+1;
     y=(double)j+data[i][j][k].y;iy0=(int)floor(y);iy1=iy0+1;
     z=(double)k+data[i][j][k].z;iz0=(int)floor(z);iz1=iz0+1;

     /*on regarde d'abord si le point deplace est dans l'image*/
     if (ix0>=0 && ix1<wdthDeb && iy0>=0 && iy1<hghtDeb && iz0>=0 && iz1<dpthDeb )
     {
      /*cas ou le point n'est entoure que par des zero*/
       v1=imdebMRI[ix0][iy0][iz0];v2=imdebMRI[ix1][iy0][iz0];
       v3=imdebMRI[ix0][iy1][iz0];v4=imdebMRI[ix0][iy0][iz1];
       v5=imdebMRI[ix1][iy1][iz0];v6=imdebMRI[ix1][iy0][iz1];
       v7=imdebMRI[ix0][iy1][iz1];v8=imdebMRI[ix1][iy1][iz1];
       if (v1==0 && v2==0 && v3==0 && v4==0 && v5==0 && v6==0 && v7==0 && v8==0) val=0.0;
       else
       {/*calcul des indices dans les tableaux de filtres precalcules*/
        ix=(int)floor(((x-(double)ix0)/pas+pas/2));
        iy=(int)floor(((y-(double)iy0)/pas+pas/2));
        iz=(int)floor(((z-(double)iz0)/pas+pas/2));

        /*calcul des indices du cadre d'application du filtre*/
        il=MAXI(0,ix0-tf);ir=MINI(wdthDeb,ix1+tf);
        jl=MAXI(0,iy0-tf);jr=MINI(hghtDeb,iy1+tf);
        kl=MAXI(0,iz0-tf);kr=MINI(dpthDeb,iz1+tf);
        
        val=0.0;

        for (ic=ix1;ic<ir;ic++) for (jc=iy1;jc<jr;jc++) for (kc=iz1;kc<kr;kc++)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fr[ix][ic-ix1]*fr[iy][jc-iy1]*fr[iz][kc-iz1]*pix;}
        for (ic=ix0;ic>il;ic--) for (jc=iy1;jc<jr;jc++) for (kc=iz1;kc<kr;kc++)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fl[ix][ix0-ic]*fr[iy][jc-iy1]*fr[iz][kc-iz1]*pix;}
        for (ic=ix1;ic<ir;ic++) for (jc=iy0;jc>jl;jc--) for (kc=iz1;kc<kr;kc++)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fr[ix][ic-ix1]*fl[iy][iy0-jc]*fr[iz][kc-iz1]*pix;}
        for (ic=ix1;ic<ir;ic++) for (jc=iy1;jc<jr;jc++) for (kc=iz0;kc>kl;kc--)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fr[ix][ic-ix1]*fr[iy][jc-iy1]*fl[iz][iz0-kc]*pix;}
        for (ic=ix0;ic>il;ic--) for (jc=iy0;jc>jl;jc--) for (kc=iz1;kc<kr;kc++)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fl[ix][ix0-ic]*fl[iy][iy0-jc]*fr[iz][kc-iz1]*pix;}
        for (ic=ix0;ic>il;ic--) for (jc=iy1;jc<jr;jc++) for (kc=iz0;kc>kl;kc--)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fl[ix][ix0-ic]*fr[iy][jc-iy1]*fl[iz][iz0-kc]*pix;}
        for (ic=ix1;ic<ir;ic++) for (jc=iy0;jc>jl;jc--) for (kc=iz0;kc>kl;kc--)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fr[ix][ic-ix1]*fl[iy][iy0-jc]*fl[iz][iz0-kc]*pix;}
        for (ic=ix0;ic>il;ic--) for (jc=iy0;jc>jl;jc--) for (kc=iz0;kc>kl;kc--)
        {pix=(double)(imdebMRI[ic][jc][kc]);
        val+=fl[ix][ix0-ic]*fl[iy][iy0-jc]*fl[iz][iz0-kc]*pix;}

       }
     } else val=0.0;
     tmp_res[i][j][k]=val;
    }

 //normalisation
 normaliser_3d(imdeb, imres, tmp_res, wdth, hght, dpth);

 if(tmp_res) free_dmatrix_3d(tmp_res);
 
 return(0);
}

/*******************************************************************************
**    inter_Bspline_3d(imdeb,champ,imres,degre)                                   
*/
/*!                                                                       
**    fonction d'interpolation par une Bspline de degre n.     
**    imres est l'image de im1 transformee par le champ  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
**   
*******************************************************************************/
int inter_Bspline_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres,int degre)
{
    int wdthDeb,hghtDeb,dpthDeb,nbpoles;
    int wdthField,hghtField,dpthField;
    int i,j,k,l,m,n;
    double x,y,z,w,w2,w4,t,t0,t1,vint,sumx,sumxy;
    double ***coeff,*row1D, ***tmp_res;
    double epsilon=2.2204460492503131E-16;
    double poles[2],xWeight[6],yWeight[6],zWeight[6]; 
    int xIndex[6],yIndex[6],zIndex[6];
   
  wdthDeb=imdeb->width;hghtDeb=imdeb->height;dpthDeb=imdeb->depth;
  wdthField=champ->width;hghtField=champ->height;dpthField=champ->depth;
  if (wdthField!=(int)imres->width || hghtField!=(int)imres->height || dpthField!=(int)imres->depth)
  {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_bspline_3d\n");return(1);}
    
  imx_copie_image_params_3d_p(imdeb,imres);

    
    // Calcul des coefficient d'interpolation a partir des echantillons de l'image de depart
    coeff=alloc_dmatrix_3d(wdthDeb, hghtDeb, dpthDeb);
    
    // image temporaire en double
    tmp_res=alloc_dmatrix_3d(wdthField, hghtField, dpthField);
    
    // Recuperation de la valeur des poles en z
    switch (degre) 
    {
    case 2:
        poles[0]=(sqrt(8.0) - 3.0);
        nbpoles=1;
        break;
    case 3:
        poles[0]=(sqrt(3.0) - 2.0);
        nbpoles=1;
        break;
    case 4:
        poles[0]=(sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0);
        poles[1]=(sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0);
        nbpoles=2;
        break;
    case 5:
        poles[0]=(sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)- 13.0 / 2.0);
        poles[1]=(sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)- 13.0 / 2.0);
        nbpoles=2;
        break;
    default:
        {fprintf(stderr,"inter_Bspline_3d: Degre de spline invalide, degres admis: {2,3,4,5}\n");return(2);}
    }

        
        // Calcul des coeeficients suivant x
        row1D=(double *)malloc(wdthDeb*sizeof(double));
        
        for (k=0; k<dpthDeb; k++)
        {
            for (j=0; j<hghtDeb; j++)
            {
                for (i=0; i<wdthDeb; i++){row1D[i]=imdeb->mri[i][j][k];}
                GetInterpolationCoefficients(row1D,wdthDeb,poles,nbpoles,epsilon);
                for (i=0; i<wdthDeb; i++){coeff[i][j][k]=row1D[i];}
            }
        }
        free(row1D);
    
    
        // Calcul des coefficients suivant y
        row1D=(double *)malloc(hghtDeb*sizeof(double));
        
        for (k=0; k<dpthDeb; k++)
        {
            for (i=0; i<wdthDeb; i++)
            {
                for (j=0; j<hghtDeb; j++){row1D[j]=coeff[i][j][k];}
                GetInterpolationCoefficients(row1D,hghtDeb,poles,nbpoles,epsilon);
                for (j=0; j<hghtDeb; j++){coeff[i][j][k]=row1D[j];}
            }
        }
        free(row1D);
    
        // Calcul des coefficients suivant z
        row1D=(double *)malloc(dpthDeb*sizeof(double));
        
        for (i=0; i<wdthDeb; i++)
        {
            for (j=0; j<hghtDeb; j++)
            {
                for (k=0; k<dpthDeb; k++){row1D[k]=coeff[i][j][k];}
                GetInterpolationCoefficients(row1D,dpthDeb,poles,nbpoles,epsilon);
                for (k=0; k<dpthDeb; k++){coeff[i][j][k]=row1D[k];}
            }
        }
        free(row1D);
        
        for (i=0;i<wdthField;i++)
            for (j=0;j<hghtField;j++)
                for (k=0;k<dpthField;k++)
                {
                x=i+champ->raw[i][j][k].x;
                y=j+champ->raw[i][j][k].y;
                z=k+champ->raw[i][j][k].z;
                
                if (x<0 || x>=wdthDeb || y<0 || y>=hghtDeb || z<0 || z>=dpthDeb)  
                    {
          tmp_res[i][j][k]=0.0;
                    }
                else
                    {
                    
                    //***********************************************
                    // Calcul de la valeur interpolee vint en (x,y,z)
                    //**********************************************/
                    
                    
                    xIndex[0]=(int)floor(x+0.5)-degre/2;
                    yIndex[0]=(int)floor(y+0.5)-degre/2;
                    zIndex[0]=(int)floor(z+0.5)-degre/2;
                    
                    for (l=0;l<degre;l++)
                        {
                        xIndex[l+1]=xIndex[l]+1;
                        yIndex[l+1]=yIndex[l]+1;
                        zIndex[l+1]=zIndex[l]+1;
                        }
                        
                    
                    // Calcul des poids de la fonction d'interpolation
                    
                    switch (degre) 
                    {
                    case 2:
                        /* x */
                        w = x - (double)xIndex[1];
                        xWeight[1] = 3.0 / 4.0 - w * w;
                        xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
                        xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
                        /* y */
                        w = y - (double)yIndex[1];
                        yWeight[1] = 3.0 / 4.0 - w * w;
                        yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
                        yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
                        /* z */
                        w = z - (double)zIndex[1];
                        zWeight[1] = 3.0 / 4.0 - w * w;
                        zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
                        zWeight[0] = 1.0 - zWeight[1] - zWeight[2];
                        break;
                    case 3:
                        /* x */
                        w = x - (double)xIndex[1];
                        xWeight[3] = (1.0 / 6.0) * w * w * w;
                        xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
                        xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
                        xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
                        /* y */
                        w = y - (double)yIndex[1];
                        yWeight[3] = (1.0 / 6.0) * w * w * w;
                        yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
                        yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
                        yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
                        /* z */
                        w = z - (double)zIndex[1];
                        zWeight[3] = (1.0 / 6.0) * w * w * w;
                        zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
                        zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
                        zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
                        break;
                    case 4:
                        /* x */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        t = (1.0 / 6.0) * w2;
                        xWeight[0] = 1.0 / 2.0 - w;
                        xWeight[0] *= xWeight[0];
                        xWeight[0] *= (1.0 / 24.0) * xWeight[0];
                        t0 = w * (t - 11.0 / 24.0);
                        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
                        xWeight[1] = t1 + t0;
                        xWeight[3] = t1 - t0;
                        xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
                        xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
                        /* y */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        t = (1.0 / 6.0) * w2;
                        yWeight[0] = 1.0 / 2.0 - w;
                        yWeight[0] *= yWeight[0];
                        yWeight[0] *= (1.0 / 24.0) * yWeight[0];
                        t0 = w * (t - 11.0 / 24.0);
                        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
                        yWeight[1] = t1 + t0;
                        yWeight[3] = t1 - t0;
                        yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
                        yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
                        /* z */
                        w = z - (double)zIndex[2];
                        w2 = w * w;
                        t = (1.0 / 6.0) * w2;
                        zWeight[0] = 1.0 / 2.0 - w;
                        zWeight[0] *= zWeight[0];
                        zWeight[0] *= (1.0 / 24.0) * zWeight[0];
                        t0 = w * (t - 11.0 / 24.0);
                        t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
                        zWeight[1] = t1 + t0;
                        zWeight[3] = t1 - t0;
                        zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
                        zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
                        break;
                    case 5:
                        /* x */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
                        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
                        t1 = (-1.0 / 12.0) * w * (t + 4.0);
                        xWeight[2] = t0 + t1;
                        xWeight[3] = t0 - t1;
                        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
                        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
                        xWeight[1] = t0 + t1;
                        xWeight[4] = t0 - t1;
                        /* y */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
                        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
                        t1 = (-1.0 / 12.0) * w * (t + 4.0);
                        yWeight[2] = t0 + t1;
                        yWeight[3] = t0 - t1;
                        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
                        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
                        yWeight[1] = t0 + t1;
                        yWeight[4] = t0 - t1;
                        /* z */
                        w = z - (double)zIndex[2];
                        w2 = w * w;
                        zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
                        t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
                        t1 = (-1.0 / 12.0) * w * (t + 4.0);
                        zWeight[2] = t0 + t1;
                        zWeight[3] = t0 - t1;
                        t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
                        t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
                        zWeight[1] = t0 + t1;
                        zWeight[4] = t0 - t1;
                        break;
                    default:
                    fprintf(stderr,"inter_Bspline_3d: Degre de spline invalide, degres admis: {2,3,4,5}\n");
          return(2);
                    }
                    // Condition limite (mirror boundaries)
                    for (l=0;l<=degre;l++)
                        {
                        if (xIndex[l]<0)
                            xIndex[l]=-xIndex[l];
                        if (xIndex[l]>=wdthDeb)
                            xIndex[l]=2*(wdthDeb-1)-xIndex[l];
                        if (yIndex[l]<0)
                            yIndex[l]=-yIndex[l];
                        if (yIndex[l]>=hghtDeb)
                            yIndex[l]=2*(hghtDeb-1)-yIndex[l];
                        if (zIndex[l]<0)
                            zIndex[l]=-zIndex[l];
                        if (zIndex[l]>=dpthDeb)
                            zIndex[l]=2*(dpthDeb-1)-zIndex[l];
                        }
                                    
                    vint=0;
                    for (n=0;n<=degre;n++) 
                        {
                        sumxy=0;
                        for (m=0;m<=degre;m++)
                            {
                            sumx=0;
                            for (l=0;l<=degre;l++)
                                {
                                sumx+=xWeight[l]*coeff[xIndex[l]][yIndex[m]][zIndex[n]];
                                }
                            sumxy+=yWeight[m]*sumx;
                            }
                        vint+=zWeight[n]*sumxy;
                        } 

                    tmp_res[i][j][k]=vint;
                    }
                }

 //normalisation       
 normaliser_3d(imdeb, imres, tmp_res, wdthField, hghtField,dpthField);
 
free_dmatrix_3d(coeff);
free_dmatrix_3d(tmp_res);    
return(0);
}


/*******************************************************************************
**              inter_Bspline_3d_2
*/
/*!
**      Fonction qui permet l'interface avec les fonctions de type interpol                                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval res : la valeur retournee par la fonction d'interpolation bspline(c-a-d 1)
*******************************************************************************/
int inter_Bspline_3d_2(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
    int res;
    res=inter_Bspline_3d(imdeb,champ,imres,2);
    return(res);
}

/*******************************************************************************
**              inter_Bspline_3d_3
*/
/*!
**      Fonction qui permet l'interface avec les fonctions de type interpol                                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval res : la valeur retournee par la fonction d'interpolation bspline(c-a-d 1)
*******************************************************************************/
int inter_Bspline_3d_3(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
    int res;
    res=inter_Bspline_3d(imdeb,champ,imres,3);
    return(res);
}

/*******************************************************************************
**              inter_Bspline_3d_4
*/
/*!
**      Fonction qui permet l'interface avec les fonctions de type interpol                                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval res : la valeur retournee par la fonction d'interpolation bspline(c-a-d 1)
*******************************************************************************/
int inter_Bspline_3d_4(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
    int res;
    res=inter_Bspline_3d(imdeb,champ,imres,4);
    return(res);
}

/*******************************************************************************
**              inter_Bspline_3d_5
*/
/*!
**      Fonction qui permet l'interface avec les fonctions de type interpol                                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval res : la valeur retournee par la fonction d'interpolation bspline(c-a-d 1)
*******************************************************************************/
int inter_Bspline_3d_5(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
    int res;
    res=inter_Bspline_3d(imdeb,champ,imres,5);
    return(res);
}

/*******************************************************************************
**    inter_VP_3d(imdeb,champ,imres,degre)
*/
/*!
**    fonction d'interpolation utilisee pour le volume partiel
**    avec les distances basees sur l'entropie conjointe
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
**
*******************************************************************************/
int inter_VP_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
  imx_copie_3d_p(imdeb, imres);

  return 1;
}  

/*******************************************************************************
**    int normaliser_3d(imdeb, imres, temp_res, wdth, hght, dpth)
*/
/*!
**    cree une image au format ipb a partir d'un tableau 3D de double
**  \param imdeb : image de depart, pour l'obtention des parametres
**  \param imres : image resultat (E/S)
**  \param temp_res : donnees de l'image a creer (en type double)
**  \param wdth, hght, dpth : la taille de l'image
**  \retval 1
**
*******************************************************************************/
int normaliser_3d(grphic3d *imdeb, grphic3d *imres, double ***tmp_res, int wdth, int hght, int dpth)
{
 int min, max;
 int i, j, k;
 double rcoeff;
 int iVal;
  
 max=0;min=0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    iVal = (int)floor(tmp_res[i][j][k]+0.5); 
    if (max<iVal)
     max=iVal;
    if (min>iVal)
     min=iVal;
   }

 //si on depasse la valeur max d'un type TYPEMRI3D, on remet tout sur une echelle ne depassant pas MAXMRI3D
 //sinon on se contente d'un arrondi
 //les valeurs negatives sont gardees dans l'image  
 if (max>(double)MAXMRI3D)
 {
  max=(int)floor(max*imdeb->rcoeff);
  min=(int)floor(min*imdeb->rcoeff);
  imx_brukermax_3d((float)max,(float)min,imres);
  rcoeff=imres->rcoeff;

  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     imres->mri[i][j][k]=(TYPEMRI3D)floor(imdeb->rcoeff*tmp_res[i][j][k]/rcoeff+0.5);
 }
 else
 {
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     imres->mri[i][j][k]=(TYPEMRI3D)floor(tmp_res[i][j][k]+0.5);

  imres->max_pixel=(TYPEMRI3D)floor(max+0.5);
  imres->min_pixel=(TYPEMRI3D)floor(min+0.5);
  imres->rcoeff=imdeb->rcoeff;
    imres->icomp=imdeb->icomp;
 }
 
 return 0;    
}

int normaliser_float_3d(grphic3d *imdeb, grphic3d *imres, float ***tmp_res, int wdth, int hght, int dpth)
{
 int min, max;
 int i, j, k;
 double rcoeff;
 int iVal;

 max=0;min=0;
 for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    iVal = (int)floor(tmp_res[i][j][k]+0.5);
    if (max<iVal)
     max=iVal;
    if (min>iVal)
     min=iVal;
   }

 //si on depasse la valeur max d'un type TYPEMRI3D, on remet tout sur une echelle ne depassant pas MAXMRI3D
 //sinon on se contente d'un arrondi
 //les valeurs negatives sont gardees dans l'image
 if (max>(double)MAXMRI3D)
 {
  max=(int)floor(max*imdeb->rcoeff);
  min=(int)floor(min*imdeb->rcoeff);
  imx_brukermax_3d((float)max,(float)min,imres);
  rcoeff=imres->rcoeff;

  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     imres->mri[i][j][k]=(TYPEMRI3D)floor(imdeb->rcoeff*tmp_res[i][j][k]/rcoeff+0.5);
 }
 else
 {
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     imres->mri[i][j][k]=(TYPEMRI3D)floor(tmp_res[i][j][k]+0.5);

  imres->max_pixel=(TYPEMRI3D)floor(max+0.5);
  imres->min_pixel=(TYPEMRI3D)floor(min+0.5);
  imres->rcoeff=imdeb->rcoeff;
    imres->icomp=imdeb->icomp;

 }

 return 0;
}

/*! @} */
/*******************************************************************************
**    GetInterpolationCoefficients ()
*/
/*!     \param  c :             input samples --> output coefficients 
**      \param  DataLength :    number of samples or coefficients
**      \param  z[] :           poles 
**      \param  NbPoles :       number of poles 
**      \param  Tolerance :     admissible relative error 
**                                       
**                                                                       
**      \brief Fonction calculant les coefficients pour l'interpolation Bspline
**              a partir de la valeur des echantillons
**   
*******************************************************************************/
void GetInterpolationCoefficients (double c[], int DataLength, double z[], int NbPoles, double Tolerance)
{
    double Lambda = 1.0;
    int n, k;
    /* special case required by mirror boundaries */
    if (DataLength == 1)
        return;
    /* compute the overall gain */
    for (k = 0; k < NbPoles; k = k + 1)
        Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
    /* apply the gain */
    for (n = 0; n < DataLength; n = n + 1)
        c[n] = c[n] * Lambda;
    /* loop over all poles */
    for (k = 0; k < NbPoles; k = k + 1) 
    {
        /* causal initialization */
        c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
        /* causal recursion */
        for (n = 1; n < DataLength; n = n + 1)
            c[n] = c[n] + z[k] * c[n-1];
        /* anticausal initialization */
        c[DataLength-1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
        /* anticausal recursion */
        for (n = DataLength-2; n >= 0; n = n - 1)
            c[n] = z[k] * (c[n+1] - c[n]);
    }
}

/*******************************************************************************
**   InitialCausalCoefficient  ()
*/
/*!     \param  c :             input samples --> output coefficients 
**      \param  DataLength :    number of samples or coefficients
**      \param  z[] :           poles 
**      \param  Tolerance :     admissible relative error 
**                                       
**      
**   
*******************************************************************************/
double InitialCausalCoefficient (double c[], int DataLength, double z, double Tolerance)
{
    double Sum, zn, z2n, iz;
    int n, Horizon=0;
    int TruncatedSum;
    /* this initialization corresponds to mirror boundaries */
    TruncatedSum = 0;
    if (Tolerance > 0.0) 
    {
        Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
        TruncatedSum = (Horizon < DataLength);
    }
    
    if (TruncatedSum) 
    {
        /* accelerated loop */
        zn = z;
        Sum = c[0];
        for (n = 1; n < Horizon; n = n + 1) 
        {
            Sum = Sum + zn * c[n];
            zn = zn * z;
        }
    return(Sum);
    }
    else 
    {
        /* full loop */
        zn = z;
        iz = 1.0 / z;
        z2n = pow(z, DataLength-1);
        Sum = c[0] + z2n * c[DataLength-1];
        z2n = z2n * z2n * iz;
        for (n = 1; n <= DataLength - 2; n = n + 1) 
        {
            Sum = Sum + (zn + z2n) * c[n];
            zn = zn * z;
            z2n = z2n * iz;
        }
    return(Sum / (1.0 - zn * zn));
    }
}

/*******************************************************************************
**   InitialAntiCausalCoefficient  ()
*/
/*!     \param  c :             input samples --> output coefficients 
**      \param  DataLength :    number of samples or coefficients
**      \param  z[] :           poles 
**                                       
**      
**   
*******************************************************************************/
double InitialAntiCausalCoefficient (double c[], int DataLength, double z)
{
    /* this initialization corresponds to mirror boundaries */
    return((z / (z * z - 1.0)) * (z * c[DataLength-2] + c[DataLength-1]));
}



/*******************************************************************************
**    inter_linear_chps_3d(data1,x,y,z,p)                                   
**                                                                       
**    fonction d'interpolation par la methode lineaire      
**   pour un champ de vecteur.
**   data1 est le champ.
**   x,y,z le point 
**   p le vecteur champ resultat pour les coordonnees x,y,zi                  
*******************************************************************************/
void inter_linear_chps_3d(field3d *champ1,double x,double y,double z,vector3d *p)
{
 int wdth,hght,dpth;
 int i2,j2,k2,i3,j3,k3;
 double xrel,yrel,zrel;
 double va,vb,vc,vd,ve,vf,vg,vh,val;
 vector3d ***data1;
 
 wdth=champ1->width;hght=champ1->height;dpth=champ1->depth;
 data1=champ1->raw;
 
/*  debut interpolation  pour x */
 i2=(int)x;j2=(int)y;k2=(int)z;
 i3=i2+1;j3=j2+1;k3=k2+1;
 
 if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
   {/* tout le voisinage est dans l'image*/
     xrel=x-(double)i2;yrel=y-(double)j2;zrel=z-(double)k2;

     va=(double)data1[i2][j2][k2].x;
     vb=(double)data1[i2][j2][k3].x;
     vc=(double)data1[i2][j3][k2].x;
     vd=(double)data1[i2][j3][k3].x;
     ve=(double)data1[i3][j2][k2].x;
     vf=(double)data1[i3][j2][k3].x;
     vg=(double)data1[i3][j3][k2].x;
     vh=(double)data1[i3][j3][k3].x;

     val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
     +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
     +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
     (*p).x=(float)floor(val+0.5);
   }
   else
   {
    if (i3>=0 && i2<wdth && j3>=0 && j2<hght && k3>=0 && k2<dpth)
    {/*une partie du voisinage est dans l'image*/

     xrel=x-(double)i2;yrel=y-(double)j2;zrel=z-(double)k2;
 if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)data1[i2][j2][k2].x;
     if (i2==-1 || j2==-1 || k3==dpth) vb=0.0; else vb=(double)data1[i2][j2][k3].x;
     if (i2==-1 || j3==hght || k2==-1) vc=0.0; else vc=(double)data1[i2][j3][k2].x;
     if (i2==-1 || j3==hght || k3==dpth) vd=0.0; else vd=(double)data1[i2][j3][k3].x;
     if (i3==wdth || j2==-1 || k2==-1) ve=0.0; else ve=(double)data1[i3][j2][k2].x;
     if (i3==wdth || j2==-1 || k3==dpth) vf=0.0; else vf=(double)data1[i3][j2][k3].x;
     if (i3==wdth || j3==hght || k2==-1) vg=0.0; else vg=(double)data1[i3][j3][k2].x;
     if (i3==wdth || j3==hght || k3==dpth) vh=0.0; else vh=(double)data1[i3][j3][k3].x;

     val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
     +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
     +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
     (*p).x=(float)floor(val+0.5);     
    }
    else
    {/*le voisinage est totalement en dehors de l'image*/
     (*p).x=0;
    }
   }

/*  debut interpolation  pour y */
 i2=(int)x;j2=(int)y;k2=(int)z;
 i3=i2+1;j3=j2+1;k3=k2+1;
 
 if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
   {/* tout le voisinage est dans l'image*/
     xrel=x-(double)i2;yrel=y-(double)j2;zrel=z-(double)k2;

     va=(double)data1[i2][j2][k2].y;
     vb=(double)data1[i2][j2][k3].y;
     vc=(double)data1[i2][j3][k2].y;
     vd=(double)data1[i2][j3][k3].y;
     ve=(double)data1[i3][j2][k2].y;
     vf=(double)data1[i3][j2][k3].y;
     vg=(double)data1[i3][j3][k2].y;
     vh=(double)data1[i3][j3][k3].y;

     val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
     +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel

     +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
     (*p).y=(float)floor(val+0.5);
   }
   else
   {
    if (i3>=0 && i2<wdth && j3>=0 && j2<hght && k3>=0 && k2<dpth)
    {/*une partie du voisinage est dans l'image*/

     xrel=x-(double)i2;yrel=y-(double)j2;zrel=z-(double)k2;
 if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)data1[i2][j2][k2].y;
     if (i2==-1 || j2==-1 || k3==dpth) vb=0.0; else vb=(double)data1[i2][j2][k3].y;
     if (i2==-1 || j3==hght || k2==-1) vc=0.0; else vc=(double)data1[i2][j3][k2].y;
     if (i2==-1 || j3==hght || k3==dpth) vd=0.0; else vd=(double)data1[i2][j3][k3].y;
     if (i3==wdth || j2==-1 || k2==-1) ve=0.0; else ve=(double)data1[i3][j2][k2].y;
     if (i3==wdth || j2==-1 || k3==dpth) vf=0.0; else vf=(double)data1[i3][j2][k3].y;
     if (i3==wdth || j3==hght || k2==-1) vg=0.0; else vg=(double)data1[i3][j3][k2].y;
     if (i3==wdth || j3==hght || k3==dpth) vh=0.0; else vh=(double)data1[i3][j3][k3].y;

     val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
     +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
     +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
     (*p).y=(float)floor(val+0.5);     
    }
    else
    {/*le voisinage est totalement en dehors de l'image*/
     (*p).y=0;
    }
   }

/*  debut interpolation  pour z */
 i2=(int)x;j2=(int)y;k2=(int)z;
 i3=i2+1;j3=j2+1;k3=k2+1;
 
 if (i2>=0 && i3<wdth && j2>=0 && j3<hght && k2>=0 && k3<dpth)
   {/* tout le voisinage est dans l'image*/
     xrel=x-(double)i2;yrel=y-(double)j2;zrel=z-(double)k2;

     va=(double)data1[i2][j2][k2].z;
     vb=(double)data1[i2][j2][k3].z;
     vc=(double)data1[i2][j3][k2].z;
     vd=(double)data1[i2][j3][k3].z;
     ve=(double)data1[i3][j2][k2].z;
     vf=(double)data1[i3][j2][k3].z;
     vg=(double)data1[i3][j3][k2].z;
     vh=(double)data1[i3][j3][k3].z;

     val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
     +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
     +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
     (*p).z=(float)floor(val+0.5);
   }
   else
   {
    if (i3>=0 && i2<wdth && j3>=0 && j2<hght && k3>=0 && k2<dpth)
    {/*une partie du voisinage est dans l'image*/
     xrel=x-(double)i2;yrel=y-(double)j2;zrel=z-(double)k2;
 if (i2==-1 || j2==-1 || k2==-1) va=0.0; else va=(double)data1[i2][j2][k2].z;
     if (i2==-1 || j2==-1 || k3==dpth) vb=0.0; else vb=(double)data1[i2][j2][k3].z;
     if (i2==-1 || j3==hght || k2==-1) vc=0.0; else vc=(double)data1[i2][j3][k2].z;
     if (i2==-1 || j3==hght || k3==dpth) vd=0.0; else vd=(double)data1[i2][j3][k3].z;
     if (i3==wdth || j2==-1 || k2==-1) ve=0.0; else ve=(double)data1[i3][j2][k2].z;
     if (i3==wdth || j2==-1 || k3==dpth) vf=0.0; else vf=(double)data1[i3][j2][k3].z;
     if (i3==wdth || j3==hght || k2==-1) vg=0.0; else vg=(double)data1[i3][j3][k2].z;
     if (i3==wdth || j3==hght || k3==dpth) vh=0.0; else vh=(double)data1[i3][j3][k3].z;

     val=(ve-va)*xrel+(vc-va)*yrel+(vb-va)*zrel+(vg-ve-vc+va)*xrel*yrel
     +(vd-vc-vb+va)*yrel*zrel+(vf-ve-vb+va)*xrel*zrel
     +(vh+ve+vc+vb-vg-vf-vd-va)*xrel*yrel*zrel+va;
     (*p).z=(float)floor(val+0.5);     
    }
    else
    {/*le voisinage est totalement en dehors de l'image*/
     (*p).z=0;
    }
   }
}

/*******************************************************************************
**    inter_labeled_3d(imdeb,champ,imres)                                   
*/                                                                       
/*!    fonction d'interpolation pour les images labelisees      
**    imres est l'image de im1 transformee par le champ                  
**  \param imdeb : image de depart
**  \param imres : image resultat (E/S)
**  \param champ : le champ
**  \retval 1
*******************************************************************************/
int inter_labeled_3d(grphic3d *imdeb, field3d *champ, grphic3d *imres)
{
 int wdthDeb,hghtDeb,dpthDeb,nb_label;
 int wdthField,hghtField,dpthField;
 int i,j,k,i1,j1,k1,l,auxi,auxj,auxk,maxlabel=0;
 vector3d   ***data;
 double *label,maxval,tmp;
 double x,y,z,dx,dy,dz;
 
 data=champ->raw;

 wdthDeb=imdeb->width;hghtDeb=imdeb->height;dpthDeb=imdeb->depth;
 wdthField=champ->width;hghtField=champ->height;dpthField=champ->depth;
 if (champ->width!=(int)imres->width || champ->height!=(int)imres->height || champ->depth!=(int)imres->depth)
 {fprintf(stderr,"le champ n'est pas de la bonne taille dans inter_bspline_3d\n");return(1);}
 imx_copie_param_3d_p(imdeb,imres);
 // imx_copie... ecrase les champs de taille de l'image
 imres->width = champ->width;
 imres->height= champ->height;
 imres->depth = champ->depth; 
 
 nb_label=imdeb->max_pixel+1;
 label=(double *)malloc(nb_label*sizeof(double));
 if (label==NULL) {fprintf(stderr,"Pb allocation memoire dans inter_labeled_3d \n");return(1);}
 
 for (i=0;i<wdthField;i++)
  for (j=0;j<hghtField;j++)
   for (k=0;k<dpthField;k++)
    {
     x=(double)i+data[i][j][k].x;
     y=(double)j+data[i][j][k].y;
     z=(double)k+data[i][j][k].z;
     i1=(int)floor(x);j1=(int)floor(y);k1=(int)floor(z);
     dx=x-i1;
     dy=y-j1;
     dz=z-k1;
     if (x>=0 && x<(wdthDeb-1) && y>=0 && y<(hghtDeb-1) && z>=0 && z<(dpthDeb-1))
     {
          for (l=0;l<nb_label;l++) label[l]=0.0;
          for (auxi=0;auxi<2;auxi++)
          for (auxj=0;auxj<2;auxj++)
          for (auxk=0;auxk<2;auxk++)
          {
                    tmp=((fabs(1.0-dx-(double)auxi))*(fabs(1.0-dy-(double)auxj))*(fabs(1.0-dz-(double)auxk)));
                    label[imdeb->mri[i1+auxi][j1+auxj][k1+auxk]]+=tmp;
          }
          maxval=0;
          for (l=0;l<nb_label;l++)
                    if(label[l]>maxval) {maxval=label[l];maxlabel=l;}
          imres->mri[i][j][k]=maxlabel;
     }
     else
     {imres->mri[i][j][k]=0;}
 }  
 
 free(label);  
 return(0);  
}
/*******************************************************************************
********************************************************************************
************** ENREGISTREMENT ET LECTURE SUR FICHIER ***************************
********************************************************************************
*******************************************************************************/


/*******************************************************************************
**   test_fichier_transf_3d(nomfichier)                                  
**                                                                        
**      Test si le fichier dont le nom est passe en parametre, contient bien une
**      transformation 3D
**  Ce test est fait en verifiant que la premiere ligne contient bien le
**      mot clef TRANSF3D         
**  renvoie 1 si c'est bon 0 sinon
*******************************************************************************/
int test_fichier_transf_3d(char *nomfichier)
{
 FILE *fichier;
 char  ligne[PATH_LEN];
 unsigned int retVal=0;
 
 fichier=fopen(nomfichier,"r");
 if (fichier)
 {
    fgets(ligne,PATH_LEN,fichier);
    fclose(fichier);
 
    if ((strstr(ligne,"TRANSF3D"))!=NULL) retVal=1;
 }

 return retVal; 
}   

/*******************************************************************************
**   quest_save_result_3d(save_type)                                  
*/                                                                        
/*!     Demande si on veut enregistrer la transforamtion dans un fichier
**  Si oui, demande sous quel format: champ ou parametres
**      et renvoi le nom du fichier (sans extension)
**  Si non, retourne NULL
**  La fonction gere l'existence d'un des deux fichiers        
**  \param save_type : 
**      - 0 sous format champ 
**      - 1 sous format parametre 
*******************************************************************************/
char *quest_save_result_3d(int *save_type)
{ 
 char   s[250],str[256],nomchp[250],nomraw[250],nomtrf[250],*nomfichres,*quest[4];
 FILE   *filechp,*fileraw,*filetrf;
 int    ans,err=0,s_type,i;
 
 s_type=1;
 sprintf(s,"Save resulting field?\n");
 ans=GET_EXIST(s,&err);
 if (ans==1)
 {/*on veut enregistrer le resultats*/
  
  /*on demande sous quel format: champ ou parametre*/
   for(i=0;i<4;i++)
    quest[i]=CALLOC(80,char);
   strcpy(quest[0],"champ");
   strcpy(quest[1],"parametres");
   strcpy(quest[2],"\0");
   s_type=GETV_QCM("Format d'enregistrement",(char **)quest);
   for(i=0;i<4;i++)
     free(quest[i]);
     
  
   nomfichres=CALLOC(250,char);
   sprintf(str,"File ? in %s",_CurrentPath);
  
   switch (s_type)
   {
    case 0: /*sous format champ*/ 
     do

     {
      strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
      sprintf(nomchp,"%s.chp",nomfichres);
      sprintf(nomraw,"%s.rw",nomfichres);
      filechp=fopen(nomchp,"r");
      fileraw=fopen(nomraw,"r");
      if (filechp!=NULL || fileraw!=NULL)
       { /*one of the two files already exist*/
        sprintf(s,"Delete old file ?\n");
        ans=GET_EXIST(s,&err);
        if (ans==1)
         {
          if (filechp!=NULL) {fclose(filechp);remove(nomchp);filechp=NULL;}
          if (fileraw!=NULL) {fclose(fileraw);remove(nomraw);fileraw=NULL;}
         } 
       }
      } while (filechp!=NULL || fileraw!=NULL);
    break;
    case 1: /*sous format parametre*/ 
     do
     {    
      strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
      put_file_extension(nomfichres,".trf",nomfichres);

      filetrf=fopen(nomtrf,"r");
      if (filetrf!=NULL)
       { /*the file already exist*/
        sprintf(s,"Delete old file ?\n");
        ans=GET_EXIST(s,&err);
        if (ans==1)
         {
          if (filetrf!=NULL) {fclose(filetrf);remove(nomtrf);filetrf=NULL;}
         } 
       }
      } while (filetrf!=NULL);
    break;
   }
 } else nomfichres=NULL;

 *save_type=s_type;
 return(nomfichres);
}

/*******************************************************************************
**   transf_to_field_3d()                                  
*/
/*!                                                                       
**   alloue et retourne un champ de deformation correspondant a la transformation
**   passee en parametre             
**   \param transfo : la transformation 
**   \retval field3d : le champs de deformation, NULL si echec
*******************************************************************************/
field3d *transf_to_field_3d(transf3d *transfo,grphic3d* imref,grphic3d *imreca)
{
  field3d *champ;
  int wdth, hght, dpth;
  int err = 0;

  if (imref)
  {
    transfo->width = imref->width;
    transfo->height = imref->height;
    transfo->depth = imref->depth;
  }

  wdth = transfo->width;
  hght = transfo->height;
  dpth = transfo->depth;

  /*allocation du champ*/
  champ = cr_field3d(wdth,hght,dpth);
  champ->dx=transfo->dx;
  champ->dy=transfo->dy;
  champ->dz=transfo->dz;
  
  /*remplissage du champ*/
  err = transf_to_field_3d_noalloc(transfo, champ,imref,imreca);
  if (err)
  {
    free_field3d(champ);
    champ = NULL;
  }

  return(champ);
}


/*******************************************************************************
**   transf_to_field_3d_noalloc()                                  
*/
/*!                                                                       
**   similaire a transf_to_field_3d(), sauf que le champ n'est pas alloue
**   
**   \param transfo : la transformation 
**   \param champ : le field3d deja alloue (E/S)
**   \retval 0
*******************************************************************************/
int transf_to_field_3d_noalloc(transf3d *transfo, field3d *champ,grphic3d* imref, grphic3d* imreca)
{
  field3d *champ_ini;
  transf3d *transfo_ini;
  scal_func_t scal_func;
  TDimension wdth, hght, dpth;
  int i, j, k, l, nb_param;
  double *p;
  vector3d *ch, ***data;
    int res;

 
  wdth=transfo->width;hght=transfo->height;dpth=transfo->depth;
  ch = transfo->ch;

  /*verification du champ*/
  if ((champ->width != (int)wdth) || (champ->height != (int)hght) || (champ->depth != (int)dpth))
    return 1;
  /*remplissage du champ*/
  switch (transfo->typetrans)
  {
//  case RIGID3D:
//    rigid_to_field_3d(transfo->nb_param,transfo->param,champ,imref,imreca);
//    break;
//  case RIGIDZOOM3D:
//    rigidz_to_field_3d(transfo->nb_param,transfo->param,champ,imref,imreca);
//    break;
//  case AFFINE3D:
//    affine_to_field_3d(transfo->nb_param,transfo->param,champ,imref,imreca);
//    break;
//  case RIGIDGLOBALZOOM3D:
//    rigid_global_zoom_to_field_3d(transfo->nb_param,transfo->param,champ,imref,imreca);
//    break;
//  case AFFINEDECOUPLE3D:
//    affine_decouple_to_field_3d(transfo->nb_param,transfo->param,champ,imref,imreca);
//    break;
  case RIGID3D:
  case RIGIDZOOM3D:
  case AFFINE3D:
  case RIGIDGLOBALZOOM3D:
  case AFFINEDECOUPLE3D:
    transf_geom_to_field3d(transfo, champ, imref, imreca, NULL);
    break;    
  case BSPLINE3D:
    switch (transfo->degre)
    {
    case 0: scal_func=Bsplined0;break;
    case 1: scal_func=Bsplined1;break;
    case 2: scal_func=Bsplined2;break;
    default: scal_func=Bsplined0;break;
    }
    p = CALLOC(transfo->nb_param,double);  
    nb_param = init_base_3d(wdth, hght, dpth, scal_func);
    for (i = BASE3D.resol; i<transfo->resol; i++)
      nb_param = base_resol_up_3d(p,nb_param);
    if (nb_param != transfo->nb_param)
      printf("ERREUR dans transf_to_field_3d\n");
        res=base_to_field_3d(nb_param, transfo->param, champ,imref,imreca);
    if(res==-1)
            return(-1);
    end_base_3d();
    free(p);
    break;
  case CHAMPVOXEL3D:
    data = champ->raw;
    l=0;
    for (i=0; i<(int)wdth; i++)
     for (j=0; j<(int)hght; j++)
      for (k=0; k<(int)dpth; k++)
      {
        data[i][j][k] = ch[l];
        l++;
      }
    break;
  case CHAMP3D:
    champ_to_field_3d(transfo,champ,imref,imreca);
    break;
  }

  /* on ajoute la transfo initiale si elle existe */
  if (transfo->trans_ini && (strlen(transfo->trans_ini) > 1))
  {
    (transfo->trans_ini)[strlen(transfo->trans_ini)-1] = 0; /* bidouille ignoble */
    transfo_ini = load_transf_3d(transfo->trans_ini);
    if ((transfo_ini->width != wdth) || (transfo_ini->height != hght)
          || (transfo_ini->depth != dpth))
    {
      resize_transf_3d(transfo_ini, wdth, hght, dpth);
    }
    champ_ini = transf_to_field_3d(transfo_ini,imref,imreca);
    add_field_3d(champ, champ_ini, champ);
    free_field3d(champ_ini); 
    free_transf3d(transfo_ini);
  }

  return(0);
}

int transf_to_field_champ_3d_noalloc(transf3d *transfo, field3d* points_choisis, field3d *champ,grphic3d *imref, grphic3d *imreca)
{
  TDimension wdth, hght, dpth;

  wdth=transfo->width;hght=transfo->height;dpth=transfo->depth;

  //verification du champ
  if ((champ->width != (int)wdth) || (champ->height != (int)hght) || (champ->depth != (int)dpth))
    return 1;
  //remplissage du champ
  switch (transfo->typetrans)
  {
  case RIGID3D:
    rigid_to_field_champ_3d(transfo->nb_param,transfo->param,points_choisis,champ,imref,imreca);
    break;
  case RIGIDZOOM3D:
    rigidz_to_field_champ_3d(transfo->nb_param,transfo->param,points_choisis,champ,imref,imreca);
    break;
  case AFFINE3D:
    affine_to_field_champ_3d(transfo->nb_param,transfo->param,points_choisis,champ,imref,imreca);
    break;
  case RIGIDGLOBALZOOM3D:
    rigid_global_zoom_to_field_champ_3d(transfo->nb_param,transfo->param,points_choisis,champ,imref,imreca);
    break;
  case AFFINEDECOUPLE3D:
    affine_decouple_to_field_champ_3d(transfo->nb_param,transfo->param,points_choisis,champ,imref,imreca);
    break;
  default :
    fprintf (stderr, "CASE NON PRIS EN CHARGE DANS transf_to_field_champ_3d_noalloc\n");
    return 1;
 }     

 return 0;
}

/*******************************************************************************
**   field_to_transf_3d()                                  
*/
/*!                                                                       
**   alloue et retourne une transfo champ millimetrique correspondant au champ voxelique     
**  \param champ : le champ
**  \param imref : useless (on ne sait jamais)
**  \param imreca: image sur laquelle s'applique le champ
**  \retval transf3d : la transformation 3d correspondante
*******************************************************************************/
//transf3d *field_to_transf_3d(field3d *champ,grphic3d *imref, grphic3d *imreca)
//{
// unsigned int i,j,k,l,wdth,hght,dpth;
// transf3d *transfo;
// vector3d *ch,***data;
// float dxreca,dyreca,dzreca;
// float dxref,dyref,dzref;
//
// wdth=champ->width;hght=champ->height;dpth=champ->depth;
// data=champ->raw;
//
// if (imref && imreca)
// {
//   if (imref->width != wdth || imref->height != hght || imref->depth != dpth )
//   {   fprintf(stderr,"pas la bonne taille\n"); return NULL;}
//   dxreca = imreca->dx; dyreca = imreca->dy; dzreca = imreca->dz;
//   dxref = imref->dx; dyref = imref->dy; dzref = imref->dz;
// }
// else
// {
//   dxreca=1.0 ; dyreca=1.0 ; dzreca=1.0 ;
//   dxref=1.0 ; dyref=1.0 ; dzref=1.0 ;
// }
//
//
// transfo=cr_transf3d(wdth,hght,dpth,NULL);
// transfo->typetrans=CHAMP3D;
//
// ch=CALLOC(wdth*hght*dpth,vector3d);
// transfo->ch=ch;
//
//
//  l=0;
//  for (i=0;i<wdth;i++)
//   for (j=0;j<hght;j++)
//    for (k=0;k<dpth;k++)
//    {
//      ch[l].x=(i*dxref+data[i][j][k].x)/dxreca-i;
//      ch[l].y=(j*dyref+data[i][j][k].y)/dyreca-j;
//      ch[l].z=(k*dzref+data[i][j][k].z)/dzreca-k;
//      l++;
//    }
//  if (imref)
//  {
//    transfo->dx = imref->dx; transfo->dy = imref->dy; transfo->dz = imref->dz;
//  }
//  else
//  {
//    transfo->dx = 1.0; transfo->dy = 1.0; transfo->dz = 1.0;
//  }
//
//
// return(transfo);
//}

transf3d *field_to_transf_3d(field3d *champ,grphic3d *imref, grphic3d *imreca)
{
 unsigned int   i,j,k,l,wdth,hght,dpth;
 transf3d *transfo;
 vector3d *ch,***data;
 float dxreca,dyreca,dzreca;
 float dxref,dyref,dzref;

 wdth=champ->width;hght=champ->height;dpth=champ->depth;
 data=champ->raw;

 if (imref && imreca)
 {
     if (imref->width != wdth || imref->height != hght || imref->depth != dpth )
     {   fprintf(stderr,"pas la bonne taille\n"); return NULL;}
   dxreca = imreca->dx; dyreca = imreca->dy; dzreca = imreca->dz;
   dxref = imref->dx; dyref = imref->dy; dzref = imref->dz;
 }
 else
 {
   dxreca=1.0 ; dyreca=1.0 ; dzreca=1.0 ;
   dxref=1.0 ; dyref=1.0 ; dzref=1.0 ;
 }


 transfo=cr_transf3d(wdth,hght,dpth,NULL);
 transfo->typetrans=CHAMP3D;

 ch=CALLOC(wdth*hght*dpth,vector3d);
 transfo->ch=ch;

  l=0;
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
      ch[l].x=(i+data[i][j][k].x)*dxreca-i*dxref;
      ch[l].y=(j+data[i][j][k].y)*dyreca-j*dyref;
      ch[l].z=(k+data[i][j][k].z)*dzreca-k*dzref;
      l++;
    }
  if (imref)
  {
    transfo->dx = imref->dx; transfo->dy = imref->dy; transfo->dz = imref->dz;
  }
  else
  {
    transfo->dx = 1.0; transfo->dy = 1.0; transfo->dz = 1.0;
  }


 return(transfo);
}


void field_to_transf_3d_noalloc(field3d *champ,transf3d* transfo,grphic3d *imref, grphic3d *imreca)
{
 unsigned int   i,j,k,l,wdth,hght,dpth;
 vector3d *ch,***data;
 float dxreca,dyreca,dzreca;
 float dxref,dyref,dzref;

 wdth=champ->width;hght=champ->height;dpth=champ->depth;
 data=champ->raw;

 if (imref && imreca)
 {
     if (imref->width != wdth || imref->height != hght || imref->depth != dpth )
     {   fprintf(stderr,"pas la bonne taille\n"); 
     
    #ifdef COMPILE_FOR_MEDIPY
    return;
    #else
    return NULL;
    #endif
 
     }
   dxreca = imreca->dx; dyreca = imreca->dy; dzreca = imreca->dz;
   dxref = imref->dx; dyref = imref->dy; dzref = imref->dz;
 }
 else
 {
   dxreca=1.0 ; dyreca=1.0 ; dzreca=1.0 ;
   dxref=1.0 ; dyref=1.0 ; dzref=1.0 ;
 }


 transfo->typetrans=CHAMP3D;

 ch=CALLOC(wdth*hght*dpth,vector3d);
 transfo->ch=ch;

  l=0;
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
      ch[l].x=(i+data[i][j][k].x)*dxreca-i*dxref;
      ch[l].y=(j+data[i][j][k].y)*dyreca-j*dyref;
      ch[l].z=(k+data[i][j][k].z)*dzreca-k*dzref;
      l++;
    }
  if (imref)
  {
    transfo->dx = imref->dx; transfo->dy = imref->dy; transfo->dz = imref->dz;
  }
  else
  {
    transfo->dx = 1.0; transfo->dy = 1.0; transfo->dz = 1.0;
  }


}


/*******************************************************************************
**   save_transf_3d(transfo,nomfichier)                                  
*/                                                                       
/*!     Enregistre la transformation dans un fichier texte
**  Le nom de fichier donne est considere sans extension
**  l'extension .trf est ajoutee au nom 
**  Si le fichier existe, on l'ecrase           
**  \param  transfo  : la transformation a sauver
**  \param nomfichier : le nom du fichier a sauvegarder

**  \retval 1
*******************************************************************************/
int save_transf_3d(transf3d *transfo, char *nomfichier)
{
 char     nom[256],nomsext[256],*ext;
 FILE     *fichier;
 double   *p,reel,tmp;
 int      entier,i;
 vector3d *ch;
 HEADER header;
  
 /*si l'extension .trf existe on la supprime*/
  strcpy(nom,nomfichier);
  strcpy(nomsext,nomfichier);
  ext=strstr(nomsext,".trf");
  if (ext!=NULL) *ext='\0';
 /*rajout de l'extension .trf*/
  sprintf(nom,"%s.trf",nomsext);
  aff_log("enregistrement de la transformation dans %s\n",nom);
 /*si un fichier existe on le supprime*/
  fichier=fopen(nom,"r");
  if (fichier!=NULL) {fclose(fichier); remove(nom);}
 
 /*ecriture du fichier*/
  /*carateristiques du fichier*/   
   header.exists=0;
   header.maxlines = 2000;
   header.numlines = 0;
   header.lines = (char **)malloc(2000 * sizeof(char *));;
   put_header(&header,"TRANSF3D",STRING,"",(int)NULL); 
   if ((transfo->trans_ini)!=NULL)
    put_header(&header,"trans_ini=",STRING,transfo->trans_ini,(int)NULL);

   
   
  /*type de transformation*/
   switch (transfo->typetrans) 
   {
    case RIGID3D:
     p=transfo->param;
     put_header(&header,"type=RIGID3D",STRING,"",(int)NULL); 
     entier = transfo->nb_param;put_header(&header, "nb_param=", ENTIER, &entier,(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     reel = p[0];   put_header(&header, "Theta_x=", REEL, &reel, (int)NULL);
     reel = p[1];   put_header(&header, "Theta_y=", REEL, &reel, (int)NULL);
     reel = p[2];   put_header(&header, "Theta_z=", REEL, &reel, (int)NULL);
     reel = p[3];   put_header(&header, "Tx=", REEL, &reel, (int)NULL);
     reel = p[4];   put_header(&header, "Ty=", REEL, &reel, (int)NULL);
     reel = p[5];   put_header(&header, "Tz=", REEL, &reel, (int)NULL);
     reel = p[6];   put_header(&header, "Xg=", REEL, &reel, (int)NULL);
     reel = p[7];   put_header(&header, "Yg=", REEL, &reel, (int)NULL);
     reel = p[8];   put_header(&header, "Zg=", REEL, &reel, (int)NULL);
     break;
    case RIGIDZOOM3D:
     p=transfo->param;
     put_header(&header,"type=RIGIDZOOM3D",STRING,"",(int)NULL); 
     entier = transfo->nb_param;put_header(&header, "nb_param=", ENTIER, &entier,(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     reel = p[0];   put_header(&header, "Theta_x=", REEL, &reel, (int)NULL);
     reel = p[1];   put_header(&header, "Theta_y=", REEL, &reel, (int)NULL);
     reel = p[2];   put_header(&header, "Theta_z=", REEL, &reel, (int)NULL);
     reel = p[3];   put_header(&header, "Sx=", REEL, &reel, (int)NULL);
     reel = p[4];   put_header(&header, "Sy=", REEL, &reel, (int)NULL);
     reel = p[5];   put_header(&header, "Sz=", REEL, &reel, (int)NULL);
     reel = p[6];   put_header(&header, "Tx=", REEL, &reel, (int)NULL);
     reel = p[7];   put_header(&header, "Ty=", REEL, &reel, (int)NULL);
     reel = p[8];   put_header(&header, "Tz=", REEL, &reel, (int)NULL);
     reel = p[9];   put_header(&header, "Xg=", REEL, &reel, (int)NULL);
     reel = p[10];  put_header(&header, "Yg=", REEL, &reel, (int)NULL);
     reel = p[11];  put_header(&header, "Zg=", REEL, &reel, (int)NULL);
     break;
    case AFFINE3D:
     p=transfo->param;
     put_header(&header,"type=AFFINE3D",STRING,"",(int)NULL); 
     entier = transfo->nb_param;put_header(&header, "nb_param=", ENTIER, &entier,(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     reel = p[0];   put_header(&header, "a00=", REEL, &reel, (int)NULL);
     reel = p[1];   put_header(&header, "a01=", REEL, &reel, (int)NULL);
     reel = p[2];   put_header(&header, "a02=", REEL, &reel, (int)NULL);
     reel = p[3];   put_header(&header, "a10=", REEL, &reel, (int)NULL);
     reel = p[4];   put_header(&header, "a11=", REEL, &reel, (int)NULL);
     reel = p[5];   put_header(&header, "a12=", REEL, &reel, (int)NULL);
     reel = p[6];   put_header(&header, "a20=", REEL, &reel, (int)NULL);
     reel = p[7];   put_header(&header, "a21=", REEL, &reel, (int)NULL);
     reel = p[8];   put_header(&header, "a22=", REEL, &reel, (int)NULL);
     reel = p[9];   put_header(&header, "Tx=", REEL, &reel, (int)NULL);
     reel = p[10];  put_header(&header, "Ty=", REEL, &reel, (int)NULL);
     reel = p[11];  put_header(&header, "Tz=", REEL, &reel, (int)NULL);
     reel = p[12];  put_header(&header, "Xg=", REEL, &reel, (int)NULL);
     reel = p[13];  put_header(&header, "Yg=", REEL, &reel, (int)NULL);
     reel = p[14];  put_header(&header, "Zg=", REEL, &reel, (int)NULL);
     break;
    case RIGIDGLOBALZOOM3D:
     p=transfo->param;
     put_header(&header,"type=RIGIDGLOBALZOOM3D",STRING,"",(int)NULL);
     entier = transfo->nb_param;put_header(&header, "nb_param=", ENTIER, &entier,(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     reel = p[0];   put_header(&header, "Theta_x=", REEL, &reel, (int)NULL);
     reel = p[1];   put_header(&header, "Theta_y=", REEL, &reel, (int)NULL);
     reel = p[2];   put_header(&header, "Theta_z=", REEL, &reel, (int)NULL);
     reel = p[3];   put_header(&header, "S=", REEL, &reel, (int)NULL);
     reel = p[4];   put_header(&header, "Tx=", REEL, &reel, (int)NULL);
     reel = p[5];   put_header(&header, "Ty=", REEL, &reel, (int)NULL);
     reel = p[6];   put_header(&header, "Tz=", REEL, &reel, (int)NULL);
     reel = p[7];   put_header(&header, "Xg=", REEL, &reel, (int)NULL);
     reel = p[8];   put_header(&header, "Yg=", REEL, &reel, (int)NULL);
     reel = p[9];   put_header(&header, "Zg=", REEL, &reel, (int)NULL);
     break;
    case AFFINEDECOUPLE3D:
     p=transfo->param;
     put_header(&header,"type=AFFINEDECOUPLE3D",STRING,"",(int)NULL);
     entier = transfo->nb_param;put_header(&header, "nb_param=", ENTIER, &entier,(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     reel = p[0];   put_header(&header, "Theta_x=", REEL, &reel, (int)NULL);
     reel = p[1];   put_header(&header, "Theta_y=", REEL, &reel, (int)NULL);
     reel = p[2];   put_header(&header, "Theta_z=", REEL, &reel, (int)NULL);
     reel = p[3];   put_header(&header, "Tx=", REEL, &reel, (int)NULL);
     reel = p[4];   put_header(&header, "Ty=", REEL, &reel, (int)NULL);
     reel = p[5];   put_header(&header, "Tz=", REEL, &reel, (int)NULL);
     reel = p[6];   put_header(&header, "Sx=", REEL, &reel, (int)NULL);
     reel = p[7];   put_header(&header, "Sy=", REEL, &reel, (int)NULL);
     reel = p[8];   put_header(&header, "Sz=", REEL, &reel, (int)NULL);
     reel = p[9];   put_header(&header, "Skewyx=", REEL, &reel, (int)NULL);
     reel = p[10];  put_header(&header, "Skewzx=", REEL, &reel, (int)NULL);
     reel = p[11];  put_header(&header, "Skewzy=", REEL, &reel, (int)NULL);
     reel = p[12];  put_header(&header, "Xg=", REEL, &reel, (int)NULL);
     reel = p[13];  put_header(&header, "Yg=", REEL, &reel, (int)NULL);
     reel = p[14];  put_header(&header, "Zg=", REEL, &reel, (int)NULL);
     break;
    case BSPLINE3D:
     p=transfo->param;
     put_header(&header,"type=BSPLINE3D",STRING,"",(int)NULL); 
     entier = transfo->nb_param;put_header(&header, "nb_param=", ENTIER, &entier,(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     entier = transfo->resol;put_header(&header, "resol=", ENTIER, &entier,(int)NULL);
     entier = transfo->degre;put_header(&header, "degre=", ENTIER, &entier,(int)NULL);
     /*enregistrement du fichier binaire contenant les parametres*/
     {
      char nomprm[250],*nomtmp;
      
      sprintf(nomprm,"%s.prm",nomsext);
      fichier=fopen(nomprm,"wb");
      //fwrite(p,sizeof(double),transfo->nb_param,fichier);
      imx_fwrite(p,transfo->nb_param,sizeof(double),fichier);
      fclose(fichier);
      
      /*suppresion du chemin si precise*/
      nomtmp=strrchr(nomprm,'/');
      if (nomtmp!=NULL) strcpy(nomprm,nomtmp+1);
      put_header(&header,"fichier_param=",STRING,nomprm,(int)NULL);
     }
     break;
    case CHAMP3D:
     ch=transfo->ch;
     put_header(&header,"type=CHAMP3D",STRING,"",(int)NULL); 
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     reel = transfo->dx;put_header(&header, "dx=", REEL, &reel,(int)NULL);
     reel = transfo->dy;put_header(&header, "dy=", REEL, &reel,(int)NULL);
     reel = transfo->dz;put_header(&header, "dz=", REEL, &reel,(int)NULL);
     /*enregistrement du fichier binaire contenant le champ*/
     {
      char nomchp[250],*nomtmp;
      
      sprintf(nomchp,"%s.chp",nomsext);
      fichier=fopen(nomchp,"w");
      //fwrite(ch,sizeof(vector3d),transfo->width*transfo->height*transfo->depth,fichier);
      
      /* pour assurer la compatibilite entre HP et linux */
      if(_is_bigEndian==0)
      {
       for(i=0;i<(int)(transfo->width*transfo->height*transfo->depth);i++)
       {
         tmp=ch[i].x;
         ch[i].x=ch[i].z;
         ch[i].z=(float)tmp;
       }
      }
      imx_fwrite(ch,transfo->width*transfo->height*transfo->depth,sizeof(vector3d),fichier);
      fclose(fichier);
      
      /*suppresion du chemin si precise*/
      nomtmp=strrchr(nomchp,'/');
      if (nomtmp!=NULL) strcpy(nomchp,nomtmp+1);
      put_header(&header,"fichier_champ=",STRING,nomchp,(int)NULL);
     }
     break;
    case CHAMPVOXEL3D:
     ch=transfo->ch;
     put_header(&header,"type=CHAMPVOXEL3D",STRING,"",(int)NULL);
     entier = transfo->width;put_header(&header, "width=", ENTIER, &entier,(int)NULL);
     entier = transfo->height;put_header(&header, "height=", ENTIER, &entier,(int)NULL);
     entier = transfo->depth;put_header(&header, "depth=", ENTIER, &entier,(int)NULL);
     /*enregistrement du fichier binaire contenant le champ*/
     {
      char nomchp[250],*nomtmp;

      sprintf(nomchp,"%s.chp",nomsext);
      fichier=fopen(nomchp,"w");
      //fwrite(ch,sizeof(vector3d),transfo->width*transfo->height*transfo->depth,fichier);

      /* pour assurer la compatibilite entre HP et linux */
      if(_is_bigEndian==0)
      {
       for(i=0;i<(int)(transfo->width*transfo->height*transfo->depth);i++)
       {
         tmp=ch[i].x;
         ch[i].x=ch[i].z;
         ch[i].z=(float)tmp;
       }
      }
      imx_fwrite(ch,(int)(transfo->width*transfo->height*transfo->depth),sizeof(vector3d),fichier);
      fclose(fichier);

      /*suppresion du chemin si precise*/
      nomtmp=strrchr(nomchp,'/');
      if (nomtmp!=NULL) strcpy(nomchp,nomtmp+1);
      put_header(&header,"fichier_champ=",STRING,nomchp,(int)NULL);
     }
     break;
   }
/*protoize:???*/ /* ok now?*/      
   save_header(&header/*???*/,nom);
 
 return(1); 
}

/*******************************************************************************
**    load_transf_3d(nomfichier)                                  
*/                                                                       
/*!     charge une transfo depuis un fichier .trf
**  alloue et retourne un  pointeur sur une transfo
**
**  \param nomfichier : le nom du fichier contenant la transfo
**  \retval transf3d : un pointeur sur la transfo, NULL en cas d'echec
*******************************************************************************/
transf3d *load_transf_3d(char *nomfichier)
{
 char   nom[256],*ext,typet[50],chemin[250],nomsext[255];
 FILE   *fichier;
 double *p,tmp;
 vector3d *ch;
 int    nb_param,wdth,hght,dpth,i;
 double dx,dy,dz;
 transf3d *transfo;
 TRANSFO3D typetrans;
 
 nb_param=0;p=NULL;
 /*test de l'extension*/
  strcpy(nomsext,nomfichier);
  ext=strstr(nomsext,".trf");
  if (ext!=NULL) *ext=0;
  sprintf(nom,"%s.trf",nomsext);
  /*extraction du chemin*/
  strcpy(chemin,nomsext);
  ext=strrchr(chemin,'/');
  if (ext!=NULL) *(ext+1)=0; else strcpy(chemin,"");
 
  /*lecture */
  fichier = fopen(nomfichier, "rb");
  if (!fichier)
  {
    printf("File does not exist: %s \n", nomfichier);
    return NULL;
  }
  fclose (fichier);
  
 wdth=atoi(getheader_interfile(nom,"width=",ENTIER,(int)NULL));
 hght=atoi(getheader_interfile(nom,"height=",ENTIER,(int)NULL));
 dpth=atoi(getheader_interfile(nom,"depth=",ENTIER,(int)NULL));
 dx=atof(getheader_interfile(nom,"dx=",REEL,(int)NULL));
 dy=atof(getheader_interfile(nom,"dy=",REEL,(int)NULL));
 dz=atof(getheader_interfile(nom,"dz=",REEL,(int)NULL));
 typetrans=RIGID3D;
 strcpy(typet,getheader_interfile(nom,"type=",STRING,(int)NULL));
 if ((strncmp(typet,"CHAMP3D",7))==0) typetrans=CHAMP3D;
 if ((strncmp(typet,"RIGID3D",7))==0) typetrans=RIGID3D;
 if ((strncmp(typet,"RIGIDZOOM3D",11))==0) typetrans=RIGIDZOOM3D;
 if ((strncmp(typet,"AFFINE3D",8))==0) typetrans=AFFINE3D;
 if ((strncmp(typet,"RIGIDGLOBALZOOM3D",17))==0) typetrans=RIGIDZOOM3D;
 if ((strncmp(typet,"AFFINEDECOUPLE3D",16))==0) typetrans=RIGIDZOOM3D;
 if ((strncmp(typet,"BSPLINE3D",9))==0) typetrans=BSPLINE3D;
 /*allocation de transfo*/
 transfo=cr_transf3d(wdth,hght,dpth,getheader_interfile(nom,"trans_ini=",STRING,(int)NULL));
 /*lecture de p*/
 switch (typetrans)
 {
  case RIGID3D:
   nb_param=atoi(getheader_interfile(nom,"nb_param=",ENTIER,(int)NULL));
   p=CALLOC(nb_param,double);
   transfo->nb_param=nb_param;transfo->param=p;
   p[0]=atof(getheader_interfile(nom,"Theta_x=",REEL,(int)NULL));
   p[1]=atof(getheader_interfile(nom,"Theta_y=",REEL,(int)NULL));
   p[2]=atof(getheader_interfile(nom,"Theta_z=",REEL,(int)NULL));
   p[3]=atof(getheader_interfile(nom,"Tx=",REEL,(int)NULL));
   p[4]=atof(getheader_interfile(nom,"Ty=",REEL,(int)NULL));
   p[5]=atof(getheader_interfile(nom,"Tz=",REEL,(int)NULL));
   p[6]=atof(getheader_interfile(nom,"Xg=",REEL,(int)NULL));
   p[7]=atof(getheader_interfile(nom,"Yg=",REEL,(int)NULL));
   p[8]=atof(getheader_interfile(nom,"Zg=",REEL,(int)NULL));
   break;
  case RIGIDZOOM3D:
   nb_param=atoi(getheader_interfile(nom,"nb_param=",ENTIER,(int)NULL));
   p=CALLOC(nb_param,double);
   transfo->nb_param=nb_param;transfo->param=p;
   p[0]=atof(getheader_interfile(nom,"Theta_x=",REEL,(int)NULL));
   p[1]=atof(getheader_interfile(nom,"Theta_y=",REEL,(int)NULL));
   p[2]=atof(getheader_interfile(nom,"Theta_z=",REEL,(int)NULL));
   p[3]=atof(getheader_interfile(nom,"Sx=",REEL,(int)NULL));
   p[4]=atof(getheader_interfile(nom,"Sy=",REEL,(int)NULL));
   p[5]=atof(getheader_interfile(nom,"Sz=",REEL,(int)NULL));
   p[6]=atof(getheader_interfile(nom,"Tx=",REEL,(int)NULL));
   p[7]=atof(getheader_interfile(nom,"Ty=",REEL,(int)NULL));
   p[8]=atof(getheader_interfile(nom,"Tz=",REEL,(int)NULL));
   p[9]=atof(getheader_interfile(nom,"Xg=",REEL,(int)NULL));
   p[10]=atof(getheader_interfile(nom,"Yg=",REEL,(int)NULL));
   p[11]=atof(getheader_interfile(nom,"Zg=",REEL,(int)NULL));
   break;
  case AFFINE3D:
   nb_param=atoi(getheader_interfile(nom,"nb_param=",ENTIER,(int)NULL));
   p=CALLOC(nb_param,double);
   transfo->nb_param=nb_param;transfo->param=p;
   p[0]=atof(getheader_interfile(nom,"a00=",REEL,(int)NULL));
   p[1]=atof(getheader_interfile(nom,"a01=",REEL,(int)NULL));
   p[2]=atof(getheader_interfile(nom,"a02=",REEL,(int)NULL));
   p[3]=atof(getheader_interfile(nom,"a10=",REEL,(int)NULL));
   p[4]=atof(getheader_interfile(nom,"a11=",REEL,(int)NULL));
   p[5]=atof(getheader_interfile(nom,"a12=",REEL,(int)NULL));
   p[6]=atof(getheader_interfile(nom,"a20=",REEL,(int)NULL));
   p[7]=atof(getheader_interfile(nom,"a21=",REEL,(int)NULL));
   p[8]=atof(getheader_interfile(nom,"a22=",REEL,(int)NULL));
   p[9]=atof(getheader_interfile(nom,"Tx=",REEL,(int)NULL));
   p[10]=atof(getheader_interfile(nom,"Ty=",REEL,(int)NULL));
   p[11]=atof(getheader_interfile(nom,"Tz=",REEL,(int)NULL));
   p[12]=atof(getheader_interfile(nom,"Xg=",REEL,(int)NULL));
   p[13]=atof(getheader_interfile(nom,"Yg=",REEL,(int)NULL));
   p[14]=atof(getheader_interfile(nom,"Zg=",REEL,(int)NULL));
   break;
  case RIGIDGLOBALZOOM3D:
   nb_param=atoi(getheader_interfile(nom,"nb_param=",ENTIER,(int)NULL));
   p=CALLOC(nb_param,double);
   transfo->nb_param=nb_param;transfo->param=p;
   p[0]=atof(getheader_interfile(nom,"Theta_x=",REEL,(int)NULL));
   p[1]=atof(getheader_interfile(nom,"Theta_y=",REEL,(int)NULL));
   p[2]=atof(getheader_interfile(nom,"Theta_z=",REEL,(int)NULL));
   p[3]=atof(getheader_interfile(nom,"S=",REEL,(int)NULL));
   p[6]=atof(getheader_interfile(nom,"Tx=",REEL,(int)NULL));
   p[7]=atof(getheader_interfile(nom,"Ty=",REEL,(int)NULL));
   p[8]=atof(getheader_interfile(nom,"Tz=",REEL,(int)NULL));
   p[9]=atof(getheader_interfile(nom,"Xg=",REEL,(int)NULL));
   p[10]=atof(getheader_interfile(nom,"Yg=",REEL,(int)NULL));
   p[11]=atof(getheader_interfile(nom,"Zg=",REEL,(int)NULL));
   break;
  case AFFINEDECOUPLE3D:
   nb_param=atoi(getheader_interfile(nom,"nb_param=",ENTIER,(int)NULL));
   p=CALLOC(nb_param,double);
   transfo->nb_param=nb_param;transfo->param=p;
   p[0]=atof(getheader_interfile(nom,"Theta_x=",REEL,(int)NULL));
   p[1]=atof(getheader_interfile(nom,"Theta_y=",REEL,(int)NULL));
   p[2]=atof(getheader_interfile(nom,"Theta_z=",REEL,(int)NULL));
   p[3]=atof(getheader_interfile(nom,"Tx=",REEL,(int)NULL));
   p[4]=atof(getheader_interfile(nom,"Ty=",REEL,(int)NULL));
   p[5]=atof(getheader_interfile(nom,"Tz=",REEL,(int)NULL));
   p[6]=atof(getheader_interfile(nom,"Sx=",REEL,(int)NULL));
   p[7]=atof(getheader_interfile(nom,"Sy=",REEL,(int)NULL));
   p[8]=atof(getheader_interfile(nom,"Sz=",REEL,(int)NULL));
   p[9]=atof(getheader_interfile(nom,"Skewyx=",REEL,(int)NULL));
   p[10]=atof(getheader_interfile(nom,"Skewzx=",REEL,(int)NULL));
   p[11]=atof(getheader_interfile(nom,"Skewzy=",REEL,(int)NULL));
   p[12]=atof(getheader_interfile(nom,"Xg=",REEL,(int)NULL));
   p[13]=atof(getheader_interfile(nom,"Yg=",REEL,(int)NULL));
   p[14]=atof(getheader_interfile(nom,"Zg=",REEL,(int)NULL));
   break;
  case BSPLINE3D:
   nb_param=atoi(getheader_interfile(nom,"nb_param=",ENTIER,(int)NULL));
   p=CALLOC(nb_param,double);
   transfo->nb_param=nb_param;transfo->param=p;
   transfo->resol=atoi(getheader_interfile(nom,"resol=",ENTIER,(int)NULL));
   transfo->degre=atoi(getheader_interfile(nom,"degre=",ENTIER,(int)NULL));
   /*lecture du fichier binaire de parametres*/
   {
   char nomprm[256],nomtmp[256];
   
   strcpy(nomtmp,getheader_interfile(nom,"fichier_param=",STRING,(int)NULL));
   nomtmp[strlen(nomtmp)-1]=0;
   sprintf(nomprm,"%s%s",chemin,nomtmp);
   fichier=fopen(nomprm,"rb");
   if (fichier==NULL) {printf("erreur a l'ouverture du fichier %s \n",nomprm); return NULL;}
   else
    {
     //fread(p,sizeof(double),nb_param,fichier);
     imx_fread(p,nb_param,sizeof(double),fichier);
     fclose(fichier);
    } 
   }
  break;
  case CHAMP3D:
   ch=CALLOC(wdth*hght*dpth,vector3d);
   transfo->ch=ch;
   /*lecture du fichier binaire contenant le champ*/
   {
   char nomprm[256],nomtmp[256];
   
   strcpy(nomtmp,getheader_interfile(nom,"fichier_champ=",STRING,(int)NULL));
   nomtmp[strlen(nomtmp)-1]=0;
   sprintf(nomprm,"%s%s",chemin,nomtmp);
   fichier=fopen(nomprm,"rb");
   if (fichier==NULL) {printf("erreur a l'ouverture du fichier %s \n",nomprm); return NULL;}
   else
    {
        //fread(ch,sizeof(vector3d),wdth*hght*dpth,fichier);
        imx_fread(ch,wdth*hght*dpth,sizeof(vector3d),fichier);
    fclose(fichier);

        if(_is_bigEndian==0)
        {
        for(i=0;i<(int)(transfo->width*transfo->height*transfo->depth);i++)
            {
        tmp=ch[i].x;
        ch[i].x=ch[i].z;
        ch[i].z=(float)tmp;
        }
        }
    } 
   }
  break;
  case CHAMPVOXEL3D:
   ch=CALLOC(wdth*hght*dpth,vector3d);
   transfo->ch=ch;
   /*lecture du fichier binaire contenant le champ*/
   {
   char nomprm[256],nomtmp[256];

   strcpy(nomtmp,getheader_interfile(nom,"fichier_champ=",STRING,(int)NULL));
   nomtmp[strlen(nomtmp)-1]=0;
   sprintf(nomprm,"%s%s",chemin,nomtmp);
   fichier=fopen(nomprm,"rb");
   if (fichier==NULL) {printf("erreur a l'ouverture du fichier %s \n",nomprm); return NULL;}
   else
    {
    //fread(ch,sizeof(vector3d),wdth*hght*dpth,fichier);
      imx_fread(ch,wdth*hght*dpth,sizeof(vector3d),fichier);
      fclose(fichier);

      if(_is_bigEndian==0)
      {
        for(i=0;i<(int)(transfo->width*transfo->height*transfo->depth);i++)
        {
          tmp=ch[i].x;
          ch[i].x=ch[i].z;
          ch[i].z=(float)tmp;
        }
      }
    }
   }
   break;
 }
  
 transfo->typetrans=typetrans;
 transfo->dx = (float)dx;
 transfo->dy = (float)dy;
 transfo->dz = (float)dz; 
 return(transfo); 
}

/*******************************************************************************
********************************************************************************
******************************* DIVERS  ****************************************
********************************************************************************
*******************************************************************************/

/*!

!!!!! NE PAS UTILISER   !!!!!!!!!!

  redimensionne une transformation 3D
  \param wdth,hght,dpth : les nouveaux parametres de taille
  \param transfo : la transformation a redimensionner (E/S)
  \retval 0 si succes, 1 si echec
!!!!! NE PAS UTILISER   !!!!!!!!!!

*/
int resize_transf_3d(transf3d *transfo, int wdth, int hght, int dpth)
{
  int err = 0;
  int i;
  double *param = transfo->param;

fprintf(stderr," PROBLEME : UTILISATION DE LA FONCTION resize_transf_3d ???\n");
  switch (transfo->typetrans)
  {
  case RIGID3D:
    /* tx, ty, tz */
    param[3] *= (double) wdth / transfo->width;
    param[4] *= (double) hght / transfo->height;
    param[5] *= (double) dpth / transfo->depth; 
    /* xc, yc, zc */
    param[6] *= (double) wdth / transfo->width;
    param[7] *= (double) hght / transfo->height;
    param[8] *= (double) dpth / transfo->depth; 
    /* et pour finir width, height et depth */
    transfo->width = wdth;
    transfo->height = hght;
    transfo->depth = dpth;
    break; 
  case RIGIDZOOM3D:
    /* tx, ty, tz */
    param[6] *= (double) wdth / transfo->width;
    param[7] *= (double) hght / transfo->height;
    param[8] *= (double) dpth / transfo->depth; 
    /* xc, yc, zc */
    param[9] *= (double) wdth / transfo->width;
    param[10] *= (double) hght / transfo->height;
    param[11] *= (double) dpth / transfo->depth; 
    /* et pour finir width, height et depth */
    transfo->width = wdth;
    transfo->height = hght;
    transfo->depth = dpth;
    break; 
  case AFFINE3D:
    /* tx, ty, tz */
    param[9] *= (double) wdth / transfo->width;
    param[10] *= (double) hght / transfo->height;
    param[11] *= (double) dpth / transfo->depth; 
    /* xc, yc, zc */
    param[12] *= (double) wdth / transfo->width;
    param[13] *= (double) hght / transfo->height;
    param[14] *= (double) dpth / transfo->depth; 
    /* et pour finir width, height et depth */
    transfo->width = wdth;
    transfo->height = hght;
    transfo->depth = dpth;
    break;
  case RIGIDGLOBALZOOM3D:
    /* tx, ty, tz */
    param[4] *= (double) wdth / transfo->width;
    param[5] *= (double) hght / transfo->height;
    param[6] *= (double) dpth / transfo->depth;
    /* xc, yc, zc */
    param[7] *= (double) wdth / transfo->width;
    param[8] *= (double) hght / transfo->height;
    param[9] *= (double) dpth / transfo->depth;
    /* et pour finir width, height et depth */
    transfo->width = wdth;
    transfo->height = hght;
    transfo->depth = dpth;
    break;
  case AFFINEDECOUPLE3D:
    /* tx, ty, tz */
    param[3] *= (double) wdth / transfo->width;
    param[4] *= (double) hght / transfo->height;
    param[5] *= (double) dpth / transfo->depth;
    /* xc, yc, zc */
    param[12] *= (double) wdth / transfo->width;
    param[13] *= (double) hght / transfo->height;
    param[14] *= (double) dpth / transfo->depth;
    /* et pour finir width, height et depth */
    transfo->width = wdth;
    transfo->height = hght;
    transfo->depth = dpth;
    break; 
  case BSPLINE3D:
    /* tous les parametres sont des vecteurs a mettre a l'echelle */
    for (i=0; i<transfo->nb_param/3; i++)
    {
      param[3*i] *= (double) wdth / transfo->width;
      param[3*i+1] *= (double) hght / transfo->height;
      param[3*i+2] *= (double) dpth / transfo->depth;
    }
    /* et pour finir width, height et depth */
    transfo->width = wdth;
    transfo->height = hght;
    transfo->depth = dpth;
    break;
  case CHAMP3D:
    /* la, on nage dans l'inconnu, il faudrait carrement interpoler ou
    sous-echantillonner le champ */
    err = 1;
    break;
  case CHAMPVOXEL3D:
    /* la, on nage dans l'inconnu, il faudrait carrement interpoler ou
    sous-echantillonner le champ */
    err = 1;
    break;
  }

  return err;
}

/*!
  copie une transformation 3D
  \param src : transformation source
  \param dest : trasformation destination (doit etre allouee)
  \retval 0 si succes, 1 si echec
*/
int copie_transf_3d(transf3d *src, transf3d *dest)
{
  int err = 0;
  int i,j,k,l;
  int wdth,hght,dpth;

  if (dest->param)
    FREE(dest->param);
  dest->param = NULL;
  if (dest->ch)
    FREE(dest->ch);
  dest->ch = NULL;

  //copie de parametres de dimension
  wdth = src->width;
  hght = src->height;
  dpth = src->depth;
  dest->width = wdth; dest->height = hght; dest->depth = dpth;
  
  dest->nb_param = src->nb_param;
  dest->dx = src->dx;
  dest->dy = src->dy;
  dest->dz = src->dz;
  if (src->nb_param)
    dest->param = MALLOC(src->nb_param,double);
  for (i=0;i<src->nb_param;i++)
    dest->param[i]=src->param[i];

  dest->typetrans = src->typetrans;

  if (src->typetrans == CHAMP3D)
  {
    dest->ch=CALLOC(wdth*hght*dpth,vector3d);
    l=0;
    for (i=0;i<wdth;i++)
     for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
      {dest->ch[l]=src->ch[l];l++;} 
  }

  return err;
}

int realloc_transf_param_3d(transf3d *transf,int nb_param)
{
  int i;
  int err=0;
  double *new_param;
  
  new_param = MALLOC(nb_param,double);
  if (transf->param)
  {
    for (i=0;i<transf->nb_param;i++)
    {
      new_param[i]=transf->param[i];
    } 
    free(transf->param);
  }
  transf->param = new_param;
  return err;
}


/*******************************************************************************
**   apply_transf_3d()                                  
*/
/*!                                                                       
**      appliquer a une image une transformation contenu dans un fichier          
*******************************************************************************/

void apply_transf_3d(void)
{
  int im_deb, im_res, err=0, inter_type;
  char nomfichier[PATH_LEN]/*, *quest[6]*/;


  /*question sur les images*/
  im_deb = GET_PLACE3D(TEXT0030);
  im_res = GET_PLACE3D(TEXT0006);

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

  /*methode d'interpolation*/
  inter_type = imx_query_interpolation_type_3D(0);

  imx_apply_transf_3d(im_res, im_deb, nomfichier, inter_type);

  show_picture_3d(im_res);

}

/***********************************************************/
/*!  \ingroup Recalage
**
**
**  \brief Application d'une transformation a l'aide d'une fonction d'interpolation
**          a une image 3d 
**
**
**  \param imdeb : image source
**  \param imres : image resultat (E/S)
**  \param tranfso : la transformation 
**  \param inter_fct : fct d'interpolation
**
***********************************************************/
int imx_apply_transf_3d_p(grphic3d *imres, grphic3d *imdeb, transf3d *transfo, InterpolationFct inter_fct)
{
  field3d *champ=NULL;
  int err = 0;
  grphic3d *imTemp=NULL;

  // Creation de imtemp si necessaire
  if (imdeb==imres)
  {
   imTemp=cr_grphic3d(imdeb);
   if (!imTemp) { fprintf (stderr, "memory allocation failed in imx_apply_transf_3d_p\n"); return 1; }
   imx_copie_3d_p(imdeb, imTemp);
  }
  else imTemp=imdeb;

  //generation du champ
  if (transfo->dx && transfo->dy && transfo->dz)
  { // la transfo est en millimetre, mais on connait l'espace de depart (imdeb->dx,dy,dz) et l'espace d'arrive (transof->dx,dy,dz)
    err=resize_grphic3d_buffer(imres, transfo->width, transfo->height, transfo->depth);
    if (err) goto end_func;
    imres->dx=transfo->dx ; imres->dy=transfo->dy ; imres->dz=transfo->dz;
    imres->width=transfo->width ; imres->height=transfo->height ; imres->depth=transfo->depth ;
    champ = transf_to_field_3d(transfo,imres,imTemp);  
  }
  else
  { // on ne connait aucune info...
    err=resize_grphic3d_buffer(imres, imTemp->width, imTemp->height, imTemp->depth);
    if (err) goto end_func;
    imres->dx=imTemp->dx ; imres->dy=imTemp->dy ; imres->dz=imTemp->dz;
    imres->width=imTemp->width ; imres->height=imTemp->height ; imres->depth=imTemp->depth ;
    champ = transf_to_field_3d(transfo,NULL,NULL);  
  }

  if (! champ) { err = 1; goto end_func; }
 
  // calcul de l'image resultat avec la methode d'interpolation specifiee 
  inter_fct(imTemp,champ,imres);
  // mise a jour de imres 
  imx_inimaxminpixel_3d_p(imres);
  
end_func:
  // liberation memoire 
  if (champ) free_field3d(champ);
  if (imTemp!=imdeb) free_grphic3d(imTemp);

  return err;
}

/***********************************************************/
/*!  \ingroup Recalage
**
**
**  \brief Application d'une transformation definie dans un fichier de type *.trf
**          a une image 3d 
**
**
**  \param imdeb : image source
**  \param imres : image resultat (E/S)
**  \param nomfichier : le nom du fichier de la transformation (*.trf)
**  \param inter_type : choix de l'interpolation
**
***********************************************************/
int imx_apply_transf_3d(int im_res, int im_deb, char *nomfichier, int inter_type)
{
  transf3d *transfo;
  int err = 0;
  char tmp[PATH_LEN];
  InterpolationFct inter_fct;
  grphic3d * imdeb = NULL;
  grphic3d * imres = NULL;
  
  imdeb = ptr_img_3d(im_deb);
  imres = ptr_img_3d(im_res);
    
  put_file_extension(nomfichier,".trf",tmp);
  
  /*chargement de la transformation*/
  transfo = load_transf_3d(tmp);
  if (! transfo)
  {
    err = 1;
    goto error1;
  }
  
  inter_fct = imx_choose_interpolation_fct(inter_type);
  
  err = imx_apply_transf_3d_p(imres, imdeb, transfo, inter_fct);

  imx_reinitialise_visual_params_3d_p(imres);

  free_transf3d(transfo); 
error1:
  return err;
}

/***********************************************************/
/*!  \ingroup Recalage
**
**
**  \brief Application d'une transformation anisotropique a l'aide d'une fonction d'interpolation
**          a une image 3d
**
**
**  \param imref    : image reference
**  \param imtransf : image sur laquelle appliquer une transfo
**  \param imres    : image resultat (E/S)
**  \param tranfso  : la transformation
**  \param inter_fct: fct d'interpolation
**
***********************************************************/
//int imx_apply_transf_anisotrope_3d_p(grphic3d *imres,grphic3d *imref, grphic3d *imtransf, transf3d *transfo, InterpolationFct inter_fct)
//{
//  field3d *champ;
//  grphic3d *imtres;
//
//
//  /* Creation de imtemp, copie de imlabeled->pos dans imtemp, et modif taille */
//  /*generation du champ*/
//
//  champ = transf_to_field_3d(transfo,imref,imtransf);
//
//  if (!champ)   /* some error creating the field */
//  {
//    return 1;
//  }
//
//  imtres = cr_grphic3d(imref);
//  if (!imtres) // error in creating temp image
//  {
//    free_field3d(champ);
//    return 1;
//  }
//
//  /* calcul de l'image resultat avec la methode d'interpolation specifiee */
//  /* on ignore le code d'erreur, parce que la plupart des fonctions ne le
//  retournent pas (les fonctions CPP ont meme ete declarees comme retournant un
//  void - ce qui soit est une erreur, soit implique qu'il faut changer les protos
//  des fonctions C et du type InterpolFct) */
//  /*err =*/ inter_fct(imtransf,champ,imtres);
//
//  /* mise a jour de imres */
//  imx_inimaxminpixel_3d_p(imtres);
//  imx_copie_3d_p(imtres,imres);
//  imx_copie_dimension_params_3d_p(imref,imres);
//    /* liberation memoire */
//  free_field3d(champ);
//  free_grphic3d(imtres);
//
//  return 0;
//}


/*******************************************************************************
**   apply_transf_seuil_3d()                                  
*/
/*!                                                                       
**      appliquer a une image seuille une transformation contenu dans un fichier          
*******************************************************************************/
void apply_transf_seuil_3d(void)
{
  int im_deb, im_res, err=0, inter_type;
  char nomfichier[PATH_LEN]/*, *quest[6]*/;
    
  /*question sur les images*/
  im_deb = GET_PLACE3D(TEXT0030);
  im_res = GET_PLACE3D(TEXT0006);

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

  /*methode d'interpolation*/
  inter_type = imx_query_interpolation_type_3D(0);
 
 
        
  imx_apply_transf_seuil_3d(im_res, im_deb, nomfichier, inter_type);

  show_picture_3d(im_res);

}

/***********************************************************/
/*!  \ingroup Recalage
**
**
**  \brief Application d'une transformation definie dans un fichier de type *.trf
**          a une image 3d 
**
**
**  \param imdeb : image source
**  \param imres : image resultat (E/S)
**  \param nomfichier : le nom du fichier de la transformation (*.trf)
**  \param inter_type : choix de l'interpolation
**
***********************************************************/
int imx_apply_transf_seuil_3d(int im_res, int im_deb, char *nomfichier, int inter_type)
{
  transf3d *transfo;
  int err = 0;
  unsigned int i,j,k;
  char tmp[PATH_LEN];
  InterpolationFct inter_fct,inter_label;
  grphic3d * imdeb = NULL;
  grphic3d * imres = NULL;
    grphic3d *maskdeb,*maskres;
  
  imdeb = ptr_img_3d(im_deb);
  imres = ptr_img_3d(im_res);
  
    maskdeb=cr_grphic3d(imdeb);imx_copie_3d_p(imdeb,maskdeb);
  maskres=cr_grphic3d(imres);imx_copie_3d_p(imres,maskres);  
  
    maskdeb->max_pixel=1;
    maskdeb->min_pixel=0;
    
    for (i=0;i<imdeb->width;i++)
    for (j=0;j<imdeb->height;j++)
    for (k=0;k<imdeb->depth;k++)
        {
        if (imdeb->mri[i][j][k]!=0) maskdeb->mri[i][j][k]=1;
        }
        
    put_file_extension(nomfichier,".trf",tmp);
  
  /*chargement de la transformation*/
  transfo = load_transf_3d(tmp);
  if (! transfo)
  {
    err = 1;
    goto error1;
  }
  
  inter_fct = imx_choose_interpolation_fct(inter_type);
  inter_label = imx_choose_interpolation_fct(9);
  
  err = imx_apply_transf_3d_p(imres, imdeb, transfo, inter_fct);
    err = imx_apply_transf_3d_p(maskres, maskdeb, transfo,inter_label);
    
    for (i=0;i<imres->width;i++)
    for (j=0;j<imres->height;j++)
    for (k=0;k<imres->depth;k++)
        {
        if (maskres->mri[i][j][k]==0) imres->mri[i][j][k]=0;
        }
        
  imx_reinitialise_visual_params_3d_p(imres);

  free_transf3d(transfo); 
    free_grphic3d(maskdeb);free_grphic3d(maskres);
error1:
  return err;
}

/*******************************************************************************
**   visu_transf_3d()                                  
**                                                                       
**      visualise la transfo ou d'un champ sous forme d'un quadrillage deforme          
**  ATTENTION : on visualise une transformation en mm.!!!
**
*******************************************************************************/
void visu_transf_3d(void)
{
 field3d *champ;
 transf3d *transfo;
 int    im_res,e=0,i,j,k,type;
 char   nomfichier[256],*quest[9];
 grphic3d *imres;
 
 
 champ=NULL;transfo=NULL;
 

 /*fichier contenant le champ*/
 strcpy(nomfichier,GET_FILE("*.trf",&e));

 if(e != 0)
     return ;
 
  put_file_extension(nomfichier,".trf",nomfichier);

  if ((test_fichier_transf_3d(nomfichier))!=1) 
  {PUT_ERROR("This file contains no 3D transformation");return;}

 
 
 /*question sur les images*/
 im_res=GET_PLACE3D(TEXT0031);
 imres=ptr_img_3d(im_res);

 /*type de visualisation*/
 for(i=0;i<9;i++)
    quest[i]=CALLOC(80,char);
  strcpy(quest[0],"module");
  strcpy(quest[1],"quadrillage");
  strcpy(quest[2],"Jacobien");
  strcpy(quest[3],"Jacobien(bspline1)");
  strcpy(quest[4],"Integral Jacobien(bspline1)");
  strcpy(quest[5],"Composante suivant x");
  strcpy(quest[6],"Composante suivant y");
  strcpy(quest[7],"Composante suivant z");
  strcpy(quest[8],"\0");
  
  type=GETV_QCM("Methode",(char **)quest);
 for(i=0;i<9;i++)
    free(quest[i]);
 
  transfo=load_transf_3d(nomfichier);
   
 imres->width=transfo->width;imres->height=transfo->height;imres->depth=transfo->depth;
 
 
 
 
 
 /*calcul du quadrillage*/
 if (type==2 && transfo->typetrans==BSPLINE3D)
  jacobien_bspline1_3d(transfo, ptr_img_3d(im_res));
 else 
        {   
        if (type==3)
            jacobien_champ_interpolation_bspline1(transfo,imres);
        else
            {
            if (type==4)
                integrale_jacobien_bspline1(transfo,imres);
            else
            {
             /*generation du champ*/         
            champ=transf_to_field_3d(transfo,NULL,NULL);
            
            //if (type==2)
                {
                for(i=0;i<transfo->width;i++)
                for(j=0;j<transfo->height;j++)
                for(k=0;k<transfo->depth;k++)
                 {
                 champ->raw[i][j][k].x=champ->raw[i][j][k].x/transfo->dx;
                 champ->raw[i][j][k].y=champ->raw[i][j][k].y/transfo->dy;
                 champ->raw[i][j][k].z=champ->raw[i][j][k].z/transfo->dz;    
                 }
                }
            
            printf("Attention, le resultat est donnee en voxel et non en millimetre ! \n");         
            
            imx_visu_field_3d(champ,imres,type);
            free_field3d(champ);
            }
            }
        }
 // Mise a jour des dx,dy,dz de l'image resultat
 imres->dx=transfo->dx;
 imres->dy=transfo->dy;
 imres->dz=transfo->dz;
 

 show_picture_3d(im_res);
 free_transf3d(transfo);
 
}



/*******************************************************************************
**      jacobien_champ_interpolation_bspline1(ch,imres)
**  
**      jacobien d'un champ de deformation en interpoalant par une spline de degre 1
**      On se ram�ne au cas du modele de d�formation bspline qui � chaque voxel
**  associe un jeu de param�tre 
*******************************************************************************/
int jacobien_champ_interpolation_bspline1(transf3d *transfo_src, grphic3d *imres)
{
transf3d *transfo;
int wdth,hght,dpth,maxsize,i,j,k,l,resol;
int nb_param;
double *param;
field3d *champ;
wdth=transfo_src->width;
hght=transfo_src->height;
dpth=transfo_src->depth;

maxsize=MAXI(wdth,hght);
maxsize=MAXI(maxsize,dpth);

transfo=cr_copy_transf3d(transfo_src);
transfo->typetrans=BSPLINE3D;
transfo->degre=1;
transfo->width=maxsize;
transfo->height=maxsize;
transfo->depth=maxsize;
  
resol = (int)floor(log(maxsize)/log(2)+0.5);

nb_param=3*pow((pow(2,resol)-1),3);;
param=(double *) malloc(nb_param*sizeof(double));


champ=transf_to_field_3d(transfo_src,NULL,NULL);

for (i=0;i<nb_param;i++)
param[i]=0;
  
for (i=1;i<maxsize;i++)
for (j=1;j<maxsize;j++)
for (k=1;k<maxsize;k++)
    if ((i<wdth)&&(j<hght)&&(k<dpth))
        {
        l = TOP_conv_ind(i-1,j-1,k-1,nb_param);
        param[l]    =champ->raw[i][j][k].x;
        param[l+1]=champ->raw[i][j][k].y;
        param[l+2]=champ->raw[i][j][k].z;
        }
        
        
transfo->nb_param=nb_param;
transfo->param=param;
transfo->resol=resol;

jacobien_bspline1_3d(transfo,imres);

free(param);
free_field3d(champ);
transfo->typetrans=CHAMP3D;
free_transf3d(transfo);
return(1);

}


/*******************************************************************************
**      integrale_jacobien_bspline1(ch,imres)
**  
**      integrale du jacobien d'un champ de deformation en interpoalant par une spline de degre 1
**      On se ram�ne au cas du modele de d�formation bspline qui � chaque voxel
**  associe un jeu de param�tre 
*******************************************************************************/
int integrale_jacobien_bspline1(transf3d *transfo_src, grphic3d *imres)
{


if (transfo_src->typetrans!=BSPLINE3D) // on cree une transfo Bspline 3D a partir de la transfo donnee
    {
    transf3d *transfo;
    int wdth,hght,dpth,maxsize,i,j,k,l,resol;
    int nb_param;
    double *param;
    field3d *champ;

    wdth=transfo_src->width;
    hght=transfo_src->height;
    dpth=transfo_src->depth;

    maxsize=MAXI(wdth,hght);
    maxsize=MAXI(maxsize,dpth);

    transfo=cr_copy_transf3d(transfo_src);
    transfo->typetrans=BSPLINE3D;
    transfo->degre=1;
    transfo->width=maxsize;
    transfo->height=maxsize;
    transfo->depth=maxsize;
  
    resol = (int)floor(log(maxsize)/log(2)+0.5);

    nb_param=3*pow((pow(2,resol)-1),3);;
    param=(double *) malloc(nb_param*sizeof(double));


    champ=transf_to_field_3d(transfo_src,NULL,NULL);

    for (i=0;i<nb_param;i++)
    param[i]=0;
      
    for (i=1;i<maxsize;i++)
    for (j=1;j<maxsize;j++)
    for (k=1;k<maxsize;k++)
        if ((i<wdth)&&(j<hght)&&(k<dpth))
        {
        l = TOP_conv_ind(i-1,j-1,k-1,nb_param);
        param[l]    =champ->raw[i][j][k].x;
        param[l+1]=champ->raw[i][j][k].y;
        param[l+2]=champ->raw[i][j][k].z;
        }
        
        transfo->nb_param=nb_param;
    transfo->param=param;
    transfo->resol=resol;

    integrale_jacobien_bspline1_compute(transfo,imres);

    free(param);
    free_field3d(champ);
    transfo->typetrans=CHAMP3D;
    free_transf3d(transfo);


    }
else
    {
    integrale_jacobien_bspline1_compute(transfo_src,imres);
    }



return(1);

}


/*******************************************************************************
**      jacobien_transf3d_to_image_3d(transfo,imres)
**  
**      jacobien d'une deformation vers une images
*******************************************************************************/
int jacobien_transf3d_to_image_3d(transf3d *transfo, grphic3d *imres)
{
 int     (*scal_func)(int pas, double **f, double **df);
 double  *p;
 int     nb_param;
 double *fx,*fy,*fz,*dfx,*dfy,*dfz;
 int    *x0,*y0,*z0,*x1,*y1,*z1;
 int    wdth,hght,dpth,i,j,k,l,i0,i1,j0,j1,k0,k1;
 double aa1,bb1,cc1,aa2,bb2,cc2,aa3,bb3,cc3,d23,d31,d12,det,rcoeff,ax,ay,az     ;
 double  ***J,minJ,maxJ;
 
 
 wdth=transfo->width;hght=transfo->height;dpth=transfo->depth;
 J=alloc_dmatrix_3d(wdth,hght,dpth);
 /*pour l'instant, seulement pour la transfo en Bspline*/
 

 switch (transfo->degre)
  {
   case 0: scal_func=Bsplined0;break;
   case 1: scal_func=Bsplined1;break;
   case 2: scal_func=Bsplined2;break;
   default: scal_func=Bsplined0;break;
  }
 p=CALLOC(transfo->nb_param,double);  
 nb_param=init_base_3d(wdth,hght,dpth,scal_func);
 for (i=BASE3D.resol;i<transfo->resol;i++) nb_param=base_resol_up_3d(p,nb_param);
 if (nb_param!=transfo->nb_param) printf("ERREUR dans transf_to_field_3d\n");
 for (i=0;i<nb_param;i++) p[i]=(transfo->param)[i];
   
 x0=BASE3D.x0;y0=BASE3D.y0;z0=BASE3D.z0;x1=BASE3D.x1;y1=BASE3D.y1;z1=BASE3D.z1;
 fx=BASE3D.fx;fy=BASE3D.fy;fz=BASE3D.fz;dfx=BASE3D.dfx;dfy=BASE3D.dfy;dfz=BASE3D.dfz;
  
 minJ=100000.0;maxJ=-100000000.0;   
 for (i=0;i<wdth;i++) 
  for (j=0;j<hght;j++)
   for (k=0;k<dpth;k++)
   {
    /*initialisation de la jacobienne a zero en ce point*/
    aa1=1.0;aa2=aa3=bb1=0.0;bb2=1.0;bb3=cc1=cc2=0.0;cc3=1.0;
    
    for (l=0;l<nb_param/3;l++)
    {
     i0=x0[l];i1=x1[l];j0=y0[l];j1=y1[l];k0=z0[l];k1=z1[l];
     if (i>=i0 && j>=j0 && k>=k0 && i<i1 && j<j1 && k<k1)
     {
      ax=p[3*l];ay=p[3*l+1];az=p[3*l+2];
      aa1=aa1+ax*dfx[i-i0]*fy[j-j0]*fz[k-k0];
      bb1=bb1+ax*fx[i-i0]*dfy[j-j0]*fz[k-k0];
      cc1=cc1+ax*fx[i-i0]*fy[j-j0]*dfz[k-k0];
      aa2=aa2+ay*dfx[i-i0]*fy[j-j0]*fz[k-k0];
      bb2=bb2+ay*fx[i-i0]*dfy[j-j0]*fz[k-k0];
      cc2=cc2+ay*fx[i-i0]*fy[j-j0]*dfz[k-k0];
      aa3=aa3+az*dfx[i-i0]*fy[j-j0]*fz[k-k0];
      bb3=bb3+az*fx[i-i0]*dfy[j-j0]*fz[k-k0];
      cc3=cc3+az*fx[i-i0]*fy[j-j0]*dfz[k-k0];
     }
    }
    /*calcul du jacobien*/
    d23=bb2*cc3-bb3*cc2;d31=bb3*cc1-bb1*cc3;d12=bb1*cc2-bb2*cc1;
    det=aa1*d23+aa2*d31+aa3*d12;
    J[i][j][k]=det;
    if (det<minJ) minJ=det;
    if (det>maxJ) maxJ=det;
   }

 end_base_3d();
 free(p);
 
 /*calcul de l'image*/
 imx_brukermax_3d((float)maxJ,(float)minJ,imres);
 rcoeff=imres->rcoeff;
  
  imres->width=wdth;imres->height=hght;imres->depth=dpth;
  
  for (i=0;i<wdth;i++)
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
    {
     imres->mri[i][j][k]=(int)(J[i][j][k]/rcoeff);
    }

 printf("Jacobien minimum (essai): %.2f \n",minJ);
 imx_inimaxminpixel_3d_p(imres);
 free_dmatrix_3d(J);
 
 return(1);
}


/*******************************************************************************
**      imx_visu_field_3d(champ,imres,type)
**  
**      visualisation d'un champ sous forme d'une image
*******************************************************************************/
int imx_visu_field_3d(field3d *ch, grphic3d *imres, int type)
{
 grphic3d *imtemp;
 int i,j,k,wdth,hght,dpth,pas;
 
 wdth=ch->width;hght=ch->height;dpth=ch->depth;
 imres->width=wdth;imres->height=hght;imres->depth=dpth;
 
 
 switch (type)
 {
  case 0 :/*module*/
    module_field_to_image_3d(ch,imres);
  break;
  case 1 : /*quadrillage*/
    imtemp=cr_grphic3d(imres);
    imtemp->rcoeff=1.0;imtemp->icomp=0.0;imtemp->max_pixel=imtemp->cutoff_max=1;
    pas=8;
    for (i=0;i<wdth;i++) for (j=0;j<hght;j++) for (k=0;k<dpth;k++) imtemp->mri[i][j][k]=0;
    /*for (i=0;i<wdth;i=i+pas) for (j=0;j<hght;j++) for (k=0;k<dpth;k++) imtemp->mri[i][j][k]=1;
    for (i=0;i<wdth;i++) for (j=0;j<hght;j=j+pas) for (k=0;k<dpth;k++) imtemp->mri[i][j][k]=1;
    for (i=0;i<wdth;i++) for (j=0;j<hght;j++) for (k=0;k<dpth;k=k+pas) imtemp->mri[i][j][k]=1;*/
    for (i=0;i<wdth;i=i+pas) for (j=0;j<hght;j++) imtemp->mri[i][j][(int)(0.5*dpth)]=1;
    for (i=0;i<wdth;i++) for (j=0;j<hght;j=j+pas) imtemp->mri[i][j][(int)(0.5*dpth)]=1;
    
    for (i=0;i<wdth;i++) imtemp->mri[i][hght-1][(int)(0.5*dpth)]=1;
    for (j=0;j<hght;j++) imtemp->mri[wdth-1][j][(int)(0.5*dpth)]=1;
   
    inter_nearest_3d(imtemp,ch,imres);
    free_grphic3d(imtemp);
  break;
  case 2 : /*jaocibien*/
    jacobien_field_to_image_3d(ch,imres);
  break;
  case 5 : /* composante suivant x*/
        composante_field_to_image_3d(ch,imres,0);
  break;
    
  case 6 : /* composante suivant y*/
        composante_field_to_image_3d(ch,imres,1);
  break;
    
  case 7 : /* composante suivant z*/
        composante_field_to_image_3d(ch,imres,2);
  break;
        
 }
   
 return(1);
}




/*******************************************************************************
********************************************************************************
************************ DEFORMATION DE SYNTHESE *******************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        defor_synth_3d                                                
**                                                                        
**       genere un fichier chp contenant un champ de synthese                          
*******************************************************************************/
void defor_synth_3d_point(void)
{
 char nomfichres[500],str[500],*ext;
 int i,err=0;
 int npoints;
 vector3d *vecteur;
 vector3d *pts;
 

 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 
 /*test de l'extension*/
  ext=strstr(nomfichres,".chp");
  if (ext!=NULL) *ext=0;
  
  npoints = (int)GET_FLOAT("nbre de points a modifier", 1, &err);
  if (npoints<=0)
    npoints=1;
  vecteur = malloc(npoints*sizeof(vector3d));
  pts = malloc(npoints*sizeof(vector3d));
    
 for (i=0;i<npoints;i++)
    {
    pts[i].x=(int)GET_FLOAT("coord pts en X", 1, &err);
    pts[i].y=(int)GET_FLOAT("coord pts en Y", 1, &err);
    pts[i].z=(int)GET_FLOAT("coord pts en Z", 1, &err);     
    vecteur[i].x = GET_FLOAT("vecteur deplacement en X", 1, &err);
    vecteur[i].y = GET_FLOAT("vecteur deplacement en Y", 1, &err);
    vecteur[i].z = GET_FLOAT("vecteur deplacement en Z", 1, &err);      
    }
    
 imx_defor_synth_3d_point(nomfichres,256,256,256,npoints,pts,vecteur);

 free(vecteur);
 free(pts);
 
}


/*******************************************************************************
**        imx_defor_synth_3d                                                
**                                                                        
**       genere un fichier chp contenant un champ de synthese                          
*******************************************************************************/
int imx_defor_synth_3d_point(char *nomfichres, int wdth, int hght, int dpth,int npts,vector3d *Pts,vector3d *Dep)
{
 void Convolve(float *y, long int nbpts, int wnd,double * filter);
 field3d  *ch;
 transf3d *transfo;
 int    i,j,k;
 vector3d ***data;
 int taille=7;
 double *triangle;
 float *timecrs;
 triangle = malloc((2*taille+1)*sizeof(double));
 for (i=0;i<taille;i++)
    triangle[i] = triangle[2*taille+-i] = (double)(i+1)/(double)(taille+1);
 triangle[taille] = 1;
 //Le triangle n est pas normalise sinon, ca changerait la valeur entree par user...
    
 
 ch=cr_field3d(wdth,hght,dpth);
 data=ch->raw;


for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
    data[i][j][k].x=data[i][j][k].y=data[i][j][k].z=0;

for (i=0;i<npts;i++)
    {
    data[(int)Pts[i].x][(int)Pts[i].y][(int)Pts[i].z].x = Dep[i].x;     
    data[(int)Pts[i].x][(int)Pts[i].y][(int)Pts[i].z].y = Dep[i].y;
    data[(int)Pts[i].x][(int)Pts[i].y][(int)Pts[i].z].z = Dep[i].z;
    }
        
 /* Filtrage dans le sens des X */
  if( (timecrs = (float*)malloc(sizeof(float)*wdth)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return 0;
  }
    
  for(j=0 ; j<hght ; j++)
    for(k=0 ; k<dpth ; k++)
    {
    
      for (i=0;i<wdth;i++)
        timecrs[i]=data[i][j][k].x;
    
      Convolve(timecrs,wdth,2*taille+1,triangle);
    
      for (i=0;i<wdth;i++)
        data[i][j][k].x = timecrs[i];
    
      for (i=0;i<wdth;i++)
        timecrs[i]=data[i][j][k].y;
    
    if (j)                      
      Convolve(timecrs,wdth,2*taille+1,triangle);
    
      for (i=0;i<wdth;i++)
        data[i][j][k].y = timecrs[i];
    
      for (i=0;i<wdth;i++)
        timecrs[i]=data[i][j][k].z;
    
      Convolve(timecrs,wdth,2*taille+1,triangle);
    
      for (i=0;i<wdth;i++)
        data[i][j][k].z = timecrs[i];
    

    }

  free(timecrs);

  /* Filtrage dans le sens des Y */
  if( (timecrs = (float*)malloc(sizeof(float)*hght)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return 0;
  }
  for(i=0 ; i<wdth ; i++)
    for(k=0 ; k<dpth ; k++)
    {
      for (j=0;j<hght;j++)
        timecrs[j] = data[i][j][k].x;
      Convolve(timecrs,wdth,2*taille+1,triangle);
      for (j=0;j<hght;j++)
        data[i][j][k].x = timecrs[j];
      for (j=0;j<hght;j++)
        timecrs[j] = data[i][j][k].y;
      Convolve(timecrs,wdth,2*taille+1,triangle);
      for (j=0;j<hght;j++)
        data[i][j][k].y = timecrs[j];
      for (j=0;j<hght;j++)
        timecrs[j] = data[i][j][k].z;
      Convolve(timecrs,wdth,2*taille+1,triangle);
      for (j=0;j<hght;j++)
        data[i][j][k].z = timecrs[j];
    }

  free(timecrs);

    /*Filtrage  Dans le sens des Z */
  if( (timecrs = (float*)malloc(sizeof(float)*dpth)) == NULL)
  {
    PUT_WNDMSG( ERR_MEM_ALLOC);
    return 0;
  }
  for(i=0 ; i<wdth ; i++)
    for(j=0 ; j<hght ; j++)
    {
      for (k=0;k<dpth;k++)
        timecrs[k]=data[i][j][k].x;
      Convolve(timecrs,wdth,2*taille+1,triangle);
      for (k=0;k<dpth;k++)
        data[i][j][k].x=timecrs[k];
      for (k=0;k<dpth;k++)
        timecrs[k]=data[i][j][k].y;
      Convolve(timecrs,wdth,2*taille+1,triangle);
      for (k=0;k<dpth;k++)
        data[i][j][k].y=timecrs[k];
      for (k=0;k<dpth;k++)
        timecrs[k]=data[i][j][k].z;
      Convolve(timecrs,wdth,2*taille+1,triangle);
      for (k=0;k<dpth;k++)
        data[i][j][k].z=timecrs[k];   
    }

 free(timecrs);
 free(triangle);
 transfo=field_to_transf_3d(ch,NULL,NULL);
 
 save_transf_3d(transfo,nomfichres);
  
 free_transf3d(transfo);
 free_field3d(ch);

 
 return(1);
}

void Convolve(float *y, long int nbpts, int wnd,double * filter)
{
  int l,p,swap;
  float *yres;

  swap = (wnd-1)/2; 
  
  if( (yres = (float*)malloc(sizeof(float)*nbpts)) == NULL)
  {
    printf("pb allocation memmoire\n");    
//    return 0;
  }


  for (l=0;l<nbpts;l++)
    yres[l]=0.;


  for (l=swap;l<nbpts-swap;l++)
  {
      for(p=-swap;p<=swap;p++)
          yres[l]+=(double)y[l+p]*filter[p+swap];
  }

  for (l=0;l<nbpts;l++)
    y[l]=yres[l];

  free(yres);
}


/*******************************************************************************
**        defor_synth_3d                                                
**                                                                        
**       genere un fichier chp contenant un champ de synthese                          
*******************************************************************************/
void defor_synth_3d(void)
{
 char nomfichres[500],str[500],*ext;
 int err=0;

 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 
 /*test de l'extension*/
  ext=strstr(nomfichres,".chp");
  if (ext!=NULL) *ext=0;
  

   imx_defor_synth_3d(nomfichres,128,128,128);
 
 return;
}


/*******************************************************************************
**        imx_defor_synth_3d                                                
**                                                                        
**       genere un fichier chp contenant un champ de synthese                          
*******************************************************************************/
int imx_defor_synth_3d(char *nomfichres, int wdth, int hght, int dpth)
{
 field3d  *ch;
 transf3d *transfo;
 vector3d *points,*vecteurs;
 int    nbp,i,j,k,l,pas;
 
 
 
 ch=cr_field3d(wdth,hght,dpth);

 pas=30;
 

 /*creation d'un ensemble de vecteurs*/
 nbp=(wdth/pas)*(hght/pas)*(dpth/pas);
 points=CALLOC(nbp,vector3d);
 vecteurs=CALLOC(nbp,vector3d);
 l=0;
 for (i=pas/2;i<wdth;i=i+pas)
  for (j=pas/2;j<hght;j=j+pas)
   for (k=pas/2;k<dpth;k=k+pas)
   { 
    points[l].x=(float)i;
    points[l].y=(float)j;
    points[l].z=(float)k;
    vecteurs[l].x=(float)((rand()-16384)*10/16384);
    vecteurs[l].y=(float)((rand()-16384)*10/16384);
    vecteurs[l].z=(float)((rand()-16384)*10/16384);
    /*printf("%f %f %f %f %f
    %f\n",points[i].x,points[i].y,points[i].z,vecteurs[i].x,vecteurs[i].y,vecteurs[i].z);*/
    l++;
   }
 
 printf("nbp:%d l:%d \n",nbp,l);
 inter_field_thinplate_3d(ch,points,vecteurs,nbp);

 transfo=field_to_transf_3d(ch,NULL,NULL);
 
 save_transf_3d(transfo,nomfichres);
  
 free_transf3d(transfo);
 free_field3d(ch);

 
 return(1);
}



/*******************************************************************************
********************************************************************************
************************ INTERPOLATION D'UN CHAMP ******************************
*********************A PARTIR DE VECTEURS SPECIFIE *****************************
********************************************************************************
*******************************************************************************/

/*******************************************************************************
**        inter_field_thinplate_3d()                                              
**                                                                        
**       calcul le champ interpole (en thin plate spline a partir des donnees 
**  contenues dans les talbeau points et vecteur
**  Le champ est suppos�alloue �la bonne taille avant l'appel
**  ch:champ a interpoler
**  points:tableau contenant les points ou le champ est connu
**  vecteurs:tableau contenant la valeur du champ au points definis dans points
**  nbp: nombre de points                        
*******************************************************************************/
int inter_field_thinplate_3d(field3d *ch, vector3d *points, vector3d *vecteurs, int nb_points)
{
  int    i,j,k,m,iu,ju,ku,wdth,hght,dpth;
  double *xp,*yp,*zp,*vx,*vy,*vz;
  double **L,**Linv; /*meme notation que dans l'article de reference Bookstein89*/
  double *Yx,*Yy,*Yz,*Wx,*Wy,*Wz;
  double r,***Ur;
  int    xpc,ypc,zpc;  
  vector3d ***data;
  

  data=ch->raw;

  
  wdth=ch->width;hght=ch->height;dpth=ch->depth;
  
  /*allocation memoire des matrices et vecteurs*/
  L=alloc_dmatrix(nb_points+4,nb_points+4);
  Yx=alloc_dvector(nb_points+4);
  Yy=alloc_dvector(nb_points+4);
  Yz=alloc_dvector(nb_points+4);
  xp=CALLOC(nb_points,double);
  yp=CALLOC(nb_points,double);
  zp=CALLOC(nb_points,double);
  vx=CALLOC(nb_points,double);
  vy=CALLOC(nb_points,double);
  vz=CALLOC(nb_points,double);
  Ur=alloc_dmatrix_3d(wdth,hght,dpth);
  
  /*precalcul de la matrice Ur*/
  for(i=0;i<wdth;i++) 
   for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
     {
      r=sqrt((double)i*(double)i+(double)j*(double)j+(double)k*(double)k);
      Ur[i][j][k]=r;
     }

  /*remplissage de xp,yp,vx et vy */
  for (i=0;i<nb_points;i++)
  {
   xp[i]=(double)(points[i].x);yp[i]=(double)(points[i].y);zp[i]=(double)(points[i].z);
   vx[i]=(double)(vecteurs[i].x);vy[i]=(double)(vecteurs[i].y);vz[i]=(double)(vecteurs[i].z);
  }
 
  /*remplissage de la matrice L*/
  /*on utilise le fait que L est symetrique=> on ne rempli que la moitie*/
  /*remplissage d la diagonale*/
  for (i=0;i<nb_points+4;i++) L[i][i]=0.0;
  /*remplissage de la moitie en haut a droite*/
  for (i=0;i<nb_points-1;i++)
   for (j=i+1;j<nb_points;j++)
    {
     r=sqrt((xp[i]-xp[j])*(xp[i]-xp[j])+(yp[i]-yp[j])*(yp[i]-yp[j])+(zp[i]-zp[j])*(zp[i]-zp[j]));
     L[i][j]=r;
    }
  for (i=0;i<nb_points;i++)
    {
     L[i][nb_points]=1.0;
     L[i][nb_points+1]=xp[i];
     L[i][nb_points+2]=yp[i];
     L[i][nb_points+3]=zp[i];
    }
  for (i=nb_points;i<nb_points+4;i++)
   for (j=i+1;j<nb_points+4;j++)
    {
     L[i][j]=0.0;
    }
  /*maintenant on exploite la symetrie*/
  for (i=1;i<nb_points+4;i++)
   for (j=0;j<i;j++)
     L[i][j]=L[j][i];  
  
  /*inversion de la matrice L*/
   Linv=matrix_inversion(L,nb_points+4);

   
  /*remplissage des vecteurs Yx et Yy et Yz*/
  for(i=0;i<nb_points;i++)
   {Yx[i]=xp[i]+vx[i];Yy[i]=yp[i]+vy[i];Yz[i]=zp[i]+vz[i];}
  for (i=nb_points;i<nb_points+4;i++)
   {Yx[i]=0.0;Yy[i]=0.0;Yz[i]=0.0;}
   
  /*calcul des vecteurs Wx et Wy*/
  Wx=matrix_multiply_vector(Linv,nb_points+4,nb_points+4,Yx);
  Wy=matrix_multiply_vector(Linv,nb_points+4,nb_points+4,Yy);
  Wz=matrix_multiply_vector(Linv,nb_points+4,nb_points+4,Yz);

  /*Calcul du champ resultat*/
   /*calcul de la composante lineaire*/
   for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
     for (k=0;k<dpth;k++)
      {
       data[i][j][k].x=(float)(Wx[nb_points]+(Wx[nb_points+1]-1.0)*(double)i+Wx[nb_points+2]*(double)j+Wx[nb_points+3]*(double)k);
       data[i][j][k].y=(float)(Wy[nb_points]+Wy[nb_points+1]*(double)i+(Wy[nb_points+2]-1.0)*(double)j+Wy[nb_points+3]*(double)k);    
       data[i][j][k].z=(float)(Wz[nb_points]+Wz[nb_points+1]*(double)i+Wz[nb_points+2]*(double)j+(Wz[nb_points+3]-1.0)*(double)k);    
      }
   /*calcul des composantes en UR*/
   for (m=0;m<nb_points;m++)
   {
    xpc=(int)(xp[m]);ypc=(int)(yp[m]);zpc=(int)(zp[m]);
    for (i=0;i<wdth;i++)
     for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++)
       {
        if (i<=xpc) iu=xpc-i; else iu=i-xpc;
        if (j<=ypc) ju=ypc-j; else ju=j-ypc;
        if (k<=zpc) ku=zpc-k; else ku=k-zpc;
        data[i][j][k].x+=(float)(Wx[m]*Ur[iu][ju][ku]);
        data[i][j][k].y+=(float)(Wy[m]*Ur[iu][ju][ku]);
        data[i][j][k].z+=(float)(Wz[m]*Ur[iu][ju][ku]);
       }
   }

  /*liberation memoire des matrices et vecteurs*/
  free_dmatrix(L,nb_points+4,nb_points+4);
  free_dmatrix(Linv,nb_points+4,nb_points+4);
  free_dvector(Yx,nb_points+4);
  free_dvector(Yy,nb_points+4);
  free_dvector(Yz,nb_points+4);
  free_dvector(Wx,nb_points+4);
  free_dvector(Wy,nb_points+4);
  free_dvector(Wz,nb_points+4);
  free_dmatrix_3d(Ur);
  free(xp);free(yp);free(zp);free(vx);free(vy);free(vz);  
  
  return(1);

}

/*fonction de test - ne sert a rien - a detruire*/
void affine_to_affinesscg()
{
/*
    char * src;
    char * dst;
    int e;
    transf3d *transfo1;
    
    src = strdup(GET_FILE(POSx_SAVE,POSy_SAVE,"transfo affine",&e));
    dst = strdup(GET_FILE(POSx_SAVE,POSy_SAVE,"transfo a sauvegarder",&e));
    
    transfo1=load_transf_3d(src);
    transfo1->nb_param  = affine_to_affinesscg_3d(transfo1->param);
    transfo1->typetrans = AFFINE3D;
    save_transf_3d(transfo1, dst);
    free(src);
    free(dst);
    free_transf3d(transfo1);
*/
}

//----------------------------------//
//                                  //
//   INVERSION DE TRANSFORMATION    //
//                                  //
//----------------------------------//

/*******************************************************************************
**     inverser_transformation_3d()
*/
/*!
**     inverse une transformation 3D
**
*******************************************************************************/
void inverser_transformation_3d()
{
 char nomfichdeb[300]="\0", nomfichres[300]="\0";
 int err=0;

 //fichier contenant la transformation a inverser
 {
  sprintf(nomfichdeb,"%s",GET_FILE("fichier de transformation a inverser",&err));

    if (err)
      return ;

  put_file_extension(nomfichdeb,".trf",nomfichdeb);

    if (test_fichier_transf_3d(nomfichdeb) != 1)
    {
      PUT_ERR("This file contains no 3D transformation");
      return ;
    }   

 }

 //fichier contenant la transformation resultat
 {
  sprintf(nomfichres,"%s",GET_FILE("fichier de la transformatino resultat",&err));

  if (err) return ;

  put_file_extension(nomfichres,".trf",nomfichres);

 }

 //inversion de la transformation
 imx_inverser_transformation_3d(nomfichdeb, nomfichres);

}

/*******************************************************************************
**     imx_inverser_transformation_3d(nom_fich_transfo_deb, nom_fich_transfo_res)
*/
/*!
**     inverse une transformation 3D
**
**  \param  nom_fich_transfo_deb : nom du fichier contenant la transformation en entree
**  \param  nom_fich_transfo_res : nom du fichier pour enregistrer la transformation en sortie
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_inverser_transformation_3d(char *nom_fich_transfo_deb, char *nom_fich_transfo_res)
{
 transf3d *transfodeb=NULL, *transfores=NULL;
 int res=0;

 //chargement de la transformation a inverser
 transfodeb = load_transf_3d(nom_fich_transfo_deb);
 
 if (!transfodeb)
 {
  fprintf(stderr,"Impossible de charger la transfo a inverser dans imx_inverser_transformation_3d\n");
  res=1; goto end_func;
 }

 transfores=cr_copy_transf3d(transfodeb);
 if (!transfores)
 {
  fprintf(stderr,"Erreur d'allocation memoire dans imx_inverser_transformation_3d\n");
  res=2; goto end_func;
 }

 //calcul de la transformation inverse
 res=imx_inverser_transformation_3d_p(transfodeb, transfores);
 if (res) goto end_func;

 //enregistrement de la transfo resultat
 if (nom_fich_transfo_res!=NULL)
 {
  save_transf_3d(transfores,nom_fich_transfo_res);
 }

end_func :

 //desallocations memoire
 if (transfores) free_transf3d(transfores);
 if (transfodeb) free_transf3d(transfodeb);

 return res;
}

/*******************************************************************************
**     imx_inverser_transformation_3d_p(transfo_deb, transfo_res)
*/
/*!
**     inverse une transformation 3D
**
**  \param  transfos_deb : la transformation en entree
**  \param  transfo_res : la transformation en sortie
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_inverser_transformation_3d_p(transf3d *transfo_deb, transf3d *transfo_res)
{
 int res=0;
 
 //verifications d'usage
 if ((!transfo_deb)||(!transfo_res))
 {
  fprintf(stderr,"parametres d'entree null dans imx_inverser_transformation_3d_p\n");
  return 1;
 }

 switch (transfo_deb->typetrans)
 {
  case RIGID3D:      res=imx_inverser_transformation_rigide_3d_p(transfo_deb, transfo_res);
                     break;
  case RIGIDZOOM3D: res=imx_inverser_transformation_rigidezoom_3d_p(transfo_deb, transfo_res);
                     break;
  case AFFINE3D:     res=imx_inverser_transformation_affine_3d_p(transfo_deb, transfo_res);
                     break;
  default: fprintf(stderr,"on ne sait pas inverser ce type de transformation  dans imx_inverser_transformation_3d_p\n");
           res=2; break;  
 }
  
 return res;
}

/*******************************************************************************
**     imx_inverser_transformation_rigide_3d_p(transfo_deb, transfo_res)
*/
/*!
**     inverse une transformation rigide 3D
**
**  \param  transfos_deb : la transformation en entree
**  \param  transfo_res : la transformation en sortie
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_inverser_transformation_rigide_3d_p(transf3d *transfo_deb, transf3d *transfo_res)
{
 int res=0, i, nb_param=0;
 double **R=NULL, **RI=NULL;
 double tetax,tetay,tetaz;
 double cxx,sxx,cyy,syy,czz,szz;
 double *paramdeb=NULL, *paramres=NULL;
 double *trans_res=NULL;
  
 //verifications d'usage
 if ((!transfo_deb)||(!transfo_res))
 {
  fprintf(stderr,"parametres d'entree null dans imx_inverser_transformation_rigide_3d_p\n");
  res=1; goto end_func;
 }

 if (transfo_deb->typetrans!=RIGID3D)
 {
  fprintf(stderr,"mauvais type de transformation en entree de imx_inverser_transformation_rigide_3d_p\n");
  res=1; goto end_func;
 }

 //init
 res=copie_transf_3d(transfo_deb, transfo_res);
 if (res)
 {
  fprintf(stderr,"erreur dans la copie de la transformation dans imx_inverser_transformation_rigide_3d_p\n");
  res=2; goto end_func;
 }  
 paramdeb=transfo_deb->param; paramres=transfo_res->param; nb_param=transfo_deb->nb_param;
 for (i=0;i<nb_param;i++) paramres[i]=0.0;

 R=alloc_dmatrix(3,3);
 if (!R)
 {
  fprintf(stderr,"Erreur d'allocation memoire dans imx_inverser_transformation_rigide_3d_p\n");
  res=2; goto end_func;
 }

 //calcul de la matrice de rotation
 tetax=paramdeb[0]; tetay=paramdeb[1]; tetaz=paramdeb[2];
 cxx=cos(tetax);sxx=sin(tetax);cyy=cos(tetay);syy=sin(tetay);czz=cos(tetaz);szz=sin(tetaz);
 R[0][0]=cyy*czz;              R[0][1]=-cyy*szz;             R[0][2]=-syy;
 R[1][0]=cxx*szz-sxx*syy*czz;  R[1][1]=cxx*czz+sxx*syy*szz;  R[1][2]=-sxx*cyy;
 R[2][0]=sxx*szz+cxx*syy*czz;  R[2][1]=sxx*czz-cxx*syy*szz;  R[2][2]=cxx*cyy;
 //calcul de la matrice de rotation inverse
 RI=matrix_inversion (R, 3);
 if (!RI)
 {
  fprintf(stderr,"Erreur d'inversion de matrice de rotation dans imx_inverser_transformation_rigide_3d_p\n");
  res=3; goto end_func;
 }
 //recuperation de la translation resultat
 trans_res=matrix_multiply_vector(RI, 3, 3, &(paramdeb[3]));
 if (!trans_res)
 {
  fprintf(stderr,"Erreur de multiplication de matrice dans imx_inverser_transformation_rigide_3d_p\n");
  res=3; goto end_func;
 }
 for (i=0;i<3;i++) paramres[i+3]=-trans_res[i];

 //recuperation des angles de la rotation resultat
 res=decomposer_rotation(RI, paramres);
 if (res)
 {
  fprintf(stderr,"Erreur dedecomposition de rotation dans imx_inverser_transformation_rigide_3d_p\n");
  res=3; goto end_func;
 }

 //recuperation du centre
 for (i=6;i<9;i++) paramres[i]=paramdeb[i];

end_func:

 if (R) free_dmatrix(R,3,3);
 if (RI) free_dmatrix(RI,3,3);
 if (trans_res) free_dvector(trans_res, 3);
 
 return res; 
}

/*******************************************************************************
**     imx_inverser_transformation_affine_3d_p(transfo_deb, transfo_res)
*/
/*!
**     inverse une transformation affine 3D
**
**  \param  transfos_deb : la transformation en entree
**  \param  transfo_res : la transformation en sortie
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_inverser_transformation_affine_3d_p(transf3d *transfo_deb, transf3d *transfo_res)
{
 int res=0, i, nb_param=0;
 double **A=NULL, **AI=NULL;
 double *paramdeb=NULL, *paramres=NULL;
 double *trans_res=NULL;
  
 //verifications d'usage
 if ((!transfo_deb)||(!transfo_res))
 {
  fprintf(stderr,"parametres d'entree null dans imx_inverser_transformation_affine_3d_p\n");
  res=1; goto end_func;
 }
 
 if (transfo_deb->typetrans!=AFFINE3D)
 {
  fprintf(stderr,"mauvais type de transformation en entree de imx_inverser_transformation_affine_3d_p\n");
  res=1; goto end_func;
 }

 //init
 res=copie_transf_3d(transfo_deb, transfo_res);
 if (res)
 {
  fprintf(stderr,"erreur dans la copie de la transformation dans imx_inverser_transformation_affine_3d_p\n");
  res=2; goto end_func;
 }
 paramdeb=transfo_deb->param; paramres=transfo_res->param; nb_param=transfo_deb->nb_param;
 for (i=0;i<nb_param;i++) paramres[i]=0.0;

 A=alloc_dmatrix(3,3);
 if (!A)
 {
  fprintf(stderr,"Erreur d'allocation memoire dans imx_inverser_transformation_affine_3d_p\n");
  res=2; goto end_func;
 }

 //calcul de la matrice affine
 A[0][0]=paramdeb[0];  A[0][1]=paramdeb[1];  A[0][2]=paramdeb[2];
 A[1][0]=paramdeb[3];  A[1][1]=paramdeb[4];  A[1][2]=paramdeb[5];
 A[2][0]=paramdeb[6];  A[2][1]=paramdeb[7];  A[2][2]=paramdeb[8];
 //calcul de la matrice de rotation inverse
 AI=matrix_inversion (A, 3);
 if (!AI)
 {
  fprintf(stderr,"Erreur d'inversion de matrice de rotation dans imx_inverser_transformation_affine_3d_p\n");
  res=3; goto end_func;
 }
 
 //recuperation des parametres affine
 paramres[0]=AI[0][0];  paramres[1]=AI[0][1];  paramres[2]=AI[0][2];
 paramres[3]=AI[1][0];  paramres[4]=AI[1][1];  paramres[5]=AI[1][2];
 paramres[6]=AI[2][0];  paramres[7]=AI[2][1];  paramres[8]=AI[2][2];
 
 //recuperation de la translation resultat
 trans_res=matrix_multiply_vector(AI, 3, 3, &(paramdeb[9]));
 if (!trans_res)
 {
  fprintf(stderr,"Erreur de multiplication de matrice dans imx_inverser_transformation_affine_3d_p\n");
  res=3; goto end_func;
 }
 for (i=0;i<3;i++) paramres[i+9]=-trans_res[i];

 //recuperation du centre
 for (i=12;i<15;i++) paramres[i]=paramdeb[i];

end_func:

 if (A) free_dmatrix(A,3,3);
 if (AI) free_dmatrix(AI,3,3);
 if (trans_res) free_dvector(trans_res, 3);

 return res;

}

/*******************************************************************************
**     imx_inverser_transformation_rigidezoom_3d_p(transfo_deb, transfo_res)
*/
/*!
**     inverse une transformation rigide + zoom 3D
**
**  \param  transfos_deb : la transformation en entree
**  \param  transfo_res : la transformation en sortie
**  \retval : 0 en cas de succes
*******************************************************************************/
int imx_inverser_transformation_rigidezoom_3d_p(transf3d *transfo_deb, transf3d *transfo_res)
{
 int res=0, i;
 transf3d *transfo_tmp=NULL;
 double *tmp_param=NULL;

  //verifications d'usage
 if ((!transfo_deb)||(!transfo_res))
 {
  fprintf(stderr,"parametres d'entree null dans imx_inverser_transformation_rigidezoom_3d_p\n");
  res=1; goto end_func;
 }
 
 if (transfo_deb->typetrans!=RIGIDZOOM3D)
 {
  fprintf(stderr,"mauvais type de transformation en entree de imx_inverser_transformation_rigide_3d_p\n");
  res=1; goto end_func;
 }

 transfo_tmp=cr_copy_transf3d(transfo_deb);
 if (!transfo_tmp)
 {
  fprintf(stderr,"Erreur d'allocation memoire dans imx_inverser_transformation_rigidezoom_3d_p\n");
  res=2; goto end_func;
 }

 //on se met en affine
 transfo_tmp->typetrans=AFFINE3D;
 transfo_tmp->nb_param=15;
 tmp_param=CALLOC(15, double);
 if (!tmp_param)
 {
  fprintf(stderr,"Erreur d'allocation memoire a echoue dans imx_inverser_transformation_rigidezoom_3d_p\n");
  res=2; goto end_func;
 }
 for (i=0;i<12;i++) tmp_param[i]=transfo_tmp->param[i];
 rigidz_to_affine_3d(tmp_param);
 FREE(transfo_tmp->param); transfo_tmp->param=tmp_param;

 //inversion de la transformation
 res = imx_inverser_transformation_affine_3d_p(transfo_tmp, transfo_res);
 
end_func:

 if (transfo_tmp) free_transf3d(transfo_tmp);

 return res;  
}

/*******************************************************************************
**        imx_query_interpolation_type_3D
*/
/*!
**       \brief choix de l'interpolation
**
**      ouvre une boite de dialogue pour le choix du type de l'interpolation
**      \retval int : le numero de l'interpolation
**          - nearest -> 0
**          - linear -> 1
**          - sin card -> 2
**          - quick sin card2 -> 3
**          - etc ...
**

*******************************************************************************/
int imx_query_interpolation_type_3D(int dist_type)
{
  int inter_type=0;
  /*on propose tous les choix mais on precise que ce choix n'aura pas
  d'influence si on a choisi une distance se basant sur l'interpolation volume partiel*/
  char *query[100];

    query[0]="nearest";
    query[1]="linear";
    query[2]="sin card.";
    query[3]="quick sin card2";
    query[4]="quick sin card3";
    query[5]="bspline 2";
    query[6]="bspline 3";
    query[7]="bspline 4";
    query[8]="bspline 5";
    query[9]="label";
    query[10]="\0";
    query[11]=NULL;

  /*on affiche quand meme les possibilites dans le cas d'une distance se basant sur l'interpolation volume partiel
  (cf fonctionnement des automates)
  mais on signale que le choix n'est pas pris en compte*/
  if (((dist_type>=11)&&(dist_type<=15))||(dist_type==19))
  {
    GETV_QCM("!!Interpolation non prise en compte!!",(char **)query);
    inter_type=10;
  }
  else
   inter_type=GETV_QCM("Interpolation",(char **)query);

   return inter_type;
}

/*------------------------------------------------*/
/*!
**  \brief determine la fonction d'interpolation
**
**  \param  inter_type : le type d'interpolation choisie
**  \retval InterpolationFct : un pointeur sur la fonction d'interpolatio
*/
/*------------------------------------------------*/
InterpolationFct
imx_choose_interpolation_fct(int inter_type)
{
    InterpolationFct interpol;
    switch (inter_type)
    {
        case 0: interpol=inter_nearest_3d;aff_log("NEAREST ");break;
        case 1: interpol=inter_linear_3d;aff_log("LINEAR ");break;
        case 2: interpol=inter_sinc_3d;aff_log("SINC ");break;
        case 3: interpol=inter_qsinc2_3d;aff_log("QSINC2 ");break;
        case 4: interpol=inter_qsinc3_3d;aff_log("QSINC3 ");break;
        case 5: interpol=inter_Bspline_3d_2;aff_log("BSPLINE2 ");break;
        case 6: interpol=inter_Bspline_3d_3;aff_log("BSPLINE3 ");break;
        case 7: interpol=inter_Bspline_3d_4;aff_log("BSPLINE4 ");break;
        case 8: interpol=inter_Bspline_3d_5;aff_log("BSPLINE5 ");break;
        case 9: interpol=inter_labeled_3d;aff_log("LABEL ");break;
    case 10: interpol=inter_VP_3d;aff_log("Volume Partiel ");break;


/*#ifdef HAS_IMLIB3D
    case 5: interpol=CPP_Imx_SplineResample_2;aff_log("BSPLINE2 ");break;
    case 6: interpol=CPP_Imx_SplineResample_3;aff_log("BSPLINE3 ");break;
    case 7: interpol=CPP_Imx_SplineResample_4;aff_log("BSPLINE4 ");break;
    case 8: interpol=CPP_Imx_SplineResample_5;aff_log("BSPLINE5 ");break;
#endif // HAS_IMLIB3D */
    default:PUT_WARN("Use by default of linear interpolation");
            interpol=inter_linear_3d;aff_log("LINEAR ");break;
    }
    return interpol;
}

/*******************************************************************************
**   apply_inverse_primitives_transf_3d()                                  
*/
/*!                                                                       
**      appliquer a une image une transformation contenu dans un fichier          
*******************************************************************************/

void apply_inverse_primitives_transf_3d(void)
{
  int im_deb, im_res, err=0;
  char nomfichier[PATH_LEN]/*, *quest[6]*/;


  /*question sur les images*/
  im_deb = GET_PLACE3D(TEXT0030);
  im_res = GET_PLACE3D(TEXT0006);

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

  
  imx_apply_inverse_primitives_transf_3d(im_res, im_deb, nomfichier);

  show_picture_3d(im_res);

}

/***********************************************************/
/*!  \ingroup Recalage
**
**
**  \brief Application d'une transformation definie dans un fichier de type *.trf
**          a une image 3d 
**
**
**  \param imdeb : image source
**  \param imres : image resultat (E/S)
**  \param nomfichier : le nom du fichier de la transformation (*.trf)
**  \param inter_type : choix de l'interpolation
**
***********************************************************/
int imx_apply_inverse_primitives_transf_3d(int im_res, int im_deb, char *nomfichier)
{
  transf3d *transfo;
  int err = 0;
  char tmp[PATH_LEN];
  grphic3d * imdeb = NULL;
  grphic3d * imres = NULL;
  
  imdeb = ptr_img_3d(im_deb);
  imres = ptr_img_3d(im_res);
    
  put_file_extension(nomfichier,".trf",tmp);
  
  /*chargement de la transformation*/
  transfo = load_transf_3d(tmp);
  if (! transfo)
  {
    err = 1;
    goto error1;
  }
  
  
  err = imx_apply_inverse_primitives_transf_3d_p(imres, imdeb, transfo);

  imx_reinitialise_visual_params_3d_p(imres);

  free_transf3d(transfo); 
error1:
  return err;
}

/***********************************************************/
/*!  \ingroup Recalage
**
**
**  \brief Application d'une transformation a l'aide d'une fonction d'interpolation
**          a une image 3d 
**
**
**  \param imdeb : image source
**  \param imres : image resultat (E/S)
**  \param tranfso : la transformation 
**  \param inter_fct : fct d'interpolation
**
***********************************************************/
int imx_apply_inverse_primitives_transf_3d_p(grphic3d *imres, grphic3d *imdeb, transf3d *transfo)
{
  field3d *champ=NULL;
  int err = 0,i,j,k,wdth,hght,dpth,wdthr,hghtr,dpthr,x,y,z;
  grphic3d *imTemp=NULL;

    wdth=imdeb->width;
    hght=imdeb->height;
    dpth=imdeb->depth;

    wdthr=imres->width;
    hghtr=imres->height;
    dpthr=imres->depth;

  // Creation de imtemp si necessaire
  if (imdeb==imres)
  {
   imTemp=cr_grphic3d(imdeb);
   if (!imTemp) { fprintf (stderr, "memory allocation failed in imx_apply_transf_3d_p\n"); return 1; }
   imx_copie_3d_p(imdeb, imTemp);
  }
  else imTemp=imdeb;

  //generation du champ
  if (transfo->dx && transfo->dy && transfo->dz)
  { // la transfo est en millimetre, mais on connait l'espace de depart (imdeb->dx,dy,dz) et l'espace d'arrive (transof->dx,dy,dz)
    err=resize_grphic3d_buffer(imres, transfo->width, transfo->height, transfo->depth);
    if (err) goto end_func;
    imres->dx=transfo->dx ; imres->dy=transfo->dy ; imres->dz=transfo->dz;
    imres->width=transfo->width ; imres->height=transfo->height ; imres->depth=transfo->depth ;
    champ = transf_to_field_3d(transfo,imres,imTemp);  
  }
  else
  { // on ne connait aucune info...
    err=resize_grphic3d_buffer(imres, imTemp->width, imTemp->height, imTemp->depth);
    if (err) goto end_func;
    imres->dx=imTemp->dx ; imres->dy=imTemp->dy ; imres->dz=imTemp->dz;
    imres->width=imTemp->width ; imres->height=imTemp->height ; imres->depth=imTemp->depth ;
    champ = transf_to_field_3d(transfo,NULL,NULL);  
  }

  if (! champ) { err = 1; goto end_func; }
 
    for (i=0;i<wdthr;i++)
    for (j=0;j<hghtr;j++)
    for (k=0;k<dpthr;k++)
        imres->mri[i][j][k]=0;
            
    
    for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
    for (k=0;k<dpth;k++)
        if (imdeb->mri[i][j][k]>0)
            {
            x=floor(i+champ->raw[i][j][k].x+0.5);
            y=floor(j+champ->raw[i][j][k].y+0.5);
            z=floor(k+champ->raw[i][j][k].z+0.5);
            
            if ((x>=0)&&(x<wdthr)&&(y>=0)&&(y<hghtr)&&(z>=0)&&(z<dpthr))
                imres->mri[x][y][z]=imdeb->mri[i][j][k];
            
            
            }
    
 
 
 
  // mise a jour de imres 
  imx_inimaxminpixel_3d_p(imres);
  
end_func:
  // liberation memoire 
  if (champ) free_field3d(champ);
  if (imTemp!=imdeb) free_grphic3d(imTemp);

  return err;
}


/*******************************************************************************
**   convert_trf_to_chp()                                  
*/
/*!                                                                       
**      Convertir une transfo en champ          
*******************************************************************************/

void convert_trf_to_chp(void)
{
int err=0;
char nomfichier[PATH_LEN],nomfichres[PATH_LEN],str[PATH_LEN];
transf3d* transfo,*transfores;
field3d* ch;
    
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
  
 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 put_file_extension(nomfichres,".trf",nomfichres);


  transfo=load_transf_3d(nomfichier);
  ch=transf_to_field_3d(transfo,NULL,NULL);
  transfores=field_to_transf_3d(ch,NULL,NULL);
  transfores->dx=transfo->dx;
  transfores->dy=transfo->dy;
  transfores->dz=transfo->dz;
         
  save_transf_3d(transfores,nomfichres);
  free_transf3d(transfo);
  free_transf3d(transfores);
  free_field3d(ch);
  
}



/*******************************************************************************
**   subsample_chp()                                  
*/
/*!                                                                       
**      sous echantillonne un champ          
*******************************************************************************/

void subsample_chp(void)
{
int err=0,zoom,puissance, i, j, k;
int wdth, hght,dpth;
char nomfichier[PATH_LEN],nomfichres[PATH_LEN],str[PATH_LEN];
transf3d* transfo,*transfores;
field3d* ch, *chres;
    
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
  
 /*nom du fichier resultat*/
 sprintf(str,"File ? in %s",_CurrentPath);
 strcpy(nomfichres,SAVE_FILE(str,NULL,&err));
 put_file_extension(nomfichres,".trf",nomfichres);

 puissance=GET_INT("puissance de 2 pour reduction ", 1, &err);

 zoom=pow(2.0,puissance);

  transfo=load_transf_3d(nomfichier);
  ch=transf_to_field_3d(transfo,NULL,NULL);
  
  wdth=transfo->width/zoom;
  hght=transfo->height/zoom;
  dpth=transfo->depth/zoom;
  
  chres=cr_field3d(wdth,  hght,  dpth);
  
   
  for (i=0; i<wdth; i++)
  for (j=0; j<hght; j++)
  for (k=0; k<dpth; k++)
    {
    chres->raw[i][j][k]=ch->raw[i*zoom][j*zoom][k*zoom];
    }
  
  transfores=field_to_transf_3d(chres,NULL,NULL);
  transfores->dx=transfo->dx*zoom;
  transfores->dy=transfo->dy*zoom;
  transfores->dz=transfo->dz*zoom;
         
  save_transf_3d(transfores,nomfichres);
  free_transf3d(transfo);
  free_transf3d(transfores);
  free_field3d(ch);
  free_field3d(chres);

}



/*******************************************************************************
**   warpDeformationField_menu()                                  
*/
/*!                                                                       
**      Déforme un champ de déformation          
*******************************************************************************/

void warpDeformationField_menu(void)
{
int err=0, inter_type;
int wdth, hght,dpth;
char nomfichierRef[PATH_LEN],nomfichierTrf[PATH_LEN], nomfichRes[PATH_LEN],str[PATH_LEN];
transf3d* transfoRef,*transfo, *transfoRes;
field3d* chRef,  *chRes;


    
/*fichier contenant la transformation à déformer */
{
    sprintf(nomfichierRef,"%s",GET_FILE("Champ a deformer (*.trf) ",&err));
    
    if (err)
        return ;
    
    put_file_extension(nomfichierRef,".trf",nomfichierRef);

    if (test_fichier_transf_3d(nomfichierRef) != 1)
    {
        PUT_ERROR("This file contains no 3D transformation");
        return ;
    }
}
    

/*fichier contenant la transformation*/
{
    sprintf(nomfichierTrf,"%s",GET_FILE("Transformation (*.trf)",&err));
    
    if (err)
    return ;
    
    put_file_extension(nomfichierTrf,".trf",nomfichierTrf);

    if (test_fichier_transf_3d(nomfichierTrf) != 1)
    {
    PUT_ERROR("This file contains no 3D transformation");
    return ;
    }
}

/*methode d'interpolation*/
inter_type = imx_query_interpolation_type_3D(0);

/*nom du fichier resultat*/
sprintf(str,"File ? in %s",_CurrentPath);
strcpy(nomfichRes,SAVE_FILE(str,NULL,&err));
put_file_extension(nomfichRes,".trf",nomfichRes);


transfoRef=load_transf_3d(nomfichierRef);
chRef=transf_to_field_3d(transfoRef,NULL,NULL);

transfo=load_transf_3d(nomfichierTrf);


wdth=transfo->width;
hght=transfo->height;
dpth=transfo->depth;

chRes=cr_field3d(wdth,  hght,  dpth);


warpDeformationField_p(chRef, transfo, chRes, inter_type);
                       
transfoRes=field_to_transf_3d(chRes,NULL,NULL);
transfoRes->dx=chRes->dx;
transfoRes->dy=chRes->dy;
transfoRes->dz=chRes->dz;

save_transf_3d(transfoRes,nomfichRes);
free_transf3d(transfo);
free_transf3d(transfoRef);
free_transf3d(transfoRes);
free_field3d(chRef);
free_field3d(chRes);
}

/*******************************************************************************
**   warpDeformationField_p()                                  
*/
/*!                                                                       
 **     Déforme un champ de déformation          
 *******************************************************************************/
void warpDeformationField_p(field3d* chRef, transf3d* transfo, field3d* chRes, int inter_type)
{
InterpolationFct inter_fct;
grphic3d* ux, *uy, *uz, *imtemp;
int wdth, hght, dpth;

wdth=chRes->width;
hght=chRes->height;
dpth=chRes->depth;

chRes->dx=transfo->dx;
chRes->dy=transfo->dy;
chRes->dz=transfo->dz;


inter_fct = imx_choose_interpolation_fct(inter_type);

imtemp=cr_grphic3d_modif(chRef->width,chRef->height,chRef->depth,0.0,1.0,0);

ux=cr_grphic3d_modif(wdth,hght,dpth,0.0,1.0,0);
uy=cr_grphic3d_modif(wdth,hght,dpth,0.0,1.0,0);
uz=cr_grphic3d_modif(wdth,hght,dpth,0.0,1.0,0);

composante_field_to_image_3d(chRef, imtemp, 0);
imx_apply_transf_3d_p(ux,imtemp, transfo, inter_fct);

composante_field_to_image_3d(chRef, imtemp, 1);
imx_apply_transf_3d_p(uy, imtemp, transfo, inter_fct);

composante_field_to_image_3d(chRef, imtemp, 2);
imx_apply_transf_3d_p(uz, imtemp, transfo, inter_fct);


images_to_field_3d(ux,uy,uz,chRes);

reorientDeformationField_p(chRes, transfo, chRes);

free_grphic3d(imtemp);
free_grphic3d(ux);
free_grphic3d(uy);
free_grphic3d(uz);
}

/*******************************************************************************
**   reorientDeformationField_p()                                  
*/
/*!                                                                       
 **     Réoriente un champ de déformation          
 *******************************************************************************/
void reorientDeformationField_p(field3d* chRef, transf3d* transfo, field3d* chRes)
{
int wdth, hght, dpth;
int i,j,k,l;
field3d* chTemp, *chTransfo;
double** J, **invJ,  *vec, *res, normvec, normres;



J= alloc_dmatrix(3,3);
vec=alloc_dvector(3);

wdth=chRef->width;
hght=chRef->height;
dpth=chRef->depth;

if (chRes==chRef)
    chTemp=cr_field3d(wdth,  hght,  dpth);
else
    chTemp=chRes;

chTransfo=transf_to_field_3d(transfo,NULL,NULL);


for (i=0; i<wdth; i++)
for (j=0; j<hght; j++)
for (k=0; k<dpth; k++)
    {
    
    if ((i==22)&&(j==45)&&(k==16))
        printf("coucou\n");
        
        
    eval_matrice_jacobienne_3d(chTransfo, i, j, k,  J);
    invJ=matrix_inversion(J,3);

    vec[0]=chRef->raw[i][j][k].x;
    vec[1]=chRef->raw[i][j][k].y;
    vec[2]=chRef->raw[i][j][k].z;
    
    normvec=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    

    res=matrix_multiply_vector(invJ,3,3,vec);
    normres=sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
    
    if (normres>0)
        for (l=0; l<3; l++) res[l]=res[l]*normvec/normres;
    else 
        for (l=0; l<3; l++) res[l]=0;
    
    chTemp->raw[i][j][k].x=res[0];
    chTemp->raw[i][j][k].y=res[1];
    chTemp->raw[i][j][k].z=res[2];

    free_dvector(res,3);
    free_dmatrix(invJ,3,3);
    }

if (chRes==chRef)
    {
    for (i=0; i<wdth; i++)
    for (j=0; j<hght; j++)
    for (k=0; k<dpth; k++)
        {
        chRes->raw[i][j][k].x=chTemp->raw[i][j][k].x;
        chRes->raw[i][j][k].y=chTemp->raw[i][j][k].y;
        chRes->raw[i][j][k].z=chTemp->raw[i][j][k].z;
        }
    free_field3d(chTemp);
    }

free_field3d(chTransfo);
free_dmatrix(J,3,3);
free_dvector(vec,3);
}



/*******************************************************************************
**     eval_matrice_jacobienne_3d(ch,i,j,k,mat)             
**                                                                  
**     jacobien d'un champ dans une image                             
*******************************************************************************/
void eval_matrice_jacobienne_3d(field3d *ch, int i, int j, int k, double** J)
{
double J1,J2,J3,J4,J5,J6,J7,J8,J9;
double filt[5];
int l;
vector3d ***data; 

data= ch->raw;
 
filt[1]=0.3326;filt[2]=0.0612;filt[3]=0.015;
J1=J5=J9=1.0;J2=J3=J4=J6=J7=J8=0.0;

if ((i>=3)&&(i<(ch->width-3))&&(j>=3)&&(j<(ch->height-3))&&(k>=3)&&(k<(ch->depth-3)))
     for (l=1;l<=3;l++)
     {
     J1+=filt[l]*(data[i+l][j][k].x-data[i-l][j][k].x);
     J2+=filt[l]*(data[i][j+l][k].x-data[i][j-l][k].x);
     J3+=filt[l]*(data[i][j][k+l].x-data[i][j][k-l].x);
     J4+=filt[l]*(data[i+l][j][k].y-data[i-l][j][k].y);
     J5+=filt[l]*(data[i][j+l][k].y-data[i][j-l][k].y);
     J6+=filt[l]*(data[i][j][k+l].y-data[i][j][k-l].y);
     J7+=filt[l]*(data[i+l][j][k].z-data[i-l][j][k].z);
     J8+=filt[l]*(data[i][j+l][k].z-data[i][j-l][k].z);
     J9+=filt[l]*(data[i][j][k+l].z-data[i][j][k-l].z); 
     }

J[0][0]=J1;
J[0][1]=J2;
J[0][2]=J3;
J[1][0]=J4;
J[1][1]=J5;
J[1][2]=J6;
J[2][0]=J7;
J[2][1]=J8;
J[2][2]=J9;

 }

