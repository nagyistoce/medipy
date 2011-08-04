/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h> 	
#include <config.h>
#include <time.h>

#include "recalage/chps_3d.h"
#include "recalage/topo/mtch_topo_3d.h"
#include "recalage/topo/divers_topo.h"
#include "math/imx_matrix.h"
#include "noyau/gui/imx_picture3d.h"
#include "math/imx_bspline.h"	/* pour les fonctions d'interpolation */
#include "noyau/mani_3d.h"
#ifdef __GTK_GUI
#include "noyau/gui/imx_3dmo.h"
#endif /* __GTK_GUI */
#include "noyau/io/imx_log.h"
#include "math/oper_3d.h"
#include "math/operation.h"
#include "traitement/trai_3d.h"
#include "math/ana_3d.h"
#include "noyau/io/imx_export_file.h"
#include "outils/imx_misc.h"
#include "recalage/topo/normalisation_topo.h"
#include "outils/imx_sort.h"
#include "noyau/imx_3d.h"
#include "noyau/imagix.h"
#include "noyau/gtkimagix.h" 
#include "combinaisons_sujetcpp.h"

#include "appel_combinaisons_sujet.h"

 
/******************************************************************************
** -- cout_combinaisons_sujet ---------------------------------------------------------------
** 
**   demande au user les entrées:
**  -k : on va étudier toutes les combinaison a k éléments (ou encore a k moments)
**  -dossier_sujet : le répertoire du sujet contenant les fichier d'erreur trf des 16 moments pour ce sujet la
**  -im_rec : c'est l'atlas de référence, chaque voxel nul de cet atlas ne sera pas inclu pour le calcul d erreur
**  -chemin_cout_sujet : le chemin du fichier qui contiendra le cout de cette combinaison
**
******************************************************************************/

void cout_combinaisons_sujet()
{
  int result,err,k,taille_img,im_rec;
  char* dossier_sujet =CALLOC(80,char);
  char* chemin_cout_sujet =CALLOC(80,char);

  k= GET_INT("combinaison de k? ", 4, &err);
  sprintf(dossier_sujet,"File ? in %s",_CurrentPath);
  strcpy(dossier_sujet,SAVE_FILE(dossier_sujet,NULL,&err));
  im_rec=GET_PLACE3D("Image entre 0 et 1 l image référence");///aprés ca sera l image référence
  sprintf(chemin_cout_sujet,"File ? in %s",_CurrentPath);
  strcpy(chemin_cout_sujet,SAVE_FILE(chemin_cout_sujet,NULL,&err));

  
  result = imx_cout_combinaisons_sujet(dossier_sujet,k,im_rec,chemin_cout_sujet);

}


 
/******************************************************************************
** -- add_fich_cout ---------------------------------------------------------------
** 
**   demande au user les entrées:
**  -chemin_cout1 : le fichier contient dans chaque ligne une combinaison de moment et le cout correspondant pour un 1er sujet donné
**  -chemin_cout2 : le fichier contient les meme combinaison de moment que chemin_cout1 mais les cout d'un 2eme sujet
**  -chemin_cout_add : contient les combinaison de moment existant dans chemin_cout1 et chemin_cout2 et fait la somme des deux cout
**
******************************************************************************/

void add_fich_cout()
{
  int err;
  char* chemin_cout1 =CALLOC(80,char);
  char* chemin_cout2 =CALLOC(80,char);
  char* chemin_cout_add =CALLOC(80,char);



  sprintf(chemin_cout1,"%s",GET_FILE("*.txt",&err));
  if (err)   return ;
  sprintf(chemin_cout2,"%s",GET_FILE("*.txt",&err));
  if (err)   return ;
  sprintf(chemin_cout_add,"File ? in %s",_CurrentPath);
  strcpy(chemin_cout_add,SAVE_FILE(chemin_cout_add,NULL,&err));

  imx_add_fich_cout(chemin_cout1,chemin_cout2,chemin_cout_add);
}

 
 
/******************************************************************************
** -- imx_add_fich_cout ---------------------------------------------------------------
** 
**  -chemin_cout1 : le fichier contient dans chaque ligne une combinaison de moment et le cout correspondant pour un 1er sujet donné
**  -chemin_cout2 : le fichier contient les meme combinaison de moment que chemin_cout1 mais les cout d'un 2eme sujet
**  -chemin_cout_add : contient les combinaison de moment existant dans chemin_cout1 et chemin_cout2 et fait la somme des deux cout
**
******************************************************************************/

void imx_add_fich_cout(char* chemin_cout1,char* chemin_cout2,char* chemin_cout_add)
{
  FILE * fic1, * fic2, * fic3;
  char moment1[80],moment2[80];

/****La somme dans un fichier temporaire*/
  fic1 = fopen(chemin_cout1, "r"); 	
  fic2 = fopen(chemin_cout2, "r");
  fic3 = fopen("tmpSC.txt", "w"); 

  while( fscanf(fic1,"%s",moment1) !=EOF )
  {
	fscanf(fic2,"%s",moment2);
	if(atof(moment1)==0 && atof(moment2)==0)	fprintf(fic3, "%s \t",moment1);			
	else fprintf(fic3, "%f \n",atof(moment1) + atof(moment2));   				
  }
fflush(fic3);
fclose(fic1);
fclose(fic2);
fclose(fic3);
/*La copie de La somme dans le fichier réel*/
rename("tmpSC.txt",chemin_cout_add);
}

 
/******************************************************************************
** -- apprentissage_erreur_moment-----------------------------------------------------
** 
**   demande au user les entrées:
**  -chemin_erreur_sujets : le fichier contenant la liste des fichier erreur(trf)(des différent sujets) sur lesquelle on fera l'apprentissage
**  -chemin_apprentissage_erreur : le fichier dans lequel sera mis les 9 statisristique d erreur pour chaque voxel
**  -im_ref : c'est l'atlas de référence, pour chaque voxel non nul de cet atlas sera calculé les 9 statistique
**
******************************************************************************/

void apprentissage_erreur_moment()
{
  int im_ref,err;
  char* chemin_erreur_sujets =CALLOC(80,char);
  char* chemin_apprentissage_erreur =CALLOC(80,char);

  im_ref=GET_PLACE3D("Image entre 0 et 1 l image  référence");///aprés ca sera l image référence

  sprintf(chemin_erreur_sujets,"%s",GET_FILE("*.txt",&err));
  if (err)   return ;

  sprintf(chemin_apprentissage_erreur,"File ? in %s",_CurrentPath);
  strcpy(chemin_apprentissage_erreur,SAVE_FILE(chemin_apprentissage_erreur,NULL,&err));
 
  imx_apprentissage_erreur_moment(im_ref,chemin_erreur_sujets,chemin_apprentissage_erreur);
}
 
/******************************************************************************
** -- imx_apprentissage_erreur_moment-----------------------------------------------------
** 
**  -chemin_erreur_sujets : le fichier contenant les chemins des fichier erreur(trf)(des différent sujets) sur lesquelle on fera l'apprentissage
**  -chemin_apprentissage_erreur : le fichier dans lequel sera mis les 9 statisristique d erreur pour chaque voxel
**  -im_ref : c'est l'atlas de référence, pour chaque voxel non nul de cet atlas sera calculé les 9 statistique
**
******************************************************************************/
void imx_apprentissage_erreur_moment(int im_ref,char* chemin_erreur_sujets,char* chemin_apprentissage_erreur)
{
  grphic3d *imref;

  FILE * fic1, * fic2, * fic3;
  char* chemin_erreur =CALLOC(80,char);
  field3d *champ;
  transf3d *transfo;
  vector3d ***data;
  TDimension wdth, hght,dpth;
  int i,j,k;
  int nbr_suj = 0;
  float Mx,My,Mz,Vx,Vy,Vz,Vxy,Vxz,Vyz;     		
int n;
/*********Initialisation du fichier d'apprentissage d'erreur chemin_apprentissage_erreur a ZERO pour tte les statistiques***********/
imref=ptr_img_3d(im_ref);
fic2 = fopen(chemin_apprentissage_erreur, "w");
  wdth=imref->width;
  hght=imref->height;
  dpth=imref->depth;
  for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
  for (k=0;k<dpth;k++)
  {
  	if(imref->mri[i][j][k] != 0)
  	{
	fprintf(fic2, "%f %f %f %f %f %f %f %f %f \n",0,0,0,0,0,0,0,0,0);	
	}
  }
  close(fic2);
/*********ouverture et parcours du fichier contenant la liste des fichiers d'erreur***********/
  fic1 = fopen(chemin_erreur_sujets, "r"); 	
  while( fscanf(fic1,"%s",chemin_erreur) !=EOF )
  {
	nbr_suj = nbr_suj+1;
	/*******Ouvrir le fichier contenant l'apprentissage d'erreur fic2 et un autre fichier temporaire fic3*******/
        fic2 = fopen(chemin_apprentissage_erreur, "r"); 
	fic3 = fopen("tmpAP.txt", "w"); 
	/*******Ouvrir le fichier trf dont le chemin est chemin_erreur du nouveau sujet *******/
	put_file_extension(chemin_erreur,".trf",chemin_erreur);
	if ((test_fichier_transf_3d(chemin_erreur))!=1)
	  {PUT_ERROR("This file contains no 3D transformation");return;}
	transfo=load_transf_3d(chemin_erreur);
	champ=transf_to_field_3d(transfo,NULL,NULL);

	data=champ->raw;
n = 0;
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
	{
		if(imref->mri[i][j][k] != 0)
		{

			/***Lecture des 9 statistique de fic2, leur mise a jour puis leur enregistrement dans fic3*****/
			//Moy de x, Moy de y,Moy de z, Var de x,Var de y,Var de z,Var de x y,Var de x z,Var de y z
			fscanf(fic2,"%f %f %f %f %f %f %f %f %f \n",&Mx,&My,&Mz,&Vx,&Vy,&Vz,&Vxy,&Vxz,&Vyz);

			fprintf(fic3,"%f %f %f %f %f %f %f %f %f \n",Mx+(data[i][j][k]).x,My+(data[i][j][k]).y,Mz+(data[i][j][k]).z,Vx+pow((data[i][j]	
			[k]).x,2),Vy+pow((data[i][j][k]).y,2),Vz+pow((data[i][j][k]).z,2),Vxy+ (data[i][j][k]).y * (data[i][j][k]).x, Vxz+ (data[i][j]
			[k]).z * (data[i][j][k]).x, Vyz+ (data[i][j][k]).y * (data[i][j][k]).z);
n=n+1;

		}
	}

	/*******fermer le fichier contenant l'apprentissage d'erreur et le fichier temporaire*******/
	fflush(fic3);
	close(fic2);
	close(fic3);
	/*La copie des nouvelles val statistique dans le fichier réel et suppression du fichier temporaire*/
	rename("tmpAP.txt",chemin_apprentissage_erreur);

  }
  
free_field3d(champ);
free_transf3d(transfo);  

fflush(stdout);
close(fic1);

  /**************Finir le calcul des statistiques**********/
  fic2 = fopen(chemin_apprentissage_erreur, "r"); 
  fic3 = fopen("tmpAP.txt", "w"); 

n=0;

  while( fscanf(fic2,"%f %f %f %f %f %f %f %f %f \n",&Mx,&My,&Mz,&Vx,&Vy,&Vz,&Vxy,&Vxz,&Vyz) !=EOF )
  {

n=n+1;

  fprintf(fic3, "%f %f %f %f %f %f %f %f %f \n",Mx/nbr_suj,My/nbr_suj,Mz/nbr_suj,(Vx/nbr_suj)-pow((Mx/nbr_suj),2),(Vy/nbr_suj)-pow((My/nbr_suj),2),(Vz/nbr_suj)-pow((Mz/nbr_suj),2), (Vxy/nbr_suj)-(Mx/nbr_suj)*(My/nbr_suj), (Vxz/nbr_suj)-(Mx/nbr_suj)*(Mz/nbr_suj), (Vyz/nbr_suj)-(Mz/nbr_suj)*(My/nbr_suj));
  }
  fflush(fic3);
  close(fic2);
  close(fic3);
  /*La copie des nouvelles val statistique dans le fichier réel*/ 
  rename("tmpAP.txt",chemin_apprentissage_erreur);


}


/******************************************************************************
** -- EMV_champs_moments ---------------------------------------------------------------
** 
**   demande au user les entrées:
**  -chemin_champ_covr : le fichier contient dans chaque ligne le recalage trf d un moment et le fichier des matrice de covariance de ce moment
**  -chemin_champ_combine : le fichier contient champs trf résultat de combinaison de la liste des moments retenu dans fich chemin_champ_covr
**  -im_ref : c'est l'atlas de référence, pour chaque voxel non nul de cet atlas sera calculé les 9 statistique
**
******************************************************************************/

void EMV_champs_moments()
{
  int err, im_ref;
  char* chemin_champ_covr =CALLOC(80,char);
  char* chemin_champ_combine =CALLOC(80,char);

  im_ref=GET_PLACE3D("Image entre 0 et 1 l image  référence");

  sprintf(chemin_champ_covr,"%s",GET_FILE("*.txt",&err));
  if (err)   return ;

  /*nom du fichier combiné resultat*/
  sprintf(chemin_champ_combine,"File ? in %s",_CurrentPath);
  strcpy(chemin_champ_combine,SAVE_FILE(chemin_champ_combine,NULL,&err));
  put_file_extension(chemin_champ_combine,".trf",chemin_champ_combine);

  imx_EMV_champs_moments(im_ref,chemin_champ_covr,chemin_champ_combine);

}
/******************************************************************************
** -- imx_EMV_champs_moments ---------------------------------------------------------------
** 
**  -chemin_champ_covr : le fichier contient dans chaque ligne le recalage trf d un moment et le fichier des matrice de covariance de ce moment
**  -chemin_champ_combine : la combinaison est au sens de maximum de vraisemblance
**  -im_ref : c'est l'atlas de référence, pour chaque voxel non nul de cet atlas sera calculé les 9 statistique
**
******************************************************************************/

void imx_EMV_champs_moments(int im_ref,char* chemin_champ_covr,char* chemin_champ_combine)
{
  grphic3d *imref;

  int err,i,j,k,indice;
  FILE * fic_cov, *fic_normalisation, *fic_tmp;					
  FILE * fic;
  transf3d *transf,*transfres;
  field3d *ch=NULL,*chres=NULL;
  TDimension wdth, hght,dpth;
  vector3d ***data,***datares;
  double dx,dy,dz;
  char chemin_champ[80],chemin_cov[80];
  double Mx,My,Mz,Vx,Vy,Vz,Vxy,Vxz,Vyz;
  double *ch_vox=NULL;						
  double *tmp_vox=NULL;
  double **cov_vox=NULL;
  double **inv_cov_vox;					

  imref=ptr_img_3d(im_ref);

/********************espace mémoire du champ résultat, son intitialisation a 0 et initialisation des somme de cov inverse a 0*************************/
  fic_normalisation = fopen("normalisation.txt", "w"); 

  wdth=imref->width;
  hght=imref->height;
  dpth=imref->depth;
  chres=cr_field3d(wdth,hght,dpth);

  for(i=0;i<wdth;i++)
  for(j=0;j<hght;j++)
  for(k=0;k<dpth;k++)
  {
  chres->raw[i][j][k].x=0;					///////////////ja i le droit de donner à comme valeur de champs pour les voxel nul 
  chres->raw[i][j][k].y=0;
  chres->raw[i][j][k].z=0;
  if(imref->mri[i][j][k] != 0) fprintf(fic_normalisation, "%f %f %f %f %f %f \n",0,0,0,0,0,0);
		
  }
  fflush(fic_normalisation);
  close(fic_normalisation);

/********************Parcour du fichier contenant les mesure (cad les champs trf des moment retenus et leur fichier de covariance)*************************/
  fic = fopen(chemin_champ_covr, "r"); 
  while( fscanf(fic,"%s %s \n",chemin_champ,chemin_cov) != EOF )
  {

	/***lecture du champs trf de la mesure et chargement de son fichier de covariance****/
	put_file_extension(chemin_champ,".trf",chemin_champ);
	if ((test_fichier_transf_3d(chemin_champ))!=1)
	  {PUT_ERROR("This file contains no 3D transformation");return;}
	transf=load_transf_3d(chemin_champ);
	dx = transf->dx;
	dy = transf->dy;
	dz = transf->dz;
	ch=transf_to_field_3d(transf,NULL,NULL);
	free_transf3d(transf);
        fic_cov = fopen(chemin_cov, "r"); 

        fic_normalisation = fopen("normalisation.txt", "r"); 
        fic_tmp = fopen("tmpNRM.txt", "w"); 
	/****Parcours de l image ipb et recherche des voxels non nuls***/
	for (i=0;i<wdth;i++)
	for (j=0;j<hght;j++)
	for (k=0;k<dpth;k++)
	{
		if(imref->mri[i][j][k] != 0)
		{
			/**** si voxel non null mettre matrice de covariance en mémoire, mise en mémoire du champs en ce voxel aussi*****/
			fscanf(fic_cov,"%lf %lf %lf %lf %lf %lf %lf %lf %lf \n",&Mx,&My,&Mz,&Vx,&Vy,&Vz,&Vxy,&Vxz,&Vyz);
			cov_vox=alloc_dmatrix(3,3);
			cov_vox[0][0] = Vx; cov_vox[1][1] = Vy;cov_vox[2][2] = Vz;         
			cov_vox[0][1] = Vxy;cov_vox[1][0] = Vxy;
			cov_vox[0][2] = Vxz;cov_vox[2][0] = Vxz;
			cov_vox[1][2] = Vyz;cov_vox[2][1] = Vyz;
			ch_vox=alloc_dvector(3);
			ch_vox[0] = ch->raw[i][j][k].x;ch_vox[1] = ch->raw[i][j][k].y;ch_vox[2] = ch->raw[i][j][k].z;
			/*****calcul d un nouveau champs a sommer : inverser la matrice cov, son produit avec le champ puis son ajout au champs réel*****/
			inv_cov_vox=matrix_inversion(cov_vox, 3);	
/***le tableau normalisation contient la somme des inverse des matrice de covariance en chaque voxel******/
fscanf(fic_normalisation,"%lf %lf %lf %lf %lf %lf \n",&Vx,&Vy,&Vz,&Vxy,&Vxz,&Vyz);
fprintf(fic_tmp,"",Vx+inv_cov_vox[0][0],Vy+inv_cov_vox[1][1],Vz+inv_cov_vox[2][2],Vxy+inv_cov_vox[0][1],Vxz+inv_cov_vox[0][2],Vyz+inv_cov_vox[1][2]);


			tmp_vox = matrix_multiply_vector(inv_cov_vox,3,3,ch_vox);
			chres->raw[i][j][k].x = chres->raw[i][j][k].x + tmp_vox[0];
			chres->raw[i][j][k].y = chres->raw[i][j][k].y + tmp_vox[1];
			chres->raw[i][j][k].z = chres->raw[i][j][k].z + tmp_vox[2];

			free_dmatrix(cov_vox,3,3);
			free_dmatrix(inv_cov_vox,3,3);
			free_dvector(ch_vox,3);
			free_dvector(tmp_vox,3);
		}
	}


  	fflush(fic_tmp);
	close(fic_tmp);
	close(fic_normalisation);
	rename("tmpNRM.txt","normalisation.txt");

  	free_field3d(ch);
	close(fic_cov);
  }
  close(fic);
/********************Parcour du champs combiné résultat et sa normalisation (division par l inverse de la somme des inverse des covariance****************/
  fic_normalisation = fopen("normalisation.txt", "r"); 

  for (i=0;i<wdth;i++)
  for (j=0;j<hght;j++)
  for (k=0;k<dpth;k++)
  {
	if(imref->mri[i][j][k] != 0)
	{
			fscanf(fic_normalisation,"%lf %lf %lf %lf %lf %lf \n",&Vx,&Vy,&Vz,&Vxy,&Vxz,&Vyz);
			/**** chargement du champs somme et la somme de matrice de covariance inverse ce voxel(dans tableau normalisation)**********/
			ch_vox=alloc_dvector(3);
			ch_vox[0] = chres->raw[i][j][k].x;ch_vox[1] = chres->raw[i][j][k].y;ch_vox[2] = chres->raw[i][j][k].z;
			cov_vox=alloc_dmatrix(3,3);
			cov_vox[0][0] = Vx; cov_vox[1][1] =  Vy;cov_vox[2][2] = Vz;         
			cov_vox[0][1] = Vxy;cov_vox[1][0] = Vxy;	
			cov_vox[0][2] = Vxz;cov_vox[2][0] = Vxz;
			cov_vox[1][2] = Vyz;cov_vox[2][1] = Vyz;

			/**** inverser la matrice (de la somme de matrice de covariance inverse) et son produit avec le champ somme en ce voxel********/
			inv_cov_vox=matrix_inversion(cov_vox, 3);
			tmp_vox = matrix_multiply_vector(inv_cov_vox,3,3,ch_vox);
			chres->raw[i][j][k].x = tmp_vox[0];
			chres->raw[i][j][k].y = tmp_vox[1];
			chres->raw[i][j][k].z = tmp_vox[2];
			free_dmatrix(cov_vox,3,3);
			free_dmatrix(inv_cov_vox,3,3);
			free_dvector(ch_vox,3);
			free_dvector(tmp_vox,3);

	}
  }

  close(fic_normalisation);
/********************Enregistrement du champs combiné résultat*************************/
  transfres=field_to_transf_3d(chres,NULL,NULL);
  transfres->dx=dx;					
  transfres->dy=dy;
  transfres->dz=dz;
  save_transf_3d(transfres,chemin_champ_combine);
  free_field3d(chres);
  free_transf3d(transfres);

  
}




/******************************************************************************
** -- calcul_statistiques-----------------------------------------------------
** 
**   demande au user les entrées:
**  -im_ref : c'est l'image modélisant le trf d'erreur pour lequel on calcul les stat
**  -chemin_statistiques : le fichier qui contiendra les statistiques 
**
******************************************************************************/

void calcul_statistiques()
{
  int im_ref,err;
  char* chemin_statistiques = CALLOC(80,char);

  im_ref=GET_PLACE3D("Image entre 0 et 1 l image  référence");


  sprintf(chemin_statistiques,"File ? in %s",_CurrentPath);
  strcpy(chemin_statistiques,SAVE_FILE(chemin_statistiques,NULL,&err));
 
  imx_calcul_statistiques(im_ref,chemin_statistiques);
}
 
/******************************************************************************
** -- imx_calcul_statistiques-----------------------------------------------------
** 
**   demande au user les entrées:
**  -im_ref : c'est l'atlas de référence, pour chaque voxel non nul de cet atlas sera calculé les 9 statistique
**  -chemin_erreur : le fichier contenant le fichier erreur(trf) d'un sujet
**  -chemin_statistiques : le fichier qui contiendra les statistiques 
**
******************************************************************************/

void imx_calcul_statistiques(int im_ref,char* chemin_statistiques)
{
  int 	i,j,k,Nbpas;
  int wdth,hght,dpth,inb;
  float moy,ecty;
  float rcoeff;
  double a1;
  float tmp;
  double std_tmp, moy_tmp;
  grphic3d *imref;
  FILE * fic;


  imref=ptr_img_3d(im_ref);
  wdth=imref->width;
  hght=imref->height;
  dpth=imref->depth;
  rcoeff=imref->rcoeff;



/*   Calcul du nombre de points, de la moyenne et de l'ecart type  */
  moy=0.0;
  ecty=0.0;
  std_tmp=0.0;
  moy_tmp=0.0;
  inb=0;
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
     for (k=0;k<dpth;k++)
      {

        a1=imref->mri[i][j][k]*rcoeff;
        if (a1 != 0.0)
	  {
	  moy_tmp+=a1;
	  std_tmp+=a1*a1;
	  inb=inb+1;
	  }
      }
   moy_tmp/=(double)(inb);
   std_tmp/=(double)(inb);
   std_tmp = std_tmp - moy_tmp*moy_tmp;
   std_tmp=(float)pow((double)std_tmp,(double) 0.5);

   moy = (float) moy_tmp;
   ecty=(float) std_tmp;


  fic = fopen(chemin_statistiques, "a"); 
  fprintf(fic,"%f %f \n",moy,ecty);
  close(fic);


}






