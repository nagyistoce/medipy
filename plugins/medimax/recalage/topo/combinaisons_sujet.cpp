/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

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







#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h> 	
#include "noyau/imagix.h"
#include<map>
#include<list>
#include<iostream>
#include<vector>
#include<set>
#include<algorithm>

using namespace std;
		    
class Combin
{
public :
 long int n ;
 long int k ;
 long int* data;


 ~Combin()
 	{delete data;}
 Combin(long n, long k)
	{
	this->n = n;
	this->k = k;
	this->data = new long[k];
	for (long i = 0; i < k; ++i)
		this->data[i] = i;
	} 

	bool Successor();


}; // Combination class

bool Combin::Successor()
{
if (data[0] == n - k)
	return false;

long* ansdata = new long[k];

long i;
for (i = 0; i < this->k; ++i)
	ansdata[i] = this->data[i];

for (i = this->k - 1; i > 0 && ansdata[i] == this->n - this->k + i; --i)
	;

++ansdata[i];

for (long j = i; j < this->k - 1; ++j)
ansdata[j+1] = ansdata[j] + 1;

for (i = 0; i < this->k; ++i)
	this->data[i] = ansdata[i];
delete ansdata;
return true;
}


//Extern C car ces fonctions sont utilis�es � partir d un code C
extern "C"
{ 
 
/******************************************************************************
** -- imx_cout_combinaisons_sujet ---------------------------------------------------------------
** 
**   demande au user les entr�es:
**  -k : on va �tudier toutes les combinaison a k �l�ments (ou encore a k moments)
**  -dossier_sujet : le r�pertoire du sujet contenant les 16 fichier d'erreur trf des 16 moments pour ce sujet la
**  -im_rec : c'est l'atlas de r�f�rence, chaque voxel nul de cet atlas ne sera pas inclu pour le calcul d erreur
**  ----->Pour une combinaison de k moments donn�s : -pour chaque voxel non nul de l atlas je retient son erreur min sur les k erreur trf des k moments.
** 						     -Le cout de la combinaison est alors la somme sur l image des min retenu
**  -chemin_cout_sujet : le fichier mis dans le r�pertoire du sujet : il contiendra toute les combinaison possible d k moment et pour chaque **										   combinaison son cout
**
******************************************************************************/

int imx_cout_combinaisons_sujet(char* dossier_sujet,int k,int im_rec, char* chemin_cout_sujet)
{
int n=16;
char* chemin_err_trf=CALLOC(250,char);
FILE * fic;


int l,m,o,r;
Combin* c = new Combin(n,k);
int tutu  = 0;
field3d *champ;
transf3d *transfo;
vector3d ***data;
TDimension wdth, hght,dpth;

double cout_combinaison,cout;
double ***image_cout;

grphic3d *imrec;
imrec=ptr_img_3d(im_rec);
int taille=imrec->width;

char* moment[16];
moment[0]=CALLOC(50,char);sprintf(moment[0],"errN1_R3_%d.trf",taille);
moment[1]=CALLOC(50,char);sprintf(moment[1],"errN1_R5_%d.trf",taille);
moment[2]=CALLOC(50,char);sprintf(moment[2],"errN1_R7_%d.trf",taille);
moment[3]=CALLOC(50,char);sprintf(moment[3],"errN1_R9_%d.trf",taille);
moment[4]=CALLOC(50,char);sprintf(moment[4],"errN2_R3_%d.trf",taille);
moment[5]=CALLOC(50,char);sprintf(moment[5],"errN2_R5_%d.trf",taille);
moment[6]=CALLOC(50,char);sprintf(moment[6],"errN2_R7_%d.trf",taille);
moment[7]=CALLOC(50,char);sprintf(moment[7],"errN2_R9_%d.trf",taille);
moment[8]=CALLOC(50,char);sprintf(moment[8],"errN3_R3_%d.trf",taille);
moment[9]=CALLOC(50,char);sprintf(moment[9],"errN3_R5_%d.trf",taille);
moment[10]=CALLOC(50,char);sprintf(moment[10],"errN3_R7_%d.trf",taille);
moment[11]=CALLOC(50,char);sprintf(moment[11],"errN3_R9_%d.trf",taille);
moment[12]=CALLOC(50,char);sprintf(moment[12],"errN4_R3_%d.trf",taille);
moment[13]=CALLOC(50,char);sprintf(moment[13],"errN4_R5_%d.trf",taille);
moment[14]=CALLOC(50,char);sprintf(moment[14],"errN4_R7_%d.trf",taille);
moment[15]=CALLOC(50,char);sprintf(moment[15],"errN4_R9_%d.trf",taille);

printf("%s %d %d \n",dossier_sujet,k,im_rec);
/*************************Ouvrir le fichier dans lequel on va enregistrer les cout des combinaison Puis boucler sur les combinaisons*******************/
//sprintf(chemin_cout_sujet,"%scout_%d_%d.txt",dossier_sujet,taille,k); g�r� par l'automate
fic = fopen(chemin_cout_sujet, "w"); 						
///boucler sur les combinaisons
while (true)
	{
	/*************************rester sur la premi�re combinaison ou avancer a la suivante en cours*******************/
	bool a = true;
	if (tutu==0)	tutu=1;
	else 		a = c->Successor();

	/*************************Traitement d une combinaison en cours si elle existe*******************/	
	if (a==true)
		{

		/*************************Boucler sur les moments de cette combinaison*******************/
		for (r=0;r<k;r++)
		{
		/*************************ouvrir le fichier trf erreur trf de la composante num�ro r de la combinaison*****************/
		sprintf(chemin_err_trf,"%s%s",dossier_sujet,moment[c->data[r]]);		
		put_file_extension(chemin_err_trf,".trf",chemin_err_trf);
		if ((test_fichier_transf_3d(chemin_err_trf))!=1) 
		  {PUT_ERROR("This file contains no 3D transformation");return(0);}
		transfo=load_transf_3d(chemin_err_trf);
		champ=transf_to_field_3d(transfo,NULL,NULL);
		free_transf3d(transfo);
		wdth=champ->width;hght=champ->height;dpth=champ->depth;
		data=champ->raw;								
		if(r == 0)/***************initialisation de l image cout au premier champ trf*****/
		{

			image_cout = alloc_dmatrix_3d(wdth, hght, dpth);			
	 		for (l=0;l<wdth;l++)
	  		for (m=0;m<hght;m++)
	   		for (o=0;o<dpth;o++)
	   			{
		        	if (imrec->mri[l][m][o] != 0)
				image_cout[l][m][o] = pow((data[l][m][o]).y ,2)+pow((data[l][m][o]).x ,2)+pow((data[l][m][o]).z ,2);
	   			}
			free_field3d(champ);
		}
		else /***************comparaison de ce nouveeau trf correspondant a un nouveau moment de la combinaison a l image image_cout*****/
		{
	 		for (l=0;l<wdth;l++)
	  		for (m=0;m<hght;m++)
	   		for (o=0;o<dpth;o++)
	   			{
		        	if (imrec->mri[l][m][o] != 0)
				   {
					cout = pow((data[l][m][o]).y ,2)+pow((data[l][m][o]).x ,2)+pow((data[l][m][o]).z ,2);
					if(cout < image_cout[l][m][o]) image_cout[l][m][o] = cout;
				   }
	   			}
			free_field3d(champ);
		}
		}
		/*************************Parcourir l image de cout et faire la somme de cout******/
		/// Initialiser 
		cout_combinaison = 0;
		/// Parcourir les voxel et mettre a jour le cout
	 	for (l=0;l<wdth;l++)
	  	for (m=0;m<hght;m++)
	   	for (o=0;o<dpth;o++)
		{
        	if (imrec->mri[l][m][o] != 0)	cout_combinaison = cout_combinaison + image_cout[l][m][o];
		}

		free_dmatrix_3d(image_cout);
		/*************************Enregister le cout et la combibnaison correspondante dans le fichier******/	   			
		for (r=0;r<k;r++)
		fprintf(fic, "%s \t",moment[c->data[r]]);
		fprintf(fic, "%f \t \n",cout_combinaison);		
		}
	else
		{
		delete c;
		break;
		}

	}
/* Fermeture du fichier des cout des combinaisons pour ce sujet */
fflush(fic);
fclose(fic);

return 1;
}

}
