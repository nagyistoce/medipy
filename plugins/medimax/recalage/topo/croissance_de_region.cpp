/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include <config.h>
#include "croissance_de_region.hpp"
#include "normalisation_topo.h"
#include "math/imx_matrix.h"


extern "C"
{
/*********************************************************************/
/*      						CRvoxel_compare_distance	 							         */
/*********************************************************************/

bool CRvoxel_compare_distance(CRvoxel& a, CRvoxel& b)
{
  return (a.distance < b.distance);
}

/*********************************************************************/
/*      		insertion_suppression_liste_trie	 							         */
/*********************************************************************/

int insertion_suppression_liste_trie(list<CRvoxel> &liste_voxel,CRvoxel &tmp)
{
std::list<CRvoxel>::iterator it = liste_voxel.begin();
int stop=1,res=0;


if (it->distance>tmp.distance)
	{
	stop=0;
	if ((it->i!=tmp.i)||(it->j!=tmp.j)||(it->k!=tmp.k))	
	 {liste_voxel.push_front(tmp);res=1;}
	}

while ((stop>0)&&(it != liste_voxel.end()))
{
	if ((it->i==tmp.i)&&(it->j==tmp.j)&&(it->k==tmp.k))
		{
		stop=0;
		}	
	else
	{
		if (it->distance>tmp.distance)
		{
		stop=0;
		it--;
		liste_voxel.insert(it,tmp);
		res=1;
		}
	}
it++;
}

if (stop>0)
	{liste_voxel.push_back(tmp);res=1;}


if ((stop==0)&&(res==1))
	while ((res>0)&&(it != liste_voxel.end()))
	{
	if ((it->i==tmp.i)&&(it->j==tmp.j)&&(it->k==tmp.k))
		{
		liste_voxel.erase(it);
		res=0;
		}	
		it++;
	}
	
return(res);
}

/*********************************************************************/
/*      						insertion_liste_trie	  								         */
/*********************************************************************/

int insertion_liste_trie(list<CRvoxel> &liste_voxel,CRvoxel &tmp)
{
std::list<CRvoxel>::iterator it = liste_voxel.begin();
int stop=1,res=0;


if (it->distance>tmp.distance)
	{
	stop=0;
	liste_voxel.push_front(tmp);
	}

while ((stop>0)&&(it != liste_voxel.end()))
{
		if (it->distance>tmp.distance)
		{
		stop=0;
		it--;
		liste_voxel.insert(it,tmp);
		}
	
it++;
}

if (stop>0)
	{liste_voxel.push_back(tmp);}

	
return(res);
}

/*********************************************************************/
/*      				raffine_segment_croissance_de_region	  		         */
/*********************************************************************/

void raffine_segment_croissance_de_region(grphic3d * im,grphic3d * im_cl,int nb_classe)
{
int i,j,k,wdth,hght,dpth,voisinage,label,indice,size,Ntot;
int wdth1,hght1,dpth1,*nb,***carte_label;
double  min,*moyenne,***distance;

std::list<CRvoxel>  liste_voxel;
std::list<CRvoxel>  liste_voxel_compl;
CRvoxel tmp;
CRvoxel tmp_compl;

tmp.distance=HUGE_VAL;    
insertion_liste_trie(liste_voxel_compl,tmp);
		
nb=(int *) malloc (nb_classe * sizeof(int));
moyenne=(double *) malloc (nb_classe * sizeof(double));


wdth=im->width;
hght=im->height;
dpth=im->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;


distance=alloc_dmatrix_3d(wdth,hght,dpth);
carte_label=alloc_imatrix_3d(wdth,hght,dpth);

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	distance[i][j][k]=HUGE_VAL;
	
/* initialisation des moyennes de chaque + du nombre d'element par classe */
for (i=0; i<nb_classe; i++)
	{
	nb[i]=0;
	moyenne[i]=0;
	}

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im_cl->mri[i][j][k]>0)
		{
		label=im_cl->mri[i][j][k]-1;
		moyenne[label]+=im->mri[i][j][k];
		nb[label]++;
		}

for (i=0; i<nb_classe; i++)
	if (nb[i]>0)
	moyenne[i]=(double)moyenne[i]/nb[i];
	else
	moyenne[i]=HUGE_VAL;

Ntot=0;
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im_cl->mri[i][j][k]==-1)
	Ntot++;

/* initialisation de la liste des voxels de la frontière */
for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im_cl->mri[i][j][k]==-1)
		{
		
		 
		voisinage=0;
			if (i>0) if (im_cl->mri[i-1][j][k]>0) voisinage++;
			if (j>0) if (im_cl->mri[i][j-1][k]>0) voisinage++;
			if (k>0) if (im_cl->mri[i][j][k-1]>0) voisinage++;
			if (i<wdth1) if (im_cl->mri[i+1][j][k]>0) voisinage++;
			if (j<hght1) if (im_cl->mri[i][j+1][k]>0) voisinage++;
			if (k<dpth1) if (im_cl->mri[i][j][k+1]>0) voisinage++;
			
				if (voisinage>=1)
						{
						min=HUGE_VAL;
						label=-1;
						if (i>0) 			if (im_cl->mri[i-1][j][k]>0) 	if (fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i-1][j][k]-1])<min) {min=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i-1][j][k]-1]); label=im_cl->mri[i-1][j][k];}
						if (j>0) 			if (im_cl->mri[i][j-1][k]>0) 	if (fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j-1][k]-1])<min) {min=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j-1][k]-1]); label=im_cl->mri[i][j-1][k];}
						if (k>0) 			if (im_cl->mri[i][j][k-1]>0) 	if (fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j][k-1]-1])<min) {min=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j][k-1]-1]); label=im_cl->mri[i][j][k-1];}
						if (i<wdth1) 	if (im_cl->mri[i+1][j][k]>0)	if (fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i+1][j][k]-1])<min) {min=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i+1][j][k]-1]); label=im_cl->mri[i+1][j][k];}
						if (j<hght1) 	if (im_cl->mri[i][j+1][k]>0)	if (fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j+1][k]-1])<min) {min=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j+1][k]-1]); label=im_cl->mri[i][j+1][k];}
						if (k<dpth1) 	if (im_cl->mri[i][j][k+1]>0)	if (fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j][k+1]-1])<min) {min=fabs(im->mri[i][j][k]-moyenne[im_cl->mri[i][j][k+1]-1]); label=im_cl->mri[i][j][k+1];}
						
						tmp.i=i;
						tmp.j=j;
						tmp.k=k;
						tmp.label=label;
						tmp.distance=min;
			
      			liste_voxel.push_back(tmp); 
						im_cl->mri[i][j][k]=-10; // c'est pour dire que ce point figure dans la liste 
						distance[i][j][k]=min;  
						carte_label[i][j][k]=label;
						}
		
		  //liste_voxel.merge(liste_tmp, CRvoxel_compare_distance);			
		}
	
 liste_voxel.sort(CRvoxel_compare_distance);	

size=liste_voxel.size();
while (size>0)
	{
	
	
	tmp=liste_voxel.front();
	tmp_compl=liste_voxel_compl.front();
	
	if (tmp.distance<tmp_compl.distance)
		{
		liste_voxel.pop_front();
		i=tmp.i;
		j=tmp.j;
		k=tmp.k;
		label=tmp.label;
		indice=label-1;
		}
	else
		{
		liste_voxel_compl.pop_front();
		i=tmp_compl.i;
		j=tmp_compl.j;
		k=tmp_compl.k;
		label=tmp_compl.label;
		indice=label-1;
		}
	
	if (im_cl->mri[i][j][k]==-10)
		{
		im_cl->mri[i][j][k]=label;
		size--;
		Ntot--;
		/* mise à jour de la moyenne et du nombre */
		moyenne[indice]=moyenne[indice]*nb[indice]+im->mri[i][j][k];
		nb[indice]++;
		moyenne[indice]=moyenne[indice]/nb[indice];
	
	
		/* on complete la liste */	
		if (i>0) 			if 	(im_cl->mri[i-1][j][k]==-10)	
			{tmp.i=i-1; tmp.j=j; tmp.k=k; tmp.distance=fabs(im->mri[i-1][j][k]-moyenne[im_cl->mri[i][j][k]-1]);	tmp.label=im_cl->mri[i][j][k];
			if ((tmp.distance<distance[tmp.i][tmp.j][tmp.k])&&(tmp.label!=carte_label[tmp.i][tmp.j][tmp.k]))	
				{insertion_suppression_liste_trie(liste_voxel_compl,tmp);carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
				distance[tmp.i][tmp.j][tmp.k]=MINI(tmp.distance,distance[tmp.i][tmp.j][tmp.k]);}

		if (j>0) 			if 	(im_cl->mri[i][j-1][k]==-10)	
			{tmp.i=i; tmp.j=j-1; tmp.k=k; tmp.distance=fabs(im->mri[i][j-1][k]-moyenne[im_cl->mri[i][j][k]-1]);	tmp.label=im_cl->mri[i][j][k];
				if ((tmp.distance<distance[tmp.i][tmp.j][tmp.k])&&(tmp.label!=carte_label[tmp.i][tmp.j][tmp.k])) 
				{insertion_suppression_liste_trie(liste_voxel_compl,tmp);carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
				distance[tmp.i][tmp.j][tmp.k]=MINI(tmp.distance,distance[tmp.i][tmp.j][tmp.k]);}
				
				
		if (k>0) 			if	(im_cl->mri[i][j][k-1]==-10)	
			{tmp.i=i; tmp.j=j; tmp.k=k-1; tmp.distance=fabs(im->mri[i][j][k-1]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];	
			if ((tmp.distance<distance[tmp.i][tmp.j][tmp.k])&&(tmp.label!=carte_label[tmp.i][tmp.j][tmp.k]))	
			{insertion_suppression_liste_trie(liste_voxel_compl,tmp);carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
			distance[tmp.i][tmp.j][tmp.k]=MINI(tmp.distance,distance[tmp.i][tmp.j][tmp.k]);}
			
			
		if (i<wdth1)	if	(im_cl->mri[i+1][j][k]==-10) 	
			{tmp.i=i+1; tmp.j=j; tmp.k=k; tmp.distance=fabs(im->mri[i+1][j][k]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];	
			if ((tmp.distance<distance[tmp.i][tmp.j][tmp.k])&&(tmp.label!=carte_label[tmp.i][tmp.j][tmp.k]))	
			{insertion_suppression_liste_trie(liste_voxel_compl,tmp);carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
			distance[tmp.i][tmp.j][tmp.k]=MINI(tmp.distance,distance[tmp.i][tmp.j][tmp.k]);}
			
			
		if (j<hght1) 	if	(im_cl->mri[i][j+1][k]==-10)	
			{tmp.i=i; tmp.j=j+1; tmp.k=k; tmp.distance=fabs(im->mri[i][j+1][k]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];			
				if ((tmp.distance<distance[tmp.i][tmp.j][tmp.k])&&(tmp.label!=carte_label[tmp.i][tmp.j][tmp.k]))	
				{insertion_suppression_liste_trie(liste_voxel_compl,tmp);carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
				distance[tmp.i][tmp.j][tmp.k]=MINI(tmp.distance,distance[tmp.i][tmp.j][tmp.k]);}
				
		if (k<dpth1) 	if	(im_cl->mri[i][j][k+1]==-10)	
		{tmp.i=i; tmp.j=j; tmp.k=k+1; tmp.distance=fabs(im->mri[i][j][k+1]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
			if ((tmp.distance<distance[tmp.i][tmp.j][tmp.k])&&(tmp.label!=carte_label[tmp.i][tmp.j][tmp.k]))	
			{insertion_suppression_liste_trie(liste_voxel_compl,tmp);carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
			distance[tmp.i][tmp.j][tmp.k]=MINI(tmp.distance,distance[tmp.i][tmp.j][tmp.k]);}
		
		
		
		
		if (i>0) 			if 	(im_cl->mri[i-1][j][k]==-1)	
			{tmp.i=i-1; tmp.j=j; tmp.k=k;  tmp.distance=fabs(im->mri[i-1][j][k]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
				insertion_liste_trie(liste_voxel,tmp);size++; distance[tmp.i][tmp.j][tmp.k]=tmp.distance;
				im_cl->mri[tmp.i][tmp.j][tmp.k]=-10; carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
				
		if (j>0) 			if 	(im_cl->mri[i][j-1][k]==-1)
			{tmp.i=i; tmp.j=j-1; tmp.k=k;  tmp.distance=fabs(im->mri[i][j-1][k]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
				insertion_liste_trie(liste_voxel,tmp);size++;distance[tmp.i][tmp.j][tmp.k]=tmp.distance;
				im_cl->mri[tmp.i][tmp.j][tmp.k]=-10; carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
				
		if (k>0) 			if	(im_cl->mri[i][j][k-1]==-1)	
			{tmp.i=i; tmp.j=j; tmp.k=k-1;  tmp.distance=fabs(im->mri[i][j][k-1]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
			insertion_liste_trie(liste_voxel,tmp);size++;distance[tmp.i][tmp.j][tmp.k]=tmp.distance;
			im_cl->mri[tmp.i][tmp.j][tmp.k]=-10; carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
			
		if (i<wdth1)	if	(im_cl->mri[i+1][j][k]==-1) 
			{tmp.i=i+1; tmp.j=j; tmp.k=k;  tmp.distance=fabs(im->mri[i+1][j][k]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
				insertion_liste_trie(liste_voxel,tmp);size++;distance[tmp.i][tmp.j][tmp.k]=tmp.distance;
				im_cl->mri[tmp.i][tmp.j][tmp.k]=-10;carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
				
		if (j<hght1) 	if	(im_cl->mri[i][j+1][k]==-1)	
			{tmp.i=i; tmp.j=j+1; tmp.k=k;  tmp.distance=fabs(im->mri[i][j+1][k]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
				insertion_liste_trie(liste_voxel,tmp);size++;distance[tmp.i][tmp.j][tmp.k]=tmp.distance;
				im_cl->mri[tmp.i][tmp.j][tmp.k]=-10;carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}
		
		if (k<dpth1) 	if	(im_cl->mri[i][j][k+1]==-1)	
			{tmp.i=i; tmp.j=j; tmp.k=k+1;  tmp.distance=fabs(im->mri[i][j][k+1]-moyenne[im_cl->mri[i][j][k]-1]); tmp.label=im_cl->mri[i][j][k];
				insertion_liste_trie(liste_voxel,tmp);size++;distance[tmp.i][tmp.j][tmp.k]=tmp.distance;
				im_cl->mri[tmp.i][tmp.j][tmp.k]=-10;carte_label[tmp.i][tmp.j][tmp.k]=tmp.label;}

		}
		printf("nombre de voxel dans la liste %d   Nb voxel restant à traiter %d  \n",size, Ntot);
	}

free (nb);
free(moyenne);
free_dmatrix_3d(distance);
free_imatrix_3d(carte_label);
}

/*********************************************************************/
/*      							swendsen_wang	  											         */
/*********************************************************************/
void swendsen_wang(grphic3d * im,double phi,int nb_classe)
{
int i,j,k,i0,j0,k0,wdth,hght,dpth,wdth1,hght1,dpth1,label,new_label,longueur;
double tirage;
grphic3d *im_visite;
std::list<CRvoxel>  liste_voxel;
CRvoxel tmp;

wdth=im->width;
hght=im->height;
dpth=im->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;

im_visite = cr_grphic3d_modif(wdth,hght,dpth,0.0,1,1);

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	im_visite->mri[i][j][k]=0;

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
	if (im_visite->mri[i][j][k]==0)
		{
		//printf("%d %d %d\n",i,j,k);
		label=im->mri[i][j][k];
		
		tirage=drand48();
		new_label=(int)floor(tirage*nb_classe);
		
		im->mri[i][j][k]=new_label;
		im_visite->mri[i][j][k]=1;
		
		tmp.i=i; tmp.j=j; tmp.k=k;
		liste_voxel.push_back(tmp);
		longueur=1;
		while (longueur>0)
			{
			//printf("taille de la liste : %d \n",longueur);
			tmp=liste_voxel.front();
			liste_voxel.pop_front();
			longueur--;
			
			/* on cherche les voisins et on empile dans la liste */
			i0=tmp.i;j0=tmp.j;k0=tmp.k;
			
			if (i0>0) if ((im_visite->mri[i0-1][j0][k0]==0)&&(im->mri[i0-1][j0][k0]==label))
					{tirage=drand48(); 
						if (tirage<(1-exp(-phi))) 
							{im->mri[i0-1][j0][k0]=new_label;im_visite->mri[i0-1][j0][k0]=1;
							tmp.i=i0-1;tmp.j=j0; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;}}
			
			if (j0>0) if ((im_visite->mri[i0][j0-1][k0]==0)&&(im->mri[i0][j0-1][k0]==label))
					{tirage=drand48(); 
						if (tirage<(1-exp(-phi))) 
							{im->mri[i0][j0-1][k0]=new_label;im_visite->mri[i0][j0-1][k0]=1;
							tmp.i=i0;tmp.j=j0-1; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;}}
			
			if (k0>0) if ((im_visite->mri[i0][j0][k0-1]==0)&&(im->mri[i0][j0][k0-1]==label))
					{tirage=drand48(); 
						if (tirage<(1-exp(-phi))) 
							{im->mri[i0][j0][k0-1]=new_label;im_visite->mri[i0][j0][k0-1]=1;
							tmp.i=i0;tmp.j=j0; tmp.k=k0-1; liste_voxel.push_back(tmp);longueur++;}}
			
			if (i0<wdth1) if ((im_visite->mri[i0+1][j0][k0]==0)&&(im->mri[i0+1][j0][k0]==label))
					{tirage=drand48(); 
						if (tirage<(1-exp(-phi))) 
							{im->mri[i0+1][j0][k0]=new_label;im_visite->mri[i0+1][j0][k0]=1;
							tmp.i=i0+1;tmp.j=j0; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;}}
			
			if (j0<hght1) if ((im_visite->mri[i0][j0+1][k0]==0)&&(im->mri[i0][j0+1][k0]==label))
					{tirage=drand48(); 
						if (tirage<(1-exp(-phi))) 
							{im->mri[i0][j0+1][k0]=new_label;im_visite->mri[i0][j0+1][k0]=1;
							tmp.i=i0;tmp.j=j0+1; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;}}
			
			if (k0<dpth1) if ((im_visite->mri[i0][j0][k0+1]==0)&&(im->mri[i0][j0][k0+1]==label))
					{tirage=drand48(); 
						if (tirage<(1-exp(-phi))) 
							{im->mri[i0][j0][k0+1]=new_label;im_visite->mri[i0][j0][k0+1]=1;
							tmp.i=i0;tmp.j=j0; tmp.k=k0+1; liste_voxel.push_back(tmp);longueur++;}}
							
					
					
			}
		
		}


free_grphic3d(im_visite);
}


/*********************************************************************/
/*        		 apply_norm_mean_std_classe_proba_reg				           */
/*********************************************************************/
/*void apply_norm_mean_std_classe_proba_reg(grphic3d *im1,grphic3d *im1_cl,grphic3d *imres,double *moyenne1,double *moyenne2,double *std1,double *std2,double *covariance, Param_Gaussienne_2d *param, int nb_classe)
{
int i,j,k,i0,j0,k0,wdth,hght,dpth,wdth1,hght1,dpth1,label,longueur,longueur_prec,nb;
double tirage;
grphic3d *im_visite;
std::list<CRvoxel>  liste_voxel;
CRvoxel tmp;

wdth=im1->width;
hght=im1->height;
dpth=im1->depth;

wdth1=wdth-1;
hght1=hght-1;
dpth1=dpth-1;

im_visite = cr_grphic3d_modif(wdth,hght,dpth,0.0,1,1);

for (i=0;i<wdth;i++)
for (j=0;j<hght;j++)
for (k=0;k<dpth;k++)
		if ((im1->mri[i][j][k]>0)&&(im_visite->mri[i][j][k]==0))
		{
		label=im1_cl->mri[i][j][k];
		im_visite->mri[i][j][k]=1;
		
		tmp.i=i; tmp.j=j; tmp.k=k;
		liste_voxel.push_back(tmp);
		longueur=1;longueur_prec=0;
		while (longueur>longueur_prec)
			{
			longueur_prec=longueur;
			tmp=liste_voxel.back();
			liste_voxel.pop_back();
			liste_voxel.push_front(tmp);
	
			/// on cherche les voisins et on empile dans la liste 
			i0=tmp.i;j0=tmp.j;k0=tmp.k;
			
			if (i0>0) if ((im_visite->mri[i0-1][j0][k0]==0)&&(im1_cl->mri[i0-1][j0][k0]==label))
					{tmp.i=i0-1;tmp.j=j0; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;im_visite->mri[tmp.i][tmp.j][tmp.k]=1;}
			
			if (j0>0) if ((im_visite->mri[i0][j0-1][k0]==0)&&(im1_cl->mri[i0][j0-1][k0]==label))
					{tmp.i=i0;tmp.j=j0-1; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;im_visite->mri[tmp.i][tmp.j][tmp.k]=1;}
			
			if (k0>0) if ((im_visite->mri[i0][j0][k0-1]==0)&&(im1_cl->mri[i0][j0][k0-1]==label))
					{tmp.i=i0;tmp.j=j0; tmp.k=k0-1; liste_voxel.push_back(tmp);longueur++;im_visite->mri[tmp.i][tmp.j][tmp.k]=1;}
			
			if (i0<wdth1) if ((im_visite->mri[i0+1][j0][k0]==0)&&(im1_cl->mri[i0+1][j0][k0]==label))
					{tmp.i=i0+1;tmp.j=j0; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;im_visite->mri[tmp.i][tmp.j][tmp.k]=1;}
			
			if (j0<hght1) if ((im_visite->mri[i0][j0+1][k0]==0)&&(im1_cl->mri[i0][j0+1][k0]==label))
					{tmp.i=i0;tmp.j=j0+1; tmp.k=k0; liste_voxel.push_back(tmp);longueur++;im_visite->mri[tmp.i][tmp.j][tmp.k]=1;}
			
			if (k0<dpth1) if ((im_visite->mri[i0][j0][k0+1]==0)&&(im1_cl->mri[i0][j0][k0+1]==label))
					{tmp.i=i0;tmp.j=j0; tmp.k=k0+1; liste_voxel.push_back(tmp);longueur++;im_visite->mri[tmp.i][tmp.j][tmp.k]=1;}
			
			}
		
		nb=0;
		/// on compte le pourcentage de voxel bien recale 
		 for (std::list<CRvoxel>::iterator it = liste_voxel.begin(); it != liste_voxel.end(); ++it)
    {
      i0=it->i;j0=it->j;k0=it->k;
			
			if (im1_cl->mask->mri[i0][j0][k0]<0)
				nb++;
    }
			
		tirage=drand48();
		
		if (tirage<1.0*nb/longueur)
			{
			/// on tire un label au sort 
			label=(int)floor(drand48()*nb_classe)+1;
			 for (std::list<CRvoxel>::iterator it = liste_voxel.begin(); it != liste_voxel.end(); ++it)
    		{
      		i0=it->i;j0=it->j;k0=it->k;
					im1_cl->mri[i0][j0][k0]=label;
			
    		}
			}	
		}
		
apply_norm_mean_std_classe_reg(im1,im1_cl,imres,moyenne1,moyenne2,std1,std2,covariance,param, nb_classe);

free_grphic3d(im_visite);

}
*/
}
