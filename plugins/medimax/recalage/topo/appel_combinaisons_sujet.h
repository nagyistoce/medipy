/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef appel_combinaisons_sujet_H
#define appel_combinaisons_sujet_H

void cout_combinaisons_sujet();
void add_fich_cout();
void imx_add_fich_cout(char* chemin_cout1,char* chemin_cout2,char* chemin_cout_add);
void apprentissage_erreur_moment();
void imx_apprentissage_erreur_moment(int im_ref,char* chemin_erreur_sujets,char* chemin_apprentissage_erreur);
void EMV_champs_moments();
void imx_EMV_champs_moments(int im_ref,char* chemin_champ_covr,char* chemin_champ_combine);
void calcul_statistiques();
void imx_calcul_statistiques(int im_ref, char* chemin_statistiques);
#endif
