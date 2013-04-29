/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef CROISSANCE_DE_REGION_H
#define CROISSANCE_DE_REGION_H


#include <iostream>

#include <algorithm>
#include <vector>
#include <list>



using namespace std;


extern "C"
{

#include "noyau/imagix.h"
#include "noyau/gtkimagix.h"
#include "noyau/imx_3d.h"
}

struct CRvoxel
{
  int i,j,k;
  int label; /* codé entre 1 et nb_classe */
  double distance;
};

//void raffine_segment_croissance_de_region(grphic3d *im,grphic3d *im_cl,int nb_classe);


#endif
