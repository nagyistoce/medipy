/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "vtkMisc.h"

#include <vtkCellLocator.h>


int
vtkMisc
::intersect_with_line(vtkCellLocator* locator,
                       double a0[3], double a1[3], double tol,
                       double t[1], double x[3], double pcoords[3],
                       int subId[1])
{
    return locator->IntersectWithLine(a0, a1, tol, t[0], x, pcoords, subId[0]);
}

int
vtkMisc
::intersect_with_line(vtkCellLocator* locator,
                       double a0[3], double a1[3], double tol,
                       double t[1], double x[3], double pcoords[3],
                       int subId[1], vtkIdType cellId[1])
{
    return locator->IntersectWithLine(a0, a1, tol, t[0], x, pcoords, subId[0],
                                      cellId[0]);
}
