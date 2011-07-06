/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef medipy_components_visualization_3d_vtknormalsbyplanefitting_h
#define medipy_components_visualization_3d_vtknormalsbyplanefitting_h

#include <vtkObject.h>

class vtkFloatArray;
class vtkPointLocator;
class vtkPointSet;

class vtkNormalsByPlaneFitting : public vtkObject
{
public :
    vtkTypeRevisionMacro(vtkNormalsByPlaneFitting,vtkObject);
    static vtkFloatArray* from_plane(vtkPointSet* pointSet, vtkPointLocator* locator,
                                     float radius);

    static vtkFloatArray* by_projection(vtkPointSet* inner, vtkPointSet* outer,
                                        vtkPointLocator* innerLocator);

    static vtkFloatArray* choose_best_normal(vtkFloatArray* planeNormals,
                                             vtkFloatArray* projectionNormals,
                                             float angle);
    static vtkFloatArray* project_normals(vtkPointSet* source, vtkPointSet* destination, vtkPointLocator* sourceLocator);
//public :
//    vtkIdList* Wave;
//    vtkIdList* Wave2;
};

#endif // medipy_components_visualization_3d_vtknormalsbyplanefitting_h
