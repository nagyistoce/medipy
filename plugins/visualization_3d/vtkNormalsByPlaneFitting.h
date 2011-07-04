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
