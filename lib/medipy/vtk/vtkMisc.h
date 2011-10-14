#ifndef medipy_vtk_addons_misc_h
#define medipy_vtk_addons_misc_h

#include <vtkObject.h>

class vtkCellLocator;

class vtkMisc : public vtkObject
{
public :
    vtkTypeRevisionMacro(vtkMisc,vtkObject);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Function not wrapped in Python
    static int intersect_with_line(vtkCellLocator* locator,
                                    double a0[3], double a1[3], double tol,
                                    double t[1], double x[3], double pcoords[3],
                                    int subId[1]);

    // Function not wrapped in Python
    static int intersect_with_line(vtkCellLocator* locator,
                                    double a0[3], double a1[3], double tol,
                                    double t[1], double x[3], double pcoords[3],
                                    int subId[1], vtkIdType cellId[1]);
};

#endif // medipy_vtk_addons_misc_h
