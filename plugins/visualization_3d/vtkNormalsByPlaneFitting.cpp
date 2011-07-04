#include "vtkNormalsByPlaneFitting.h"

#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkPriorityQueue.h>
#include <vtkSmartPointer.h>

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector_fixed_ref.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

vnl_vector_fixed<float, 3>
computeCoG(vtkPoints* points, vtkIdList* ids)
{
    float const * const dataArray = reinterpret_cast<float*>(points->GetData()->GetVoidPointer(0));

    vnl_vector_fixed<float, 3> cog(0.);
    for(int i=0; i<ids->GetNumberOfIds(); ++i)
    {
        float const * const p = dataArray+3*ids->GetId(i);
        cog[0] += p[0];
        cog[1] += p[1];
        cog[2] += p[2];
    }
    cog /= points->GetNumberOfPoints();

    return cog;
}

float meanDistance(vtkPoints* points, vtkIdList* ids, vnl_vector_fixed<float, 3> const & cog)
{
    float const * const dataArray = reinterpret_cast<float*>(points->GetData()->GetVoidPointer(0));

    float distance=0;

    for(int i=0; i<ids->GetNumberOfIds(); ++i)
    {
        float const * const p = dataArray+3*ids->GetId(i);
        float const x = p[0]-cog[0];
        float const y = p[1]-cog[1];
        float const z = p[2]-cog[2];
        distance += std::sqrt(x*x+y*y+z*z);
    }

    return distance;
}

void
normal_from_plane(vtkPoints* points, vtkPointLocator* locator, double const * center,
                  float radius, float* normal)
{
    // http://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
    // sections 3 and 4

    // Find all points within radius
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    locator->FindPointsWithinRadius(radius, center, ids);

    // FIXME
    if(ids->GetNumberOfIds() < 3)
    {
        std::cout << center[0] << " " << center[1] << " " << center[2] << " failed : "
                  << ids->GetNumberOfIds() << " points" << std::endl;
        return;
    }

    float const * const dataArray = reinterpret_cast<float*>(points->GetData()->GetVoidPointer(0));

    // Find center of gravity
    vnl_vector_fixed<float, 3> cog = computeCoG(points, ids);

    float const distance = meanDistance(points, ids, cog);
    // FIXME if distance<tolerance (1e-6)
    float const scale = 1./distance;
    vnl_vector_fixed<float, 3> const translation(-cog[0]*scale, -cog[1]*scale, -cog[2]*scale);

    vnl_matrix_fixed<float, 4, 4> m(0.);
    for(int i=0; i<ids->GetNumberOfIds(); ++i)
    {
        float const * const p = dataArray+3*ids->GetId(i);
        // Normalize point
        float const x = p[0]*scale+translation[0];
        float const y = p[1]*scale+translation[1];
        float const z = p[2]*scale+translation[2];

        m[0][3] += x;
        m[1][3] += y;
        m[2][3] += z;

        m[0][0] += x*x;
        m[1][1] += y*y;
        m[2][2] += z*z;

        m[0][1] += x*y;
        m[1][2] += y*z;
        m[0][2] += x*z;
    }
    m[3][3] = ids->GetNumberOfIds();

    m[1][0] = m[0][1];
    m[2][1] = m[1][2];
    m[2][0] = m[0][2];

    m[3][0] = m[0][3];
    m[3][1] = m[1][3];
    m[3][2] = m[2][3];

    vnl_svd<float> const svd(m);
    vnl_vector<float> x = svd.left_nullvector();

    // Re-transform the point back to the real world
    vnl_matrix_fixed <float, 4, 4> transposed_transformation(0.);
    transposed_transformation[0][0] = scale;
    transposed_transformation[1][1] = scale;
    transposed_transformation[2][2] = scale;
    transposed_transformation[3][3] = 1.;
    transposed_transformation[3][0] = translation[0];
    transposed_transformation[3][1] = translation[1];
    transposed_transformation[3][2] = translation[2];
    x = transposed_transformation*x;

    // Plane coefficients are (a,b,c,d), normal is (a,b,c)
    normal[0] = x[0];
    normal[1] = x[1];
    normal[2] = x[2];
    // Normalize the vector
    vnl_vector_fixed_ref<float, 3> n(normal);
    n /= n.magnitude();

}

vtkFloatArray*
vtkNormalsByPlaneFitting
::from_plane(vtkPointSet* pointSet, vtkPointLocator* locator, float radius)
{
    vtkFloatArray* normals = vtkFloatArray::New();
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(pointSet->GetNumberOfPoints());
    normals->FillComponent(0,0);
    normals->FillComponent(1,0);
    normals->FillComponent(2,0);

    pointSet->ComputeBounds();
    double* bounds = pointSet->GetBounds();
    vnl_vector_fixed<float, 3> center(0.5*(bounds[0]+bounds[1]),
                                      0.5*(bounds[2]+bounds[3]),
                                      0.5*(bounds[4]+bounds[5]));

    float const * const pointsArray = reinterpret_cast<float*>(pointSet->GetPoints()->GetData()->GetVoidPointer(0));
    for(unsigned int i=0; i<pointSet->GetNumberOfPoints(); ++i)
    {
        normal_from_plane(pointSet->GetPoints(), locator, pointSet->GetPoint(i),
                          radius, normals->GetPointer(3*i));
    }

    for(unsigned int i=0; i<pointSet->GetNumberOfPoints(); ++i)
    {
        vnl_vector_fixed_ref_const<float, 3> p(pointsArray+3*i);
        vnl_vector_fixed_ref<float, 3> n(normals->GetPointer(3*i));
        vnl_vector_fixed<float, 3> const d(p[0]-center[0], p[1]-center[1], p[2]-center[2]);

        float const dot = n[0]*d[0]+n[1]*d[1]+n[2]*d[2];
        float const angle=std::acos(dot/(n.magnitude()*d.magnitude()));

        if(angle > M_PI/2.)
        {
            n *= -1;
        }
    }

    return normals;

#if 0
    // Orient the normals, code taken from vtkPolyDataNormals
    {
        // No need to check this->Consistency. It's implied.

        // Ok, here's the basic idea: the "left-most" polygon should
        // have its outward pointing normal facing left. If it doesn't,
        // reverse the vertex order. Then use it as the seed for other
        // connected polys. To find left-most polygon, first find left-most
        // point, and examine neighboring polys and see which one
        // has a normal that's "most aligned" with the X-axis. This process
        // will need to be repeated to handle all connected components in
        // the mesh. Report bugs/issues to cvolpe@ara.com.
        int foundLeftmostCell;
        vtkIdType leftmostCellID = -1, currentPointID, currentCellID;
        vtkIdType *leftmostCells;
        unsigned short nleftmostCells;
        vtkIdType *cellPts;
        vtkIdType nCellPts;
        int cIdx;
        double bestNormalAbsXComponent;
        int bestReverseFlag;
        vtkPriorityQueue *leftmostPoints = vtkPriorityQueue::New();
        this->Wave = vtkIdList::New();
        this->Wave->Allocate(numPolys / 4 + 1, numPolys);
        this->Wave2 = vtkIdList::New();
        this->Wave2->Allocate(numPolys / 4 + 1, numPolys);

        // Put all the points in the priority queue, based on x coord
        // So that we can find leftmost point
        leftmostPoints->Allocate(numPts);
        for (ptId = 0; ptId < numPts; ptId++)
        {
            leftmostPoints->Insert(inPts->GetPoint(ptId)[0], ptId);
        }

        // Repeat this while loop as long as the queue is not empty,
        // because there may be multiple connected components, each of
        // which needs to be seeded independently with a correctly
        // oriented polygon.
        while (leftmostPoints->GetNumberOfItems())
        {
            foundLeftmostCell = 0;
            // Keep iterating through leftmost points and cells located at
            // those points until I've got a leftmost point with
            // unvisited cells attached and I've found the best cell
            // at that point
            do
            {
                currentPointID = leftmostPoints->Pop();
                this->OldMesh->GetPointCells(currentPointID, nleftmostCells,
                                             leftmostCells);
                bestNormalAbsXComponent = 0.0;
                bestReverseFlag = 0;
                for (cIdx = 0; cIdx < nleftmostCells; cIdx++)
                {
                    currentCellID = leftmostCells[cIdx];
                    if (this->Visited[currentCellID] == VTK_CELL_VISITED)
                    {
                        continue;
                    }
                    this->OldMesh->GetCellPoints(currentCellID, nCellPts,
                                                 cellPts);
                    vtkPolygon::ComputeNormal(inPts, nCellPts, cellPts, n);
                    // Ok, see if this leftmost cell candidate is the best
                    // so far
                    if (fabs(n[0]) > bestNormalAbsXComponent)
                    {
                        bestNormalAbsXComponent = fabs(n[0]);
                        leftmostCellID = currentCellID;
                        // If the current leftmost cell's normal is pointing to the
                        // right, then the vertex ordering is wrong
                        bestReverseFlag = (n[0] > 0);
                        foundLeftmostCell = 1;
                    } // if this normal is most x-aligned so far
                } // for each cell at current leftmost point
            }
            while (leftmostPoints->GetNumberOfItems() && !foundLeftmostCell);

            if (foundLeftmostCell)
            {
                // We've got the seed for a connected component! But do
                // we need to flip it first? We do, if it was pointed the wrong
                // way to begin with, or if the user requested flipping all
                // normals, but if both are true, then we leave it as it is.
                if (bestReverseFlag ^ this->FlipNormals)
                {
                    this->NewMesh->ReverseCell(leftmostCellID);
                    this->NumFlips++;
                }
                this->Wave->InsertNextId(leftmostCellID);
                this->Visited[leftmostCellID] = VTK_CELL_VISITED;
                this->TraverseAndOrder();
                this->Wave->Reset();
                this->Wave2->Reset();
            } // if found leftmost cell
        } // Still some points in the queue
        this->Wave->Delete();
        this->Wave2->Delete();
        vtkDebugMacro(<<"Reversed ordering of " << this->NumFlips << " polygons");
    } // automatically orient normals
#endif
}

vtkFloatArray*
vtkNormalsByPlaneFitting
::by_projection(vtkPointSet* inner, vtkPointSet* outer, vtkPointLocator* innerLocator)
{
    vtkFloatArray* normals = vtkFloatArray::New();
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(outer->GetNumberOfPoints());

    for(unsigned int i=0; i<outer->GetNumberOfPoints(); ++i)
    {
        double* outerPoint = outer->GetPoint(i);
        vtkIdType const innerId = innerLocator->FindClosestPoint(outerPoint);
        double* innerPoint = inner->GetPoint(innerId);
        vnl_vector_fixed<double, 3> n(outerPoint[0]-innerPoint[0],
                                      outerPoint[1]-innerPoint[1],
                                      outerPoint[2]-innerPoint[2]);
        n /= n.magnitude();
        normals->GetPointer(3*i)[0] = n[0];
        normals->GetPointer(3*i)[1] = n[1];
        normals->GetPointer(3*i)[2] = n[2];
    }

    return normals;
}

vtkFloatArray*
vtkNormalsByPlaneFitting
::choose_best_normal(vtkFloatArray* planeNormals,
                     vtkFloatArray* projectionNormals, float maxAngle)
{
    vtkFloatArray* normals = vtkFloatArray::New();
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(planeNormals->GetNumberOfTuples());

    for(unsigned int i=0; i<planeNormals->GetNumberOfTuples(); ++i)
    {
        vnl_vector_fixed_ref<float, 3> n1(planeNormals->GetPointer(3*i));
        vnl_vector_fixed_ref<float, 3> n2(projectionNormals->GetPointer(3*i));

        float const dot = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
        float const angle=std::acos(dot/(n1.magnitude()*n2.magnitude()));

        if(angle > maxAngle)
        {
            normals->GetPointer(3*i)[0] = n2[0];
            normals->GetPointer(3*i)[1] = n2[1];
            normals->GetPointer(3*i)[2] = n2[2];
        }
        else
        {
            normals->GetPointer(3*i)[0] = n1[0];
            normals->GetPointer(3*i)[1] = n1[1];
            normals->GetPointer(3*i)[2] = n1[2];
        }
    }

    return normals;
}

vtkFloatArray*
vtkNormalsByPlaneFitting
::project_normals(vtkPointSet* source, vtkPointSet* destination, vtkPointLocator* sourceLocator)
{
    vtkFloatArray* normals = vtkFloatArray::New();
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(destination->GetNumberOfPoints());

    float* sourceNormals = reinterpret_cast<float*>(
        source->GetPointData()->GetArray("probe_normals")->GetVoidPointer(0));

    for(vtkIdType destinationId=0; destinationId<destination->GetNumberOfPoints(); ++destinationId)
    {
        double* destinationPoint = destination->GetPoint(destinationId);
        vtkIdType const sourceId = sourceLocator->FindClosestPoint(destinationPoint);
        float* sourceNormal = sourceNormals+3*sourceId;
        normals->SetTuple(destinationId, sourceNormal);
    }

    return normals;
}
