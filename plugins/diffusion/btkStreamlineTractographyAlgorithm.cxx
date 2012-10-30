/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 22/08/2012
  Author(s): Julien Pontabry (pontabry@unistra.fr)
  
  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.
  
==========================================================================*/

#include "btkStreamlineTractographyAlgorithm.h"


// VTK includes
/*#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyLine.h"*/


namespace btk
{

StreamlineTractographyAlgorithm::StreamlineTractographyAlgorithm() : m_StepSize(0.5), m_UseRungeKuttaOrder4(false), m_ThresholdAngle(M_PI/3.0f), Superclass()
{
    // ----
    m_Calculator.SetDimension(3);
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PropagateSeed(Self::PhysicalPoint point)
{
    //
    // Estimate Fiber
    //

    std::vector< Self::PhysicalPoint > points;
    points.push_back(point);

    if(m_UseRungeKuttaOrder4)
    {
        Self::PropagateSeedRK4(points);
    }
    else // m_UseRungeKuttaOrder4 = false
    {
        Self::PropagateSeedRK1(points);
    }

    //
    // Build graphical fiber
    //

    if(points.size() > 1)
    {
        /*m_CurrentFiber = vtkSmartPointer< vtkPolyData >::New();

        // Graphical representation structures
        vtkSmartPointer< vtkPoints >  vpoints = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer< vtkCellArray > lines = vtkSmartPointer< vtkCellArray >::New();
        vtkSmartPointer< vtkPolyLine >   line = vtkSmartPointer< vtkPolyLine >::New();

        line->GetPointIds()->SetNumberOfIds(points.size());

        for(unsigned int i = 0; i < points.size(); i++)
        {
            vpoints->InsertNextPoint(-points[i][0], -points[i][1], points[i][2]);
            line->GetPointIds()->SetId(i,i);
        }

        lines->InsertNextCell(line);
        m_CurrentFiber->SetPoints(vpoints);
        m_CurrentFiber->SetLines(lines);*/

        m_CurrentFiber.SetSize(points.size());
        for(unsigned int i = 0; i < points.size(); i++) {
            itk::Vector<float,3> v;
            v[0] = -points[i][0];
            v[1] = -points[i][1];
            v[2] = points[i][2];
            m_CurrentFiber[i] = v;
        }
    }
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PropagateSeedRK4(std::vector< Self::PhysicalPoint > &points)
{
    // Usefull constants
    float stepSize_2 = m_StepSize / 2.f;
    float stepSize_6 = m_StepSize / 6.f;

    bool stop = false;
    PhysicalPoint nextDirection = MeanDirectionsAt(points.back())[0];

    do
    {
        // This use the RK4 to compute the next point (Runge-Kutta, order 4)
        Self::PhysicalPoint lastPoint = points.back();

        PhysicalPoint k1 = nextDirection;
        Self::PhysicalPoint nextPoint;
        for (unsigned int i=0; i<3; i++) { nextPoint[i] = lastPoint[i] + (k1[i]*stepSize_2); }
        std::vector< PhysicalPoint > directions = MeanDirectionsAt(nextPoint);
        PhysicalPoint k2 = Self::SelectClosestDirection(directions, k1);

        if(k2 != PhysicalPoint())
        {
            for (unsigned int i=0; i<3; i++) { nextPoint[i] = lastPoint[i] + (k2[i]*stepSize_2); }
            directions = MeanDirectionsAt(nextPoint);
            PhysicalPoint k3 = Self::SelectClosestDirection(directions, k2);

            if(k3 != PhysicalPoint())
            {
                for (unsigned int i=0; i<3; i++) { nextPoint[i] = lastPoint[i] + (k3[i]*m_StepSize); }
                directions = MeanDirectionsAt(nextPoint);
                PhysicalPoint k4 = Self::SelectClosestDirection(directions, k3);

                if(k4 != PhysicalPoint())
                {
                    for (unsigned int i=0; i<3; i++) { nextPoint[i] = lastPoint[i] + (k1[i] + (k2[i]*2.f) + (k3[i]*2.f) + k4[i]) * stepSize_6; }

                    // Check if the physical point is in the mask
                    Self::MaskImage::IndexType maskIndex;
                    m_Mask->TransformPhysicalPointToIndex(nextPoint, maskIndex);

                    if(m_Mask->GetPixel(maskIndex) == 0)
                    {
                        stop = true;
                    }
                    else // m_Mask->GetPixel(maskIndex) != 0
                    {
                        // Add the new point
                        points.push_back(nextPoint);

                        // Search next direction
                        std::vector< PhysicalPoint > meanDirections = MeanDirectionsAt(nextPoint);
                        nextDirection = Self::SelectClosestDirection(meanDirections, k1);

                        if(nextDirection == PhysicalPoint())
                        {
                            stop = true;
                        }
                    }
                }
                else
                {
                    stop = true;
                }
            }
            else
            {
                stop = true;
            }
        }
        else
        {
            stop = true;
        }
    } while(!stop);
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PropagateSeedRK1(std::vector< Self::PhysicalPoint > &points)
{
    bool stop = false;
    PhysicalPoint nextDirection = MeanDirectionsAt(points.back())[0];

    do
    {
        // This use the RK0 (Euler) to compute the next point (Runge-Kutta, order 0, Euler method)
        PhysicalPoint k1 = nextDirection;
        Self::PhysicalPoint nextPoint;
        for (unsigned int i=0; i<3; i++) { nextPoint[i] = points.back()[i] + (k1[i]*m_StepSize); }

        // Check if the physical point is in the mask
        Self::MaskImage::IndexType maskIndex;
        m_Mask->TransformPhysicalPointToIndex(nextPoint, maskIndex);

        if(m_Mask->GetPixel(maskIndex) == 0)
        {
            stop = true;
        }
        else // m_Mask->GetPixel(maskIndex) != 0
        {
            // Add the new point
            points.push_back(nextPoint);

            // Search next direction
            std::vector< PhysicalPoint > meanDirections = MeanDirectionsAt(nextPoint);
            nextDirection = Self::SelectClosestDirection(meanDirections, k1); // critère d'arrêt ici

            if(nextDirection==PhysicalPoint())
            {
                stop = true;
            }
        }
    } while(!stop);
}

//----------------------------------------------------------------------------------------

StreamlineTractographyAlgorithm::PhysicalPoint StreamlineTractographyAlgorithm::SelectClosestDirection(std::vector< PhysicalPoint > &meanDirections, 
                                                                                                        PhysicalPoint &previousVector)
{
    unsigned int meanDirectionsSize = meanDirections.size();
    PhysicalPoint nextDirection;

    if(meanDirectionsSize == 1)
    {
        nextDirection = meanDirections[0];
    }
    else if(meanDirectionsSize > 1)
    {
        // The next direction is choosen to be the closer to the previous one.
        float   minDotProduct = std::abs(meanDirections[0][0]*previousVector[0]+meanDirections[0][1]*previousVector[1]+meanDirections[0][2]*previousVector[2]);
        unsigned int minIndex = 0;

        for(unsigned int i = 1; i < meanDirections.size(); i++)
        {
            float dotProduct = std::abs(meanDirections[i][0]*previousVector[0]+meanDirections[i][1]*previousVector[1]+meanDirections[i][2]*previousVector[2]);

            if(dotProduct < minDotProduct)
            {
                minDotProduct = dotProduct;
                minIndex      = i;
            }
        } // for each mean direction

        nextDirection = meanDirections[minIndex];
    }

    return nextDirection;
}



std::vector< StreamlineTractographyAlgorithm::PhysicalPoint >  StreamlineTractographyAlgorithm::MeanDirectionsAt(PhysicalPoint vector)
{
    itk::Vector<float,6> dt6;    
    for( unsigned int i=0; i<6; i++) {
        m_Adaptor->SetExtractComponentIndex( i );
        m_Adaptor->SetImage( m_InputModelImage );
        m_Adaptor->Update();
        m_ModelImageFunction->SetInputImage( m_Adaptor );
        dt6[i] = m_ModelImageFunction->Evaluate(vector);
    } 
    //ModelImage::PixelType tensor = m_ModelImageFunction->EvaluateAtContinuousIndex(vector);

    InputMatrixType dt33(3,3);
    dt33[0][0] =  (double) dt6[0];
    dt33[1][1] =  (double) dt6[3];
    dt33[2][2] =  (double) dt6[5];
    dt33[0][1] =  (double) dt6[1];
    dt33[0][2] =  (double) dt6[2];
    dt33[1][2] =  (double) dt6[4];
    dt33[1][0] =  dt33[0][1];
    dt33[2][0] =  dt33[0][2];
    dt33[2][1] =  dt33[1][2];

    EigenValuesArrayType eigenvalues;
    EigenVectorMatrixType eigenvectors;

    m_Calculator.ComputeEigenValuesAndVectors(dt33,eigenvalues,eigenvectors);

    /*ModelImage::PixelType::EigenValuesArrayType   eigenValues;
    ModelImage::PixelType::EigenVectorsMatrixType eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);*/

    std::vector< PhysicalPoint > meanDirections;
    PhysicalPoint p;
    p[0] = eigenvectors(2,0);
    p[1] = eigenvectors(2,1);
    p[2] = eigenvectors(2,2);
    meanDirections.push_back(p);

    return meanDirections;
}




} // namespace btk
