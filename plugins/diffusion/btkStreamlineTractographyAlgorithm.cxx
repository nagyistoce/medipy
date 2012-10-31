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

#ifndef BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_cxx
#define BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_cxx

#include "btkStreamlineTractographyAlgorithm.h"


namespace btk
{

template<typename ModelType, typename MaskType>
StreamlineTractographyAlgorithm<ModelType, MaskType>
::StreamlineTractographyAlgorithm() : m_StepSize(0.5), m_UseRungeKuttaOrder4(false), m_ThresholdAngle(M_PI/3.0f), m_ThresholdFA(0.2), Superclass()
{
    m_Calculator.SetDimension(3);
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void
StreamlineTractographyAlgorithm<ModelType, MaskType>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------


template<typename ModelType, typename MaskType>
typename StreamlineTractographyAlgorithm<ModelType, MaskType>::FiberType
StreamlineTractographyAlgorithm<ModelType, MaskType>
::PropagateSeed(PointType const &seed)
{
    FiberType currentFiber;

    // Propagation
    if(m_UseRungeKuttaOrder4) {
        currentFiber = PropagateSeedRK4(seed);
    }
    else {
        currentFiber = PropagateSeedRK0(seed);
    }

    return currentFiber;
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
typename StreamlineTractographyAlgorithm<ModelType, MaskType>::FiberType
StreamlineTractographyAlgorithm<ModelType, MaskType>
::PropagateSeedRK0(PointType const &seed)
{
    bool stop = false;
    VectorType d1 = PropagationDirectionT2At(seed,stop);
    VectorType d2 = -d1;

    std::vector< PointType > points1;  
    std::vector< PointType > points2;

    // Track towards one direction
    PointType nextPoint = seed;
    while (!stop) {
        // Insert new point 
        points1.push_back(nextPoint);
        // This use the RK0 (Euler) to compute the next point (Runge-Kutta, order 0, Euler method)
        nextPoint = points1.back() + (d1*m_StepSize);
        // Propagation
        VectorType d = PropagationDirectionT2At(nextPoint,stop);
        // Select the best propagation direction
        d1 = SelectClosestDirectionT2(d,d1,stop);
    }

    // Track towards the opposite direction
    stop = false;
    nextPoint = seed;
    while (!stop) {
        points2.push_back(nextPoint);
        nextPoint = points2.back() + (d2*m_StepSize);
        VectorType d = PropagationDirectionT2At(nextPoint,stop);
        d2 = SelectClosestDirectionT2(d,d2,stop);
    }

    // Concatenate
    std::vector< PointType > points;
    for (unsigned int i=0; i<points2.size()-1; i++) { points.push_back(points2[points2.size()-1-i]); }
    for (unsigned int i=0; i<points1.size(); i++) { points.push_back(points1[i]); }

    return points;
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
typename StreamlineTractographyAlgorithm<ModelType, MaskType>::VectorType
StreamlineTractographyAlgorithm<ModelType, MaskType>
::PropagationDirectionT2At(PointType const &point, bool &stop)
{
    VectorType pd;

    // Check if the point is in the mask
    if (this->m_Mask.IsNotNull()) {
        typename MaskType::IndexType maskIndex;
        this->m_Mask->TransformPhysicalPointToIndex(point, maskIndex);
        if(this->m_Mask->GetPixel(maskIndex)==0) { stop = true; }
    }
    else {
        typename ModelType::IndexType modelIndex;
        this->m_InputModel->TransformPhysicalPointToIndex(point, modelIndex);
        if ( (modelIndex[0]>=(this->m_size[0]-1)) || (modelIndex[1]>=(this->m_size[1]-1)) || (modelIndex[2]>=(this->m_size[2]-1)) ||
             (modelIndex[0]<=0) || (modelIndex[1]<=0) || (modelIndex[2]<=0) ) { stop = true; }
    }

    if (!stop) {
        itk::Vector<float,6> dt6;    
        typename InterpolateModelType::Pointer m_InterpolateModelFunction = InterpolateModelType::New();
        typename VectorImageToImageAdaptorType::Pointer m_Adaptor = VectorImageToImageAdaptorType::New(); 
        for( unsigned int i=0; i<6; i++) {
            m_Adaptor->SetExtractComponentIndex( i );
            m_Adaptor->SetImage( this->m_InputModel );
            m_Adaptor->Update();
            m_InterpolateModelFunction->SetInputImage( m_Adaptor );
            dt6[i] = m_InterpolateModelFunction->Evaluate(point);
        } 

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

        pd[0] = eigenvectors(2,0);
        pd[1] = eigenvectors(2,1);
        pd[2] = eigenvectors(2,2);

        // Check if the point is within the white matter
        double ev1 = eigenvalues[2];
        double ev2 = eigenvalues[1];
        double ev3 = eigenvalues[0];
        double fa = 0.5 * ((ev1-ev2)*(ev1-ev2) + (ev2-ev3)*(ev2-ev3) + (ev3-ev1)*(ev3-ev1)) / (ev1*ev1 + ev2*ev2 + ev3*ev3);
        if (fa>0) { fa = sqrt(fa); }
        else { fa = 0; }
        if ((float)fa<m_ThresholdFA) {
            stop = true;
        }
    }

    return pd;
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
typename StreamlineTractographyAlgorithm<ModelType, MaskType>::VectorType
StreamlineTractographyAlgorithm<ModelType, MaskType>
::SelectClosestDirectionT2(VectorType const &currentDirection, VectorType const &previousDirection, bool &stop)
{
    // The new direction is choosen to be the closer to the previous one.
    float dotProduct = currentDirection*previousDirection;
    VectorType nextDirection;
    if (dotProduct<0) { nextDirection = -currentDirection; }
    else { nextDirection = currentDirection; }

    // Check if the propagation angle is allowed
    float angle = std::acos( std::abs(dotProduct) );
    if (angle > m_ThresholdAngle) { stop = true; }

    return nextDirection;
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
typename StreamlineTractographyAlgorithm<ModelType, MaskType>::FiberType
StreamlineTractographyAlgorithm<ModelType, MaskType>
::PropagateSeedRK4(PointType const &seed)
{
    // Usefull constants
    /*float stepSize_2 = m_StepSize / 2.f;
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
    } while(!stop);*/
}


} // namespace btk

#endif
