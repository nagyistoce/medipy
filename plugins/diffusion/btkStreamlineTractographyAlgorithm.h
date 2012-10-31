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

#ifndef BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_H
#define BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_H

// ITK includes
#include "itkMacro.h"
#include "itkSmartPointer.h"
#include "itkDiffusionTensor3DReconstructionImageFilter.h"
#include <itkSymmetricEigenAnalysis.h>

// Local includes
#include "btkTractographyAlgorithm.h"

namespace btk
{

/**
 * @brief Define a deterministic streamline propagation algorithm.
 * @author Julien Pontabry
 * @ingroup Tractography
 */

template<typename ModelType, typename MaskType>
class StreamlineTractographyAlgorithm : public btk::TractographyAlgorithm<ModelType, MaskType>
{
    public:
        typedef StreamlineTractographyAlgorithm                         Self;
        typedef btk::TractographyAlgorithm<ModelType, MaskType>         Superclass;
        typedef itk::SmartPointer< Self >                               Pointer;
        typedef itk::SmartPointer< const Self >                         ConstPointer;

        typedef typename Superclass::PointType               PointType;
        typedef typename Superclass::VectorType              VectorType;
        typedef typename Superclass::FiberType               FiberType;
        typedef typename Superclass::InterpolateModelType             InterpolateModelType;
        typedef typename Superclass::VectorImageToImageAdaptorType    VectorImageToImageAdaptorType;

        itkNewMacro(Self);
        itkTypeMacro(StreamlineTractographyAlgorithm, btk::TractographyAlgorithm);

        itkSetMacro(StepSize,float);
        itkGetConstMacro(StepSize,float);

        itkSetMacro(ThresholdAngle,float);
        itkGetConstMacro(ThresholdAngle,float);

        itkSetMacro(ThresholdFA,float);
        itkGetConstMacro(ThresholdFA,float);

        itkSetMacro(UseRungeKuttaOrder4,bool);
        itkGetConstMacro(UseRungeKuttaOrder4,bool);

    protected:

        typedef vnl_matrix< double >                InputMatrixType;
        typedef itk::FixedArray< double,3 >         EigenValuesArrayType;
        typedef itk::Matrix< double,3,3 >           EigenVectorMatrixType;
        typedef itk::SymmetricEigenAnalysis< InputMatrixType,EigenValuesArrayType,EigenVectorMatrixType > CalculatorType;

        StreamlineTractographyAlgorithm();
        // Print a message on output stream.
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;
        // Propagation using the streamline tractography algorithm at a seed point.
        virtual FiberType PropagateSeed(PointType const &seed);
        // Compute the propagation direction at a given point for a second order diffusion tensor model (ie, principal direction).
        VectorType PropagationDirectionT2At(PointType const &point, bool &stop);
        // Propagate a seed using the Runge-Kutta method at order 4.
        FiberType PropagateSeedRK4(PointType const &seed);
        // Propagate a seed using the Euler method (Runge-Kutta method at order 0).
        FiberType PropagateSeedRK0(PointType const &seed);
        // Select the best direction (+ or -), depending on the previous direction.
        VectorType SelectClosestDirectionT2(VectorType const &currentDirection, VectorType const &previousDirection, bool &stop);

    private:
        //Step size between two points of output.
        float m_StepSize;
        //True if the RK4 method is used or false if the RK0 (Euler) method is used.
        bool m_UseRungeKuttaOrder4;
        // Allowed angle for propagation.
        float m_ThresholdAngle;
        // FA criterion for propagation.
        float m_ThresholdFA;
        // Calculator to extract principal direction
        CalculatorType m_Calculator;
};

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkStreamlineTractographyAlgorithm.cxx"
#endif

#endif // BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_H
