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

#ifndef BTK_TRACTOGRAPHY_ALGORITHM_H
#define BTK_TRACTOGRAPHY_ALGORITHM_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkPoint.h"
#include "itkContinuousIndex.h"
#include "itkVector.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorImageToImageAdaptor.h"

namespace btk
{

/**
 * @brief Representation of an abstract tractography algorithm.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
template<typename ModelType, typename MaskType>
class ITK_EXPORT TractographyAlgorithm : public itk::ProcessObject
{
    public:
        typedef TractographyAlgorithm                       Self;
        typedef itk::ProcessObject                          Superclass;
        typedef itk::SmartPointer< Self >                   Pointer;
        typedef itk::SmartPointer< const Self >             ConstPointer;

        // TODO : template parameters    
        /*typedef itk::VectorImage< float,3 >                 ModelType;        
        typedef itk::Image< short,3 >                       MaskType;*/

        typedef typename ModelType::PointType               PointType;
        typedef typename PointType::VectorType              VectorType;   
        typedef typename std::vector< PointType >           FiberType;    


        itkTypeMacro(TractographyAlgorithm, itk::ProcessObject);

        itkSetMacro(Mask, typename MaskType::Pointer);
        itkGetConstMacro(Mask, typename MaskType::Pointer);

        itkSetMacro(InputModel, typename ModelType::Pointer);
        itkGetConstMacro(InputModel, typename ModelType::Pointer);

        //itkSetMacro(Seeds, itk::Vector< VectorType >);
        //itkGetConstMacro(Seeds, itk::Vector< VectorType >);

        // Run the algorithm.
        virtual void ThreadedUpdate(PointType seed, unsigned int index, int threadId);
        virtual void Update();
        // Static function used as a "callback" by the MultiThreader.  The threading library will call this routine for each thread, which will delegate the
        // control to ThreadedUpdate().
        static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );
        // Accessor to output estimated fibers
        FiberType GetOutputFiber(unsigned int i) const;
        // Get number of output fibers
        unsigned int GetNumberOfSeeds() const;
        unsigned int GetNumberOfFibers() const;
        // Get fiber length
        static typename VectorType::ValueType GetLength(FiberType const &fiber);
        // Set seeds (index must be in range of m_Seeds size)
        void SetSeed(unsigned int index, PointType seed);
        void AppendSeed(PointType seed);

    protected:
        TractographyAlgorithm();
        // Print a message on output stream.
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;
        // Propagate using the tractography algorithm at a seed point.
        virtual FiberType PropagateSeed(PointType const &point) = 0;
        // Internal structure used for passing image data into the threading library */
        struct ThreadStruct {
            Pointer Filter;
        };

    protected:

        typedef itk::VectorImageToImageAdaptor< float,3 >                                        VectorImageToImageAdaptorType; 
        typedef itk::LinearInterpolateImageFunction< VectorImageToImageAdaptorType,float >       InterpolateModelType;

        // Mask image determining where the tractography algorithm can process.
        typename MaskType::Pointer m_Mask;
        // Tensor image.
        typename ModelType::Pointer m_InputModel;
        typename ModelType::SizeType m_size;
        // Estimated fibers.
        std::vector< FiberType > m_OutputFibers;
        // Seeds.
        std::vector< PointType > m_Seeds;
        // Interpolation tools.
        //InterpolateModelType::Pointer m_InterpolateModelFunction;
        //VectorImageToImageAdaptorType::Pointer m_Adaptor;

};

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkTractographyAlgorithm.cxx"
#endif

#endif // BTK_TRACTOGRAPHY_ALGORITHM_H
