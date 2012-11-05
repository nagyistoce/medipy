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


#ifndef _BTK_TRACTOGRAPHY_ALGORITHM_cxx
#define _BTK_TRACTOGRAPHY_ALGORITHM_cxx


#include "itkTractographyAlgorithm.h"


// STL includes
#include "algorithm"

namespace itk
{

template<typename ModelType, typename MaskType>
TractographyAlgorithm<ModelType, MaskType>
::TractographyAlgorithm(): m_Mask(NULL), m_InputModel(NULL)
{
    //---
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::PrintSelf(std::ostream &os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::ThreadedUpdate(PointType seed, unsigned int index, int threadId)
{
    // Start tractography from seed point
    typename InterpolateModelType::Pointer interpolate = InterpolateModelType::New();
    typename VectorImageToImageAdaptorType::Pointer adaptor = VectorImageToImageAdaptorType::New(); 
    adaptor->SetImage( m_InputModel );
    interpolate->SetInputImage( adaptor );
    FiberType currentFiber = PropagateSeed(seed,interpolate,adaptor);
    // Save fiber
    m_OutputFibers[index] = currentFiber;
}

template<typename ModelType, typename MaskType>
ITK_THREAD_RETURN_TYPE 
TractographyAlgorithm<ModelType, MaskType>
::ThreaderCallback( void *arg )
{
    ThreadStruct *str;
    unsigned int threadId, threadCount;

    threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
    threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

    str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

    // execute the actual method with appropriate parameters
    int seedsPerThread = str->Filter->m_Seeds.size()/threadCount;
    std::vector<PointType> seeds;
    std::vector<unsigned int> indexs;
    for(unsigned int i=seedsPerThread*threadId; i<seedsPerThread*(threadId+1); i++) {
        seeds.push_back(str->Filter->m_Seeds[i]);
        indexs.push_back(i);
    }
    if ((seedsPerThread*threadCount+threadId)<str->Filter->m_Seeds.size()) {
        seeds.push_back(str->Filter->m_Seeds[seedsPerThread*threadCount+threadId]);
        indexs.push_back(seedsPerThread*threadCount+threadId);
    }

    //std::cout << "thread " << threadId << "; size " << seeds.size() << "; nb " << seedsPerThread << std::endl;
    for(unsigned int i=0; i<seeds.size(); i++) {
        PointType seed = seeds[i];
        unsigned int index = indexs[i];
        str->Filter->ThreadedUpdate(seed, index, threadId);
    }

    return ITK_THREAD_RETURN_VALUE;
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::Update()
{
    if (m_InputModel.IsNotNull()) {

         m_size = m_InputModel->GetLargestPossibleRegion().GetSize();

        // Allocate memory to store fibers
        m_OutputFibers.resize(m_Seeds.size());

        // Set up the multithreaded processing
        ThreadStruct str;
        str.Filter = this;
        this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
        this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

        // multithread the execution
        this->GetMultiThreader()->SingleMethodExecute();

        /* typename std::vector< PointType >::iterator it;

        for (it=m_Seeds.begin(); it<m_Seeds.end(); it++) {

            // Get the seed point to process
            PointType seed = *it;
            // Start tractography from seed point
            FiberType currentFiber = PropagateSeed(seed);
            // Save fiber
            m_OutputFibers.push_back(currentFiber);

        } // for each seed*/
    }
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::SetSeed(unsigned int index, PointType seed)
{
    m_Seeds[index] = seed;
}

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::AppendSeed(PointType seed)
{
    m_Seeds.push_back(seed);
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
unsigned int
TractographyAlgorithm<ModelType, MaskType>
::GetNumberOfSeeds() const
{
    return m_Seeds.size();
}

template<typename ModelType, typename MaskType>
unsigned int
TractographyAlgorithm<ModelType, MaskType>
::GetNumberOfFibers() const
{
    return m_OutputFibers.size();
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
typename TractographyAlgorithm<ModelType, MaskType>::FiberType
TractographyAlgorithm<ModelType, MaskType>
::GetOutputFiber(unsigned int i) const
{
    if (i<m_OutputFibers.size()) {
        return m_OutputFibers[i];
    }
    else {
        throw "Wrong index to access fibers!";
    }
}

} // namespace itk

#endif
