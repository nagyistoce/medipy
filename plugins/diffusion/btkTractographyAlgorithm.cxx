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


#include "btkTractographyAlgorithm.h"


// STL includes
#include "algorithm"

namespace btk
{

template<typename ModelType, typename MaskType>
TractographyAlgorithm<ModelType, MaskType>
::TractographyAlgorithm(): m_Mask(NULL), m_InputModel(NULL)
{
    // Create interpolate function on model image
    m_InterpolateModelFunction = InterpolateModelType::New();
    m_Adaptor = VectorImageToImageAdaptorType::New(); 
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

template<typename ModelType, typename MaskType>
void 
TractographyAlgorithm<ModelType, MaskType>
::Update()
{
    if (m_InputModel.IsNotNull()) {

        m_size = m_InputModel->GetLargestPossibleRegion().GetSize();
        typename std::vector< PointType >::iterator it;

        for (it=m_Seeds.begin(); it<m_Seeds.end(); it++) {

            // Get the seed point to process
            PointType seed = *it;
            // Start tractography from seed point
            FiberType currentFiber = PropagateSeed(seed);
            // Save fiber
            m_OutputFibers.push_back(currentFiber);

        } // for each seed
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

} // namespace btk

#endif
