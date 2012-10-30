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

#include "btkTractographyAlgorithm.h"


// STL includes
#include "algorithm"

namespace btk
{

TractographyAlgorithm::TractographyAlgorithm()
{
    // Create interpolate function on model image
    m_InterpolateModelFunction = InterpolateModelType::New();
    m_Adaptor = VectorImageToImageAdaptorType::New(); 
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::Update()
{
    unsigned int numberOfSeeds = m_Seeds.GetNumberOfComponents();

    for (unsigned int l=0; l<numberOfSeeds; l++) {

        // Get the seed point to process
        VectorType seed_ = m_Seeds[l];
        PointType seed;
        for (unsigned int ll=0; ll<seed_.GetNumberOfComponents(); ll++) { seed[ll] = seed_[ll]; }
         // Check if the physical point is in the mask
        MaskType::IndexType maskIndex;
        m_Mask->TransformPhysicalPointToIndex(seed, maskIndex);

        if (m_Mask->GetPixel(maskIndex) != 0) {

            // Init fiber under construction.
            FiberType m_CurrentFiber;
            // Start tractography from seed point
            PropagateSeed(seed);
            // Save fiber
            if (m_CurrentFiber.GetNumberOfElements()>1) { m_OutputFibers.push_back(m_CurrentFiber); }

        }
    } // for each seed
}

//----------------------------------------------------------------------------------------

TractographyAlgorithm::FiberType TractographyAlgorithm::GetOutputFiber(unsigned int i) const
{
    return m_OutputFibers[i];
}

} // namespace btk
