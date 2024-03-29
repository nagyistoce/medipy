/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "itkHMCSegmentationImageFilter.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>

#include <itkConstantPadImageFilter.h>
#include <itkImage.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "EM.h"
#include "HilbertCurveChainGenerator.h"
#include "HMCInitializer.h"

namespace itk
{

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
typename HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>::OutputImagePointer
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::GetSegmentationImage()
{
    return dynamic_cast<OutputImageType *>(this->ProcessObject::GetOutput(0));
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
typename HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>::OutputImagePointer
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::GetOutliersImage()
{
    return dynamic_cast<OutputImageType *>(this->ProcessObject::GetOutput(1));
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::HMCSegmentationImageFilter()
: m_MaskImage(0), m_FlairImage(-1), m_Iterations(5), m_Modalities(0), 
  m_DisplayOutliers(true), m_OutliersCriterion(0), m_Threshold(0)
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(0);
 
    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::~HMCSegmentationImageFilter()
{
    // Nothing to do.
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::GenerateData()
{
    std::vector<InputImageConstPointer> images;
    for(unsigned int i=0; i < this->m_Modalities; ++i)
    {
        images.push_back(this->GetInput(i));
    }

    std::vector<InputImageConstPointer> atlas;
    for(unsigned int i=this->m_Modalities; i <this->GetNumberOfInputs(); ++i)
    {
        atlas.push_back(this->GetInput(i));
    }
    
    std::vector<InputImageConstPointer> all_images;
    for(unsigned int i=0; i <this->GetNumberOfInputs(); ++i)
    {
        all_images.push_back(this->GetInput(i));
    }
    
    //parcours Hilbert-Peano pour récuperer des vecteurs et suppression des 
    // données non masquées
    ChainGeneratorType chain_generator;
    chain_generator(all_images, this->m_MaskImage);
    
    typedef typename ChainGeneratorType::ImageChainsType ImageChainsType;
    
    ImageChainsType const & image_chains = chain_generator.GetImageChains();
    unsigned int const chain_size=image_chains.columns();
    
    ImageChainsType chain(images.size(), chain_size);
    for(unsigned int i=0; i!=images.size(); ++i)
    {
        chain.set_row(i, image_chains.get_row(i));
    }
    
    ImageChainsType chain_atlas(atlas.size(), chain_size);
    for(unsigned int i=0; i!=atlas.size(); ++i)
    {
        chain_atlas.set_row(i, image_chains.get_row(i+images.size()));
    }
    
    //Initialisation des paramètres: KMean, aij, PIi
    HMCInitializer hmc_initializer;
    hmc_initializer(chain, images.size(), atlas.size());
    
    vnl_matrix<double> Mu = hmc_initializer.GetMu();
    std::vector<vnl_matrix<double> > Sigma = hmc_initializer.GetSigma();
    vnl_vector<double> PIi = hmc_initializer.GetPIi();
    vnl_matrix<double> aij = hmc_initializer.Getaij();
    
    //Estimation des paramètres (EM)
    EM em(PIi, aij, Mu, Sigma);
    em.HMCRobusteAtlasFlair(chain, chain_atlas, this->m_Iterations,
                            this->m_FlairImage, this->m_OutliersCriterion, 
                            m_Threshold);
    
    // MPM
    vnl_vector<int> ChainSeg(chain_size);
    vnl_vector<int> ChainLesions(chain_size);
    em.SegMPMRobusteAtlas(chain, chain_atlas, ChainSeg, ChainLesions, 
                          this->m_DisplayOutliers);
    
    //Parcours d Hilbert-Peano inverse pour obtenir l image segmentée
    //mettre l'image à la bonne taille
    ScanType const & scan = chain_generator.GetScan();
    
    OutputImagePointer segmentation = this->GetSegmentationImage();
    segmentation->SetRegions(images[0]->GetLargestPossibleRegion());
    segmentation->Allocate();
    Self::_chain_to_image(ChainSeg, scan, segmentation);
    
    OutputImagePointer outliers = this->GetOutliersImage();
    outliers->SetRegions(images[0]->GetLargestPossibleRegion());
    outliers->Allocate();
    Self::_chain_to_image(ChainLesions, scan, outliers);
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
DataObject::Pointer
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::MakeOutput(unsigned int index)
{
    DataObject::Pointer output;
 
    switch(index)
    {
    case 0:
        output = (OutputImageType::New()).GetPointer();
        break;
    case 1:
        output = (OutputImageType::New()).GetPointer();
        break;
    default:
        std::cerr << "No output " << index << std::endl;
        output = NULL;
    }
    
    return output.GetPointer();
}

template<typename TInputImage, typename TMaskImage, typename TOutputImage>
void
HMCSegmentationImageFilter<TInputImage, TMaskImage, TOutputImage>
::_chain_to_image(vnl_vector<int> const & chain, ScanType const & scan, 
    OutputImagePointer image)
{
    image->FillBuffer(0);
    
    for(unsigned int i=0; i<scan.size(); ++i)
    {
        int const value = chain[i];
        typename OutputImageType::IndexType const index = image->ComputeIndex(scan[i]);
        
        image->SetPixel(index, value);
    }
}

}
