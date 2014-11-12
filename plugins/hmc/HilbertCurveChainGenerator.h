/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _520dc57b_a1d1_4849_8e56_f982399ef678
#define _520dc57b_a1d1_4849_8e56_f982399ef678

#include <stdint.h>
#include <vector>

#include <itkImage.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

/**
 * @brief Generate chains from a vector of images and an optional mask, 
 *        according to a Hilbert-Peano scan.
 */
template<typename TImage, typename TMask>
class HilbertCurveChainGenerator
{
public:
    typedef HilbertCurveChainGenerator Self;
    
    typedef TImage ImageType;
    typedef typename ImageType::ConstPointer ImageConstPointer;
    
    typedef TMask MaskType;
    typedef typename MaskType::ConstPointer MaskConstPointer;
    
    typedef std::vector<uint32_t> ScanType;
    
    typedef vnl_matrix<double> ImageChainsType;

    HilbertCurveChainGenerator();
    ~HilbertCurveChainGenerator();

    void operator()(std::vector<ImageConstPointer> const & images, 
                    MaskConstPointer mask=MaskConstPointer());
    
    ImageChainsType const & GetImageChains() const;
    ScanType const & GetScan() const;

private:
    typedef itk::Image<uint32_t, 3> Cube;
    
    ImageChainsType m_ImageChains;
    ScanType m_Scan;
    
    static unsigned int _find_cube_length(ImageType const * image);
    static Cube::Pointer _hilbert_peano_scan(unsigned int cube_size);
};

#include "HilbertCurveChainGenerator.txx"

#endif // _520dc57b_a1d1_4849_8e56_f982399ef678
