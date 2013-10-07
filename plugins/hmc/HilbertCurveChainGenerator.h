/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _520dc57b_a1d1_4849_8e56_f982399ef678
#define _520dc57b_a1d1_4849_8e56_f982399ef678

#include <vector>

#include <itkImage.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

/**
 * @brief Generate chains from a vector of images and an optional mask, 
 *        according to a Hilbert-Peano scan.
 */
class HilbertCurveChainGenerator
{
public:
    typedef HilbertCurveChainGenerator Self;
    
    typedef itk::Image<float, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;
    
    typedef itk::Image<float, 3> MaskType;
    typedef MaskType::Pointer MaskPointer;
    
    typedef itk::Image<int, 3> ScanType;
    typedef ScanType::Pointer ScanPointer;
    typedef ScanType::ConstPointer ScanConstPointer;
    
    typedef vnl_matrix<double> ImageChainsType;
    typedef vnl_vector<int> MaskChainType;

    HilbertCurveChainGenerator();
    ~HilbertCurveChainGenerator();

    void operator()(std::vector<ImagePointer> const & images, 
                    MaskPointer mask=MaskPointer());
    
    ImageChainsType const & GetImageChains() const;
    MaskChainType const & GetMaskChain() const;
    ScanPointer GetScan() const;

private:
    ImageChainsType m_ImageChains;
    MaskChainType m_MaskChain;
    ScanPointer m_Scan;
    
    static unsigned int _find_cube_length(ImagePointer image);
    static ScanPointer _hilbert_peano_scan(int cube_size);

    void _image_scan(std::vector<ImagePointer> const & images, MaskPointer mask);
};

#endif // _520dc57b_a1d1_4849_8e56_f982399ef678
