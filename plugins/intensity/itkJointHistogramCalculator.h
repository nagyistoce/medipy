/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _5ef7344d_e3a2_4baf_905e_685124938e09
#define _5ef7344d_e3a2_4baf_905e_685124938e09

#include <ostream>

#include <itkIndent.h>
#include <itkHistogram.h>
#include <itkNumericTraits.h>
#include <itkObject.h>

namespace itk
{

/**
 * @brief Compute the joint histogram of two images.
 * 
 * This class can compute the joint histogram using two methods:
 *   * the "usual" (and default) one, where each pair of intensities contributes
 *     to one pixel in the joint histgoram
 *   * a method based on linear interpolation, where for each voxel, a line
 *     segment from the current pair of intensity to the pair of neighbor 
 *     intensities is drawn on the joint histogram.
 * 
 * The two input images shall have the same region. An optional mask can be
 * specified: in that case, the mask shall have the same region as the images.
 * When a mask is specified, only pixels that are in the mask participate to the
 * joint histogram. A mask value (defaulting to one) controls more finely the
 * pixels of the mask which are considered.
 */
template<typename TImage, typename TMask>
class JointHistogramCalculator: public Object
{
public:
    /** Standard typedefs */
    typedef JointHistogramCalculator Self;
    typedef Object Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(JointHistogramCalculator, Object);

    /** standard New() method support */
    itkNewMacro(Self);

    typedef TImage ImageType;
    typedef typename ImageType::PixelType PixelType;
    
    typedef TMask MaskType;
    typedef typename MaskType::PixelType MaskPixelType;

    typedef typename NumericTraits<PixelType>::RealType ValueRealType;
    typedef Statistics::Histogram<ValueRealType, 2> HistogramType;
    
    struct Method
    {
        enum Type
        {
            NEAREST_NEIGHBOR,
            LINEAR_INTERPOLATION,
        };
    };
    
    itkGetConstObjectMacro(Image1, ImageType);
    itkSetObjectMacro(Image1, ImageType);
    
    itkGetConstObjectMacro(Image2, ImageType);
    itkSetObjectMacro(Image2, ImageType);

    itkGetMacro(BinsCount1, unsigned int);
    itkSetMacro(BinsCount1, unsigned int);
    
    itkGetMacro(BinsCount2, unsigned int);
    itkSetMacro(BinsCount2, unsigned int);
    
    itkGetConstObjectMacro(Mask, MaskType);
    itkSetObjectMacro(Mask, MaskType);
    
    itkGetMacro(MaskValue, MaskPixelType);
    itkSetMacro(MaskValue, MaskPixelType);
    
    itkGetEnumMacro(Method, typename Method::Type);
    itkSetEnumMacro(Method, typename Method::Type);
    void SetMethodToNearestNeighbor() { this->SetMethod(Method::NEAREST_NEIGHBOR); }
    void SetMethodToLinearInterpolation() { this->SetMethod(Method::LINEAR_INTERPOLATION); }

    virtual void Compute();

    itkGetObjectMacro(Histogram, HistogramType);
    itkGetConstObjectMacro(Histogram, HistogramType);

protected:
    JointHistogramCalculator();
    virtual ~JointHistogramCalculator() {};
    void PrintSelf(std::ostream & os, Indent indent) const;
    
private:
    typename ImageType::Pointer m_Image1;
    typename ImageType::Pointer m_Image2;
    
    unsigned int m_BinsCount1;
    unsigned int m_BinsCount2;
    typename MaskType::Pointer m_Mask;
    MaskPixelType m_MaskValue;
    typename Method::Type m_Method;
    
    typename HistogramType::Pointer m_Histogram;
    
    JointHistogramCalculator(Self const &); //purposely not implemented
    void operator=(Self const &); //purposely not implemented
    
    void InitializeHistogram();
    
    template<typename TIterator>
    void UpdateNearestNeighbor(TIterator const & p1, TIterator const & p2);
    
    template<typename TIterator>
    void UpdateLinearInterpolation(TIterator const & p1, TIterator const & p2);
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkJointHistogramCalculator.txx"
#endif

#endif // _5ef7344d_e3a2_4baf_905e_685124938e09
