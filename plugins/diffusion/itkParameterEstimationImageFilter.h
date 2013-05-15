/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _63915978_9892_4957_a7ab_d27e2eedbd6e
#define _63915978_9892_4957_a7ab_d27e2eedbd6e

#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>

namespace itk
{

/**
 * \class ParameterEstimationImageFilter
 * \brief Estimate the mean and variance from a second order diffusion tensor image
 * 
 */

template<typename TTensorImage, typename TMeanImage=TTensorImage, 
    typename TVarianceImage=itk::Image<
        typename TTensorImage::PixelType::ValueType, TTensorImage::ImageDimension> >
class ParameterEstimationImageFilter : 
    public ImageToImageFilter<TTensorImage, TVarianceImage>
{
public :
    /** Standard class typedefs. */
    typedef ParameterEstimationImageFilter Self;
    typedef ImageToImageFilter<TTensorImage, TVarianceImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ParameterEstimationImageFilter, ImageToImageFilter);

    typedef itk::Image<typename TTensorImage::PixelType::ValueType, 
                       TTensorImage::ImageDimension> MaskType;

    /** @brief Return the size of the neighborhood in the plane direction. */
    itkGetConstMacro(SizePlane, unsigned int);
    /** @brief Set the size of the neighborhood in the plane direction, default to 3. */
    itkSetMacro(SizePlane, unsigned int);
    
    /** @brief Return the size of the neighborhood in the normal direction. */
    itkGetConstMacro(SizeDepth, unsigned int);
    /** @brief Set the size of the neighborhood in the normal direction, default to 3. */
    itkSetMacro(SizeDepth, unsigned int);
    
    /** @brief Return the mask. */
    itkGetConstObjectMacro(Mask, MaskType);
    /** @brief Set the mask, default to NULL (no mask is used). */
    itkSetObjectMacro(Mask, MaskType);
    
    /** @brief Return the mean image. */
    TMeanImage const * GetMeanImage() const;
    
    /** @brief Return the variance image. */
    TVarianceImage const * GetVarianceImage() const;

protected :
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    
    ParameterEstimationImageFilter();
    ~ParameterEstimationImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    DataObject::Pointer MakeOutput(unsigned int index);
    void AllocateOutputs();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    unsigned int m_SizePlane;
    unsigned int m_SizeDepth;
    typename MaskType::Pointer m_Mask;
    
    ParameterEstimationImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParameterEstimationImageFilter.txx"
#endif

#endif // _63915978_9892_4957_a7ab_d27e2eedbd6e
