/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkParameterEstimationImageFilter_h
#define _itkParameterEstimationImageFilter_h


#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>


namespace itk
{

/**
 * \class ParameterEstimationImageFilter
 * \brief Estimate the mean and variance from a second order diffusion tensor image
 * 
 */

template<typename TInputImage, typename TOutputImage>
class ParameterEstimationImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef ParameterEstimationImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ParameterEstimationImageFilter, ImageToImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType     InputImageType;
    typedef typename Superclass::OutputImageType    OutputImageType;
    typedef typename TOutputImage::PixelType        OutputPixelType;
    typedef typename TInputImage::PixelType         InputPixelType;
    typedef typename InputPixelType::ValueType      InputValueType;
    typedef typename OutputPixelType::ValueType     OutputValueType;
    typedef Image<float, 3>                         MaskType;

    /** Accessors */
    itkSetMacro(SizePlane, unsigned int);
    itkGetConstMacro(SizePlane, unsigned int);
    itkSetMacro(SizeDepth, unsigned int);
    itkGetConstMacro(SizeDepth, unsigned int);
    itkGetConstMacro(NumberOfElements, unsigned int);
    void SetMask(MaskType *m);
    MaskType* GetMask() const;

protected :
    ParameterEstimationImageFilter();
    ~ParameterEstimationImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    unsigned int m_SizePlane;
    unsigned int m_SizeDepth;
    unsigned int m_NumberOfElements;
    unsigned int shift_plane;
    unsigned int shift_depth;
    typename InputImageType::SizeType size;
    bool masked;
    MaskType::Pointer mask;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParameterEstimationImageFilter.txx"
#endif

#endif 

