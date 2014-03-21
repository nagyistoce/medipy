/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _itkWeightedMeanImageFilter_h
#define _itkWeightedMeanImageFilter_h

#include <ostream>
#include <vector>

#include <itkImageToImageFilter.h>
#include <itkPoint.h>
#include <itkSmartPointer.h>

namespace itk
{

/**
 * \class WeightedMeanImageFilter
 * \brief Compute a Weighted Mean 
 *
 * @pre TInputImage must be a Image.
 *  
 */
template<typename TInputImage, typename TOutputImage, typename TTensorImage>
class WeightedMeanImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /// @brief Standard class typedefs.
    typedef WeightedMeanImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Output  */
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::OutputImagePixelType OutputImagePixelType;
    typedef typename Superclass::OutputImagePointer OutputImagePointer;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    
    /** Input */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;
    typedef typename Superclass::InputImagePointer InputImagePointer;
    typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    
    /** Run-time type information (and related methods). */
    itkTypeMacro(WeightedMeanImageFilter, ImageToImageFilter);
    
    typedef TTensorImage TensorImageType;
    typedef typename TensorImageType::Pointer TensorImagePointer;
    /// @brief Return the Tensor Image. 
    itkGetConstObjectMacro(TensorImage, TensorImageType);
    /// @brief Set the Tensor Image, default to NULL (no Tensor Image is used). 
    itkSetObjectMacro(TensorImage, TensorImageType);
    
    /// @brief Return the size of the neighborhood in the plane direction.
    itkGetConstMacro(SizePlaneX, unsigned int);
    /// @brief Set the size of the neighborhood in the plane direction, default to 3.
    itkSetMacro(SizePlaneX, unsigned int);
    
    /// @brief Return the size of the neighborhood in the plane direction.
    itkGetConstMacro(SizePlaneY, unsigned int);
    /// @brief Set the size of the neighborhood in the plane direction, default to 3.
    itkSetMacro(SizePlaneY, unsigned int);
    
    /// @brief Return the size of the neighborhood in the normal direction.
    itkGetConstMacro(SizeDepth, unsigned int);
    /// @brief Set the size of the neighborhood in the normal direction, default to 3.
    itkSetMacro(SizeDepth, unsigned int);

protected :
    WeightedMeanImageFilter();
    ~WeightedMeanImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

private :
    unsigned int m_SizePlaneX;
    unsigned int m_SizePlaneY;
    unsigned int m_SizeDepth;
    
    TensorImagePointer m_TensorImage;

    WeightedMeanImageFilter(Self const &); // purposely not implemented
    Self const & operator=(Self const &); // purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWeightedMeanImageFilter.txx"
#endif

#endif
