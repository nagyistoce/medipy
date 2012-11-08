/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkBootstrapParameterEstimationImageFilter_h
#define _itkBootstrapParameterEstimationImageFilter_h


#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>


namespace itk
{

/**
 * \class BootstrapParameterEstimationImageFilter
 * \brief Estimate the mean and variance from diffusion data using a bootstrap procedure
 * 
 */

template<typename TInputImage, typename TOutputImage>
class BootstrapParameterEstimationImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef BootstrapParameterEstimationImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BootstrapParameterEstimationImageFilter, ImageToImageFilter);

    /** Useful typedefs */
    typedef typename Superclass::InputImageType     InputImageType;
    typedef typename Superclass::OutputImageType    OutputImageType;
    typedef typename TOutputImage::PixelType        OutputPixelType;
    typedef typename TInputImage::PixelType         InputPixelType;
    typedef typename OutputPixelType::ValueType     OutputValueType;

    typedef Image<float, 3>                         MaskType;
    typedef Point<float,3>                          DirectionType; 
    typedef float                                   BvalType;
    typedef vnl_matrix<float>                       BMatrixType; 
    typedef vnl_matrix<float>                       TensorType; 

    /** Accessors */
    itkSetMacro(BVal, BvalType);
    itkGetConstMacro(BVal, BvalType);
    itkSetMacro(SizePlane, unsigned int);
    itkGetConstMacro(SizePlane, unsigned int);
    itkSetMacro(SizeDepth, unsigned int);
    itkGetConstMacro(SizeDepth, unsigned int);
    itkSetMacro(NumberOfBootstrap, unsigned int);
    itkGetConstMacro(NumberOfBootstrap, unsigned int);
    itkSetMacro(UseSpatialBootstrap, bool);
    itkGetConstMacro(UseSpatialBootstrap, bool);
    itkSetObjectMacro(Mask, MaskType);
    itkGetObjectMacro(Mask, MaskType);
    itkGetConstMacro(NumberOfElements, unsigned int);
    void SetGradientDirection(unsigned int i, DirectionType bvec);
    DirectionType GetGradientDirection(unsigned int i);

protected :
    BootstrapParameterEstimationImageFilter();
    ~BootstrapParameterEstimationImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, int);

    std::vector< vnl_matrix<float> > SpatialBootstrapInit(typename OutputImageType::IndexType idx);
    std::vector<unsigned int> RandomSpatial(unsigned int N);
    std::vector< TensorType > SpatialBootstrapGenerator(std::vector< vnl_matrix<float> > signals);

    void ComputeParameters( std::vector< TensorType > tensors, typename OutputImageType::IndexType idx,
                            typename OutputImageType::Pointer &output_mean, typename OutputImageType::Pointer &output_var);

private :
    unsigned int m_SizePlane;
    unsigned int m_SizeDepth;
    unsigned int shift_plane;
    unsigned int shift_depth;
    unsigned int m_NumberOfElements;
    unsigned int m_NumberOfBootstrap;
    MaskType::Pointer m_Mask;
    bool m_UseSpatialBootstrap;

    std::vector<DirectionType> directions;
    BvalType m_BVal;
    BMatrixType bmatrix;
    BMatrixType invbmatrix;

    typename InputImageType::SizeType size;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBootstrapParameterEstimationImageFilter.txx"
#endif

#endif 

