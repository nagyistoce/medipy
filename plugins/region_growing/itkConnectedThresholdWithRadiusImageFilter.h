/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef e5339eb8_a2c0_4cc7_9105_2f3b9613b89a
#define e5339eb8_a2c0_4cc7_9105_2f3b9613b89a

#include <itkImageToImageFilter.h>
#include <itkSimpleDataObjectDecorator.h>

#include <vector>

namespace itk {

/**
 * \class ConnectedThresholdWithRadiusImageFilter
 * \brief Label pixels that are connected to a seed, lie within a range of values,
 * and are within a distance of seeds.
 * 
 * ConnectedThresholdImageFilter labels pixels with ReplaceValue that are
 * connected to an initial Seed, lie within [Lower, Upper] threshold range, and
 * are within a Radius from one of the seeds.
 *
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT ConnectedThresholdWithRadiusImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef ConnectedThresholdWithRadiusImageFilter Self;
    typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods).  */
    itkTypeMacro(ConnectedThresholdWithRadiusImageFilter, ImageToImageFilter);

    /** Useful typedefs. */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef SimpleDataObjectDecorator<InputImagePixelType> InputPixelObjectType;
    typedef typename Superclass::OutputImagePixelType OutputImagePixelType;

    // Get/Set the lower threshold value
    InputImagePixelType GetLower() const;
    void SetLower(InputImagePixelType threshold);

    // Get/Set the upper threshold value
    InputImagePixelType GetUpper() const;
    void SetUpper(InputImagePixelType threshold);

    // Get/Set the lower threshold process object.
    InputPixelObjectType * GetLowerInput();
    void SetLowerInput(InputPixelObjectType const * input);

    // Get/Set the upper threshold process object.
    InputPixelObjectType * GetUpperInput();
    void SetUpperInput(InputPixelObjectType const * input);

    // Get/Set the maximum radius of growing. Defaults to 1. */
    itkGetMacro(Radius, float);
    itkSetMacro(Radius, float);

    void SetSeed(InputImageIndexType const & seed);
    void AddSeed(InputImageIndexType const & seed);
    void ClearSeeds();

    // Get/Set the replace value. Defaults to 1. */
    itkGetMacro(ReplaceValue, OutputImagePixelType);
    itkSetMacro(ReplaceValue, OutputImagePixelType);

protected:

    InputImagePixelType m_Lower;
    InputImagePixelType m_Upper;
    std::vector<InputImageIndexType> m_Seeds;
    float m_Radius;
    OutputImagePixelType m_ReplaceValue;

    ConnectedThresholdWithRadiusImageFilter();
    ~ConnectedThresholdWithRadiusImageFilter();
    void PrintSelf(std::ostream & os, Indent indent) const;
    void GenerateData();

private:
    ConnectedThresholdWithRadiusImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // namespace itk

#include "itkConnectedThresholdWithRadiusImageFilter.txx"

#endif // e5339eb8_a2c0_4cc7_9105_2f3b9613b89a
