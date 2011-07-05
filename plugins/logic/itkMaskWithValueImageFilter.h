/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef medipy_components_segmentation_itkmaskwithvalueimagefilter_h
#define medipy_components_segmentation_itkmaskwithvalueimagefilter_h

#include <itkBinaryFunctorImageFilter.h>
#include <itkNumericTraits.h>


namespace itk
{

/** \class MaskWithValueImageFilter
 * \brief Implements an operator for pixel-wise masking of the input
 * image with the mask.
 *
 * For each pixel, the result is as follows :
 *
 *        if pixel_from_mask_image != background_value
 *             pixel_output_image = pixel_input_image
 *        else
 *             pixel_output_image = outside_value
 *
 * The pixel from the input 1 is cast to the pixel type of the output image.
 *
 * Note that the input and the mask images must be of the same size.
 *
 * \warning Any pixel value other than 0 will not be masked out.
 *
 * \sa MaskNegatedImageFilter
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Functor {

template<typename TInput, typename TMask, typename TOutput = TInput>
class MaskWithValue
{
public:
    MaskWithValue()
    : m_BackgroundValue(NumericTraits<TMask>::Zero),
      m_OutsideValue(NumericTraits<TOutput>::Zero)
    {
        // Nothing else
    }

    ~MaskWithValue()
    {
        // Nothing else
    }

    bool operator!=(MaskWithValue const & ) const
    {
        return false;
    }

    bool operator==(MaskWithValue const & other) const
    {
        return !(*this != other);
    }

    inline TOutput operator()(TInput const & A, TMask const & B) const
    {
        if (B != this->m_BackgroundValue)
        {
            return static_cast<TOutput>(A);
        }
        else
        {
            return m_OutsideValue;
        }
    }

    /** Method to explicitly set the background value of the mask */
    void SetBackgroundValue(TMask const & backgroundValue)
    {
        m_BackgroundValue = backgroundValue;
    }

    /** Method to get the outside value of the mask */
    const TMask & GetBackgroundValue() const
    {
        return m_BackgroundValue;
    }

    /** Method to explicitly set the outside value of the mask */
    void SetOutsideValue(TOutput const & outsideValue)
    {
        m_OutsideValue = outsideValue;
    }

    /** Method to get the outside value of the mask */
    const TOutput & GetOutsideValue() const
    {
        return m_OutsideValue;
    }

private:
    TMask m_BackgroundValue;
    TOutput m_OutsideValue;
}; // class MaskWithValue

} // namespace Functor

template <typename TInputImage,
          typename TMaskImage,
          typename TOutputImage=TInputImage>
class MaskWithValueImageFilter :
    public BinaryFunctorImageFilter<
         TInputImage,TMaskImage,TOutputImage,
         Functor::MaskWithValue<typename TInputImage::PixelType,
                                typename TMaskImage::PixelType,
                                typename TOutputImage::PixelType> >
{
public:
    /** Standard class typedefs. */
    typedef MaskWithValueImageFilter Self;
    typedef BinaryFunctorImageFilter<
        TInputImage,TMaskImage,TOutputImage,
        Functor::MaskWithValue<typename TInputImage::PixelType,
                               typename TMaskImage::PixelType,
                               typename TOutputImage::PixelType> > Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const>  ConstPointer;

    typedef typename TOutputImage::PixelType OutputPixelType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Runtime information support. */
    itkTypeMacro(MaskWithValueImageFilter, BinaryFunctorImageFilter);

    /** Set the background value of the mask. Defaults to 0 */
    void SetBackgroundValue(OutputPixelType const & backgroundValue)
    {
        if(this->GetBackgroundValue() != backgroundValue)
        {
            this->Modified();
            this->GetFunctor().SetBackgroundValue(backgroundValue);
        }
    }

    OutputPixelType const & GetBackgroundValue() const
    {
        return this->GetFunctor().GetBackgroundValue();
    }

    /** Set the outside value of the mask. Defaults to 0 */
    void SetOutsideValue(OutputPixelType const & outsideValue)
    {
        if( this->GetOutsideValue() != outsideValue)
        {
            this->Modified();
            this->GetFunctor().SetOutsideValue(outsideValue);
        }
    }

    OutputPixelType const & GetOutsideValue() const
    {
        return this->GetFunctor().GetOutsideValue();
    }

#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro(MaskEqualityComparableCheck,
        (Concept::EqualityComparable<typename TMaskImage::PixelType>));
    itkConceptMacro(InputConvertibleToOutputCheck,
        (Concept::Convertible<typename TInputImage::PixelType,
                              typename TOutputImage::PixelType>));
    /** End concept checking */
#endif

protected:
    MaskWithValueImageFilter() {}
    virtual ~MaskWithValueImageFilter() {}

    void PrintSelf(std::ostream &os, Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "OutsideValue: "  << this->GetOutsideValue() << std::endl;
    }

private:
    MaskWithValueImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
}; // class MaskWithValueImageFilter

} // namespace itk


#endif // medipy_components_segmentation_itkmaskwithvalueimagefilter_h
