/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef b0f7acb3_a658_426b_9de9_f752e64f9468
#define b0f7acb3_a658_426b_9de9_f752e64f9468

#include <itkImageToImageFilter.h>
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
 */
template <typename TInputImage,
          typename TMaskImage,
          typename TOutputImage=TInputImage>
class MaskWithValueImageFilter :
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef MaskWithValueImageFilter Self;

    typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Runtime information support. */
    itkTypeMacro(MaskWithValueImageFilter, ImageToImageFilter);

    /** Convenient typedefs. */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::OutputImageType OutputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;
    typedef typename Superclass::OutputImagePixelType OutputImagePixelType;

    typedef TMaskImage MaskImageType;
    typedef typename TMaskImage::PixelType MaskImagePixelType;

    MaskImageType const * GetMask();
    void SetMask(MaskImageType const * image);

    /** Background value of the mask. Defaults to 0 */
    itkGetConstMacro(BackgroundValue, MaskImagePixelType);
    itkSetMacro(BackgroundValue, MaskImagePixelType);

    /** Outside value. Defaults to 0 */
    itkGetConstMacro(OutsideValue, OutputImagePixelType);
    itkSetMacro(OutsideValue, OutputImagePixelType);

protected:
    MaskWithValueImageFilter();
    virtual ~MaskWithValueImageFilter();
    void PrintSelf(std::ostream &os, Indent indent) const;
    void ThreadedGenerateData(typename InputImageType::RegionType const & region, int thread_id);

private:
    MaskWithValueImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    MaskImagePixelType m_BackgroundValue;
    OutputImagePixelType m_OutsideValue;

}; // class MaskWithValueImageFilter

} // namespace itk

#include "itkMaskWithValueImageFilter.txx"

#endif // b0f7acb3_a658_426b_9de9_f752e64f9468
