/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _06618cc8_2f01_4f27_83fd_a9ffbebb12c0
#define _06618cc8_2f01_4f27_83fd_a9ffbebb12c0

#include <itkUnaryFunctorImageFilter.h>

#include "itkDTIScalarCalculator.h"

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class RadialDiffusivityImageFilter: public UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    RadialDiffusivityCalculator<
        typename TInputImage::PixelType, typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef RadialDiffusivityImageFilter  Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage, 
    RadialDiffusivityCalculator< 
        typename TInputImage::PixelType, typename TOutputImage::PixelType> > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<Self const> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(RadialDiffusivityImageFilter, UnaryFunctorImageFilter);
protected:
  RadialDiffusivityImageFilter() {}
  virtual ~RadialDiffusivityImageFilter() {}

private:
  RadialDiffusivityImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

#endif // _06618cc8_2f01_4f27_83fd_a9ffbebb12c0
