/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _c9d282c1_9368_4820_ad56_9922bda995b7
#define _c9d282c1_9368_4820_ad56_9922bda995b7

#include <itkUnaryFunctorImageFilter.h>

#include "itkDTIScalarCalculator.h"

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class MeanDiffusivityImageFilter: public UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    MeanDiffusivityCalculator<
        typename TInputImage::PixelType, typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef MeanDiffusivityImageFilter  Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage, 
    MeanDiffusivityCalculator< 
        typename TInputImage::PixelType, typename TOutputImage::PixelType> > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<Self const> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(MeanDiffusivityImageFilter, UnaryFunctorImageFilter);
protected:
  MeanDiffusivityImageFilter() {}
  virtual ~MeanDiffusivityImageFilter() {}

private:
  MeanDiffusivityImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

#endif // _c9d282c1_9368_4820_ad56_9922bda995b7
