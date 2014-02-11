/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _e0569d02_a881_4c6a_ab72_ae5d057aedeb
#define _e0569d02_a881_4c6a_ab72_ae5d057aedeb

#include <itkUnaryFunctorImageFilter.h>

#include "itkDTIScalarCalculator.h"

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class AxialDiffusivityImageFilter: public UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    AxialDiffusivityCalculator<
        typename TInputImage::PixelType, typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef AxialDiffusivityImageFilter  Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage, 
    AxialDiffusivityCalculator< 
        typename TInputImage::PixelType, typename TOutputImage::PixelType> > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<Self const> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(AxialDiffusivityImageFilter, UnaryFunctorImageFilter);
protected:
  AxialDiffusivityImageFilter() {}
  virtual ~AxialDiffusivityImageFilter() {}

private:
  AxialDiffusivityImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

#endif // _e0569d02_a881_4c6a_ab72_ae5d057aedeb
