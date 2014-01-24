/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7616bbcf_5945_4091_b9fc_9c7eb753dd87
#define _7616bbcf_5945_4091_b9fc_9c7eb753dd87

#include <itkUnaryFunctorImageFilter.h>

#include "itkDTIScalarCalculator.h"

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class FractionalAnisotropyImageFilter: public UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    FractionalAnisotropyCalculator<
        typename TInputImage::PixelType, typename TOutputImage::PixelType> >
{
public:
  /** Standard class typedefs. */
  typedef FractionalAnisotropyImageFilter  Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage, 
    FractionalAnisotropyCalculator< 
        typename TInputImage::PixelType, typename TOutputImage::PixelType> > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<Self const> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(FractionalAnisotropyImageFilter, UnaryFunctorImageFilter);
protected:
  FractionalAnisotropyImageFilter() {}
  virtual ~FractionalAnisotropyImageFilter() {}

private:
  FractionalAnisotropyImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

}

#endif // _7616bbcf_5945_4091_b9fc_9c7eb753dd87
