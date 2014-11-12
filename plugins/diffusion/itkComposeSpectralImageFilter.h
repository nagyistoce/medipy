/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/


#ifndef _itkComposeSpectralImageFilter_h
#define _itkComposeSpectralImageFilter_h


#include <itkDataObject.h>
#include <itkProcessObject.h>
#include <itkSmartPointer.h>

namespace itk
{

/**
 * \class ComposeSpectralImageFilter
 * \brief Compute the tensor image from eigenvalues and corresponding eigenvectors
 *
 * This filter does not derived from ImageToImageFilter and from ImageSource (because
 * the InputImage and the OutputImage are VectorImage type)
 * Consequently, some codes are duplicated from ImageSource or ImageToImageFilter.
 * 
 */

template<typename TInputImage, typename TOutputImage>
class ComposeSpectralImageFilter : public ProcessObject
{
public :
    /** Standard class typedefs. */
    typedef ComposeSpectralImageFilter Self;
    typedef ProcessObject Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;
    
    /** Output */
    typedef TOutputImage                    OutputImageType;
    typedef typename OutputImageType::Pointer        OutputImagePointer;
    typedef typename OutputImageType::PixelType      OutputImagePixelType;
    
    /** Input */
    typedef TInputImage                    InputImageType;
    typedef typename InputImageType::Pointer        InputImagePointer;
    typedef typename InputImageType::PixelType      InputImagePixelType;

    typedef DataObject::Pointer DataObjectPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ComposeSpectralImageFilter, ProcessObject);
    
    /** Set/Get the image input of this filter. */
    void SetInput(unsigned int idx, InputImageType const * input);
    InputImageType const * GetInput();
    InputImageType const * GetInput(unsigned int idx);
    
    void SetEigenValues(InputImageType const * input);
    void SetEigenVectors(InputImageType const * input);
    virtual void Update();
    
    /** Get the images output of this filter*/
    DataObjectPointer MakeOutput(unsigned int idx);
    OutputImageType * GetOutput();
    OutputImageType const * GetOutput() const;
    /// @brief Return the eigenvalues image.
    InputImageType const * GetEigenValuesImage() const;
    /// @brief Return the eigenvectors image.
    InputImageType const * GetEigenVectorsImage() const;
    
protected :
    ComposeSpectralImageFilter();
    ~ComposeSpectralImageFilter() {}
    
    void PrintSelf(std::ostream& os, Indent indent) const;
    void AllocateOutputs();
    
    /// @brief Return the eigenvalues image.
    InputImageType * GetEigenValuesImage();
    /// @brief Return the eigenvectors image.
    InputImageType * GetEigenVectorsImage();
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkComposeSpectralImageFilter.txx"
#endif

#endif 

