/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _7444a613_2a7c_421b_af39_25b8ab1ce0e3
#define _7444a613_2a7c_421b_af39_25b8ab1ce0e3

#include <ostream>

#include <itkImage.h>
#include <itkIndent.h>
#include <itkMacro.h>
#include <itkObject.h>
#include <vtkPiecewiseFunction.h>
#include <vtkSmartPointer.h>

namespace itk
{

/**
 * @brief Compute a transfer function approximating a joint histogram.
 */
template<typename THistogram>
class JointHistogramTransferFunctionCalculator: public Object
{
public:
    /** Standard typedefs */
    typedef JointHistogramTransferFunctionCalculator Self;
    typedef Object Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(JointHistogramTransferFunctionCalculator, Object);

    /** standard New() method support */
    itkNewMacro(Self);
    
    typedef THistogram HistogramType;
    typedef typename HistogramType::Pointer HistogramPointer;
    typedef typename HistogramType::MeasurementType MeasurementType;
    
    typedef vtkPiecewiseFunction TransferFunctionType;
    typedef vtkSmartPointer<TransferFunctionType> TransferFunctionPointer;
    
    itkGetObjectMacro(Histogram, HistogramType);
    itkGetConstObjectMacro(Histogram, HistogramType);
    itkSetObjectMacro(Histogram, HistogramType);
    
    itkGetMacro(ResampleFactor, unsigned long)
    itkSetMacro(ResampleFactor, unsigned long);
    
    void Compute();
    
    TransferFunctionType const * GetTransferFunction() const;

protected:
    JointHistogramTransferFunctionCalculator();
    virtual ~JointHistogramTransferFunctionCalculator() {};
    void PrintSelf(std::ostream & os, Indent indent) const;

private:
    HistogramPointer m_Histogram;
    unsigned long m_ResampleFactor;
    TransferFunctionPointer m_TransferFunction;
    
    JointHistogramTransferFunctionCalculator(Self const &); //purposely not implemented
    void operator=(Self const &); //purposely not implemented
    
    TransferFunctionPointer MaximumProbabilityTransferFunction() const;
    std::vector<MeasurementType> ConfidenceWeights() const;
    
    std::vector<MeasurementType> Column(int const x) const;
    std::vector<MeasurementType> Resample(std::vector<MeasurementType> const & function) const;
    std::vector<MeasurementType> Smooth(std::vector<MeasurementType> const & function, int radius) const;
    TransferFunctionPointer Smooth(TransferFunctionPointer function, int factor) const;
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkJointHistogramTransferFunctionCalculator.txx"
#endif

#endif // _7444a613_2a7c_421b_af39_25b8ab1ce0e3
