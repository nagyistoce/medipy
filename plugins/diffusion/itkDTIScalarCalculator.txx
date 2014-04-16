/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/
 
#ifndef _28d89b95_c219_46b1_81bf_b09c8e7a8441
#define _28d89b95_c219_46b1_81bf_b09c8e7a8441

#include "itkDTIScalarCalculator.h"

namespace itk
{

template<typename TInput, typename TOutput>
DTIScalarCalculator<TInput, TOutput>
::DTIScalarCalculator() 
{ 
    this->m_Calculator.SetDimension(3); 
    this->m_Calculator.SetOrderEigenValues(true);
}

template<typename TInput, typename TOutput>
bool
DTIScalarCalculator<TInput, TOutput>
::operator==(DTIScalarCalculator const & other) const
{
    return (
        this->m_Calculator.GetOrder() == other.m_Calculator.GetOrder() &&
        this->m_Calculator.GetOrderEigenValues() == other.m_Calculator.GetOrderEigenValues() &&
        this->m_Calculator.GetOrderEigenMagnitudes() == other.m_Calculator.GetOrderEigenMagnitudes()
    );
}

template<typename TInput, typename TOutput>
bool
DTIScalarCalculator<TInput, TOutput>
::operator!=(DTIScalarCalculator const & other) const
{
    return !(*this == other);
}

template<typename TInput, typename TOutput>
TOutput
DTIScalarCalculator<TInput, TOutput>
::operator()(TInput const & input) const
{
    return this->ComputeScalar(input);
}

template<typename TInput, typename TOutput>
void
DTIScalarCalculator<TInput, TOutput>
::ComputeEigensystem(TInput const & input, 
                     MatrixType & eigenVectors, VectorType & eigenValues) const
{
    MatrixType matrix;
    matrix[0][0] = input[0]; matrix[0][1] = input[1]; matrix[0][2] = input[2];
    matrix[1][0] = input[1]; matrix[1][1] = input[3]; matrix[1][2] = input[4];
    matrix[2][0] = input[2]; matrix[2][1] = input[4]; matrix[2][2] = input[5];

    this->m_Calculator.ComputeEigenValuesAndVectors(
        matrix, eigenValues, eigenVectors);
}

template<typename TInput, typename TOutput>
double 
DTIScalarCalculator<TInput, TOutput>
::ComputeNorm(TInput const & input) const
{
    double const xx = input[0];
    double const xy = input[1];
    double const xz = input[2];
    double const yy = input[3];
    double const yz = input[4];
    double const zz = input[5];

    return (xx*xx + yy*yy + zz*zz + 2.0*(xy*xy + xz*xz + yz*yz));
}

template<typename TInput, typename TOutput>
double 
DTIScalarCalculator<TInput, TOutput>
::ComputeTrace(TInput const & input) const
{
    double const xx = input[0];
    double const yy = input[3];
    double const zz = input[5];
    
    return xx+yy+zz;
}

template<typename TInput, typename TOutput>
TOutput
AxialDiffusivityCalculator<TInput, TOutput>
::ComputeScalar(TInput const & input) const
{
    MatrixType eigenVectors;
    VectorType eigenValues;
    this->ComputeEigensystem(input, eigenVectors, eigenValues);
    // Eigenvalues are sorted in ascending order
    return static_cast<TOutput>(eigenValues[2]);
}

template<typename TInput, typename TOutput>
TOutput
FractionalAnisotropyCalculator<TInput, TOutput>
::ComputeScalar(TInput const & input) const
{
    // From equation 28 in
    // http://lmi.bwh.harvard.edu/papers/pdfs/2002/westinMEDIA02.pdf
    
    double result = 0.0;
    
    double const norm = this->ComputeNorm(input);
    if(norm>0.0)
    {
        double const trace = this->ComputeTrace(input);
        double const numerator = 3.0*norm-trace*trace;
        if(numerator>0.0)
        {
            result = vcl_sqrt(numerator/(2*norm));
        }
    }
    
    return static_cast<TOutput>(result);
}

template<typename TInput, typename TOutput>
TOutput
MeanDiffusivityCalculator<TInput, TOutput>
::ComputeScalar(TInput const & input) const
{
    double const trace = this->ComputeTrace(input);
    return static_cast<TOutput>(trace/3.0);
}

template<typename TInput, typename TOutput>
TOutput
RadialDiffusivityCalculator<TInput, TOutput>
::ComputeScalar(TInput const & input) const
{
    MatrixType eigenVectors;
    VectorType eigenValues;
    this->ComputeEigensystem(input, eigenVectors, eigenValues);
    // Eigenvalues are sorted in ascending order
    return 0.5*(eigenValues[0]+eigenValues[1]);
}

}

#endif // _28d89b95_c219_46b1_81bf_b09c8e7a8441
