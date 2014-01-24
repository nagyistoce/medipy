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
    MatrixType matrix;
    matrix[0][0] = input[0]; matrix[0][1] = input[1]; matrix[0][2] = input[2];
    matrix[1][0] = input[1]; matrix[1][1] = input[3]; matrix[1][2] = input[4];
    matrix[2][0] = input[2]; matrix[2][1] = input[4]; matrix[2][2] = input[5];

    MatrixType eigenVectors;
    VectorType eigenValues;
    this->m_Calculator.ComputeEigenValuesAndVectors(
        matrix, eigenValues, eigenVectors);
    
    return static_cast<TOutput>(this->ComputeScalar(eigenVectors, eigenValues));
}

template<typename TInput, typename TOutput>
double
AxialDiffusivityCalculator<TInput, TOutput>
::ComputeScalar(MatrixType const & eigenVectors, VectorType const & eigenValues) const
{
    // Eigenvalues are sorted in ascending order
    return eigenValues[2];
}

template<typename TInput, typename TOutput>
double
FractionalAnisotropyCalculator<TInput, TOutput>
::ComputeScalar(MatrixType const & eigenVectors, VectorType const & eigenValues) const
{
    double fa = 0.5 * (
        std::pow(eigenValues[0]-eigenValues[1], 2) + 
        std::pow(eigenValues[1]-eigenValues[2], 2) + 
        std::pow(eigenValues[2]-eigenValues[0], 2)
    ) / (
        std::pow(eigenValues[0], 2)+
        std::pow(eigenValues[1], 2)+
        std::pow(eigenValues[2], 2));
    if(fa>0)
    { 
        fa = sqrt(fa); 
    }
    else
    {
        fa = 0;
    }
    
    return fa;
}

template<typename TInput, typename TOutput>
double
MeanDiffusivityCalculator<TInput, TOutput>
::ComputeScalar(MatrixType const & eigenVectors, VectorType const & eigenValues) const
{
    return (eigenValues[0]+eigenValues[1]+eigenValues[2])/3.0;
}

template<typename TInput, typename TOutput>
double
RadialDiffusivityCalculator<TInput, TOutput>
::ComputeScalar(MatrixType const & eigenVectors, VectorType const & eigenValues) const
{
    // Eigenvalues are sorted in ascending order
    return 0.5*(eigenValues[0]+eigenValues[1]);
}

}

#endif // _28d89b95_c219_46b1_81bf_b09c8e7a8441
