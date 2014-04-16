/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _ecc3ccbf_a10d_4224_9fe8_2b0f9d73904c
#define _ecc3ccbf_a10d_4224_9fe8_2b0f9d73904c

#include <itkSymmetricEigenAnalysis.h>

namespace itk
{

/**
 * @brief Base class for functors that compute scalar values from diffusion tensors.
 */
template<typename TInput, typename TOutput>
class DTIScalarCalculator
{
public:
    DTIScalarCalculator();
    ~DTIScalarCalculator() {}
    
    bool operator==(DTIScalarCalculator const & other) const;
    bool operator!=(DTIScalarCalculator const & other) const;
    
    TOutput operator()(TInput const & input) const;

protected:
    typedef Matrix<double, 3, 3> MatrixType;
    typedef FixedArray<double, 3> VectorType;

    itk::SymmetricEigenAnalysis<MatrixType, VectorType> m_Calculator;
    
    /**
     * @brief Compute the eigensystem of diffusion tensor.
     */
    void ComputeEigensystem(TInput const & input, 
                            MatrixType & eigenVectors, VectorType & eigenValues) const;
    
    /**
     * @brief Compute the norm of the diffusion tensor.
     * 
     * The norm is defined as the square root of the sum of the squared elements
     * of the tensor, which equals the square root of the sum of the squared
     * eigenvalues (Processing and visualization for diffusion tensor MRI; 
     * Westin et al.; Medical Image Analysis vol. 6 pp.93-108; 2002).
     */
    double ComputeNorm(TInput const & input) const;
    
    /**
     * @brief Compute the trace of the diffusion tensor.
     * 
     * The norm is defined as the sum of the diagonal elements of the tensor, 
     * which equals the square root of the sum of the eigenvalues (Processing 
     * and visualization for diffusion tensor MRI; Westin et al.; Medical Image 
     * Analysis vol. 6 pp.93-108; 2002).
     */
    double ComputeTrace(TInput const & input) const;
    
    virtual TOutput ComputeScalar(TInput const & input) const=0;
};

template<typename TInput, typename TOutput>
class AxialDiffusivityCalculator: public DTIScalarCalculator<TInput, TOutput>
{
public:
    AxialDiffusivityCalculator() {}
    ~AxialDiffusivityCalculator() {}
protected:
    typedef typename DTIScalarCalculator<TInput, TOutput>::MatrixType MatrixType;
    typedef typename DTIScalarCalculator<TInput, TOutput>::VectorType VectorType;
    virtual TOutput ComputeScalar(TInput const & input) const;
};

template<typename TInput, typename TOutput>
class FractionalAnisotropyCalculator: public DTIScalarCalculator<TInput, TOutput>
{
public:
    FractionalAnisotropyCalculator() {}
    ~FractionalAnisotropyCalculator() {}
protected:
    typedef typename DTIScalarCalculator<TInput, TOutput>::MatrixType MatrixType;
    typedef typename DTIScalarCalculator<TInput, TOutput>::VectorType VectorType;
    virtual TOutput ComputeScalar(TInput const & input) const;
};

template<typename TInput, typename TOutput>
class MeanDiffusivityCalculator: public DTIScalarCalculator<TInput, TOutput>
{
public:
    MeanDiffusivityCalculator() {}
    ~MeanDiffusivityCalculator() {}
protected:
    typedef typename DTIScalarCalculator<TInput, TOutput>::MatrixType MatrixType;
    typedef typename DTIScalarCalculator<TInput, TOutput>::VectorType VectorType;
    virtual TOutput ComputeScalar(TInput const & input) const;
};

template<typename TInput, typename TOutput>
class RadialDiffusivityCalculator: public DTIScalarCalculator<TInput, TOutput>
{
public:
    RadialDiffusivityCalculator() {}
    ~RadialDiffusivityCalculator() {}
protected:
    typedef typename DTIScalarCalculator<TInput, TOutput>::MatrixType MatrixType;
    typedef typename DTIScalarCalculator<TInput, TOutput>::VectorType VectorType;
    virtual TOutput ComputeScalar(TInput const & input) const;
};

}

#include "itkDTIScalarCalculator.txx"

#endif // _ecc3ccbf_a10d_4224_9fe8_2b0f9d73904c
