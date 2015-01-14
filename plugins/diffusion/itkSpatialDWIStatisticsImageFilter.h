/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _cc50c230_e3f3_47e4_a1dd_e83b112a5dc7
#define _cc50c230_e3f3_47e4_a1dd_e83b112a5dc7

#include <ostream>
#include <itkSmartPointer.h>
#include "itkDWIStatisticsImageFilter.h"

namespace itk
{

/**
 * @brief Estimate the mean and standard deviation of a DTI image using a 
 *        spatial neighborhood.
 *
 * @pre TTensorImage must be a VectorImage.
 *  
 */
template<typename TTensorImage, typename TMeanImage=TTensorImage, 
    typename TStandardDeviationImage=itk::Image<
        typename TTensorImage::PixelType::ValueType, TTensorImage::ImageDimension>,
    typename TMaskImage=TStandardDeviationImage>
class SpatialDWIStatisticsImageFilter : public DWIStatisticsImageFilter<
    TTensorImage, TMeanImage, TStandardDeviationImage, TMaskImage>
{
public :
    /// @brief Standard class typedefs.
    typedef SpatialDWIStatisticsImageFilter Self;
    typedef DWIStatisticsImageFilter<
        TTensorImage, TMeanImage, TStandardDeviationImage, TMaskImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /// @brief Run-time type information (and related methods).
    itkTypeMacro(SpatialDWIStatisticsImageFilter, DWIStatisticsImageFilter);

    /// @brief Return the size of the neighborhood in the plane direction.
    itkGetConstMacro(SizePlane, unsigned int);
    
    /// @brief Set the size of the neighborhood in the plane direction, default to 3.
    itkSetMacro(SizePlane, unsigned int);
    
    /// @brief Return the size of the neighborhood in the normal direction.
    itkGetConstMacro(SizeDepth, unsigned int);
    
    /// @brief Set the size of the neighborhood in the normal direction, default to 3.
    itkSetMacro(SizeDepth, unsigned int);

protected :
    typedef typename TStandardDeviationImage::RegionType OutputImageRegionType;
    
    SpatialDWIStatisticsImageFilter();
    ~SpatialDWIStatisticsImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void BeforeThreadedGenerateData();
    void ThreadedGenerateData(OutputImageRegionType const & outputRegionForThread, ThreadIdType);

private :
    unsigned int m_SizePlane;
    unsigned int m_SizeDepth;
    
    SpatialDWIStatisticsImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatialDWIStatisticsImageFilter.txx"
#endif

#endif // _cc50c230_e3f3_47e4_a1dd_e83b112a5dc7
