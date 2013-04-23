/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/
 
#ifndef _984c06e5_e85e_4006_97b7_cf033e6f1274
#define _984c06e5_e85e_4006_97b7_cf033e6f1274

#include "itkPlane.h"

#include <ostream>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkIndent.h>
#include <vnl/vnl_cross.h>
#include <vnl/vnl_vector.h>

namespace itk
{

template<typename TCoordRep>
Plane<TCoordRep>
::Plane()
{
    this->m_Origin.Fill(0);
    this->m_Normal.Fill(0);
    this->m_Normal[2]=1;
    this->m_UnitNormal = this->m_Normal;
    
    this->Compute3Points();
}

template<typename TCoordRep>
Plane<TCoordRep>
::Plane(PointType const & P1, PointType const & P2, PointType const & P3)
: m_P1(P1), m_P2(P2), m_P3(P3)
{
    this->ComputeOriginAndNormal();
}

template<typename TCoordRep>
Plane<TCoordRep>
::Plane(PointType const & origin, VectorType const & normal)
: m_Origin(origin), m_Normal(normal)
{
    this->Compute3Points();
}

template<typename TCoordRep>
Plane<TCoordRep>
::Plane(Self const & other)
: m_P1(other.m_P1), m_P2(other.m_P2), m_P3(other.m_P3)
{
    this->ComputeOriginAndNormal();
}

template<typename TCoordRep>
typename Plane<TCoordRep>::Self const &
Plane<TCoordRep>
::operator=(Self const & other)
{
    if(this != &other)
    {
        this->m_P1 = other.GetP1();
        this->m_P2 = other.GetP2();
        this->m_P3 = other.GetP3();
        this->ComputeOriginAndNormal();
    }
    return *this;
}

template<typename TCoordRep>
Plane<TCoordRep>
::~Plane()
{
    // Nothing to do
}

template<typename TCoordRep>
typename Plane<TCoordRep>::PointType
Plane<TCoordRep>
::GetP1() const
{
    return this->m_P1;
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::SetP1(PointType const & P1)
{
    this->m_P1 = P1;
    this->ComputeOriginAndNormal();
}

template<typename TCoordRep>
typename Plane<TCoordRep>::PointType
Plane<TCoordRep>
::GetP2() const
{
    return this->m_P2;
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::SetP2(PointType const & P2)
{
    this->m_P2 = P2;
    this->ComputeOriginAndNormal();
}

template<typename TCoordRep>
typename Plane<TCoordRep>::PointType
Plane<TCoordRep>
::GetP3() const
{
    return this->m_P3;
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::SetP3(PointType const & P3)
{
    this->m_P3 = P3;
    this->ComputeOriginAndNormal();
}

template<typename TCoordRep>
typename Plane<TCoordRep>::PointType
Plane<TCoordRep>
::GetOrigin() const
{
    return this->m_Origin;
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::SetOrigin(PointType const & origin)
{
    this->m_Origin = origin;
    this->Compute3Points();
}

template<typename TCoordRep>
typename Plane<TCoordRep>::VectorType
Plane<TCoordRep>
::GetNormal() const
{
    return this->m_Normal;
}

template<typename TCoordRep>
typename Plane<TCoordRep>::VectorType
Plane<TCoordRep>
::GetUnitNormal() const
{
    return this->m_UnitNormal;
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::SetNormal(VectorType const & normal)
{
    this->m_Normal = normal;
    this->m_UnitNormal = this->m_Normal/this->m_Normal.GetNorm();
    this->Compute3Points();
}

template<typename TCoordRep>
typename Plane<TCoordRep>::PointType::CoordRepType
Plane<TCoordRep>
::GetDistance(PointType const & p) const
{
    TCoordRep const distance =
        this->m_UnitNormal[0]*(p[0]-this->m_P1[0])+
        this->m_UnitNormal[1]*(p[1]-this->m_P1[1])+
        this->m_UnitNormal[2]*(p[2]-this->m_P1[2]);

    return distance;
}

template<typename TCoordRep>
typename Plane<TCoordRep>::PointType
Plane<TCoordRep>
::Reflect(PointType const & p) const
{
    TCoordRep const distance = this->GetDistance(p);
    return p-2*distance*this->m_UnitNormal;
}

template<typename TCoordRep>
template<typename TImage>
void
Plane<TCoordRep>
::Draw(TImage* image, typename TImage::PixelType const & value, 
       double tolerance) const
{
    typedef ImageRegionIteratorWithIndex<TImage> Iterator;
    for(Iterator it(image, image->GetRequestedRegion()); !it.IsAtEnd(); ++it)
    {
        typename TImage::PointType point;
        image->TransformIndexToPhysicalPoint(it.GetIndex(), point);
        typename PointType::CoordRepType const distance = this->GetDistance(point);
        if(std::abs(distance) <= tolerance)
        {
            it.Set(value);
        }
    }
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::Print(std::ostream & os, Indent indent) const
{
    this->PrintSelf(os, indent);
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::PrintSelf(std::ostream & os, Indent indent) const
{
    os << indent << "P1: " << this->m_P1 << "\n";
    os << indent << "P2: " << this->m_P2 << "\n";
    os << indent << "P3: " << this->m_P3 << "\n";
    os << indent << "Origin: " << this->m_Origin << "\n";
    os << indent << "Unit normal: " << this->m_UnitNormal << "\n";
    os << indent << "Normal: " << this->m_Normal << "\n";
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::ComputeOriginAndNormal()
{
    this->m_Origin = this->m_P1;
    VectorType const v1 = this->m_P2-this->m_P1;
    VectorType const v2 = this->m_P3-this->m_P1;
    // Cross-product v1*v2
    this->m_Normal[0] = v1[1]*v2[2]-v1[2]*v2[1];
    this->m_Normal[1] = v1[2]*v2[0]-v1[0]*v2[2];
    this->m_Normal[2] = v1[0]*v2[1]-v1[1]*v2[0];
    this->m_UnitNormal = this->m_Normal/this->m_Normal.GetNorm();
}

template<typename TCoordRep>
void
Plane<TCoordRep>
::Compute3Points()
{
    typedef vnl_vector_fixed<TCoordRep, 3> VnlVector;
    
    // Compute the cross products between each axis (X, Y, Z) and the unit normal
    // to determine the longest one : this will be the first vector of the plane.
    // Using any axis at random might lead to degenerate plane vectors.
    VnlVector const unit_normal = this->m_UnitNormal.GetVnlVector();
    VnlVector v1;
    TCoordRep squared_magnitude=0;
    for(unsigned int d=0; d<3; ++d)
    {
        VnlVector axis(0.);
        axis[d] = 1;
        
        VnlVector cross_product = vnl_cross_3d(unit_normal, axis);
        TCoordRep const cross_product_squared_magnitude = cross_product.squared_magnitude();
        if(cross_product_squared_magnitude > squared_magnitude)
        {
            v1 = cross_product;
            squared_magnitude = cross_product_squared_magnitude;
        }
    }
    
    // Get the second plane axis : cross product of the unit normal and v1;
    VnlVector const v2 = vnl_cross_3d(unit_normal, v1);
    
    this->m_P1 = this->m_Origin;
    std::transform(this->m_Origin.Begin(), this->m_Origin.End(), v1.begin(), 
        this->m_P2.Begin(), std::plus<TCoordRep>());
    std::transform(this->m_Origin.Begin(), this->m_Origin.End(), v2.begin(), 
        this->m_P3.Begin(), std::plus<TCoordRep>());
}

template<typename TCoordRep>
std::ostream & operator<<(std::ostream & os, Plane<TCoordRep> const & plane)
{
    plane.Print(os);
    return os;
}

} // namespace itk

#endif // _984c06e5_e85e_4006_97b7_cf033e6f1274

