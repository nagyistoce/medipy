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

namespace itk
{

template<typename TCoordRep>
Plane<TCoordRep>
::Plane()
{
    this->ComputeOriginAndNormal();
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
    this->m_P1 = other.GetP1();
    this->m_P2 = other.GetP2();
    this->m_P3 = other.GetP3();
    this->ComputeOriginAndNormal();
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
    // Get an axis from the plane : compute cross products between X, Y, Z axes
    // and the normal, find the "best" axis using the maximum absolute value of
    // the dot products.
    VectorType x_cross; 
    x_cross[0] = 0; x_cross[1] = -this->m_Normal[2]; x_cross[2] = this->m_Normal[1];
    VectorType y_cross; 
    y_cross[0] = this->m_Normal[2]; y_cross[1] =  0; y_cross[2] = -this->m_Normal[0];
    VectorType z_cross; 
    z_cross[0] = -this->m_Normal[1]; z_cross[1] = this->m_Normal[0]; z_cross[2] = 0;

    VectorType x; x[0] = 1; x[1] = 0; x[2] = 0;
    VectorType y; y[0] = 0; y[1] = 1; y[2] = 0;
    VectorType z; z[0] = 0; z[1] = 0; z[2] = 1;

    TCoordRep const x_dot = std::abs(x*this->m_Normal);
    TCoordRep const y_dot = std::abs(y*this->m_Normal);
    TCoordRep const z_dot = std::abs(z*this->m_Normal);
    VectorType v1;
    if(x_dot > y_dot)
    {
        if(x_dot > z_dot)
        {
            v1 = x_cross;
        }
        else
        {
            v1 = z_cross;
        }
    }
    else
    {
        if(y_dot > z_dot)
        {
            v1 = y_cross;
        }
        else
        {
            v1 = z_cross;
        }
    }
    
    // Get a second axis : cross product of normal and v1
    VectorType v2;
    v2[0] = this->m_Normal[1]*v1[2]-this->m_Normal[2]*v1[1];
    v2[1] = this->m_Normal[2]*v1[0]-this->m_Normal[0]*v1[2];
    v2[2] = this->m_Normal[0]*v1[1]-this->m_Normal[1]*v1[0];
    
    this->m_P1 = this->m_Origin;
    this->m_P2 = this->m_Origin+v1;
    this->m_P3 = this->m_Origin+v2;
}

template<typename TCoordRep>
std::ostream & operator<<(std::ostream & os, Plane<TCoordRep> const & plane)
{
    plane.Print(os);
}

} // namespace itk

#endif // _984c06e5_e85e_4006_97b7_cf033e6f1274

