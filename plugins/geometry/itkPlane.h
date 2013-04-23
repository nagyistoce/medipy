/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/
 
#ifndef _e0bbd6a0_294b_46f1_bba8_950e538493f8
#define _e0bbd6a0_294b_46f1_bba8_950e538493f8

#include <itkPoint.h>

namespace itk
{

/**
 * @brief A 3D plane.
 * 
 * A 3D plane can be specified either by three points or by an origin and a 
 * normal.
 */
template<typename TCoordRep=float>
class Plane
{
public :

    /** Standard typedefs. */
    typedef Plane Self;

    typedef Point<TCoordRep, 3> PointType;
    typedef typename PointType::VectorType VectorType;

    /** @brief Create an un-initialized plane. */
    Plane();
    
    /** @brief Create a plane specified by three points. */
    Plane(PointType const & P1, PointType const & P2, PointType const & P3);
    
    /** @brief Create a plane specified by its origin and normal. */
    Plane(PointType const & origin, VectorType const & normal);
    
    /** @brief Copy constructor. */
    Plane(Self const & other);
    
    /** @brief Assignment operator. */
    Self const & operator=(Self const & other);
    
    /** @brief Destructor. */
    virtual ~Plane();

    /** @brief Return the first point defining the plane. */
    PointType const & GetP1() const;
    
    /** @brief Set the first point defining the plane. */
    void SetP1(PointType const & P1);

    /** @brief Return the second point defining the plane. */
    PointType const & GetP2() const;
    
    /** @brief Set the second point defining the plane. */
    void SetP2(PointType const & P2);
    
    /** @brief Return the third point defining the plane. */
    PointType const & GetP3() const;
    
    /** @brief Set the third point defining the plane. */
    void SetP3(PointType const & P3);

    /** @brief Return the origin of the plane. */
    PointType const & GetOrigin() const;
    
    /** @brief Set the origin of the plane. */
    void SetOrigin(PointType const & origin);

    /** @brief Return the normal to the plane. */
    VectorType const & GetNormal() const;
    
    /** @brief Return the unit normal to the plane. */
    VectorType const & GetUnitNormal() const;
    
    /** @brief Set the normal to the plane. */
    void SetNormal(VectorType const & normal);

    /** @brief Return the signed distance from a point to the plane. */
    typename PointType::CoordRepType GetDistance(PointType const & p) const;
    
    /** @brief Return the reflection of a point with respect to the plane. */
    PointType Reflect(PointType const & p) const;

    /** 
     * @brief Draw the plane on an image.
     * @param image : image to draw on.
     * @param value : value taken by voxels which lie on the plane.
     * @param tolerance : maximum distance for a voxel to be considered within the plane.
     */
    template<typename TImage>
    void Draw(TImage* image, typename TImage::PixelType const & value, 
              double tolerance=0.5) const;

    void Print(std::ostream & os, Indent indent=0) const;

protected :
    void PrintSelf(std::ostream & os, Indent indent) const;

private :
    PointType m_P1;
    PointType m_P2;
    PointType m_P3;

    PointType m_Origin;
    VectorType m_Normal;
    VectorType m_UnitNormal;

    /** @brief Update the origin and normal using P1, P2, and P3. */
    void ComputeOriginAndNormal();
    
    /** @brief Update P1, P1, and P3 usig the origin and the normal. */
    void Compute3Points();
};

template<typename TCoordRep>
std::ostream & operator<<(std::ostream & os, Plane<TCoordRep> const & plane);

} // namespace itk

#include "itkPlane.txx"

#endif // _e0bbd6a0_294b_46f1_bba8_950e538493f8
