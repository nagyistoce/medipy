/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2012
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _2b3ff1a0_c368_4090_8d5a_c924cc448314
#define _2b3ff1a0_c368_4090_8d5a_c924cc448314

#include <itkObject.h>
#include <itkObjectFactory.h>

#include "itkPlane.h"

namespace itk
{

/**
 * @brief Compute the best symmetry plane of an image.
 * 
 * Implementation of 
 *   "A new symmetry-based method for mid-sagittal plane extraction in neuroimages"
 *   Ruppert, G.; Teverovskiy, L.; Yu, C.; Falcao A.; Liu, Y.
 *   ISBI 2011, p. 285--288
 *   DOI: 10.1109/ISBI.2011.5872407
 */
template<typename TInputImage, typename TCoordRep=double>
class SymmetryPlaneCalculator : public Object
{
public:
    /** Standard class typedefs. */
    typedef SymmetryPlaneCalculator Self;
    typedef Object Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SymmetryPlaneCalculator, Object);

    /** Type definition for the input image. */
    typedef TInputImage ImageType;

    /** Pointer type for the image. */
    typedef typename TInputImage::Pointer ImagePointer;

    /** Const Pointer type for the image. */
    typedef typename TInputImage::ConstPointer ImageConstPointer;
    
    typedef typename TInputImage::PixelType PixelType;
    
    /** Plane type */
    typedef Plane<TCoordRep> PlaneType;
    
    /** Get the input image. */
    itkGetConstObjectMacro(Image, ImageType);
    
    /** Set the input image. */
    itkSetObjectMacro(Image, ImageType);
    
    /** Get the mid-sagittal plane of the input image. */
    PlaneType const & GetPlane() const;
    
    /** Compute the mid-sagittal plane of the input image. */
    void Compute();
    
protected:
    SymmetryPlaneCalculator();
    virtual ~SymmetryPlaneCalculator() {};
    void PrintSelf(std::ostream& os, Indent indent) const;

private:
    ImagePointer m_Image;
    PlaneType m_Plane;
    
    SymmetryPlaneCalculator(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
    /** 
     * @brief Compute the best symmetry plane of the input image.
     * @param sizeFactor : sub-sampling factor of the input image.
     * @param voxelStep : step for the search.
     * @param initial_plane : plane at which to start the search.
     * @param range : number of voxels on each side of each point of the plane 
     *                to limit the search to. If set to 0, compute metric for
     *                all possible planes.
     */
    PlaneType get_best_plane(int sizeFactor, float voxelStep, 
        PlaneType const & initial_plane=PlaneType(), float range=0) const;
    
    /** @brief Compute a quantile of the cummulated histogram of the image. */
    static PixelType get_histogram_quantile(TInputImage * image, float quantile);
    
    /** 
     * @brief Compute the symmetry metric of an image given a plane.
     * @param optimize : use the optimization for binary images given in the paper.
     */
    template<typename TMetricImage>
    static float get_metric(TMetricImage* image, PlaneType const & plane, bool optimize);
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSymmetryPlaneCalculator.txx"
#endif // ITK_MANUAL_INSTANTIATION

#endif // _2b3ff1a0_c368_4090_8d5a_c924cc448314
