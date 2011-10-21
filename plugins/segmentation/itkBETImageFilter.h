/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef components_segmentation_betimagefilter_h
#define components_segmentation_betimagefilter_h

#include <vector>

#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector_fixed_ref.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace itk
{

/**
 * \class BETImageFilter
 * \brief Implementation of Smith's Brain Extraction Tool
 *
 * This implementation follows the details given in this paper :
 * @TECHREPORT{SmithBET,
 *     author = Smith, Stevens M.,
 *     title = BET : Brain Extraction Tool,
 *     institution = FMRIB (Oxford Center for Functional Magnetic Resonance Imaging of the Brain),
 *     number = TR00SMS2b,
 *     url = http://www.fmrib.ox.ac.uk/analysis/research/bet/resolveuid/33c59904365e93773e533a5eb8e3b7d9
 * }
 */
template<typename TInputImage, typename TOutputImage>
class BETImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef BETImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BETImageFilter, ImageToImageFilter);

    /** Type for input image. */
    typedef TInputImage InputImageType;
    typedef typename Superclass::InputImagePixelType InputImagePixelType;

    /** Type for the output image. */
    typedef TOutputImage OutputImageType;

    /** Main BET parameter (b_t). Must be between 0 and 1. */
    itkGetMacro(BT, float);
    itkSetClampMacro(BT, float, 0., 1.);

    /** Center of gravity of the brain. Will be computed if user does not provide it. */
    itkGetMacro(CenterOfGravity, typename TInputImage::PointType);
    virtual void SetCenterOfGravity(typename TInputImage::PointType & _arg)
    {
        itkDebugMacro("setting CenterOfGravity to " << _arg);
        if(this->m_CenterOfGravity != _arg)
        {
            this->m_CenterOfGravity = _arg;
            this->userDefinedCOG_ = true;
            this->Modified();
        }
    }


protected :
    BETImageFilter();
    ~BETImageFilter() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const;
    void GenerateData();

private :

    struct ThreadStruct
    {
        Self* filter;
        vtkDataArray* pointsDataArray;
        vtkDataArray* normalsDataArray;
        float E;
        float F;
        float d1;
        float d2;
        std::vector<vnl_vector_fixed<float, 3> >* displacements;
    };

    struct MeanVertexDistanceThreadStruct
    {
        vtkPolyData* model;
        std::vector<std::vector<vtkIdType> >* neighborhood;
        float* distances;
    };

    float m_BT;

    typename TInputImage::PointType m_CenterOfGravity;

    bool userDefinedCOG_;

    typename TInputImage::PixelType t2_;
    typename TInputImage::PixelType t98_;
    typename TInputImage::PixelType t_;
    float r_;
    typename TInputImage::PixelType tm_;
    vtkSmartPointer<vtkPolyData> sphere_;
    std::vector<std::vector<vtkIdType> > neighborhood_;
    float l_;

    BETImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    void thresholds();
    void centerOfGravity();
    void radius();
    void medianIntensity();
    void buildSphere(unsigned int subdivisions);
    std::vector<vtkIdType> neighbors(vtkIdType pointId);
    void meanVertexDistance();
    void voxelize();

    static ITK_THREAD_RETURN_TYPE mean_vertex_distance_callback(void* data);
    static float mean_vertex_distance_thread(vtkPolyData* model,
        std::vector<std::vector<vtkIdType> > const & neighborhood,
        vtkIdType points_begin, vtkIdType points_end);

    static ITK_THREAD_RETURN_TYPE adjust_model_callback(void* data);
    static void adjust_model(Self* filter,
        float E, float F, float d1, float d2,
        vtkDataArray* pointsDataArray, vtkDataArray* normalsDataArray,
        vtkIdType points_begin, vtkIdType points_end,
        std::vector<vnl_vector_fixed<float, 3> > & displacements);

    void
    drawSegment(vnl_vector_fixed_ref<float, TOutputImage::ImageDimension> const & p1,
                vnl_vector_fixed_ref<float, TOutputImage::ImageDimension> const & p2,
                float resolution);

    void
    drawTriangle(vnl_vector_fixed_ref<float, TOutputImage::ImageDimension> const & p1,
                 vnl_vector_fixed_ref<float, TOutputImage::ImageDimension> const & p2,
                 vnl_vector_fixed_ref<float, TOutputImage::ImageDimension> const & p3,
                 float resolution);
    void drawMesh();

    void fill(typename TOutputImage::IndexType const & seed);
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBETImageFilter.txx"
#endif

#endif // components_segmentation_betimagefilter_h
