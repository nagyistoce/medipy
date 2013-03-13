#ifndef e5a5528b_6847_4f64_b1b7_76a442d67604
#define e5a5528b_6847_4f64_b1b7_76a442d67604

#include <map>
#include <vector>

#include <itkObject.h>

namespace itk
{

template<typename TImage>
class ClustersToAnnotationsCalculator : public Object
{
public :
    /** Standard class typedefs. */
    typedef ClustersToAnnotationsCalculator Self;
    typedef Object Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ClustersToAnnotationsCalculator, Object);

    typedef TImage ImageType;
    typedef typename TImage::Pointer ImagePointer;
    typedef typename TImage::ConstPointer ImageConstPointer;
    typedef typename TImage::PixelType PixelType;
    typedef typename TImage::IndexType IndexType;
    typedef typename TImage::RegionType RegionType;

    struct Annotation
    {
        IndexType position;
        float size;
    };

    typedef std::map<PixelType, Annotation> MapType;

    itkSetConstObjectMacro(Image, ImageType);

    void Compute();

    MapType const & GetAnnotations() const;
    std::vector<typename MapType::key_type> GetAnnotationsLabels() const;
    IndexType GetAnnotationPosition(PixelType label) const;
    float GetAnnotationSize(PixelType label) const;

protected:
    ClustersToAnnotationsCalculator();
    virtual ~ClustersToAnnotationsCalculator();

private:
    ClustersToAnnotationsCalculator(Self const &); //purposely not implemented
    void operator=(Self const &); //purposely not implemented

    ImageConstPointer m_Image;
    MapType m_Annotations;
};

}

#include "itkClustersToAnnotationsCalculator.txx"

#endif // e5a5528b_6847_4f64_b1b7_76a442d67604
