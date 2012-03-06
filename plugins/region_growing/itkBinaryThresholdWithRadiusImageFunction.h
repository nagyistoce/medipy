#ifndef ae1871ab_34e0_4e67_943e_a40c4df993a7
#define ae1871ab_34e0_4e67_943e_a40c4df993a7

#include <vector>
#include <itkImageFunction.h>

namespace itk
{

template <typename TInputImage, typename TCoordRep=float>
class ITK_EXPORT BinaryThresholdWithRadiusImageFunction :
    public ImageFunction<TInputImage, bool, TCoordRep>
{
public :
    /** Standard class typedefs. */
    typedef BinaryThresholdWithRadiusImageFunction Self;
    typedef ImageFunction<TInputImage, bool, TCoordRep> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods).  */
    itkTypeMacro(BinaryThresholdWithRadiusImageFunction, ImageFunction);

    /** Superclass definitions. */
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::InputPixelType InputPixelType;
    typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
    typedef typename Superclass::OutputType OutputType;
    typedef typename Superclass::CoordRepType CoordRepType;
    typedef typename Superclass::IndexType IndexType;
    typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
    typedef typename Superclass::PointType PointType;

    virtual bool Evaluate(PointType const & point) const;
    virtual bool EvaluateAtContinuousIndex(ContinuousIndexType const & continuous_index) const;
    virtual bool EvaluateAtIndex(IndexType const & index) const;

    itkGetConstReferenceMacro(Lower, InputPixelType);
    itkGetConstReferenceMacro(Upper, InputPixelType);

    itkGetMacro(Radius, float);
    void SetRadius(float radius);

    void SetSeed(IndexType const & seed);
    void AddSeed(IndexType const & seed);
    void ClearSeeds();

    /** Values greater than or equal to the value are inside. */
    void ThresholdAbove(InputPixelType threshold);
    /** Values less than or equal to the value are inside. */
    void ThresholdBelow(InputPixelType threshold);
    /** Values that lie between lower and upper inclusive are inside. */
    void ThresholdBetween(InputPixelType lower, InputPixelType upper);
protected:
    BinaryThresholdWithRadiusImageFunction();
    ~BinaryThresholdWithRadiusImageFunction();
    void PrintSelf(std::ostream& os, Indent indent) const;

    InputPixelType m_Lower;
    InputPixelType m_Upper;
    float m_Radius;
    float m_RadiusSquared;
    std::vector<IndexType> m_Seeds;

private:
  BinaryThresholdWithRadiusImageFunction(Self const &); //purposely not implemented
  void operator=(Self const &); //purposely not implemented
};

} // namespace itk

#include "itkBinaryThresholdWithRadiusImageFunction.txx"

#endif // ae1871ab_34e0_4e67_943e_a40c4df993a7
