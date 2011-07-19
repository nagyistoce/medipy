#ifndef itkAssembleTilesImageFilter_h
#define itkAssembleTilesImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkSmartPointer.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
class AssembleTilesImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    typedef AssembleTilesImageFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(AssembleTilesImageFilter, ImageToImageFilter);

    /** Type for input image. */
    typedef TInputImage InputImageType;

    /** Type for the output image. */
    typedef TOutputImage OutputImageType;

    enum EmptyTilesPosition
    {
        BEGIN,
        END
    };

    enum AssemblyOrder
    {
        BOTTOM_TO_TOP,
        TOP_TO_BOTTOM,
    };

    /** Size of a tile in pixels. */
    itkGetMacro(TileSize, typename InputImageType::SizeType);
    itkSetMacro(TileSize, typename InputImageType::SizeType);

    /** Number of used tiles in the flat image. */
    itkGetMacro(NumberOfTiles, unsigned int);
    itkSetMacro(NumberOfTiles, unsigned int);

    /** Position of the empty tiles, defaults to BEGIN. */
    itkGetEnumMacro(EmptyTiles, EmptyTilesPosition);
    itkSetEnumMacro(EmptyTiles, EmptyTilesPosition);
    
    /** Order of the assembly, defaults to BOTTOM_TO_TOP. */
    itkGetEnumMacro(AssemblyOrder, AssemblyOrder);
    itkSetEnumMacro(AssemblyOrder, AssemblyOrder);

    /** Spacing on the new dimension, defaults to 1. */
    itkGetMacro(Spacing, double);
    itkSetMacro(Spacing, double);

    /** Origin on the new dimension, defaults to 0. */
    itkGetMacro(Origin, double);
    itkSetMacro(Origin, double);

protected :
    AssembleTilesImageFilter();
    ~AssembleTilesImageFilter();
    void PrintSelf(std::ostream& os, Indent indent) const;
    virtual void GenerateOutputInformation();
    virtual void GenerateInputRequestedRegion();
    void GenerateData();

private :

    typedef typename InputImageType::RegionType InputImageRegionType;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename InputImageType::SizeType InputImageSizeType;

    typedef typename OutputImageType::RegionType OutputImageRegionType;
    typedef typename OutputImageType::IndexType OutputImageIndexType;
    typedef typename OutputImageType::SizeType OutputImageSizeType;

    InputImageSizeType m_TileSize;
    unsigned int m_NumberOfTiles;
    EmptyTilesPosition m_EmptyTiles;
    AssemblyOrder m_AssemblyOrder;
    double m_Spacing;
    double m_Origin;

    AssembleTilesImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAssembleTilesImageFilter.txx"
#endif // ITK_MANUAL_INSTANTIATION

#endif // itkAssembleTilesImageFilter_h
