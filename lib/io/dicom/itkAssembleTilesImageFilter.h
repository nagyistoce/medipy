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

    /** Size of a tile in pixels. */
    itkGetMacro(TileSize, typename InputImageType::SizeType);
    itkSetMacro(TileSize, typename InputImageType::SizeType);

    /** Number of used tiles in the flat image. */
    itkGetMacro(NumberOfTiles, unsigned int);
    itkSetMacro(NumberOfTiles, unsigned int);

    /** Position of the empty tiles. */
    itkGetEnumMacro(EmptyTiles, EmptyTilesPosition);
    itkSetEnumMacro(EmptyTiles, EmptyTilesPosition);
    
protected :
    AssembleTilesImageFilter();
    ~AssembleTilesImageFilter() {}
    void PrintSelf(std::ostream& os, Indent indent) const;
    void GenerateData();

private :
    typename InputImageType::SizeType m_TileSize;
    unsigned int m_NumberOfTiles;
    EmptyTilesPosition m_EmptyTiles;

    AssembleTilesImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAssembleTilesImageFilter.txx"
#endif // ITK_MANUAL_INSTANTIATION

#endif // itkAssembleTilesImageFilter_h
