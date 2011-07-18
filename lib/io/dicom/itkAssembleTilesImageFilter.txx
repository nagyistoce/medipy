#ifndef itkAssembleTilesImageFilter_txx
#define itkAssembleTilesImageFilter_txx

#include "itkAssembleTilesImageFilter.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include <itkFixedArray.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace itk
{

template<typename TInputImage, typename TOutputImage>
AssembleTilesImageFilter<TInputImage, TOutputImage>::AssembleTilesImageFilter()
{
    // TODO
}

template<typename TInputImage, typename TOutputImage>
void AssembleTilesImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Tile Size : " << this->m_TileSize << "\n";
    os << indent << "Number of Tiles : " << this->m_NumberOfTiles << "\n";
    os << indent << "Empty Tiles : " << this->m_EmptyTiles << "\n";
}


template<typename TInputImage, typename TOutputImage>
void AssembleTilesImageFilter<TInputImage, TOutputImage>::GenerateData()
{
    // Check that TileSize and NumberOfTiles is compatible with InputImage

    // Compute the size of the input image in tiles
    typename InputImageType::SizeType const size_in_pixels =
        this->GetInput()->GetRequestedRegion().GetSize();
    typename InputImageType::SizeType size_in_tiles;
    std::transform(size_in_pixels.m_Size, size_in_pixels.m_Size+size_in_pixels.GetSizeDimension(),
                   this->m_TileSize.m_Size, size_in_tiles.m_Size,
                   std::divides<typename InputImageType::SizeType::SizeValueType>());

    itkDebugMacro(<< "Size in tiles : " << size_in_tiles);

    // Compute the total number of tiles
    typename InputImageType::SizeType::SizeValueType total_number_of_tiles;
    total_number_of_tiles = std::accumulate(
        size_in_tiles.m_Size, size_in_tiles.m_Size+size_in_tiles.GetSizeDimension(),
        1, std::multiplies<typename InputImageType::SizeType::SizeValueType>());
    itkDebugMacro(<< "Total number of tiles : " << total_number_of_tiles);

    typedef itk::ImageBase<TInputImage::ImageDimension> ImageBaseType;
    typename ImageBaseType::Pointer image_base = ImageBaseType::New();
    typename ImageBaseType::RegionType region;
    region.SetSize(size_in_tiles);
    // SetRequestedRegion does not call UpdateOffsetTable
    image_base->SetBufferedRegion(region);

    // First and last tile to build (last is included)
    unsigned int first_tile;
    unsigned int last_tile;
    if(this->m_EmptyTiles == END)
    {
        // The first tile is at the beginning
        first_tile = 0;
        last_tile = this->m_NumberOfTiles;
    }
    else
    {
        first_tile = total_number_of_tiles-this->m_NumberOfTiles;
        last_tile = total_number_of_tiles-1;
    }

    for(unsigned int i=first_tile; i<=last_tile; ++i)
    {
        typename InputImageType::IndexType tile_begin;
        typename InputImageType::IndexType tile_index = image_base->ComputeIndex(i);
        std::transform(tile_index.GetIndex(), tile_index.GetIndex()+tile_index.GetIndexDimension(),
                       this->m_TileSize.m_Size, tile_begin.m_Index,
                       std::multiplies<int>());

        typename InputImageType::IndexType tile_end;
        std::transform(tile_begin.GetIndex(), tile_begin.GetIndex()+tile_begin.GetIndexDimension(),
                       this->m_TileSize.GetSize(), tile_end.m_Index,
                       std::plus<int>());

        typename InputImageType::RegionType input_region;
        input_region.SetIndex(tile_begin);
        typename InputImageType::RegionType::SizeType input_size;
        std::transform(tile_end.GetIndex(), tile_end.GetIndex()+tile_end.GetIndexDimension(),
            tile_begin.GetIndex(), input_size.m_Size,
            std::minus<typename InputImageType::IndexType::IndexValueType>());
        input_region.SetSize(input_size);

        // TODO : slice direction
        typename OutputImageType::RegionType output_region;

        typename OutputImageType::RegionType::IndexType output_index;
        output_index.Fill(0);
        output_index[InputImageType::ImageDimension] = i-first_tile;
        output_region.SetIndex(output_index);

        typename OutputImageType::RegionType::SizeType output_size;
        output_size.Fill(1);
        std::copy(this->m_TileSize.GetSize(), this->m_TileSize.GetSize()+this->m_TileSize.GetSizeDimension(),
            output_size.m_Size);
        output_region.SetSize(output_size);

        itkDebugMacro(<< "InputRegion " << input_region)
        itkDebugMacro(<< "Output region " << output_region);

        // Copy region between tile_begin and tile_end to slice i of output
        itk::ImageRegionConstIterator<InputImageType> source_iterator(
            this->GetInput(), input_region);
        itk::ImageRegionIterator<OutputImageType> destination_iterator(
            this->GetOutput(), output_region);
        while(!source_iterator.IsAtEnd())
        {
            destination_iterator.Set(source_iterator.Get());
            ++source_iterator;
            ++destination_iterator;
        }

    }
}

}

#endif // itkAssembleTilesImageFilter_txx
