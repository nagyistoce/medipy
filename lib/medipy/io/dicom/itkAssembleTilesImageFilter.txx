/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

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
AssembleTilesImageFilter<TInputImage, TOutputImage>
::AssembleTilesImageFilter()
: m_NumberOfTiles(0), m_EmptyTiles(END), m_AssemblyOrder(BOTTOM_TO_TOP),
  m_Spacing(1.), m_Origin(0.)
{
    this->m_TileSize.Fill(0);

    if(InputImageType::GetImageDimension() >= OutputImageType::GetImageDimension())
    {
        itkExceptionMacro(
            << "Output image dimension ("
            << OutputImageType::ImageDimension
            << ") must be strictly greater than input image dimension ("
            << InputImageType::ImageDimension);
    }
}

template<typename TInputImage, typename TOutputImage>
AssembleTilesImageFilter<TInputImage, TOutputImage>
::~AssembleTilesImageFilter()
{
    // Nothing to do.
}

template<typename TInputImage, typename TOutputImage>
void AssembleTilesImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Tile Size : " << this->m_TileSize << "\n";
    os << indent << "Number of Tiles : " << this->m_NumberOfTiles << "\n";
    os << indent << "Empty Tiles : " << this->m_EmptyTiles << "\n";
    os << indent << "Spacing : " << this->m_Spacing << "\n";
    os << indent << "Origin : " << this->m_Origin << "\n";
}

template<typename TInputImage, typename TOutputImage>
void AssembleTilesImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
    // Do not call the superclass' implementation of this method since
    // the input and output are of different dimensions.

    unsigned int const InputImageDimension = InputImageType::GetImageDimension();
    unsigned int const OutputImageDimension = OutputImageType::GetImageDimension();

    InputImageType const * input = this->GetInput();
    if(!input)
    {
        return;
    }

    OutputImageType* output = this->GetOutput();
    if(!output)
    {
        return;
    }

    // Set the output image largest possible region.
    typename OutputImageRegionType::SizeType output_size;
    output_size.Fill(1);
    std::copy(this->m_TileSize.GetSize(), this->m_TileSize.GetSize()+InputImageDimension,
              output_size.m_Size);
    output_size[InputImageDimension] = this->m_NumberOfTiles;

    typename OutputImageRegionType::IndexType output_index;
    output_index.Fill(0);

    OutputImageRegionType output_region(output_index, output_size);
    output->SetLargestPossibleRegion(output_region);

    // Spacing : use what we know from input, use 1 for the rest
    typename InputImageType::SpacingType const & input_spacing = input->GetSpacing();
    typename OutputImageType::SpacingType output_spacing;
    output_spacing.Fill(1);
    std::copy(input_spacing.Begin(), input_spacing.End(), output_spacing.Begin());

    // Origin : use what we know from input, use 0 for the rest
    typename InputImageType::PointType const & input_origin = input->GetOrigin();
    typename OutputImageType::PointType output_origin;
    output_origin.Fill(0);
    std::copy(input_origin.Begin(), input_origin.End(), output_origin.Begin());

    // Set origin and spacing for the new dimension
    output_spacing[InputImageDimension] = this->m_Spacing;
    output_origin[InputImageDimension] = this->m_Origin;


    // Copy the direction cosines from the input to the output.
    typedef typename InputImageType::DirectionType InputDirectionType;
    typedef typename OutputImageType::DirectionType OutputDirectionType;
    InputDirectionType input_direction = input->GetDirection();
    OutputDirectionType output_direction = output->GetDirection();
    for(unsigned int i=0; i<OutputImageDimension; ++i)
    {
        for(unsigned int j=0; j<OutputImageDimension; ++j)
        {
            if(i<InputImageDimension && j<InputImageDimension)
            {
                output_direction[i][j] = input_direction[i][j];
            }
            else
            {
                if(i==j)
                {
                    output_direction[i][j] = (this->m_AssemblyOrder==BOTTOM_TO_TOP)?+1.:-1.;
                }
                else
                {
                    output_direction[i][j] = 0.;
                }
            }
        }
    }

    output->SetSpacing(output_spacing);
    output->SetOrigin(output_origin);
    output->SetDirection(output_direction);
}

template<typename TInputImage, typename TOutputImage>
void AssembleTilesImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
    this->Superclass::GenerateInputRequestedRegion();

    if(this->GetInput() == NULL)
    {
        return;
    }

    // Request the whole image
    InputImageType* input = const_cast<InputImageType*>(this->GetInput());
    input->SetRequestedRegion(input->GetLargestPossibleRegion());
}

template<typename TInputImage, typename TOutputImage>
void AssembleTilesImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
    InputImageType const * input = this->GetInput();

    // Allocate the output buffer
    OutputImageType* output = this->GetOutput();
    output->SetBufferedRegion(output->GetRequestedRegion());
    output->Allocate();

    // Compute the size of the input image in tiles
    InputImageSizeType const size_in_pixels = input->GetRequestedRegion().GetSize();
    InputImageSizeType size_in_tiles;
    std::transform(size_in_pixels.GetSize(), size_in_pixels.GetSize()+size_in_pixels.GetSizeDimension(),
                   this->m_TileSize.GetSize(), size_in_tiles.m_Size,
                   std::divides<typename InputImageSizeType::SizeValueType>());

    // Compute the total number of tiles
    typename InputImageSizeType::SizeValueType total_number_of_tiles;
    total_number_of_tiles = std::accumulate(
        size_in_tiles.GetSize(), size_in_tiles.GetSize()+size_in_tiles.GetSizeDimension(),
        1, std::multiplies<typename InputImageSizeType::SizeValueType>());

    // First and last tile to assemble (last is included)
    unsigned int first_tile;
    unsigned int last_tile;
    if(this->m_EmptyTiles == END)
    {
        // The first tile is at the beginning
        first_tile = 0;
        last_tile = this->m_NumberOfTiles-1;
    }
    else
    {
        first_tile = total_number_of_tiles-this->m_NumberOfTiles;
        last_tile = total_number_of_tiles-1;
    }

    // Use ImageBase to linear tile offset to tile index
    typedef itk::ImageBase<InputImageType::ImageDimension> ImageBaseType;
    typename ImageBaseType::Pointer image_base = ImageBaseType::New();
    typename ImageBaseType::RegionType region(size_in_tiles);
    // SetRequestedRegion does not call UpdateOffsetTable
    image_base->SetBufferedRegion(region);

    // Copy each tile to output
    for(unsigned int i=first_tile; i<=last_tile; ++i)
    {
        // Index of current tile in tile space
        InputImageIndexType const tile_index = image_base->ComputeIndex(i);

        // First index of current tile, in voxel space
        InputImageIndexType tile_begin;
        std::transform(tile_index.GetIndex(), tile_index.GetIndex()+tile_index.GetIndexDimension(),
                       this->m_TileSize.m_Size, tile_begin.m_Index,
                       std::multiplies<typename InputImageIndexType::IndexValueType>());

        InputImageRegionType const source_region(tile_begin, this->m_TileSize);

        OutputImageIndexType destination_index;
        destination_index.Fill(0);
        if(this->m_AssemblyOrder == TOP_TO_BOTTOM)
        {
            destination_index[InputImageType::ImageDimension] = last_tile-i;
        }
        else
        {
            destination_index[InputImageType::ImageDimension] = i-first_tile;
        }

        OutputImageSizeType destination_size;
        destination_size.Fill(1);
        std::copy(this->m_TileSize.GetSize(), this->m_TileSize.GetSize()+this->m_TileSize.GetSizeDimension(),
                  destination_size.m_Size);

        OutputImageRegionType const destination_region(destination_index, destination_size);

        // Copy region between tile_begin and tile_end to slice i of output
        itk::ImageRegionConstIterator<InputImageType> source(input, source_region);
        itk::ImageRegionIterator<OutputImageType> destination(output, destination_region);
        while(!source.IsAtEnd())
        {
            destination.Set(source.Get());
            ++source;
            ++destination;
        }

    }
}

}

#endif // itkAssembleTilesImageFilter_txx
