/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#include "vtkEnhancedLookupTable.h"

#include <vtkBitArray.h>
#include <vtkLookupTable.h>
#include <vtkObjectFactory.h>

vtkCxxRevisionMacro(vtkEnhancedLookupTable, "$Revision: 1 $");
vtkStandardNewMacro(vtkEnhancedLookupTable);

vtkEnhancedLookupTable
::vtkEnhancedLookupTable(int sze, int ext)
: vtkLookupTable(sze, ext), CutLow(false), CutHigh(false),
  OpenLow(false), OpenHigh(false), ZeroTransparency(false)
{
    this->SetLowValue(0, 0, 0, 0);
    this->SetHighValue(0, 0, 0, 0);
    this->SetTransparentValue(0, 0, 0, 0);
}


vtkEnhancedLookupTable
::~vtkEnhancedLookupTable()
{
    // Nothing to do
}


void
vtkEnhancedLookupTable
::SetLowValue(double r, double g, double b, double a)
{
    vtkDebugMacro(<< this->GetClassName() << " (" << this
                  << "): setting LowValue to ("
                  << r << "," << r << "," << b << "," << a << ")");
    if((this->LowValue[0] != r) ||
       (this->LowValue[1] != g) ||
       (this->LowValue[2] != b) ||
       (this->LowValue[3] != a))
    {
        this->LowValue[0] = r;
        this->LowValue[1] = g;
        this->LowValue[2] = b;
        this->LowValue[3] = a;

        this->LowValueChar[0] = (unsigned char)(255.*r);
        this->LowValueChar[1] = (unsigned char)(255.*g);
        this->LowValueChar[2] = (unsigned char)(255.*b);
        this->LowValueChar[3] = (unsigned char)(255.*a);

        this->Modified();
    }
}


void
vtkEnhancedLookupTable
::SetLowValue(double rgba[4])
{
  this->SetLowValue(rgba[0], rgba[1], rgba[2], rgba[3]);
}


void
vtkEnhancedLookupTable
::SetHighValue(double r, double g, double b, double a)
{
    vtkDebugMacro(<< this->GetClassName() << " (" << this
                  << "): setting HighValue to ("
                  << r << "," << r << "," << b << "," << a << ")");
    if((this->HighValue[0] != r) ||
       (this->HighValue[1] != g) ||
       (this->HighValue[2] != b) ||
       (this->HighValue[3] != a))
    {
        this->HighValue[0] = r;
        this->HighValue[1] = g;
        this->HighValue[2] = b;
        this->HighValue[3] = a;

        this->HighValueChar[0] = (unsigned char)(255.*r);
        this->HighValueChar[1] = (unsigned char)(255.*g);
        this->HighValueChar[2] = (unsigned char)(255.*b);
        this->HighValueChar[3] = (unsigned char)(255.*a);

        this->Modified();
    }
}


void
vtkEnhancedLookupTable
::SetHighValue(double rgba[4])
{
  this->SetHighValue(rgba[0], rgba[1], rgba[2], rgba[3]);
}


void
vtkEnhancedLookupTable
::SetTransparentValue(double r, double g, double b, double a)
{
    vtkDebugMacro(<< this->GetClassName() << " (" << this
                  << "): setting TransparentValue to ("
                  << r << "," << r << "," << b << "," << a << ")");
    if((this->TransparentValue[0] != r) ||
       (this->TransparentValue[1] != g) ||
       (this->TransparentValue[2] != b) ||
       (this->TransparentValue[3] != a))
    {
        this->TransparentValue[0] = r;
        this->TransparentValue[1] = g;
        this->TransparentValue[2] = b;
        this->TransparentValue[3] = a;

        this->TransparentValueChar[0] = (unsigned char)(255.*r);
        this->TransparentValueChar[1] = (unsigned char)(255.*g);
        this->TransparentValueChar[2] = (unsigned char)(255.*b);
        this->TransparentValueChar[3] = (unsigned char)(255.*a);

        this->Modified();
    }
}


void
vtkEnhancedLookupTable
::SetTransparentValue(double rgba[4])
{
  this->SetTransparentValue(rgba[0], rgba[1], rgba[2], rgba[3]);
}


vtkIdType
vtkEnhancedLookupTable
::GetIndex(double v)
{
    if(v==0 && this->ZeroTransparency)
    {
        return TRANSPARENT_INDEX;
    }
    if(this->CutLow &&
       (v<this->TableRange[0] || (this->OpenLow && v==this->TableRange[0])))
    {
        return LOW_VALUE_INDEX;
    }
    if(this->CutHigh &&
       (v>this->TableRange[1] || (this->OpenHigh && v==this->TableRange[1])))
    {
        return HIGH_VALUE_INDEX;
    }
    else
    {
        return this->Superclass::GetIndex(v);
    }
}


double *
vtkEnhancedLookupTable
::GetTableValue(vtkIdType indx)
{
    this->GetTableValue(indx, this->RGBA);
    return this->RGBA;
}


void
vtkEnhancedLookupTable
::GetTableValue(vtkIdType id, double rgba[4])
{
    if(id==HIGH_VALUE_INDEX)
    {
        rgba[0] = this->HighValue[0];
        rgba[1] = this->HighValue[1];
        rgba[2] = this->HighValue[2];
        rgba[3] = this->HighValue[3];
    }
    else if(id==LOW_VALUE_INDEX)
    {
        rgba[0] = this->LowValue[0];
        rgba[1] = this->LowValue[1];
        rgba[2] = this->LowValue[2];
        rgba[3] = this->LowValue[3];
    }
    else if(id==TRANSPARENT_INDEX)
    {
        rgba[0] = this->TransparentValue[0];
        rgba[1] = this->TransparentValue[1];
        rgba[2] = this->TransparentValue[2];
        rgba[3] = this->TransparentValue[3];
    }
    else
    {
        this->Superclass::GetTableValue(id, rgba);
    }
}


unsigned char *
vtkEnhancedLookupTable
::MapValue(double v)
{
    if(v==0 && this->ZeroTransparency)
    {
        return this->TransparentValueChar;
    }
    if(this->CutHigh &&
       (v>this->TableRange[1] || (this->OpenHigh && v==this->TableRange[1])))
    {
        return this->HighValueChar;
    }
    if(this->CutLow &&
       (v<this->TableRange[0] || (this->OpenLow && v==this->TableRange[0])))
    {
        return this->LowValueChar;
    }
    else
    {
        return this->Superclass::MapValue(v);
    }
}

void
vtkEnhancedLookupTable
::DeepCopy(vtkEnhancedLookupTable *lut)
{
    if (!lut)
    {
        return;
    }

    this->Superclass::DeepCopy(lut);
    this->CutLow = lut->CutLow;
    this->CutHigh = lut->CutHigh;
    this->OpenLow = lut->OpenLow;
    this->OpenHigh = lut->OpenHigh;
    for(unsigned int i=0; i<4; ++i)
    {
        this->LowValue[i] = lut->LowValue[i];
        this->HighValue[i] = lut->HighValue[i];
        this->LowValueChar[i] = lut->LowValueChar[i];
        this->HighValueChar[i] = lut->HighValueChar[i];
        this->TransparentValue[i] = lut->TransparentValue[i];
        this->TransparentValueChar[i] = lut->TransparentValueChar[i];
    }
    this->ZeroTransparency = lut->ZeroTransparency;
}


void
vtkEnhancedLookupTable
::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
    os << indent << "Cut (Low, High): (" << bool(this->CutLow) << ", "
       << bool(this->CutHigh) << ")\n";
    os << indent << "Open (Low, High): (" << bool(this->OpenLow) << ", "
       << bool(this->OpenHigh) << ")\n";
    os << indent << "(Low, High) Value: ("
       << this->LowValue[0] << ", " << this->LowValue[1] << ", "
       << this->LowValue[2] << ", " << this->LowValue[3] << "), ("
       << this->HighValue[0] << ", " << this->HighValue[1] << ", "
       << this->HighValue[2] << ", " << this->HighValue[3] << ")\n";
    os << indent << "(Low, High) Value (char): ("
           << (int) this->LowValueChar[0] << ", " << (int) this->LowValueChar[1] << ", "
           << (int) this->LowValueChar[2] << ", " << (int) this->LowValueChar[3] << "), ("
           << (int) this->HighValueChar[0] << ", " << (int) this->HighValueChar[1] << ", "
           << (int) this->HighValueChar[2] << ", " << (int) this->HighValueChar[3] << ")\n";
    os << indent << "ZeroTransparency: " << bool(this->ZeroTransparency) << "\n";
    os << indent << "Transparent Value: ("
       << this->TransparentValue[0] << ", " << this->TransparentValue[1] << ", "
       << this->TransparentValue[2] << ", " << this->TransparentValue[3] << ")\n";
    os << indent << "Transparent Value (char): ("
           << this->TransparentValueChar[0] << ", " << this->TransparentValueChar[1] << ", "
           << this->TransparentValueChar[2] << ", " << this->TransparentValueChar[3] << ")\n";
}


//----------------------------------------------------------------------------
// There is a little more to this than simply taking the log10 of the
// two range values: we do conversion of negative ranges to positive
// ranges, and conversion of zero to a 'very small number'
void vtkLookupTableLogRange(double range[2], double logRange[2])
{
    double rmin = range[0];
    double rmax = range[1];

    if (rmin == 0)
    {
        rmin = 1.0e-6 * (rmax - rmin);
        if (rmax < 0)
        {
            rmin = -rmin;
        }
    }
    if (rmax == 0)
    {
        rmax = 1.0e-6 * (rmin - rmax);
        if (rmin < 0)
        {
            rmax = -rmax;
        }
    }
    if (rmin < 0 && rmax < 0)
    {
        logRange[0] = log10(-(double) rmin);
        logRange[1] = log10(-(double) rmax);
    }
    else if (rmin > 0 && rmax > 0)
    {
        logRange[0] = log10((double) rmin);
        logRange[1] = log10((double) rmax);
    }
}


//----------------------------------------------------------------------------
// Apply log to value, with appropriate constraints.
inline double vtkApplyLogScale(double v, double range[2], double logRange[2])
{
    // is the range set for negative numbers?
    if (range[0] < 0)
    {
        if (v < 0)
        {
            v = log10(-static_cast<double> (v));
        }
        else if (range[0] > range[1])
        {
            v = logRange[0];
        }
        else
        {
            v = logRange[1];
        }
    }
    else
    {
        if (v > 0)
        {
            v = log10(static_cast<double> (v));
        }
        else if (range[0] < range[1])
        {
            v = logRange[0];
        }
        else
        {
            v = logRange[1];
        }
    }
    return v;
}


//----------------------------------------------------------------------------
// Apply shift/scale to the scalar value v and do table lookup.
inline unsigned char *vtkLinearLookup(double v, unsigned char *table,
                                      double maxIndex, double shift,
                                      double scale, double range[2],
                                      vtkEnhancedLookupTable * lut)
{
    double findx;
    if(v==0 && lut->GetZeroTransparency())
    {
        findx=vtkEnhancedLookupTable::TRANSPARENT_INDEX;
    }
    else if(lut->GetCutLow() && (v<range[0] || (v==range[0] && lut->GetOpenLow())))
    {
        findx=vtkEnhancedLookupTable::LOW_VALUE_INDEX;
    }
    else if(lut->GetCutHigh() && (v>range[1] || (v==range[1] && lut->GetOpenHigh())))
    {
        findx=vtkEnhancedLookupTable::HIGH_VALUE_INDEX;
    }
    else
    {
        findx = (v + shift) * scale;
        if(findx < 0)
        {
            findx = 0;
        }

        if(findx > maxIndex)
        {
            findx = maxIndex;
        }
    }

    if(findx == vtkEnhancedLookupTable::TRANSPARENT_INDEX)
    {
        return lut->GetTransparentValueChar();
    }
    else if(findx == vtkEnhancedLookupTable::LOW_VALUE_INDEX)
    {
        return lut->GetLowValueChar();
    }
    else if(findx == vtkEnhancedLookupTable::HIGH_VALUE_INDEX)
    {
        return lut->GetHighValueChar();
    }
    else
    {
        return &table[4*static_cast<int>(findx)];
    }
}


//----------------------------------------------------------------------------
// accelerate the mapping by copying the data in 32-bit chunks instead
// of 8-bit chunks
template<class T>
void vtkEnhancedLookupTableMapData(vtkEnhancedLookupTable *self, T *input,
                                   unsigned char *output, int length, int inIncr,
                                   int outFormat)
{
    int i = length;
    double *range = self->GetTableRange();
    double maxIndex = self->GetNumberOfColors() - 1;
    double shift, scale;
    unsigned char *table = self->GetPointer(0);
    unsigned char *cptr;
    double alpha;

    if ((alpha = self->GetAlpha()) >= 1.0) //no blending required
    {
        if (self->GetScale() == VTK_SCALE_LOG10)
        {
            double val;
            double logRange[2];
            vtkLookupTableLogRange(range, logRange);
            shift = -logRange[0];
            if (logRange[1] <= logRange[0])
            {
                scale = VTK_DOUBLE_MAX;
            }
            else
            {
                /* while this looks like the wrong scale, it is the correct scale
                 * taking into account the truncation to int that happens below. */
                scale = (maxIndex + 1) / (logRange[1] - logRange[0]);
            }
            if (outFormat == VTK_RGBA)
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_RGB)
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_LUMINANCE_ALPHA)
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    *output++ = cptr[3];
                    input += inIncr;
                }
            }
            else // outFormat == VTK_LUMINANCE
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    input += inIncr;
                }
            }
        }//if log scale

        else //not log scale
        {
            shift = -range[0];
            if (range[1] <= range[0])
            {
                scale = VTK_DOUBLE_MAX;
            }
            else
            {
                /* while this looks like the wrong scale, it is the correct scale
                 * taking into account the truncation to int that happens below. */
                scale = (maxIndex + 1) / (range[1] - range[0]);
            }

            if (outFormat == VTK_RGBA)
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_RGB)
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_LUMINANCE_ALPHA)
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    *output++ = cptr[3];
                    input += inIncr;
                }
            }
            else // outFormat == VTK_LUMINANCE
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    input += inIncr;
                }
            }
        }//if not log lookup
    }//if blending not needed

    else //blend with the specified alpha
    {
        if (self->GetScale() == VTK_SCALE_LOG10)
        {
            double val;
            double logRange[2];
            vtkLookupTableLogRange(range, logRange);
            shift = -logRange[0];
            if (logRange[1] <= logRange[0])
            {
                scale = VTK_DOUBLE_MAX;
            }
            else
            {
                /* while this looks like the wrong scale, it is the correct scale
                 * taking into account the truncation to int that happens below. */
                scale = (maxIndex + 1) / (logRange[1] - logRange[0]);
            }
            if (outFormat == VTK_RGBA)
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = static_cast<unsigned char> ((*cptr) * alpha);
                    cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_RGB)
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_LUMINANCE_ALPHA)
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    *output++ = static_cast<unsigned char> (alpha * cptr[3]);
                    input += inIncr;
                }
            }
            else // outFormat == VTK_LUMINANCE
            {
                while (--i >= 0)
                {
                    val = vtkApplyLogScale(*input, range, logRange);
                    cptr = vtkLinearLookup(val, table, maxIndex, shift, scale,
                                           logRange, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    input += inIncr;
                }
            }
        }//log scale with blending

        else //no log scale with blending
        {
            shift = -range[0];
            if (range[1] <= range[0])
            {
                scale = VTK_DOUBLE_MAX;
            }
            else
            {
                /* while this looks like the wrong scale, it is the correct scale
                 * taking into account the truncation to int that happens below. */
                scale = (maxIndex + 1) / (range[1] - range[0]);
            }

            if (outFormat == VTK_RGBA)
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = static_cast<unsigned char> ((*cptr) * alpha);
                    cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_RGB)
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    *output++ = *cptr++;
                    input += inIncr;
                }
            }
            else if (outFormat == VTK_LUMINANCE_ALPHA)
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    *output++ = static_cast<unsigned char> (cptr[3] * alpha);
                    input += inIncr;
                }
            }
            else // outFormat == VTK_LUMINANCE
            {
                while (--i >= 0)
                {
                    cptr = vtkLinearLookup(*input, table, maxIndex, shift,
                                           scale, range, self);
                    *output++ = static_cast<unsigned char> (cptr[0] * 0.30
                        + cptr[1] * 0.59 + cptr[2] * 0.11 + 0.5);
                    input += inIncr;
                }
            }
        }//no log scale
    }//alpha blending
}


//----------------------------------------------------------------------------
// Although this is a relatively expensive calculation,
// it is only done on the first render. Colors are cached
// for subsequent renders.
template<class T>
void vtkEnhancedLookupTableMapMag(vtkEnhancedLookupTable *self, T *input,
                                  unsigned char *output, int length, int inIncr,
                                  int outFormat)
{
    double tmp, sum;
    double *mag;
    int i, j;

    mag = new double[length];
    for (i = 0; i < length; ++i)
    {
        sum = 0;
        for (j = 0; j < inIncr; ++j)
        {
            tmp = (double) (*input);
            sum += (tmp * tmp);
            ++input;
        }
        mag[i] = sqrt(sum);
    }

    vtkEnhancedLookupTableMapData(self, mag, output, length, 1, outFormat);

    delete[] mag;
}


void
vtkEnhancedLookupTable
::MapScalarsThroughTable2(void *input, unsigned char *output,
                          int inputDataType, int numberOfValues,
                          int inputIncrement, int outputFormat)
{
    if (this->UseMagnitude && inputIncrement > 1)
    {
        switch (inputDataType)
        {
        vtkTemplateMacro(
            vtkEnhancedLookupTableMapMag(this,static_cast<VTK_TT*>(input),output,
                numberOfValues,inputIncrement,outputFormat);
            return);
        case VTK_BIT:
            vtkErrorMacro("Cannot compute magnitude of bit array.");
            break;
        default:
            vtkErrorMacro(<< "MapImageThroughTable: Unknown input ScalarType");
        }
    }

    switch (inputDataType)
    {
    case VTK_BIT:
        {
            vtkIdType i, id;
            vtkBitArray *bitArray = vtkBitArray::New();
            bitArray->SetVoidArray(input, numberOfValues, 1);
            vtkUnsignedCharArray *newInput = vtkUnsignedCharArray::New();
            newInput->SetNumberOfValues(numberOfValues);
            for (id = i = 0; i < numberOfValues; i++, id += inputIncrement)
            {
                newInput->SetValue(i, bitArray->GetValue(id));
            }
            vtkEnhancedLookupTableMapData(this,
                static_cast<unsigned char*>(newInput->GetPointer(0)),
                output, numberOfValues, inputIncrement, outputFormat);
            newInput->Delete();
            bitArray->Delete();
        }
        break;

    vtkTemplateMacro(
        vtkEnhancedLookupTableMapData(this,static_cast<VTK_TT*>(input),output,
            numberOfValues,inputIncrement,outputFormat));
    default:
        vtkErrorMacro(<< "MapImageThroughTable: Unknown input ScalarType");
        return;
    }
}
