/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef medipy_vtk_addons_vtkenhancedlookuptable_h
#define medipy_vtk_addons_vtkenhancedlookuptable_h

#include <vtkLookupTable.h>

/*
 * The vtkLookupTable is not precise enough when we wish to cut (map to a
 * specific color) the values outside of table range instead of clamping them
 * to the first (or last) color in the table.
 */
class vtkEnhancedLookupTable : public vtkLookupTable
{
public :
    static vtkEnhancedLookupTable * New();
    vtkTypeRevisionMacro(vtkEnhancedLookupTable, vtkLookupTable);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description : set all values below the table range to a specific value
    vtkSetMacro(CutLow, unsigned char);
    vtkGetMacro(CutLow, unsigned char);
    vtkBooleanMacro(CutLow, unsigned char);

    // Description : set all values above the table range to a specific value
    vtkSetMacro(CutHigh, unsigned char);
    vtkGetMacro(CutHigh, unsigned char);
    vtkBooleanMacro(CutHigh, unsigned char);

    // Description : set the openness of the table range. If the table range is
    // open, the values equal to the lower value of the range will not be
    // mapped to LowValue
    vtkSetMacro(OpenLow, unsigned char);
    vtkGetMacro(OpenLow, unsigned char);
    vtkBooleanMacro(OpenLow, unsigned char);

    // Description : set the openness of the table range. If the table range is
    // open, the values equal to the higher value of the range will not be
    // mapped to HighValue
    vtkSetMacro(OpenHigh, unsigned char);
    vtkGetMacro(OpenHigh, unsigned char);
    vtkBooleanMacro(OpenHigh, unsigned char);

    // Description : value to which scalar below the the lwoer value of the
    // range will be mapped to if CutLow
    virtual void SetLowValue(double r, double g, double b, double a);
    virtual void SetLowValue(double rgba[4]);
    vtkGetVector4Macro(LowValue, double);
    vtkGetVector4Macro(LowValueChar, unsigned char);

    // Description : value to which scalar below the the lower value of the
    // range will be mapped to if CutHigh
    virtual void SetHighValue(double r, double g, double b, double a);
    virtual void SetHighValue(double rgba[4]);
    vtkGetVector4Macro(HighValue, double);
    vtkGetVector4Macro(HighValueChar, unsigned char);

    // Description : if true, 0 will always be mapped to a transparent color,
    // regardless of other parameters.
    vtkGetMacro(ZeroTransparency, unsigned char);
    vtkSetMacro(ZeroTransparency, unsigned char);
    vtkBooleanMacro(ZeroTransparency, unsigned char);

    // Description : value to which 0 will be mapped to if ZeroTransparency
    virtual void SetTransparentValue(double r, double g, double b, double a);
    virtual void SetTransparentValue(double rgba[4]);
    vtkGetVector4Macro(TransparentValue, double);
    vtkGetVector4Macro(TransparentValueChar, unsigned char);

    // Description:
    // Return the table index associated with a particular value. Values mapped
    // to LowValue (resp. HighValue) will return -2 (resp. -1).
    virtual vtkIdType GetIndex(double v);

    // Description:
    // Return a rgba color value for the given index into the lookup table. Color
    // components are expressed as [0,1] double values.
    double *GetTableValue(vtkIdType indx);

    // Description:
    // Return a rgba color value for the given index into the lookup table. Color
    // components are expressed as [0,1] double values.
    void GetTableValue(vtkIdType id, double rgba[4]);

    // Description:
    // Map one value through the lookup table.
    unsigned char * MapValue(double v);

    // Description:
    // map a set of scalars through the lookup table
    void MapScalarsThroughTable2(void *input, unsigned char *output,
                                 int inputDataType, int numberOfValues,
                                 int inputIncrement, int outputIncrement);

    // Description:
    // Copy the contents from another LookupTable
    void DeepCopy(vtkEnhancedLookupTable *lut);

    //BTX
    enum SpecialIndices
    {
        HIGH_VALUE_INDEX=-1,
        LOW_VALUE_INDEX=-2,
        TRANSPARENT_INDEX=-3
    };
    //ETX

protected :

    // Cut{Low,High}, Open{Low,High} and ZeroTransparency are set to 
    // unsigned char, as the Python wrappers for VTK 5.2 will not generate
    // the SetXXX functions for bool members
    unsigned char CutLow;
    unsigned char CutHigh;
    unsigned char OpenLow;
    unsigned char OpenHigh;
    double LowValue[4];
    double HighValue[4];
    unsigned char LowValueChar[4];
    unsigned char HighValueChar[4];
    unsigned char ZeroTransparency;
    double TransparentValue[4];
    unsigned char TransparentValueChar[4];

    vtkEnhancedLookupTable(int sze=256, int ext=256);
    virtual ~vtkEnhancedLookupTable();

private:
    vtkEnhancedLookupTable(const vtkEnhancedLookupTable&); // Not implemented.
    void operator=(const vtkEnhancedLookupTable&); // Not implemented.
};

#endif // medipy_vtk_addons_vtkenhancedlookuptable_h
