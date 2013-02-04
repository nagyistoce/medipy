/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

// Code from VTK 5.0.4 with additions for alpha transparency

#ifndef vtk_addons_vtkColorTransferFunctionWithAlpha_h
#define vtk_addons_vtkColorTransferFunctionWithAlpha_h

#include <vtkScalarsToColors.h>

#define VTK_CTF_RGB           0
#define VTK_CTF_HSV           1

class vtkColorTransferFunctionWithAlpha : public vtkScalarsToColors
{
public:
    static vtkColorTransferFunctionWithAlpha *New();
    vtkTypeRevisionMacro(vtkColorTransferFunctionWithAlpha,vtkScalarsToColors);
    void DeepCopy(vtkColorTransferFunctionWithAlpha *f );

    // Description:
    // Print method for vtkColorTransferFunction
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // How many points are there defining this function?
    int GetSize() {return this->NumberOfPoints;}

    // Description:
    // Add/Remove a point to/from the function defined in RGB or HSV
    // Return the index of the point (0 based), or -1 on error.
    int AddRGBAPoint(double x, double r, double g, double b, double a);
    int AddHSVAPoint(double x, double h, double s, double v, double a);
    int RemovePoint(double x);

    // Description:
    // Add two points to the function and remove all the points
    // between them
    void AddRGBASegment(double x1, double r1, double g1, double b1, double a1,
                        double x2, double r2, double g2, double b2, double a2);
    void AddHSVASegment(double x1, double h1, double s1, double v1, double a1,
                        double x2, double h2, double s2, double v2, double a2);

    // Description:
    // Remove all points
    void RemoveAllPoints();

    // Description:
    // Returns an RGB color for the specified scalar value
    double *GetColor(double x) { return vtkScalarsToColors::GetColor(x); }
    void GetColor(double x, double rgba[4]);

    // Description:
    // Get the color components individually.
    double GetRedValue(double x);
    double GetGreenValue(double x);
    double GetBlueValue(double x);
    double GetAlphaValue(double x);

    // Description:
    // Map one value through the lookup table.
    virtual unsigned char *MapValue(double v);

    // Description:
    // Returns min and max position of all function points.
    vtkGetVector2Macro(Range, double);

    // Description:
    // Remove all points out of the new range, and make sure there is a point
    // at each end of that range.
    // Return 1 on success, 0 otherwise.
    int AdjustRange(double range[2]);

    // Description:
    // Fills in a table of n function values between x1 and x2
    void GetTable(double x1, double x2, int n, double* table);
    void GetTable(double x1, double x2, int n, float* table);
    const unsigned char *GetTable(double x1, double x2, int n);

    // Description:
    // Construct a color transfer function from a table. Function range is
    // is set to [x1, x2], each function size is set to size, and function
    // points are regularly spaced between x1 and x2. Parameter "table" is
    // assumed to be a block of memory of size [3*size]
    void BuildFunctionFromTable(double x1, double x2, int size, double *table);

    // Description:
    // Sets and gets the clamping value for this transfer function.
    vtkSetClampMacro(Clamping, int, 0, 1);
    vtkGetMacro(Clamping, int);
    vtkBooleanMacro(Clamping, int);
  
    // Description:
    // Set/Get the color space used for interpolation: RGB, or HSV.
    // In HSV mode, if HSVWrap is on, it  will take the shortest path in Hue
    // (going back through 0 if that is the shortest way around the hue circle)
    // whereas if HSVWrap is off it will not go through 0 (in order the match
    // the current functionality of vtkLookupTable)
    vtkSetClampMacro(ColorSpace, int, VTK_CTF_RGB, VTK_CTF_HSV);
    void SetColorSpaceToRGB() { this->SetColorSpace(VTK_CTF_RGB); }
    void SetColorSpaceToHSV() { this->SetColorSpace(VTK_CTF_HSV); }
    vtkGetMacro(ColorSpace, int);

    // Description:
    // Returns a list of all nodes
    // Fills from a pointer to data stored in a similar list of nodes.
    double *GetDataPointer() { return this->Function; }
    void FillFromDataPointer(int, double*);

    // Description:
    // map a set of scalars through the lookup table
    virtual void MapScalarsThroughTable2(void *input, unsigned char *output,
                                         int inputDataType, int numberOfValues,
                                         int inputIncrement, int outputIncrement);
    
    virtual vtkIdType GetNumberOfAvailableColors();

protected:
    vtkColorTransferFunctionWithAlpha();
    ~vtkColorTransferFunctionWithAlpha();

    // Determines the function value outside of defined points
    // Zero = always return 0.0 outside of defined points
    // One  = clamp to the lowest value below defined points and
    //        highest value above defined points
    int Clamping;

    // The color space in which interpolation is performed
    int ColorSpace;

    // The color function, stores values as a sequence of x,r,g,b,a
    double     *Function;
    int         FunctionSize;
    int         NumberOfPoints;

    // An evaluated color (0 to 255 RGBA A=255)
    unsigned char UnsignedCharRGBAValue[4];

    // The min and max point locations for all three transfer functions
    double Range[2];

    // Transfer functions for each color component
    // Remove after corresponding depricated methods are removed
    vtkTimeStamp BuildTime;
    unsigned char *Table;
    int TableSize;

    // Description:
    // Set the range of scalars being mapped. The set has no functionality
    // in this subclass of vtkScalarsToColors.
    virtual void SetRange(double, double) {}
    void SetRange(double rng[2]) { this->SetRange(rng[0],rng[1]); }


private:
    vtkColorTransferFunctionWithAlpha(const vtkColorTransferFunctionWithAlpha&);  // Not implemented.
    void operator=(const vtkColorTransferFunctionWithAlpha&);  // Not implemented.
};

#endif // vtk_addons_vtkColorTransferFunctionWithAlpha_h
