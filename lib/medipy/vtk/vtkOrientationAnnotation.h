/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef __vtkOrientationAnnotation_h
#define __vtkOrientationAnnotation_h

#include <vtkCornerAnnotation.h>

class vtkOrientationAnnotation : public vtkCornerAnnotation
{
public:
  vtkTypeRevisionMacro(vtkOrientationAnnotation,vtkCornerAnnotation);
  
  // Description:
  // Instantiate object with a rectangle in normaled view coordinates
  // of (0.2,0.85, 0.8, 0.95).
  static vtkOrientationAnnotation *New();

protected:
  vtkOrientationAnnotation();
  ~vtkOrientationAnnotation();

  // Description:
  // Set text actor positions given a viewport size and justification
  virtual void SetTextActorsPosition(int vsize[2]);
  virtual void SetTextActorsJustification();

private:
  vtkOrientationAnnotation(const vtkOrientationAnnotation&);  // Not implemented.
  void operator=(const vtkOrientationAnnotation&);  // Not implemented.
};


#endif



