/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSPHCubicKernel.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSPHCubicKernel - a quintic SPH interpolation kernel

// .SECTION Description
// vtkSPHCubicKernel is an smooth particle hydrodynamics interpolation kernel as
// described by D.J. Price. This is a quintic formulation.
//
// .SECTION Caveats
// See D.J. Price, Smoothed particle hydrodynamics and magnetohydrodynamics,
// J. Comput. Phys. 231:759-794, 2012. Especially equation 49.

// .SECTION See Also
// vtkSPHKernel vtkSPHInterpolator


#ifndef vtkSPHCubicKernel_h
#define vtkSPHCubicKernel_h

#include "vtkFiltersPointsModule.h" // For export macro
#include "vtkSPHKernel.h"
#include <algorithm> // For std::min()

class vtkIdList;
class vtkDoubleArray;


class VTKFILTERSPOINTS_EXPORT vtkSPHCubicKernel : public vtkSPHKernel
{
public:
  // Description:
  // Standard methods for instantiation, obtaining type information, and printing.
  static vtkSPHCubicKernel *New();
  vtkTypeMacro(vtkSPHCubicKernel,vtkSPHKernel);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Compute weighting factor given a normalized distance from a sample point.
  virtual double ComputeFunctionWeight(const double d)
  {
    double tmp1 = 2.0 - std::min(d,2.0);
    double tmp2 = 1.0 - std::min(d,1.0);
    return (0.25*tmp1*tmp1*tmp1 - tmp2*tmp2*tmp2);
  }

  // Description:
  // Compute weighting factor for derivative quantities given a normalized
  // distance from a sample point.
  virtual double ComputeGradientWeight(const double d)
  {
    double tmp1 = 2.0 - std::min(d,2.0);
    double tmp2 = 1.0 - std::min(d,1.0);
    return (0.25*tmp1*tmp1 - tmp2*tmp2);
  }

protected:
  vtkSPHCubicKernel();
  ~vtkSPHCubicKernel();

private:
  vtkSPHCubicKernel(const vtkSPHCubicKernel&);  // Not implemented.
  void operator=(const vtkSPHCubicKernel&);  // Not implemented.
};

#endif
