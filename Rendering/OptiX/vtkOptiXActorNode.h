/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXActorNode.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOptiXActorNode - links vtkActor and vtkMapper to OptiX
// .SECTION Description
// Translates vtkActor/Mapper state into OptiX rendering calls

#ifndef vtkOptiXActorNode_h
#define vtkOptiXActorNode_h

#include "vtkRenderingOptiXModule.h" // For export macro
#include "vtkActorNode.h"

class vtkActor;
class vtkCompositeDataDisplayAttributes;
class vtkDataArray;
class vtkInformationIntegerKey;
class vtkInformationObjectBaseKey;
class vtkInformationStringKey;
class vtkPiecewiseFunction;
class vtkPolyData;

class VTKRENDERINGOPTIX_EXPORT vtkOptiXActorNode :
  public vtkActorNode
{
public:
  static vtkOptiXActorNode* New();
  vtkTypeMacro(vtkOptiXActorNode, vtkActorNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Overridden to take into account my renderables time, including
  // mapper and data into mapper inclusive of composite input
  virtual vtkMTimeType GetMTime();

  // Description:
  // When added to the mapper, enables scale array and scale function.
  static vtkInformationIntegerKey* ENABLE_SCALING();

  // Description:
  // Convenience method to set enabled scaling on my renderable.
  static void SetEnableScaling(int value, vtkActor *);
  static int GetEnableScaling(vtkActor *);

  // Description:
  // Name of a point aligned, single component wide, double valued array that,
  // when added to the mapper, will be used to scale each element in the
  // sphere and cylinder representations individually.
  // When not supplied the radius is constant across all elements and
  // is a function of the Mapper's PointSize and LineWidth.
  static vtkInformationStringKey* SCALE_ARRAY_NAME();

  // Description:
  // Convenience method to set a scale array on my renderable.
  static void SetScaleArrayName(const char *scaleArrayName, vtkActor *);

  // Description:
  // A piecewise function for values from the scale array that alters the resulting
  // radii arbitrarily
  static vtkInformationObjectBaseKey* SCALE_FUNCTION();

  // Description:
  // Convenience method to set a scale function on my renderable.
  static void SetScaleFunction(vtkPiecewiseFunction *scaleFunction, vtkActor *);

protected:
  vtkOptiXActorNode();
  ~vtkOptiXActorNode();

private:
  vtkOptiXActorNode(const vtkOptiXActorNode&) VTK_DELETE_FUNCTION;
  void operator=(const vtkOptiXActorNode&) VTK_DELETE_FUNCTION;
};
#endif
