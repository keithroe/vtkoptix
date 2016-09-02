/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXCameraNode.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOptiXCameraNode - links vtkCamera to OptiX
// .SECTION Description
// Translates vtkCamera state into OptiX rendering calls

#ifndef vtkOptiXCameraNode_h
#define vtkOptiXCameraNode_h

#include "vtkRenderingOptiXModule.h" // For export macro
#include "vtkCameraNode.h"

#include <optixu/optixpp_namespace.h>

class VTKRENDERINGOPTIX_EXPORT vtkOptiXCameraNode :
  public vtkCameraNode
{
public:
  static vtkOptiXCameraNode* New();
  vtkTypeMacro(vtkOptiXCameraNode, vtkCameraNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //Description:
  //Make ospray calls to render me.
  virtual void Render(bool prepass);

protected:
  vtkOptiXCameraNode();
  ~vtkOptiXCameraNode();

private:
  vtkOptiXCameraNode(const vtkOptiXCameraNode&) VTK_DELETE_FUNCTION;
  void operator=(const vtkOptiXCameraNode&) VTK_DELETE_FUNCTION;

  optix::Program RayGenProgram;
};

#endif
