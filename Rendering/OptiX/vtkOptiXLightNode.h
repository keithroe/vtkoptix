/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXLightNode.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOptiXLightNode - links vtkLights to OptiX
// .SECTION Description
// Translates vtkLight state into OptiX rendering calls

#ifndef vtkOptiXLightNode_h
#define vtkOptiXLightNode_h

#include "vtkRenderingOptiXModule.h" // For export macro
#include "vtkLightNode.h"

class VTKRENDERINGOPTIX_EXPORT vtkOptiXLightNode :
  public vtkLightNode
{
public:
  static vtkOptiXLightNode* New();
  vtkTypeMacro(vtkOptiXLightNode, vtkLightNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //Description:
  //Make ospray calls to render me.
  virtual void Render(bool prepass);

  //Description:
  //A global multiplier to all ospray lights.
  //default is 1.0
  static void SetLightScale(double s);
  static double GetLightScale();

protected:
  vtkOptiXLightNode();
  ~vtkOptiXLightNode();

private:
  vtkOptiXLightNode(const vtkOptiXLightNode&) VTK_DELETE_FUNCTION;
  void operator=(const vtkOptiXLightNode&) VTK_DELETE_FUNCTION;

  static double LightScale;
};

#endif
