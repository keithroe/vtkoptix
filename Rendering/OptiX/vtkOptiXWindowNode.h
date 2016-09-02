/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXWindowNode.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOptiXWindowNode - links vtkRendererWindows to OptiX
// .SECTION Description
// Translates vtkRenderWindow state into OptiX rendering calls

#ifndef vtkOptiXWindowNode_h
#define vtkOptiXWindowNode_h

#include "vtkRenderingOptiXModule.h" // For export macro
#include "vtkWindowNode.h"

class VTKRENDERINGOPTIX_EXPORT vtkOptiXWindowNode :
  public vtkWindowNode
{
public:
  static vtkOptiXWindowNode* New();
  vtkTypeMacro(vtkOptiXWindowNode, vtkWindowNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //Description:
  //Make ospray calls to render me.
  virtual void Render(bool prepass);

protected:
  vtkOptiXWindowNode();
  ~vtkOptiXWindowNode();

private:
  vtkOptiXWindowNode(const vtkOptiXWindowNode&) VTK_DELETE_FUNCTION;
  void operator=(const vtkOptiXWindowNode&) VTK_DELETE_FUNCTION;
};

#endif
