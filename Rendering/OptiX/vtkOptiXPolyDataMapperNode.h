/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXPolyDataMapperNode.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOptiXPolyDataMapperNode - links vtkActor and vtkMapper to OptiX
// .SECTION Description
// Translates vtkActor/Mapper state into OptiX rendering calls

#ifndef vtkOptiXPolyDataMapperNode_h
#define vtkOptiXPolyDataMapperNode_h

#include "vtkRenderingOptiXModule.h" // For export macro
#include "vtkPolyDataMapperNode.h"
#include "vtkProperty.h"

#include <optixu/optixpp_namespace.h>

class vtkOptiXActorNode;
class vtkOptiXRendererNode;
class vtkPolyData;


class VTKRENDERINGOPTIX_EXPORT vtkOptiXPolyDataMapperNode :
  public vtkPolyDataMapperNode
{
public:
  class Geom;

  static vtkOptiXPolyDataMapperNode* New();
  vtkTypeMacro(vtkOptiXPolyDataMapperNode, vtkPolyDataMapperNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //Description:
  //Make ospray calls to render me.
  virtual void Render(bool prepass);

protected:


  vtkOptiXPolyDataMapperNode();
  ~vtkOptiXPolyDataMapperNode();

  void RenderPoly( 
          vtkOptiXRendererNode* orn,
          vtkOptiXActorNode* aNode,
          vtkPolyData* poly
          );

  Geom* MyGeom;

  void CreateNewMeshes();

private:

  vtkOptiXPolyDataMapperNode(const vtkOptiXPolyDataMapperNode&) VTK_DELETE_FUNCTION;
  void operator=(const vtkOptiXPolyDataMapperNode&) VTK_DELETE_FUNCTION;
};
#endif
