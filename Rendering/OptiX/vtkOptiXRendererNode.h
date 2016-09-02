/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXRendererNode.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkOptiXRendererNode - links vtkRenderers to OptiX
// .SECTION Description
// Translates vtkRenderer state into OptiX rendering calls

#ifndef vtkOptiXRendererNode_h
#define vtkOptiXRendererNode_h

#include "vtkRenderingOptiXModule.h" // For export macro
#include "vtkRendererNode.h"
#include "CUDA/Light.h"

#include <optixu/optixpp_namespace.h>

#include <vector> // for ivars

class vtkRenderer;
class vtkInformationIntegerKey;


class VTKRENDERINGOPTIX_EXPORT vtkOptiXRendererNode : public vtkRendererNode
{
public:
    static vtkOptiXRendererNode* New();
    vtkTypeMacro(vtkOptiXRendererNode, vtkRendererNode);
    void PrintSelf(ostream& os, vtkIndent indent);

    //Description:
    //Builds myself.
    virtual void Build(bool prepass);

    //Description:
    //Traverse graph in ospray's prefered order and render
    virtual void Render(bool prepass);

    //Description:
    //Put my results into the correct place in the provided pixel buffer.
    virtual void WriteLayer(unsigned char *buffer, float *zbuffer,
            int buffx, int buffy, int layer);

    //state beyond rendering core...

    //Description:
    //When present on renderer, controls the number of primary rays
    //shot per pixel
    //default is 1
    static vtkInformationIntegerKey* SAMPLES_PER_PIXEL();

    //Description:
    //Convenience method to set/get SAMPLES_PER_PIXEL on a vtkRenderer.
    static void SetSamplesPerPixel(int, vtkRenderer *renderer);
    static int GetSamplesPerPixel(vtkRenderer *renderer);

    //Description:
    //When present on renderer, controls the number of ospray render calls
    //for each refresh.
    //default is 1
    static vtkInformationIntegerKey* MAX_FRAMES();
    static void SetMaxFrames(int, vtkRenderer *renderer);
    static int GetMaxFrames(vtkRenderer *renderer);

    //Description:
    //When present on renderer, controls the number of ambient occlusion
    //samples shot per hit.
    //default is 4
    static vtkInformationIntegerKey* AMBIENT_SAMPLES();
    //Description:
    //Convenience method to set/get SAMPLES_PER_PIXEL on a vtkRenderer.
    static void SetAmbientSamples(int, vtkRenderer *renderer);
    static int GetAmbientSamples(vtkRenderer *renderer);

    // Description:
    // Get the last rendered ColorBuffer
    virtual unsigned char *GetBuffer() { return this->Buffer; }

    // Description:
    // Get the last rendered ZBuffer
    virtual float *GetZBuffer() { return this->ZBuffer; }

    // Description:
    // Get the OptiX Context 
    virtual optix::Context GetOptiXContext() { return this->Context; }
    
    // Description:
    // Get the top level GeometryGroup 
    virtual optix::GeometryGroup GetOptiXGeometryGroup() { return this->GeometryGroup; }

    // Description:
    // Add a Light 
    virtual void AddLight( const vtkopt::Light& light ) { this->Lights.push_back( light ); }


    // if you want to traverse your children in a specific order
    // or way override this method
    virtual void Traverse(int operation);

protected:
    vtkOptiXRendererNode();
    ~vtkOptiXRendererNode();

    //internal structures
    unsigned char *Buffer;
    float *ZBuffer;

    optix::Context       Context;
    optix::GeometryGroup GeometryGroup;
    optix::Buffer        FrameBuffer;
    optix::Buffer        DepthBuffer;
    optix::Buffer        LightBuffer;

    std::vector<vtkopt::Light> Lights;

    int NumActors;
private:
    vtkOptiXRendererNode(const vtkOptiXRendererNode&) VTK_DELETE_FUNCTION;
    void operator=(const vtkOptiXRendererNode&) VTK_DELETE_FUNCTION;
};

#endif
