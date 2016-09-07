/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXRendererNode.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOptiXRendererNode.h"

#include "vtkCamera.h"
#include "vtkCollectionIterator.h"
#include "vtkInformation.h"
#include "vtkInformationIntegerKey.h"
#include "vtkObjectFactory.h"
#include "vtkOptiXActorNode.h"
#include "vtkOptiXCameraNode.h"
#include "vtkOptiXConfig.h"
#include "vtkOptiXLightNode.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkViewNodeCollection.h"

//#include "ospray/ospray.h"
//#include "ospray/version.h"

#include <cmath>

#if 0
//debug includes
#include "vtkDataSetWriter.h"
#include "vtkImageImport.h"
#include "vtkSmartPointer.h"
#include "vtkWindowToImageFilter.h"
#include <unistd.h>
#endif

vtkInformationKeyMacro(vtkOptiXRendererNode, SAMPLES_PER_PIXEL, Integer);
vtkInformationKeyMacro(vtkOptiXRendererNode, MAX_FRAMES, Integer);
vtkInformationKeyMacro(vtkOptiXRendererNode, AMBIENT_SAMPLES, Integer);

//============================================================================
vtkStandardNewMacro(vtkOptiXRendererNode);

//----------------------------------------------------------------------------
vtkOptiXRendererNode::vtkOptiXRendererNode()
{
    std::cerr << "==================================================================" << std::endl;
    std::cerr << "RendererNode constructor" << std::endl;

  this->Buffer = NULL;
  this->ZBuffer = NULL;
  this->NumActors = 0;
}

//----------------------------------------------------------------------------
vtkOptiXRendererNode::~vtkOptiXRendererNode()
{
    std::cerr << "==================================================================" << std::endl;
    std::cerr << "RendererNode destructor" << std::endl;
    delete[] this->Buffer;
    delete[] this->ZBuffer;

    if( this->Context )
        this->Context->destroy();

//  ospRelease((OSPModel)this->OModel);
//  ospRelease((OSPRenderer)this->ORenderer);
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::SetSamplesPerPixel(int value, vtkRenderer *renderer)
{
  if (!renderer)
    {
    return;
    }
  vtkInformation *info = renderer->GetInformation();
  info->Set(vtkOptiXRendererNode::SAMPLES_PER_PIXEL(), value);
}

//----------------------------------------------------------------------------
int vtkOptiXRendererNode::GetSamplesPerPixel(vtkRenderer *renderer)
{
  if (!renderer)
    {
    return 1;
    }
  vtkInformation *info = renderer->GetInformation();
  if (info && info->Has(vtkOptiXRendererNode::SAMPLES_PER_PIXEL()))
    {
    return (info->Get(vtkOptiXRendererNode::SAMPLES_PER_PIXEL()));
    }
  return 1;
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::SetMaxFrames(int value, vtkRenderer *renderer)
{
  if (!renderer)
    {
    return;
    }
  vtkInformation *info = renderer->GetInformation();
  info->Set(vtkOptiXRendererNode::MAX_FRAMES(), value);
}

//----------------------------------------------------------------------------
int vtkOptiXRendererNode::GetMaxFrames(vtkRenderer *renderer)
{
  if (!renderer)
    {
    return 1;
    }
  vtkInformation *info = renderer->GetInformation();
  if (info && info->Has(vtkOptiXRendererNode::MAX_FRAMES()))
    {
    return (info->Get(vtkOptiXRendererNode::MAX_FRAMES()));
    }
  return 1;
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::SetAmbientSamples(int value, vtkRenderer *renderer)
{
  if (!renderer)
    {
    return;
    }
  vtkInformation *info = renderer->GetInformation();
  info->Set(vtkOptiXRendererNode::AMBIENT_SAMPLES(), value);
}

//----------------------------------------------------------------------------
int vtkOptiXRendererNode::GetAmbientSamples(vtkRenderer *renderer)
{
  if (!renderer)
    {
    return 0;
    }
  vtkInformation *info = renderer->GetInformation();
  if (info && info->Has(vtkOptiXRendererNode::AMBIENT_SAMPLES()))
    {
    return (info->Get(vtkOptiXRendererNode::AMBIENT_SAMPLES()));
    }
  return 0;
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkOptiXRendererNode::Traverse(int operation)
{
    std::cerr << "==================================================================" << std::endl;
    std::cerr << "RendererNode traverse" << std::endl;
    // do not override other passes
    if (operation != render)
    {
        this->Superclass::Traverse(operation);
        return;
    }

    this->Apply(operation,true);

    //TODO: this repeated traversal to find things of particular types
    //is bad, find something smarter
    vtkViewNodeCollection *nodes = this->GetChildren();
    vtkCollectionIterator *it = nodes->NewIterator();
    it->InitTraversal();
    while (!it->IsDoneWithTraversal())
    {
        vtkOptiXCameraNode *child =
            vtkOptiXCameraNode::SafeDownCast(it->GetCurrentObject());
        if (child)
        {
            child->Traverse(operation);
            break;
        }
        it->GoToNextItem();
    }

    //lights
    // TODO: add caching similar to meshes below
    this->Lights.clear();
    it->InitTraversal();
    while (!it->IsDoneWithTraversal())
    {
        vtkOptiXLightNode *child =
            vtkOptiXLightNode::SafeDownCast(it->GetCurrentObject());
        if (child)
        {
            child->Traverse(operation);
        }
        it->GoToNextItem();
    }
    this->LightBuffer->setSize( this->Lights.size() );
    std::cerr << "Setting " << this->Lights.size() << " lights:" << std::endl; 
    for( int i = 0; i < this->Lights.size(); ++i )
    {
        vtkopt::Light l = this->Lights[i];
        std::cerr << "light " << i << std::endl;
        std::cerr << "pos: " << l.pos.x << ", " << l.pos.y << ", " << l.pos.z << std::endl;
        std::cerr << "dir: " << l.dir.x << ", " << l.dir.y << ", " << l.dir.z << std::endl;
    }
    if( !this->Lights.empty() )
        memcpy( this->LightBuffer->map(), &this->Lights[0], this->Lights.size()*sizeof(vtkopt::Light) );
    this->LightBuffer->unmap();
    

    //actors
    it->InitTraversal();
    //since we have to spatially sort everything
    //let's see if we can avoid that in the common case when
    //the objects have not changed. Note we also cache in actornodes
    //to reuse already created ospray meshes
    unsigned int recent = 0;
    int numAct = 0; //catches removed actors
    while (!it->IsDoneWithTraversal())
    {
        vtkOptiXActorNode *child =
            vtkOptiXActorNode::SafeDownCast(it->GetCurrentObject());
        if (child)
        {
            numAct++;
            unsigned int mtime = child->GetMTime();
            if (mtime > recent)
            {
                recent = mtime;
            }
        }
        it->GoToNextItem();
    }

    bool enable_cache = true; //turn off to force rebuilds for debugging
    if ( !enable_cache ||
         (recent > this->RenderTime) ||
         (numAct != this->NumActors))
    {
        this->NumActors = numAct;
        it->InitTraversal();
        while (!it->IsDoneWithTraversal())
        {
            vtkOptiXActorNode *child =
                vtkOptiXActorNode::SafeDownCast(it->GetCurrentObject());
            if (child)
            {
                child->Traverse(operation);
            }
            it->GoToNextItem();
        }
        this->RenderTime = recent;
    }
    else
    {
    }
    it->Delete();

    /*

    OSPRenderer oRenderer = (osp::Renderer*)this->ORenderer;

    //camera
    //TODO: this repeated traversal to find things of particular types
    //is bad, find something smarter
    vtkViewNodeCollection *nodes = this->GetChildren();
    vtkCollectionIterator *it = nodes->NewIterator();
    it->InitTraversal();
    while (!it->IsDoneWithTraversal())
    {
        vtkOptiXCameraNode *child =
            vtkOptiXCameraNode::SafeDownCast(it->GetCurrentObject());
        if (child)
        {
            child->Traverse(operation);
            break;
        }
        it->GoToNextItem();
    }

    //lights
    this->Lights.clear();
    it->InitTraversal();
    while (!it->IsDoneWithTraversal())
    {
        vtkOptiXLightNode *child =
            vtkOptiXLightNode::SafeDownCast(it->GetCurrentObject());
        if (child)
        {
            child->Traverse(operation);
        }
        it->GoToNextItem();
    }
    OSPData lightArray = ospNewData(this->Lights.size(), OSP_OBJECT,
            (this->Lights.size()?&this->Lights[0]:NULL), 0);
    ospSetData(oRenderer, "lights", lightArray);

    //actors
    OSPModel oModel=NULL;
    it->InitTraversal();
    //since we have to spatially sort everything
    //let's see if we can avoid that in the common case when
    //the objects have not changed. Note we also cache in actornodes
    //to reuse already created ospray meshes
    unsigned int recent = 0;
    int numAct = 0; //catches removed actors
    while (!it->IsDoneWithTraversal())
    {
        vtkOptiXActorNode *child =
            vtkOptiXActorNode::SafeDownCast(it->GetCurrentObject());
        if (child)
        {
            numAct++;
            unsigned int mtime = child->GetMTime();
            if (mtime > recent)
            {
                recent = mtime;
            }
        }
        it->GoToNextItem();
    }

    bool enable_cache = true; //turn off to force rebuilds for debugging
    if (!this->OModel ||
            !enable_cache ||
            (recent > this->RenderTime) ||
            (numAct != this->NumActors))
    {
        this->NumActors = numAct;
        ospRelease((OSPModel)this->OModel);
        oModel = ospNewModel();
        this->OModel = oModel;
        it->InitTraversal();
        while (!it->IsDoneWithTraversal())
        {
            vtkOptiXActorNode *child =
                vtkOptiXActorNode::SafeDownCast(it->GetCurrentObject());
            if (child)
            {
                child->Traverse(operation);
            }
            it->GoToNextItem();
        }
        this->RenderTime = recent;
        ospSetObject(oRenderer,"model", oModel);
        ospCommit(oModel);
    }
    else
    {
        oModel = (OSPModel)this->OModel;
    }
    it->Delete();
  */

  this->Apply(operation,false);
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::Build(bool prepass)
{
    if (prepass)
    {
        vtkRenderer *aren = vtkRenderer::SafeDownCast(this->Renderable);
        // make sure we have a camera
        if ( !aren->IsActiveCameraCreated() )
        {
            aren->ResetCamera();
        }
    }
    this->Superclass::Build(prepass);
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::Render(bool prepass)
{
    /*
    if (prepass)
    {
        OSPRenderer oRenderer = NULL;
        if (!this->ORenderer)
        {
            ospRelease((osp::Renderer*)this->ORenderer);
            oRenderer = (osp::Renderer*)ospNewRenderer("scivis");
            this->ORenderer = oRenderer;
        }
        else
        {
            oRenderer = (osp::Renderer*)this->ORenderer;
        }

        vtkRenderer *ren = vtkRenderer::SafeDownCast(this->GetRenderable());
        int *tmp = ren->GetSize();
        this->Size[0] = tmp[0];
        this->Size[1] = tmp[1];
        if (ren->GetUseShadows())
        {
            ospSet1i(oRenderer,"shadowsEnabled",1);
        }
        else
        {
            ospSet1i(oRenderer,"shadowsEnabled",0);
        }
        ospSet1i(oRenderer,"aoSamples",
                this->GetAmbientSamples(static_cast<vtkRenderer*>(this->Renderable)));
        ospSet1i(oRenderer,"spp",
                this->GetSamplesPerPixel(static_cast<vtkRenderer*>(this->Renderable)));

        double *bg = ren->GetBackground();
        ospSet3f(oRenderer,"bgColor", bg[0], bg[1], bg[2]);
    }
    else
    {
        OSPRenderer oRenderer = (osp::Renderer*)this->ORenderer;
        ospCommit(oRenderer);

        osp::vec2i isize = {this->Size[0], this->Size[1]};
        OSPFrameBuffer osp_framebuffer = ospNewFrameBuffer
            (isize,
#if OPTIX_VERSION_MAJOR < 1 && OPTIX_VERSION_MINOR < 10 && OPTIX_VERSION_PATCH < 2
             OSP_RGBA_I8, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
#else
        OSP_FB_RGBA8, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
#endif
        ospSet1f(osp_framebuffer, "gamma", 1.0f);
        ospCommit(osp_framebuffer);
        ospFrameBufferClear(osp_framebuffer, OSP_FB_COLOR|OSP_FB_DEPTH|OSP_FB_ACCUM);
        for (int i = 0; i < this->GetMaxFrames(static_cast<vtkRenderer*>(this->Renderable)); i++)
        {
            ospRenderFrame(osp_framebuffer, oRenderer,
                    OSP_FB_COLOR|OSP_FB_DEPTH|OSP_FB_ACCUM);
        }
        const void* rgba = ospMapFrameBuffer(osp_framebuffer, OSP_FB_COLOR);
        delete[] this->Buffer;
        this->Buffer = new unsigned char[this->Size[0]*this->Size[1]*4];
        memcpy((void*)this->Buffer, rgba, this->Size[0]*this->Size[1]*sizeof(char)*4);
        ospUnmapFrameBuffer(rgba, osp_framebuffer);

        vtkCamera *cam = vtkRenderer::SafeDownCast(this->Renderable)->
            GetActiveCamera();
        double *clipValues = cam->GetClippingRange();
        double clipMin = clipValues[0];
        double clipMax = clipValues[1];
        double clipDiv = 1.0 / (clipMax - clipMin);

        const void *Z = ospMapFrameBuffer(osp_framebuffer, OSP_FB_DEPTH);
        delete[] this->ZBuffer;
        this->ZBuffer = new float[this->Size[0]*this->Size[1]];
        float *s = (float *)Z;
        float *d = this->ZBuffer;
#if 0    
        float minS = 1000.0;
        float maxS = -1000.0;
        float minD = 1000.0;
        float maxD = -10000.0;
#endif
        for (int i = 0; i < (this->Size[0]*this->Size[1]); i++, s++, d++)
        {
            *d = (*s<clipMin? 1.0 : (*s - clipMin) * clipDiv);
#if 0    
            if (*d < minD) minD = *d;
            if (*d > maxD) maxD = *d;
            if (*s < minS) minS = *s;
            if (*s > maxS) maxS = *s;
#endif
        }
        ospUnmapFrameBuffer(Z, osp_framebuffer);

        ospRelease(osp_framebuffer);
    }
*/
    if (prepass)
    {
        std::cerr << "==================================================================" << std::endl;
        std::cerr << "RendererNode prepass" << std::endl;
        vtkRenderer *aren = vtkRenderer::SafeDownCast(this->Renderable);
        std::cerr << "\tMaxFrames: " << vtkOptiXRendererNode::GetMaxFrames( aren ) << std::endl;
        std::cerr << "\tSPP      : " << vtkOptiXRendererNode::GetSamplesPerPixel( aren ) << std::endl;
        
        vtkRenderer *ren = vtkRenderer::SafeDownCast(this->GetRenderable());
        int *tmp = ren->GetSize();
        this->Size[0] = tmp[0];
        this->Size[1] = tmp[1];

        // Initialize optix
        if( !this->Context ) 
        {
            this->Context = optix::Context::create();
            this->Context->setRayTypeCount( 2 );
            this->Context->setEntryPointCount( 1 );

            this->FrameBuffer = this->Context->createBuffer(
                    RT_BUFFER_OUTPUT,
                    RT_FORMAT_UNSIGNED_BYTE4, 
                    this->Size[0],
                    this->Size[1]
                    );
            
            this->DepthBuffer = this->Context->createBuffer(
                    RT_BUFFER_OUTPUT,
                    RT_FORMAT_FLOAT, 
                    this->Size[0],
                    this->Size[1]
                    );
            
            this->LightBuffer= this->Context->createBuffer(
                    RT_BUFFER_INPUT,
                    RT_FORMAT_USER
                    );
            this->LightBuffer->setElementSize( sizeof( vtkopt::Light ) );
            this->LightBuffer->setSize( 0 );

            this->GeometryGroup = this->Context->createGeometryGroup();
            this->GeometryGroup->setAcceleration( this->Context->createAcceleration( "Trbvh" ) );

            this->Context[ "frame_buffer" ]->setBuffer( this->FrameBuffer );
            this->Context[ "depth_buffer" ]->setBuffer( this->DepthBuffer );
            this->Context[ "lights"       ]->setBuffer( this->LightBuffer );
            this->Context[ "top_object"   ]->set( this->GeometryGroup );


            /*
            this->Context->setRayGenerationProgram(
                    0,
                    this->Context->createProgramFromPTXFile(
                        VTK_OPTIX_PTX_DIR + "/cuda_compile_ptx_generated_draw_color.cu.ptx",
                        "draw_color"
                        )
                    );
            */

            this->Context->validate();
        }

        if (ren->GetUseShadows())
        {
            //this->context[ "shadows_enabled" ]->setInt( 1 );
        }
        else
        {
            //this->context[ "shadows_enabled" ]->setInt( 1 );
        }
        // this->GetAmbientSamples( static_cast<vtkRenderer*>(this->Renderable));
        // this->GetSamplesPerPixel(static_cast<vtkRenderer*>(this->Renderable)));

        //double *bg = ren->GetBackground();
        //ospSet3f(oRenderer,"bgColor", bg[0], bg[1], bg[2]);
    }
    else
    {
        // DO RENDER
        this->Context->launch( 0, this->Size[0], this->Size[1] );
        
        // TODO: avoid realloc if size matches
        delete[] this->Buffer;
        this->Buffer = new unsigned char[this->Size[0]*this->Size[1]*4];
        memcpy( this->Buffer, this->FrameBuffer->map(),  this->Size[0]*this->Size[1]*4 );
        this->FrameBuffer->unmap();
        //memset( this->Buffer, 0,  this->Size[0]*this->Size[1]*4 );

        vtkCamera *cam = vtkRenderer::SafeDownCast(this->Renderable)->
            GetActiveCamera();
        double *clipValues = cam->GetClippingRange();
        double clipMin = clipValues[0];
        double clipMax = clipValues[1];
        double clipDiv = 1.0 / (clipMax - clipMin);

        // TODO: avoid realloc if size matches
        delete[] this->ZBuffer;
        this->ZBuffer = new float[this->Size[0]*this->Size[1]];
        memcpy( this->ZBuffer, this->DepthBuffer->map(),  this->Size[0]*this->Size[1]*4 );
        this->DepthBuffer->unmap();
        //memset( this->Buffer, 0,  this->Size[0]*this->Size[1]*4 );

        std::cerr << "================================================================" << std::endl;
        std::cerr << "Buffers initialized" << std::endl;

        /*
        float *s = (float *)Z; // osp depth buffer
        float *d = this->ZBuffer;
        for (int i = 0; i < (this->Size[0]*this->Size[1]); i++, s++, d++)
        {
            *d = (*s<clipMin? 1.0 : (*s - clipMin) * clipDiv);
        }
        */
    }
}

//----------------------------------------------------------------------------
void vtkOptiXRendererNode::WriteLayer(unsigned char *buffer, float *Z,
                                       int buffx, int buffy, int layer)
{
  if (layer == 0)
    {
    for (int j = 0; j < buffy && j < this->Size[1]; j++)
      {
      unsigned char *iptr = this->Buffer + j*this->Size[0]*4;
      float *zptr = this->ZBuffer + j*this->Size[0];
      unsigned char *optr = buffer + j*buffx*4;
      float *ozptr = Z +  j*buffx;
      for (int i = 0; i < buffx && i < this->Size[0]; i++)
        {
        *optr++ = *iptr++;
        *optr++ = *iptr++;
        *optr++ = *iptr++;
        *optr++ = *iptr++;
        *ozptr++ = *zptr;
        zptr++;
        }
      }
    }
  else
    {
    for (int j = 0; j < buffy && j < this->Size[1]; j++)
      {
      unsigned char *iptr = this->Buffer + j*this->Size[0]*4;
      float *zptr = this->ZBuffer + j*this->Size[0];
      unsigned char *optr = buffer + j*buffx*4;
      float *ozptr = Z +  j*buffx;
      for (int i = 0; i < buffx && i < this->Size[0]; i++)
        {
        if (*zptr<1.0)
          {
          *optr++ = *iptr++;
          *optr++ = *iptr++;
          *optr++ = *iptr++;
          *optr++ = *iptr++;
          *ozptr = *zptr;
          }
        else
          {
          optr+=4;
          iptr+=4;
          }
        ozptr++;
        zptr++;
        }
      }
    }
}
