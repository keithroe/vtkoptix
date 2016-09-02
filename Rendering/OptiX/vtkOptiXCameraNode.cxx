/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXCameraNode.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkOptiXCameraNode.h"

#include "vtkCamera.h"
#include "vtkCollectionIterator.h"
#include "vtkObjectFactory.h"
#include "vtkOptiXConfig.h"
#include "vtkOptiXRendererNode.h"
#include "vtkRenderer.h"
#include "vtkViewNodeCollection.h"

#include <optixu/optixu_math_namespace.h>


using optix::float3;

//#include "ospray/ospray.h"

//============================================================================
vtkStandardNewMacro(vtkOptiXCameraNode);

//----------------------------------------------------------------------------
vtkOptiXCameraNode::vtkOptiXCameraNode()
{
}

//----------------------------------------------------------------------------
vtkOptiXCameraNode::~vtkOptiXCameraNode()
{
}

//----------------------------------------------------------------------------
void vtkOptiXCameraNode::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkOptiXCameraNode::Render(bool prepass)
{

    std::cerr << "==================================================================" << std::endl;
    std::cerr << "Camera node render - prepass: " << prepass << std::endl;

    if (prepass)
    {
        vtkOptiXRendererNode *orn =
            static_cast<vtkOptiXRendererNode *>(
                    this->GetFirstAncestorOfType("vtkOptiXRendererNode"));
        
        if( !this->RayGenProgram )
        {
            optix::Context context = orn->GetOptiXContext();
            this->RayGenProgram = context->createProgramFromPTXFile(
                    VTK_OPTIX_PTX_DIR + "/cuda_compile_ptx_generated_PerspectiveCamera.cu.ptx",
                    "PerspectiveCameraRayGen"
                    );
            context->setRayGenerationProgram( 0, this->RayGenProgram );
        }

        vtkRenderer* ren = vtkRenderer::SafeDownCast(orn->GetRenderable());
        int tiledSize[2];
        int tiledOrigin[2];
        ren->GetTiledSizeAndOrigin(
                &tiledSize[0], &tiledSize[1],
                &tiledOrigin[0], &tiledOrigin[1]);

        vtkCamera* cam = static_cast<vtkCamera *>(this->Renderable);
        const float  fovy   = cam->GetViewAngle();
        const float  aspect = static_cast<float>( tiledSize[0] ) /
                              static_cast<float>( tiledSize[1] );
        const float3 pos    = optix::make_float3( cam->GetPosition()[0],
                                                  cam->GetPosition()[1],
                                                  cam->GetPosition()[2] ) ;
        const float3 dir    = optix::make_float3( cam->GetDirectionOfProjection()[0],
                                                  cam->GetDirectionOfProjection()[1],
                                                  cam->GetDirectionOfProjection()[2] ) ;
        const float3 up     = optix::make_float3( cam->GetViewUp()[0],
                                                  cam->GetViewUp()[1],
                                                  cam->GetViewUp()[2] ) ;

        const float vlen = tanf( 0.5f * fovy * M_PIf / 180.0f ); 
        const float ulen = vlen * aspect;
        const float3 cam_W   = optix::normalize( dir );
        const float3 cam_U   = optix::normalize( optix::cross( dir, up ) );
        const float3 cam_V   = optix::normalize( optix::cross( cam_U, cam_W ) );

        this->RayGenProgram[ "pos" ]->setFloat( pos );
        this->RayGenProgram[ "U"   ]->setFloat( cam_U*ulen );
        this->RayGenProgram[ "V"   ]->setFloat( cam_V*vlen );
        this->RayGenProgram[ "W"   ]->setFloat( cam_W );

        /*
        OSPCamera ospCamera = ospNewCamera("perspective");
        ospSetObject(orn->GetORenderer(),"camera", ospCamera);

        vtkCamera *cam = static_cast<vtkCamera *>(this->Renderable);
        ospSetf(ospCamera,"aspect", float(tiledSize[0])/float(tiledSize[1]));
        ospSetf(ospCamera,"fovy",cam->GetViewAngle());
        double *pos = cam->GetPosition();
        ospSet3f(ospCamera,"pos",pos[0], pos[1], pos[2]);
        ospSet3f(ospCamera,"up",
        cam->GetViewUp()[0], cam->GetViewUp()[1], cam->GetViewUp()[2]);
        double *dop = cam->GetDirectionOfProjection();
        ospSet3f(ospCamera,"dir", dop[0], dop[1], dop[2]);

        ospCommit(ospCamera);
        ospRelease(ospCamera);
        */
    }
}
