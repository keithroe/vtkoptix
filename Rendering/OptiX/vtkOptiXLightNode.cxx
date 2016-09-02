/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXLightNode.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOptiXLightNode.h"

#include "vtkCollectionIterator.h"
#include "vtkLight.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkOptiXRendererNode.h"

#include <CUDA/Light.h>

#include <vector>


//============================================================================
double vtkOptiXLightNode::LightScale = 1.0;

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkOptiXLightNode);

//----------------------------------------------------------------------------
vtkOptiXLightNode::vtkOptiXLightNode()
{
}

//----------------------------------------------------------------------------
vtkOptiXLightNode::~vtkOptiXLightNode()
{
}

//----------------------------------------------------------------------------
void vtkOptiXLightNode::SetLightScale(double s)
{
  vtkOptiXLightNode::LightScale = s;
}

//----------------------------------------------------------------------------
double vtkOptiXLightNode::GetLightScale()
{
  return vtkOptiXLightNode::LightScale;
}

//----------------------------------------------------------------------------
void vtkOptiXLightNode::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkOptiXLightNode::Render(bool prepass)
{
    std::cerr << "==================================================================" << std::endl;
    std::cerr << "Light node render -- prepass: " << prepass << std::endl;
    if (prepass)
    {
        vtkOptiXRendererNode *orn =
            static_cast<vtkOptiXRendererNode *>(
                    this->GetFirstAncestorOfType("vtkOptiXRendererNode"));

        vtkLight *vlight = vtkLight::SafeDownCast(this->GetRenderable());
        
        if( vlight->GetSwitch() == 0 || 
            vtkOptiXLightNode::LightScale <= 0.0f ||
            vlight->GetIntensity() <= 0.0 )
        {
            // Ignoring light 
            return;
        }
        
        const optix::float3 color = optix::make_float3(
                static_cast<float>( vlight->GetDiffuseColor()[0] ),
                static_cast<float>( vlight->GetDiffuseColor()[1] ),
                static_cast<float>( vlight->GetDiffuseColor()[2] )
                );

        const float intensity = static_cast<float>( 
                vtkOptiXLightNode::LightScale *
                vlight->GetIntensity()
                );

        vtkopt::Light light;
        light.color = color*intensity;
        light.pos   = optix::make_float3( 0.0f );
        light.dir   = optix::make_float3( 0.0f );

        if (vlight->GetPositional())
        {
            light.type = vtkopt::Light::POSITIONAL;
            light.pos  = optix::make_float3( 
                    static_cast<float>( vlight->GetPosition()[0] ),
                    static_cast<float>( vlight->GetPosition()[1] ),
                    static_cast<float>( vlight->GetPosition()[2] )
                    );
        }
        else
        {
            light.type = vtkopt::Light::DIRECTIONAL;
            light.dir  = optix::make_float3( 
                    static_cast<float>( vlight->GetFocalPoint()[0] - vlight->GetPosition()[0] ),
                    static_cast<float>( vlight->GetFocalPoint()[1] - vlight->GetPosition()[1] ),
                    static_cast<float>( vlight->GetFocalPoint()[2] - vlight->GetPosition()[2] )
                    );
            light.dir = optix::normalize( light.dir );
        }
        orn->AddLight( light );






        /*
        if (light->GetPositional())
        {
            OSPLight ospLight = ospNewLight(oRenderer, "PointLight");
            ospSet3f(ospLight, "color",
                    color[0],
                    color[1],
                    color[2]);
            float fI = static_cast<float>
                (vtkOptiXLightNode::LightScale*
                 light->GetIntensity()ght*
                 vtkMath::Pi() //since OSP 0.10.0
                );
            ospSet1f(ospLight, "intensity", fI);
            ospSet3f(ospLight, "position",
                    light->GetPosition()[0],
                    light->GetPosition()[1],
                    light->GetPosition()[2]);
            ospCommit(ospLight);
            orn->AddLight(ospLight);
        }
        else
        {
            double direction[3];
            direction[0] = light->GetPosition()[0] - light->GetFocalPoint()[0];
            direction[1] = light->GetPosition()[1] - light->GetFocalPoint()[1];
            direction[2] = light->GetPosition()[2] - light->GetFocalPoint()[2];
            OSPLight ospLight = ospNewLight(oRenderer, "DirectionalLight");
            ospSet3f(ospLight, "color",
                    color[0],
                    color[1],
                    color[2]);
            float fI = static_cast<float>
                (vtkOptiXLightNode::LightScale*
                 light->GetIntensity()*
                 vtkMath::Pi()); //since OSP 0.10.0
            ospSet1f(ospLight, "intensity", fI);
            vtkMath::Normalize(direction);
            ospSet3f(ospLight, "direction",
                    -direction[0],-direction[1],-direction[2]);
            ospCommit(ospLight);
            orn->AddLight(ospLight);
        }
        */
    }
}
