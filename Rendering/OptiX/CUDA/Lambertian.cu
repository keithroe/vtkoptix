
#include <optix.h>
#include <optixu/optixu_math_namespace.h>

#include "Common.h"
#include "Light.h"


#include <stdio.h>
rtDeclareVariable(uint2, launch_index, rtLaunchIndex, );

using namespace optix;



rtDeclareVariable(rtObject, top_object, , );
rtBuffer<vtkopt::Light> lights;

rtDeclareVariable( RadiancePRD, prd,   rtPayload, );
rtDeclareVariable( optix::Ray,  ray,   rtCurrentRay, );
rtDeclareVariable( float,       t_hit, rtIntersectionDistance, );

rtDeclareVariable( IntersectionRecord, irec, attribute irec, ); 

RT_PROGRAM void LambertianClosestHit()
{
    const float3 N = faceforward( irec.N, -ray.direction, irec.Ng );
    const float3 P  = ray.origin + t_hit * ray.direction;
    const float3 Kd = make_float3( 0.7f, 0.7f, 0.7f );

    // light loop
    float3 Lc = make_float3( 0.0f );
    const int num_lights = lights.size();
    for( int i =0; i < num_lights; ++i )
    {
        const vtkopt::Light light = lights[i];
        float3 L;
        float  Ldist;
        float3  Lcolor;
        if( light.type == vtkopt::Light::DIRECTIONAL )
        {
            L = -light.dir;
            Ldist = 1e8f;
            Lcolor = light.color;
        }
        else
        {
            Ldist = optix::length( light.pos - P );
            L = ( light.pos-P ) / Ldist;
            Lcolor = light.color/(Ldist*Ldist);
        }

        const float N_dot_L = optix::dot( L, N );

        //
        // Calculation occlusion
        //
        float3 light_attenuation = make_float3(static_cast<float>( N_dot_L > 0.0f ));
        if( N_dot_L > 0.0f )
        {
            OcclusionPRD shadow_prd;
            shadow_prd.occlusion = make_float3( 1.0f );
            optix::Ray shadow_ray = optix::make_Ray( P, L, OCCLUSION_RAY_TYPE, 0.001f, Ldist );
            rtTrace( top_object, shadow_ray, shadow_prd );
            light_attenuation = shadow_prd.occlusion;
        }

        //
        // Calculate local lighting 
        //
        if( fmaxf(light_attenuation) > 0.0f )
        {
            Lc += fmaxf( 0.0f, N_dot_L ) * Lcolor * light_attenuation;
        }

    }
    prd.result = Lc*Kd; 
    
}


rtDeclareVariable( OcclusionPRD, shadow_prd, rtPayload, );
RT_PROGRAM void LambertianAnyHit()
{
    shadow_prd.occlusion = make_float3( 0.0f );

}
