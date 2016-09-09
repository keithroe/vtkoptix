
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

rtDeclareVariable( float3, Ks, , );
rtDeclareVariable( float3, Kd, , );
rtDeclareVariable( float,  Ns, , );

rtDeclareVariable( IntersectionRecord, irec, attribute irec, ); 

RT_PROGRAM void LambertianClosestHit()
{
    const float3 N = faceforward( irec.N, -ray.direction, irec.Ng );
    const float3 P  = ray.origin + t_hit * ray.direction;
    //const float3 Kd = make_float3( 0.7f, 0.7f, 0.7f );

    // light loop
    float3 color = make_float3( 0.0f );
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

        float3 light_attenuation = make_float3( 0.0f );
        if( N_dot_L > 0.0f )
        {
            //
            // Calculation occlusion
            //
            OcclusionPRD shadow_prd;
            shadow_prd.occlusion = make_float3( 1.0f );
            optix::Ray shadow_ray = optix::make_Ray( P, L, OCCLUSION_RAY_TYPE, 0.001f, Ldist );
            rtTrace( top_object, shadow_ray, shadow_prd );
            light_attenuation = shadow_prd.occlusion;

            //
            // Calculate local lighting 
            //
            if( fmaxf(light_attenuation) > 0.0f )
            {
                //const float3 H = optix::normalize( L - ray.direction );
                //const float  N_dot_H = optix::dot( N, H );
                const float3 R = optix::reflect( ray.direction, N );
                const float  L_dot_R = fmaxf( optix::dot( L, R ), 0.0f );;
                color += ( Kd*N_dot_L + Ks*powf( L_dot_R, Ns ) ) * Lcolor * light_attenuation;
            }
        }


    }
    prd.result = color; 
    
}


rtDeclareVariable( OcclusionPRD, shadow_prd, rtPayload, );
RT_PROGRAM void LambertianAnyHit()
{
    shadow_prd.occlusion = make_float3( 0.0f );

}
