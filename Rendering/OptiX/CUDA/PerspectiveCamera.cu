
#include <optix.h>
#include <optixu/optixu_math_namespace.h>

#include "Common.h"

using namespace optix;


rtDeclareVariable(uint2, launch_index, rtLaunchIndex, );
rtDeclareVariable(uint2, launch_dim,   rtLaunchDim, );

rtBuffer<uchar4, 2>   frame_buffer;
rtBuffer<float, 2>    depth_buffer;

rtDeclareVariable(rtObject,      top_object, , );
rtDeclareVariable(float3,        pos, , );
rtDeclareVariable(float3,        U, , );
rtDeclareVariable(float3,        V, , );
rtDeclareVariable(float3,        W, , );

RT_PROGRAM void PerspectiveCameraRayGen()
{
  
    const float2 d = make_float2(launch_index) / make_float2(launch_dim) * 2.f - 1.f;
    const float3 ray_origin    = pos;
    const float3 ray_direction = normalize(d.x*U + d.y*V + W);

    optix::Ray ray = optix::make_Ray(
            ray_origin,
            ray_direction,
            RADIANCE_RAY_TYPE,
            0.001f,
            RT_DEFAULT_MAX
            );


    RadiancePRD prd;
    prd.result = make_float3( 0.0f ); //ray_direction*0.5f + make_float3( 0.5f );

    rtTrace( top_object, ray, prd );

    const float3 c = fminf( prd.result, make_float3( 1.0f ) ); 
    frame_buffer[launch_index] = make_uchar4( c.x*255.99f, c.y*255.99f, c.z*255.99f, 255 );
    depth_buffer[launch_index] = 0.0f; 
}
