

#pragma once

#include <optix_world.h>


#define RADIANCE_RAY_TYPE   0
#define OCCLUSION_RAY_TYPE  1

struct RadiancePRD
{
    float3 result;
};


struct OcclusionPRD
{
    float3 occlusion;
};


struct IntersectionRecord 
{
    float3  N;
    float3  Ng;
};
