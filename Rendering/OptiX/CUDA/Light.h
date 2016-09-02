

#pragma once

#include <optix_world.h>


namespace vtkopt
{

struct Light
{
#if defined(__cplusplus)                                                         
      typedef optix::float3 float3;                                                  
#endif 

    enum 
    {
        POSITIONAL=0,
        DIRECTIONAL
    };

    float3 pos;
    float3 dir;
    float3 color;
    int    type;
};

} // end namespace vtkopt
