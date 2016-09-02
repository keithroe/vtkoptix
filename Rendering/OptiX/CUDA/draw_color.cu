
#include <optix.h>
#include <optixu/optixu_math_namespace.h>

rtDeclareVariable(uint2, launch_index, rtLaunchIndex, );
rtBuffer<uchar4, 2>   frame_buffer;
rtBuffer<float, 2>    depth_buffer;

RT_PROGRAM void draw_color()
{
    frame_buffer[launch_index] = make_uchar4(0, 255, 0, 255 );
    depth_buffer[launch_index] = 0.0f; 
}
