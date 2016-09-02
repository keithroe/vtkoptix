
#include <optix.h>
#include <optixu/optixu_math_namespace.h>

#include "Common.h"

using namespace optix;


rtBuffer<float3>    vertices;
rtBuffer<int3>      triangles;

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

rtDeclareVariable( IntersectionRecord, irec, attribute irec, ); 


RT_PROGRAM void TriangleMeshIntersection( int prim_idx )
{
    const int3 v_idx = triangles[prim_idx];

    const float3 p0 = vertices[ v_idx.x ];
    const float3 p1 = vertices[ v_idx.y ];
    const float3 p2 = vertices[ v_idx.z ];

    // Intersect ray with triangle
    float3 n;
    float  t, beta, gamma;
    if( intersect_triangle( ray, p0, p1, p2, n, t, beta, gamma ) )
    {
        if(  rtPotentialIntersection( t ) )
        {
            irec.Ng = normalize( n );
            rtReportIntersection( 0 );
        }
    }
}


RT_PROGRAM void TriangleMeshBoundingBox( int prim_idx, float result[6] )
{
    const int3 v_idx = triangles[prim_idx];

    const float3 v0   = vertices[ v_idx.x ];
    const float3 v1   = vertices[ v_idx.y ];
    const float3 v2   = vertices[ v_idx.z ];
    const float  area = length(cross(v1-v0, v2-v0));
  
    optix::Aabb* aabb = (optix::Aabb*)result;
  
    if( area > 0.0f && !isinf(area) )
    {
        aabb->m_min = fminf( fminf( v0, v1), v2 );
        aabb->m_max = fmaxf( fmaxf( v0, v1), v2 );
    }
    else
    {
        aabb->invalidate();
    }
}
