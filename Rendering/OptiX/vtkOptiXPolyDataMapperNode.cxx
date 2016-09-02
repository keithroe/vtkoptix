/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOptiXPolyDataMapperNode.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkOptiXPolyDataMapperNode.h"

#include "vtkActor.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformationDoubleKey.h"
#include "vtkInformation.h"
#include "vtkInformationObjectBaseKey.h"
#include "vtkMapper.h"
#include "vtkObjectFactory.h"
#include "vtkOptiXActorNode.h"
#include "vtkOptiXConfig.h"
#include "vtkOptiXRendererNode.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkTexture.h"

#include <optixu/optixu_math_namespace.h>

//#include "ospray/ospray.h"
//#include "ospray/common/OSPCommon.h"
//#include "ospray/version.h"

#include <map>

//============================================================================


//------------------------------------------------------------------------------
class vtkOptiXPolyDataMapperNode::Geom
{
//A cache for the ospray meshes we make for this actor.
//When something else in the scene changes but this actor doesn't we
//reuse instead of recreating. RendererNode has a higher level cache
//that prevent spatial sorting when nothing changes other than camera.
public:
    std::vector<optix::GeometryInstance> gis;
    std::vector<optix::Buffer>           buffers;

    ~Geom()
    {
        /*
        for( std::vector<optix::GeometryInstance>::iterator it = gis.begin();
                it != gis.end();
                ++it )
            (*it)->destroy();

        for( std::vector<optix::Buffer>::iterator it = buffers.begin();
                it != buffers.end();
                ++it )
            (*it)->destroy();
        */
    }

    void Add( optix::GeometryInstance gi )
    {
        gis.push_back( gi );
    }
    
    void AddBuffer( optix::Buffer buffer )
    {
        buffers.push_back( buffer );
    }

    void AddMyselfTo( optix::GeometryGroup geom_group )
    {
        for( std::vector<optix::GeometryInstance>::iterator it = gis.begin();
                it != gis.end();
                ++it )
        {
            geom_group->addChild( *it );
        }
    }
};


  //------------------------------------------------------------------------------
optix::Buffer vtkDataArrayToBuffer( vtkDataArray* vdata, optix::Context ctx )
{

    const int num_tuples = vdata ? vdata->GetNumberOfTuples() : 0;
    optix::Buffer buff = ctx->createBuffer(
            RT_BUFFER_INPUT,
            RT_FORMAT_FLOAT3,
            num_tuples
            );

    optix::float3* buff_data = reinterpret_cast<optix::float3*>( buff->map() );
    for (int i = 0; i < num_tuples; i++)
    {
        const double* tuple = vdata->GetTuple(i);
        buff_data[i] = optix::make_float3 (
            static_cast<float>(tuple[0]),
            static_cast<float>(tuple[1]),
            static_cast<float>(tuple[2])
            );
    }
    buff->unmap();
    return buff;
}

  //------------------------------------------------------------------------------

  /*
  void CellColorMaterials(
    vtkUnsignedCharArray *vColors,
    OSPRenderer oRenderer,
    std::vector<OSPMaterial> &ospMaterials,
    float *specColor, float specPower,
    float opacity)
  {
    int numColors = vColors->GetNumberOfTuples();
    for (int i = 0; i < numColors; i++)
      {
      double *color = vColors->GetTuple(i);
      OSPMaterial oMaterial;
      oMaterial = ospNewMaterial(oRenderer,"RayTraceMaterial");
      float diffusef[] =
        {
          static_cast<float>(color[0])/(255.0f),
          static_cast<float>(color[1])/(255.0f),
          static_cast<float>(color[2])/(255.0f)
        };
      ospSet3fv(oMaterial,"Kd",diffusef);
      float specAdjust = 2.0f/(2.0f+specPower); //since OSP 0.10.0
      float specularf[] = {
        specColor[0]*specAdjust,
        specColor[1]*specAdjust,
        specColor[2]*specAdjust
      };
      ospSet3fv(oMaterial,"Ks",specularf);
      ospSet1f(oMaterial,"Ns",specPower);
      ospSet1f(oMaterial,"d", opacity);
      ospCommit(oMaterial);
      ospMaterials.push_back(oMaterial);
      }
  }
  */


  /*

  //----------------------------------------------------------------------------
  OSPMaterial VTKToOSPTexture(vtkImageData *vColorTextureMap,
                              OSPRenderer oRenderer)
  {
    int xsize = vColorTextureMap->GetExtent()[1];
    int ysize = vColorTextureMap->GetExtent()[3];
    unsigned char *ichars =
      (unsigned char *)vColorTextureMap->GetScalarPointer();
    unsigned char *ochars = new unsigned char[(xsize+1)*(ysize+1)*3];
    unsigned char *oc = ochars;
    int comps = vColorTextureMap->GetNumberOfScalarComponents();
    for (int i = 0; i <= xsize; i++)
      {
      for (int j = 0; j <= ysize; j++)
        {
        oc[0] = ichars[0];
        oc[1] = ichars[1];
        oc[2] = ichars[2];
        oc+=3;
        ichars+=comps;
        }
      }
    osp::Texture2D *t2d;
    t2d = (osp::Texture2D*)ospNewTexture2D
       (
        osp::vec2i{xsize+1,
          ysize+1},
          OSP_TEXTURE_RGB8,
       ochars,
       OSP_TEXTURE_FILTER_NEAREST);

    OSPMaterial ospMaterial;
    ospMaterial = ospNewMaterial(oRenderer,"RayTraceMaterial");
    ospSetObject(ospMaterial, "map_Kd", ((OSPTexture2D)(t2d)));
    ospCommit(t2d);
    ospCommit(ospMaterial);
    return ospMaterial;
  }
  */

float MapThroughPWF(double in, vtkPiecewiseFunction *scaleFunction)
{
    double out = in;
    if (!scaleFunction)
    {
        out = in;
    }
    else
    {
        out = scaleFunction->GetValue(in);
    }
    return static_cast<float>(out);
}

  //----------------------------------------------------------------------------
  /*
  OSPGeometry RenderAsSpheres(osp::vec3fa *vertices,
                              std::vector<unsigned int> &indexArray,
                              std::vector<unsigned int> &rIndexArray,
                              double pointSize,
                              vtkDataArray* scaleArray,
                              vtkPiecewiseFunction *scaleFunction,
                              OSPMaterial actorMaterial,
                              int numCellMaterials,
                              OSPData CellMaterials,
                              int numPointColors,
                              osp::vec4f *PointColors,
                              OSPModel oModel
                              )
  {
    OSPGeometry ospMesh = ospNewGeometry("spheres");
    int width=4;
    int scaleOffset = -1;
    if (scaleArray != NULL)
      {
      width = 5;
      scaleOffset = 4*sizeof(float);
      }
    float *mdata = new float[width*indexArray.size()];
    int *idata = (int*)mdata;
    for (size_t i = 0; i < indexArray.size(); i++)
      {
      mdata[i*width+0] = static_cast<float>(vertices[indexArray[i]].x);
      mdata[i*width+1] = static_cast<float>(vertices[indexArray[i]].y);
      mdata[i*width+2] = static_cast<float>(vertices[indexArray[i]].z);
      int mat = 0;
      if (numCellMaterials)
        {
        mat = rIndexArray[i];
        }
      else if (numPointColors)
        {
        mat = indexArray[i];
        }
      idata[i*width+3] = mat;
      if (scaleArray != NULL)
        {
        mdata[i*width+4] = MapThroughPWF
          (*scaleArray->GetTuple(indexArray[i]),
           scaleFunction);
        }
      }
    OSPData _PointColors = NULL;
    if (numPointColors)
      {
      _PointColors = ospNewData(numPointColors, OSP_FLOAT4, &PointColors[0]);
      }
    OSPData _mdata = ospNewData(indexArray.size()*width, OSP_FLOAT, mdata);
    ospSetObject(ospMesh, "spheres", _mdata);
    ospSet1i(ospMesh, "bytes_per_sphere", width*sizeof(float));
    ospSet1i(ospMesh, "offset_center", 0*sizeof(float));
    ospSet1f(ospMesh, "radius", pointSize);
    ospSet1i(ospMesh, "offset_radius", scaleOffset);

    if (numCellMaterials)
      {
      //per cell color
      ospSet1i(ospMesh, "offset_materialID", 3*sizeof(float));
      ospSetData(ospMesh, "materialList", CellMaterials);
      }
    else if (numPointColors)
      {
      //per point color
      ospSet1i(ospMesh, "offset_colorID", 3*sizeof(float));
      ospSetData(ospMesh, "color", _PointColors);
      }
    else
      {
      //per actor color
      ospSetMaterial(ospMesh, actorMaterial);
      }

    ospAddGeometry(oModel, ospMesh);
    ospCommit(ospMesh);
    //ospCommit(oModel); //TODO: crashes without yet others don't, why?
    ospRelease(_PointColors);
    ospRelease(_mdata);
    delete[] mdata;

    return ospMesh;
  }
*/

    /*
  //----------------------------------------------------------------------------
  OSPGeometry RenderAsCylinders(osp::vec3fa *vertices,
                                std::vector<unsigned int> &indexArray,
                                std::vector<unsigned int> &rIndexArray,
                                double lineWidth,
                                vtkDataArray *scaleArray,
                                vtkPiecewiseFunction *scaleFunction,
                                OSPMaterial actorMaterial,
                                int numCellMaterials,
                                OSPData CellMaterials,
                                int numPointColors,
                                osp::vec4f *PointColors,
                                OSPModel oModel)
  {
    OSPGeometry ospMesh = ospNewGeometry("cylinders");
    int width=7;
    int scaleOffset = -1;
    if (scaleArray != NULL)
      {
      width = 8;
      scaleOffset = 7*sizeof(float);
      }
    float *mdata = new float[indexArray.size()/2*width];
    int *idata = (int*)mdata;
    for (size_t i = 0; i < indexArray.size()/2; i++)
      {
      mdata[i*width+0] = static_cast<float>(vertices[indexArray[i*2+0]].x);
      mdata[i*width+1] = static_cast<float>(vertices[indexArray[i*2+0]].y);
      mdata[i*width+2] = static_cast<float>(vertices[indexArray[i*2+0]].z);
      mdata[i*width+3] = static_cast<float>(vertices[indexArray[i*2+1]].x);
      mdata[i*width+4] = static_cast<float>(vertices[indexArray[i*2+1]].y);
      mdata[i*width+5] = static_cast<float>(vertices[indexArray[i*2+1]].z);
      int mat = 0;
      if (numCellMaterials)
        {
        mat = rIndexArray[i*2];
        }
      else if (numPointColors)
        {
        mat = indexArray[i*2];
        }
      idata[i*width+6] = mat;
      if (scaleArray != NULL)
        {
        double avg = (*scaleArray->GetTuple(indexArray[i*2+0]) +
                      *scaleArray->GetTuple(indexArray[i*2+1]))*0.5;
        mdata[i*width+7] = MapThroughPWF
          (avg,
           scaleFunction);
        }
      }
    OSPData _PointColors = NULL;
    if (numPointColors)
      {
      _PointColors = ospNewData(numPointColors, OSP_FLOAT4, &PointColors[0]);
      }
    OSPData _mdata = ospNewData(indexArray.size()/2*width, OSP_FLOAT, mdata);
    ospSetData(ospMesh, "cylinders", _mdata);
    ospSet1i(ospMesh, "bytes_per_cylinder", width*sizeof(float));
    ospSet1i(ospMesh, "offset_v0", 0);
    ospSet1i(ospMesh, "offset_v1", 3*sizeof(float));
    ospSet1f(ospMesh, "radius", lineWidth);
    ospSet1i(ospMesh, "offset_radius", scaleOffset);

    if (numCellMaterials)
      {
      //per cell color
      ospSet1i(ospMesh, "offset_materialID", 6*sizeof(float));
      ospSetData(ospMesh, "materialList", CellMaterials);
      }
    else if (numPointColors)
      {
      //per point color
      ospSet1i(ospMesh, "offset_colorID", 6*sizeof(float));
      ospSetData(ospMesh, "color", _PointColors);
      }
    else
      {
      //per actor color
      ospSetMaterial(ospMesh, actorMaterial);
      }

    ospAddGeometry(oModel, ospMesh);
    ospCommit(ospMesh);
    //ospCommit(oModel); //TODO: crashes without yet others don't, why?
    ospRelease(_PointColors);
    ospRelease(_mdata);
    delete[] mdata;

    return ospMesh;
  }
*/


void RenderAsTriangles(
        optix::Context ctx,
        optix::Buffer vert_buffer,
        optix::Buffer norm_buffer,
        std::vector<unsigned int> &indexArray,
        std::vector<unsigned int> &rIndexArray,
        vtkOptiXPolyDataMapperNode::Geom* my_geom
        )
{
    const size_t numTriangles = indexArray.size() / 3;
    optix::Buffer tri_buffer = ctx->createBuffer(
            RT_BUFFER_INPUT, 
            RT_FORMAT_INT3,
            numTriangles
            );
    my_geom->AddBuffer( tri_buffer );

    optix::int3* tri_buffer_data = reinterpret_cast<optix::int3*>( tri_buffer->map() );
    for( size_t i = 0; i < numTriangles; i++ )
    { 
        tri_buffer_data[i] = optix::make_int3( 
                static_cast<int>(indexArray[i*3+0]),
                static_cast<int>(indexArray[i*3+1]),
                static_cast<int>(indexArray[i*3+2])
                );
    }
    tri_buffer->unmap();
    
    optix::Geometry geometry = ctx->createGeometry();
    geometry->setPrimitiveCount( numTriangles );
    geometry->setIntersectionProgram( ctx->createProgramFromPTXFile(
                VTK_OPTIX_PTX_DIR + "/cuda_compile_ptx_generated_TriangleMesh.cu.ptx",
                "TriangleMeshIntersection"
                )
            );
    geometry->setBoundingBoxProgram( ctx->createProgramFromPTXFile(
                VTK_OPTIX_PTX_DIR + "/cuda_compile_ptx_generated_TriangleMesh.cu.ptx",
                "TriangleMeshBoundingBox"
                )
            );
    geometry[ "vertices"  ]->setBuffer( vert_buffer );
    geometry[ "normals"   ]->setBuffer( norm_buffer );
    geometry[ "triangles" ]->setBuffer( tri_buffer );

    optix::Material matl = ctx->createMaterial();
    matl->setClosestHitProgram( 0, ctx->createProgramFromPTXFile(
                VTK_OPTIX_PTX_DIR + "/cuda_compile_ptx_generated_Lambertian.cu.ptx",
                "LambertianClosestHit"
                )
            );
    matl->setAnyHitProgram( 1, ctx->createProgramFromPTXFile(
                VTK_OPTIX_PTX_DIR + "/cuda_compile_ptx_generated_Lambertian.cu.ptx",
                "LambertianAnyHit"
                )
            );
    my_geom->Add( ctx->createGeometryInstance( geometry, &matl, &matl+1 ) );
}

/*
  //----------------------------------------------------------------------------
  OSPGeometry RenderAsTriangles(OSPData vertices,
                                std::vector<unsigned int> &indexArray,
                                std::vector<unsigned int> &rIndexArray,
                                OSPMaterial actorMaterial,
                                int numNormals,
                                osp::vec3f *normals,
                                vtkImageData *vColorTextureMap,
                                int numTextureCoordinates,
                                float *textureCoordinates,
                                int numCellMaterials,
                                OSPData CellMaterials,
                                int numPointColors,
                                osp::vec4f *PointColors,
                                int numPointValueTextureCoords,
                                float *pointValueTextureCoords,
                                OSPModel oModel,
                                OSPRenderer oRenderer
                                )
{
    OSPGeometry ospMesh = ospNewGeometry("trianglemesh");
    ospSetData(ospMesh, "position", vertices);

    size_t numTriangles = indexArray.size() / 3;
    osp::vec3i *triangles = new osp::vec3i[numTriangles];
    for (size_t i = 0, mi = 0; i < numTriangles; i++, mi += 3)
    {
        triangles[i] = osp::vec3i{static_cast<int>(indexArray[mi + 0]),
            static_cast<int>(indexArray[mi + 1]),
            static_cast<int>(indexArray[mi + 2])};
    }
    OSPData index = ospNewData(numTriangles, OSP_INT3, &triangles[0]);
    delete[] triangles;
    ospSetData(ospMesh, "index", index);

    OSPData _normals = NULL;
    if (numNormals)
    {
        _normals = ospNewData(numNormals, OSP_FLOAT3, &normals[0]);
        ospSetData(ospMesh, "vertex.normal", _normals);
    }

    //send the texture map and texture coordiantes over
    bool _hastm = false;
    osp::vec2f *tc = NULL;
    if (numTextureCoordinates || numPointValueTextureCoords)
    {
        _hastm = true;

        if (numPointValueTextureCoords)
        {
            //using 1D texture for point value LUT
            tc = new osp::vec2f[numPointValueTextureCoords];
            for (size_t i = 0; i < numPointValueTextureCoords; i++)
            {
                tc[i] = osp::vec2f{pointValueTextureCoords[i],0};
            }
            OSPData tcs = ospNewData(numPointValueTextureCoords, OSP_FLOAT2, &tc[0]);
            ospSetData(ospMesh, "vertex.texcoord", tcs);
        }
        else if (numTextureCoordinates)
        {
            //2d texture mapping
            tc = new osp::vec2f[numTextureCoordinates/2];
            float *itc = textureCoordinates;
            for (size_t i = 0; i < numTextureCoordinates; i+=2)
            {
                float t1,t2;
                t1 = *itc;
                itc++;
                t2 = *itc;
                itc++;
                tc[i/2] = osp::vec2f{t1,t2};
            }
            OSPData tcs = ospNewData(numTextureCoordinates/2, OSP_FLOAT2, &tc[0]);
            ospSetData(ospMesh, "vertex.texcoord", tcs);
        }
        OSPMaterial ospMaterial = vtkosp::VTKToOSPTexture(vColorTextureMap,oRenderer);
        ospSetMaterial(ospMesh, ospMaterial);
    }
    delete[] tc;

    //send over cell colors, point colors or whole actor color
    OSPData _cmats = NULL;
    OSPData _PointColors = NULL;
    int *ids = NULL;
    if (!_hastm)
    {
        if (numCellMaterials)
        {
            ids = new int[numTriangles];
            for (size_t i = 0; i < numTriangles; i++)
            {
                ids[i] = rIndexArray[i*3+0];
            }
            _cmats = ospNewData(numTriangles, OSP_INT, &ids[0]);
            ospSetData(ospMesh, "prim.materialID", _cmats);
            ospSetData(ospMesh, "materialList", CellMaterials);
        }
        else if (numPointColors)
        {
            _PointColors = ospNewData(numPointColors, OSP_FLOAT4, &PointColors[0]);
            ospSetData(ospMesh, "vertex.color", _PointColors);
        }
        else
        {
            ospSetMaterial(ospMesh, actorMaterial);
        }
    }

    ospAddGeometry(oModel, ospMesh);
    ospCommit(ospMesh);
    ospCommit(oModel); //TODO: crashes without yet others don't, why?
    ospRelease(index);
    ospRelease(_normals);
    ospRelease(_PointColors);
    ospRelease(_cmats);
    delete[] ids;
    return ospMesh;
}
*/

//============================================================================
vtkStandardNewMacro(vtkOptiXPolyDataMapperNode);

//----------------------------------------------------------------------------
vtkOptiXPolyDataMapperNode::vtkOptiXPolyDataMapperNode()
{
    this->MyGeom = 0;
}

//----------------------------------------------------------------------------
vtkOptiXPolyDataMapperNode::~vtkOptiXPolyDataMapperNode()
{
    delete this->MyGeom;
}

//----------------------------------------------------------------------------
void vtkOptiXPolyDataMapperNode::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkOptiXPolyDataMapperNode::RenderPoly(
        vtkOptiXRendererNode* orn,
        vtkOptiXActorNode*    aNode,
        vtkPolyData*          poly
        )
{
    optix::Context ctx = orn->GetOptiXContext();
    vtkActor* act = vtkActor::SafeDownCast(aNode->GetRenderable());
    vtkProperty* property = act->GetProperty();
    
    std::vector<double> vertices;
    vtkPolyDataMapperNode::TransformPoints( act, poly, vertices );
    const size_t numPositions = vertices.size() / 3;
    
    optix::Buffer vert_buffer = ctx->createBuffer(
            RT_BUFFER_INPUT, 
            RT_FORMAT_FLOAT3,
            numPositions );
    this->MyGeom->AddBuffer( vert_buffer );

    optix::float3* vert_buffer_data = reinterpret_cast<optix::float3*>( vert_buffer->map() );
    for( size_t i = 0; i < numPositions; i++ )
    { 
        vert_buffer_data[i] = optix::make_float3( 
                static_cast<float>(vertices[i*3+0]),
                static_cast<float>(vertices[i*3+1]),
                static_cast<float>(vertices[i*3+2]) 
                );
    }
    vert_buffer->unmap();
    vertices.clear();
    
    //make connectivity
    std::vector<unsigned int> vertex_index;
    std::vector<unsigned int> vertex_reverse;
    std::vector<unsigned int> line_index;
    std::vector<unsigned int> line_reverse;
    std::vector<unsigned int> triangle_index;
    std::vector<unsigned int> triangle_reverse;
    std::vector<unsigned int> strip_index;
    std::vector<unsigned int> strip_reverse;
    vtkPolyDataMapperNode::MakeConnectivity(poly, property->GetRepresentation(),
            vertex_index, vertex_reverse,
            line_index, line_reverse,
            triangle_index, triangle_reverse,
            strip_index, strip_reverse);

    //create an ospray mesh for the polygon cells
    if (triangle_index.size())
    {
        //format depends on representation style
        switch( property->GetRepresentation())
        {
            case VTK_POINTS:
            {
                std::cerr << "\tUnhandled polygon representation: POINTS" << std::endl;
                break;
            }
            case VTK_WIREFRAME:
            {
                std::cerr << "\tUnhandled polygon representation: WIREFRAME" << std::endl;
                break;
            }
            default:
            {
                std::cerr << "\tCreating poly mesh" << std::endl;

                vtkDataArray* vNormals = 0; 
                if( property->GetInterpolation() != VTK_FLAT )
                {
                    vNormals = poly->GetPointData()->GetNormals();
                }
                optix::Buffer norm_buffer = vtkDataArrayToBuffer( vNormals, ctx ); 
                RenderAsTriangles(
                        ctx,
                        vert_buffer,
                        norm_buffer,
                        triangle_index,
                        triangle_reverse,
                        this->MyGeom
                        );
                /*
                osp::vec3f *normals = NULL;
                int numNormals = 0;
                if (property->GetInterpolation() != VTK_FLAT)
                {
                    vtkDataArray *vNormals = poly->GetPointData()->GetNormals();
                    if (vNormals)
                    {
                        vtkosp::VToOPointNormals
                            (vNormals, normals);
                        numNormals = vNormals->GetNumberOfTuples();
                    }
                }
                myMeshes
                    ->Add(vtkosp::RenderAsTriangles(position,
                                triangle_index, triangle_reverse,
                                oMaterial,
                                numNormals, normals,
                                vColorTextureMap,
                                numTextureCoordinates, textureCoordinates,
                                numCellMaterials, cellMaterials,
                                numPointColors, pointColors,
                                numPointValueTextureCoords, pointValueTextureCoords,
                                oModel, oRenderer
                                ));
                delete[] normals;
                */
            }
        }
    }
    else
    {
        std::cerr << "\tNo polygons present" << std::endl;
    }

    this->MyGeom->AddMyselfTo( orn->GetOptiXGeometryGroup() );


    /*
    OSPRenderer oRenderer = static_cast<OSPRenderer>(renderer);
    OSPModel oModel = static_cast<OSPModel>(model);
    vtkActor *act = vtkActor::SafeDownCast(aNode->GetRenderable());
    vtkProperty *property = act->GetProperty();
    //  vtkosp::MyGeom *myMeshes = (vtkosp::MyGeom*)this->OSPMeshes;

    //make geometry
    std::vector<double> _vertices;
    vtkPolyDataMapperNode::TransformPoints(act, poly, _vertices);
    size_t numPositions = _vertices.size()/3;
    osp::vec3fa *vertices = new osp::vec3fa[numPositions];
    for (size_t i = 0; i < numPositions; i++)
    {
        vertices[i] =
            osp::vec3fa{static_cast<float>(_vertices[i*3+0]),
                static_cast<float>(_vertices[i*3+1]),
                static_cast<float>(_vertices[i*3+2])};
    }
    OSPData position = ospNewData(numPositions, OSP_FLOAT3A, &vertices[0]);
    ospCommit(position);
    _vertices.clear();

    //make connectivity
    std::vector<unsigned int> vertex_index;
    std::vector<unsigned int> vertex_reverse;
    std::vector<unsigned int> line_index;
    std::vector<unsigned int> line_reverse;
    std::vector<unsigned int> triangle_index;
    std::vector<unsigned int> triangle_reverse;
    std::vector<unsigned int> strip_index;
    std::vector<unsigned int> strip_reverse;
    vtkPolyDataMapperNode::MakeConnectivity(poly, property->GetRepresentation(),
            vertex_index, vertex_reverse,
            line_index, line_reverse,
            triangle_index, triangle_reverse,
            strip_index, strip_reverse);

    //choose sphere and cylinder radii (for points and lines) that approximate pointsize and linewidth
    //TODO: take into account camera so that size is roughly numpixels wide on screen?
    vtkMapper *mapper = act->GetMapper();
    double length = 1.0;
    if (mapper)
    {
        length = mapper->GetLength();
    }
    double pointSize = length/1000.0 * property->GetPointSize();
    double lineWidth = length/1000.0 * property->GetLineWidth();
    //finer control over sphere and cylinders sizes
    int enable_scaling =
        vtkOptiXActorNode::GetEnableScaling(act);
    vtkDataArray *scaleArray = NULL;
    vtkPiecewiseFunction *scaleFunction = NULL;
    if (enable_scaling && mapper)
    {
        vtkInformation *mapInfo = mapper->GetInformation();
        char *scaleArrayName = (char *)mapInfo->Get
            (vtkOptiXActorNode::SCALE_ARRAY_NAME());
        scaleArray = poly->GetPointData()->GetArray(scaleArrayName);
        scaleFunction = vtkPiecewiseFunction::SafeDownCast
            (mapInfo->Get(vtkOptiXActorNode::SCALE_FUNCTION()));
    }

    //--------------------------------------------------------------------------
    //
    //per actor material
    //
    //--------------------------------------------------------------------------

    OSPMaterial oMaterial;
    oMaterial = ospNewMaterial(oRenderer,"RayTraceMaterial");
    float ambientf[] =
    {
        static_cast<float>(ambientColor[0]*property->GetAmbient()),
        static_cast<float>(ambientColor[1]*property->GetAmbient()),
        static_cast<float>(ambientColor[2]*property->GetAmbient())
    };
    float diffusef[] =
    {
        static_cast<float>(diffuseColor[0]*property->GetDiffuse()),
        static_cast<float>(diffuseColor[1]*property->GetDiffuse()),
        static_cast<float>(diffuseColor[2]*property->GetDiffuse())
    };
    float specPower =
        static_cast<float>(property->GetSpecularPower());
    float specAdjust = 2.0f/(2.0f+specPower); //since OSP 0.10.0
    float specularf[] =
    {
        static_cast<float>(property->GetSpecularColor()[0]*
                property->GetSpecular()*
                specAdjust),
        static_cast<float>(property->GetSpecularColor()[1]*
                property->GetSpecular()*
                specAdjust),
        static_cast<float>(property->GetSpecularColor()[2]*
                property->GetSpecular()*
                specAdjust)
    };

    ospSet3fv(oMaterial,"Ka",ambientf);
    double scale = property->GetDiffuse();
    if (property->GetDiffuse()==0.0)
    {
        //a workaround for ParaView, remove when ospray supports Ka
        ospSet3fv(oMaterial,"Kd",ambientf);
    }
    else
    {
        ospSet3fv(oMaterial,"Kd",diffusef);
    }
    ospSet3fv(oMaterial,"Ks",specularf);
    ospSet1f(oMaterial,"Ns",specPower);
    ospSet1f(oMaterial,"d",float(opacity));
    ospCommit(oMaterial);

    //texture
    vtkTexture *texture = act->GetTexture();
    int numTextureCoordinates = 0;
    float *textureCoordinates = NULL;
    vtkImageData* vColorTextureMap = NULL;
    if (texture)
    {
        vtkDataArray *da = poly->GetPointData()->GetTCoords();
        numTextureCoordinates = da->GetNumberOfTuples();
        textureCoordinates = new float[numTextureCoordinates*2];
        float *tp = textureCoordinates;
        for (int i=0; i<numTextureCoordinates; i++)
        {
            *tp = static_cast<float>(da->GetTuple(i)[0]);
            tp++;
            *tp = static_cast<float>(da->GetTuple(i)[1]);
            tp++;
        }
        vColorTextureMap = vtkImageData::SafeDownCast
            (texture->GetInput());
        numTextureCoordinates = numTextureCoordinates*2;
    }

    //colors from point and cell arrays
    int numCellMaterials = 0;
    OSPData cellMaterials = NULL;
    int numPointColors = 0;
    osp::vec4f *pointColors = NULL;
    int numPointValueTextureCoords = 0;
    float *pointValueTextureCoords = NULL;
    //
    //now ask mapper to do most of the work and provide use with
    //colors per cell and colors or texture coordinates per point
    vtkUnsignedCharArray *vColors = NULL;
    vtkFloatArray *vColorCoordinates = NULL;
    vtkImageData* pColorTextureMap = NULL;
    int cellFlag = -1; //mapper tells us which
    if (mapper)
    {
        mapper->MapScalars(poly, 1.0, cellFlag);
        vColors = mapper->GetColorMapColors();
        vColorCoordinates = mapper->GetColorCoordinates();
        pColorTextureMap = mapper->GetColorTextureMap();
    }
    if (vColors)
    {
        if (cellFlag==2 && mapper->GetFieldDataTupleId() > -1)
        {
            int numComp = vColors->GetNumberOfComponents();
            unsigned char *colorPtr = vColors->GetPointer(0);
            colorPtr = colorPtr + mapper->GetFieldDataTupleId()*numComp;
            // this setting (and all the other scalar colors)
            // really depends on mapper->ScalarMaterialMode
            // but I'm not sure Ka is working currently so leaving it on Kd
            float fdiffusef[] =
            {
                static_cast<float>(colorPtr[0]*property->GetDiffuse()/255.0f),
                static_cast<float>(colorPtr[1]*property->GetDiffuse()/255.0f),
                static_cast<float>(colorPtr[2]*property->GetDiffuse()/255.0f)
            };
            ospSet3fv(oMaterial,"Kd",fdiffusef);
            ospCommit(oMaterial);
        }
        else if (cellFlag==1)
        {
            //color on cell
            std::vector<OSPMaterial> cellColors;
            vtkosp::CellColorMaterials(vColors, oRenderer, cellColors, specularf,
                    float(property->GetSpecularPower()), opacity);
            numCellMaterials = cellColors.size();
            cellMaterials = ospNewData(cellColors.size(), OSP_OBJECT, &cellColors[0]);
            ospCommit(cellMaterials);
            cellColors.clear();
        }
        else if (cellFlag==0)
        {
            //color on point interpolated RGB
            numPointColors = vColors->GetNumberOfTuples();
            pointColors = new osp::vec4f[numPointColors];
            for (int i = 0; i < numPointColors; i++)
            {
                unsigned char *color = vColors->GetPointer(4 * i);
                pointColors[i] = osp::vec4f{color[0] / 255.0f,
                    color[1] / 255.0f,
                    color[2] / 255.0f,
                    1};
            }
        }
    }
    else
    {
        if (vColorCoordinates && pColorTextureMap)
        {
            //color on point interpolated values (subsequently colormapped via 1D LUT)
            numPointValueTextureCoords = vColorCoordinates->GetNumberOfTuples();
            pointValueTextureCoords = new float[numPointValueTextureCoords];
            float *tc = vColorCoordinates->GetPointer(0);
            for (int i = 0; i < numPointValueTextureCoords; i++)
            {
                pointValueTextureCoords[i] = *tc;
                tc+=2;
            }
            vColorTextureMap = pColorTextureMap;
        }

    }

    //--------------------------------------------------------------------------
    //
    // Create geometry 
    //
    //--------------------------------------------------------------------------

    //create an ospray mesh for the vertex cells
    if (vertex_index.size())
    {
        myMeshes
            ->Add(vtkosp::RenderAsSpheres(vertices,
                        vertex_index, vertex_reverse,
                        pointSize, scaleArray, scaleFunction,
                        oMaterial,
                        numCellMaterials, cellMaterials,
                        numPointColors, pointColors,
                        oModel
                        ));
    }


    //create an ospray mesh for the line cells
    if (line_index.size())
    {
        //format depends on representation style
        if (property->GetRepresentation() == VTK_POINTS)
        {
            myMeshes
                ->Add(vtkosp::RenderAsSpheres(vertices,
                            line_index, line_reverse,
                            pointSize, scaleArray, scaleFunction,
                            oMaterial,
                            numCellMaterials, cellMaterials,
                            numPointColors, pointColors,
                            oModel
                            ));
        }
        else
        {
            myMeshes
                ->Add(vtkosp::RenderAsCylinders(vertices,
                            line_index, line_reverse,
                            lineWidth, scaleArray, scaleFunction,
                            oMaterial,
                            numCellMaterials, cellMaterials,
                            numPointColors, pointColors,
                            oModel
                            ));
        }
    }

    //create an ospray mesh for the polygon cells
    if (triangle_index.size())
    {
        //format depends on representation style
        switch (property->GetRepresentation())
        {
            case VTK_POINTS:
                {
                    myMeshes
                        ->Add(vtkosp::RenderAsSpheres(vertices,
                                    triangle_index, triangle_reverse,
                                    pointSize, scaleArray, scaleFunction,
                                    oMaterial,
                                    numCellMaterials, cellMaterials,
                                    numPointColors, pointColors,
                                    oModel
                                    ));
                    break;
                }
            case VTK_WIREFRAME:
                {
                    myMeshes
                        ->Add(vtkosp::RenderAsCylinders(vertices,
                                    triangle_index, triangle_reverse,
                                    lineWidth, scaleArray, scaleFunction,
                                    oMaterial,
                                    numCellMaterials, cellMaterials,
                                    numPointColors, pointColors,
                                    oModel
                                    ));
                    break;
                }
            default:
                {
                    osp::vec3f *normals = NULL;
                    int numNormals = 0;
                    if (property->GetInterpolation() != VTK_FLAT)
                    {
                        vtkDataArray *vNormals = poly->GetPointData()->GetNormals();
                        if (vNormals)
                        {
                            vtkosp::VToOPointNormals
                                (vNormals, normals);
                            numNormals = vNormals->GetNumberOfTuples();
                        }
                    }
                    myMeshes
                        ->Add(vtkosp::RenderAsTriangles(position,
                                    triangle_index, triangle_reverse,
                                    oMaterial,
                                    numNormals, normals,
                                    vColorTextureMap,
                                    numTextureCoordinates, textureCoordinates,
                                    numCellMaterials, cellMaterials,
                                    numPointColors, pointColors,
                                    numPointValueTextureCoords, pointValueTextureCoords,
                                    oModel, oRenderer
                                    ));
                    delete[] normals;
                }
        }
    }

    if (strip_index.size())
    {
        switch (property->GetRepresentation())
        {
            case VTK_POINTS:
                {
                    myMeshes
                        ->Add(vtkosp::RenderAsSpheres(vertices,
                                    strip_index, strip_reverse,
                                    pointSize, scaleArray, scaleFunction,
                                    oMaterial,
                                    numCellMaterials, cellMaterials,
                                    numPointColors, pointColors,
                                    oModel
                                    ));
                    break;
                }
            case VTK_WIREFRAME:
                {
                    myMeshes
                        ->Add(vtkosp::RenderAsCylinders(vertices,
                                    strip_index, strip_reverse,
                                    lineWidth, scaleArray, scaleFunction,
                                    oMaterial,
                                    numCellMaterials, cellMaterials,
                                    numPointColors, pointColors,
                                    oModel
                                    ));
                    break;
                }
            default:
                {
                    osp::vec3f *normals = NULL;
                    int numNormals = 0;
                    if (property->GetInterpolation() != VTK_FLAT)
                    {
                        vtkDataArray *vNormals = poly->GetPointData()->GetNormals();
                        if (vNormals)
                        {
                            vtkosp::VToOPointNormals
                                (vNormals, normals);
                            numNormals = strip_index.size();
                        }
                    }
                    myMeshes
                        ->Add(vtkosp::RenderAsTriangles(position,
                                    strip_index, strip_reverse,
                                    oMaterial,
                                    numNormals, normals,
                                    vColorTextureMap,
                                    numTextureCoordinates, textureCoordinates,
                                    numCellMaterials, cellMaterials,
                                    numPointColors, pointColors,
                                    numPointValueTextureCoords, pointValueTextureCoords,
                                    oModel, oRenderer
                                    ));
                    delete[] normals;
                }
        }
    }
    ospRelease(position);
    delete[] vertices;
    delete[] pointColors;
    delete[] pointValueTextureCoords;
*/
}


void vtkOptiXPolyDataMapperNode::CreateNewMeshes()
{
    delete this->MyGeom;
    this->MyGeom = new vtkOptiXPolyDataMapperNode::Geom;
}

//----------------------------------------------------------------------------
void vtkOptiXPolyDataMapperNode::Render(bool prepass)
{
    std::cerr << "==================================================================" << std::endl;
    std::cerr << "Poly Data Mapper node render -- prepass: " << prepass << std::endl;
    if (prepass)
    {
        // we use a lot of params from our parent
        vtkOptiXActorNode *aNode = vtkOptiXActorNode::SafeDownCast(this->Parent);
        vtkActor *act = vtkActor::SafeDownCast(aNode->GetRenderable());

        //
        // TODO: handle visibility
        //

        vtkOptiXRendererNode *orn =
            static_cast<vtkOptiXRendererNode *>(
                    this->GetFirstAncestorOfType("vtkOptiXRendererNode"));
        
        //if there are no changes, just reuse last result
        bool enable_cache = true; //turn off to force rebuilds for debugging
        vtkMTimeType inTime = aNode->GetMTime();
        if (enable_cache && this->RenderTime >= inTime)
        {
            std::cerr << "\tReusing polydata" << std::endl;
            return;
        }
        this->RenderTime = inTime;
        
        //something changed so make new meshes
        this->CreateNewMeshes();

        vtkPolyData* poly = NULL;
        vtkMapper* mapper = act->GetMapper();
        if (mapper)
        {
            poly = (vtkPolyData*)(mapper->GetInput());
        }
        if (poly)
        {
            this->RenderPoly( orn, aNode, poly );
        }
    }
    /*
    if (prepass)
    {
        // we use a lot of params from our parent
        vtkOptiXActorNode *aNode = vtkOptiXActorNode::SafeDownCast(this->Parent);
        vtkActor *act = vtkActor::SafeDownCast(aNode->GetRenderable());

        if (act->GetVisibility() == false)
        {
            delete (vtkosp::MyGeom*)this->OSPMeshes;
            this->OSPMeshes = NULL;
            return;
        }

        vtkOptiXRendererNode *orn =
            static_cast<vtkOptiXRendererNode *>(
                    this->GetFirstAncestorOfType("vtkOptiXRendererNode"));

        //if there are no changes, just reuse last result
        bool enable_cache = true; //turn off to force rebuilds for debugging
        vtkMTimeType inTime = aNode->GetMTime();
        if (enable_cache && this->RenderTime >= inTime)
        {
            OSPModel oModel = static_cast<OSPModel>(orn->GetOModel());
            this->AddMeshesToModel(oModel);
            return;
        }
        this->RenderTime = inTime;

        //something changed so make new meshes
        this->CreateNewMeshes();

        vtkPolyData *poly = NULL;
        vtkMapper *mapper = act->GetMapper();
        if (mapper)
        {
            poly = (vtkPolyData*)(mapper->GetInput());
        }
        if (poly)
        {
            vtkProperty * property = act->GetProperty();
            this->ORenderPoly(
                    orn->GetORenderer(),
                    orn->GetOModel(),
                    aNode,
                    poly,
                    property->GetAmbientColor(),
                    property->GetDiffuseColor(),
                    property->GetOpacity()
                    );
        }
    }
    */
}
