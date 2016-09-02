find_package( OptiX REQUIRED )
find_package( CUDA 4.0 REQUIRED )

include_directories(
    "${OptiX_INCLUDE}"
    "${CUDA_INCLUDE_DIRS}"
    )


set( PTX_OUTPUT_DIR "${CMAKE_BINARY_DIR}/lib/ptx" CACHE PATH
    "Path to where generated PTX is stored"
    )

set( CUDA_GENERATED_OUTPUT_DIR ${PTX_OUTPUT_DIR} )

# Set up build .cu -> .ptx build
set( CUDA_NVCC_FLAGS "-arch sm_30 --use_fast_math" )
CUDA_COMPILE_PTX( ptx_files
    CUDA/Common.h
    CUDA/draw_color.cu
    CUDA/Lambertian.cu
    CUDA/PerspectiveCamera.cu
    CUDA/TriangleMesh.cu
    )

configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/vtkOptiXConfig.h.in"
    "${CMAKE_BINARY_DIR}/include/vtkOptiXConfig.h"
    )
include_directories(
    "${CMAKE_BINARY_DIR}/include"
    )


set(Module_SRCS
  vtkOptiXActorNode.cxx
  vtkOptiXCameraNode.cxx
  vtkOptiXCompositePolyDataMapper2Node.cxx
  vtkOptiXLightNode.cxx
  vtkOptiXPass.cxx
  vtkOptiXPolyDataMapperNode.cxx
  vtkOptiXRendererNode.cxx
  vtkOptiXViewNodeFactory.cxx
  vtkOptiXWindowNode.cxx
  )

vtk_module_library(vtkRenderingOptiX 
    ${Module_SRCS}
    ${ptx_files}
    )

target_link_libraries(${vtk-module} LINK_PUBLIC ${optix_LIBRARY})

# OptiX_Core uses MMTime which is in it's own special library.
if(WIN32)
  target_link_libraries(${vtk-module} LINK_PRIVATE Winmm)
endif()
