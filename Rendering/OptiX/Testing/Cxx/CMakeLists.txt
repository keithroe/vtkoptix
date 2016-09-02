
find_package( OptiX REQUIRED )
find_package( CUDA 4.0 REQUIRED )                                                

include_directories(
    "${OptiX_INCLUDE}"
    "${CUDA_INCLUDE_DIRS}"                                                       
    )

vtk_add_test_cxx(${vtk-module}CxxTests tests
  TestOptiXCompositePolyDataMapper2.cxx,
  TestOptiXDynamicObject.cxx,
  TestOptiXDynamicScene.cxx,
  TestOptiXImplicits.cxx,
  TestOptiXLayers.cxx,
  TestOptiXLights.cxx,
  TestOptiXMultiBlock.cxx
  TestOptiXMultiBlockPartialArrayFieldData.cxx,
  TestOptiXPass.cxx,
  TestOptiXRenderMesh.cxx,
  TestOptiXScalarBar.cxx
  TestOptiXTime.cxx,
  TestOptiXWindow.cxx
  #TODO: interaction, pure offscreen, volumes
  )
vtk_test_cxx_executable(${vtk-module}CxxTests tests RENDERING_FACTORY
  vtkOptiXTestInteractor.cxx
  )


