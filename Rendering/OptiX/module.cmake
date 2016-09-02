
vtk_module(vtkRenderingOptiX
  DEPENDS
    vtkRenderingSceneGraph
    #todo promote compositedatadisplayattributes to rendering/core
    vtkRendering${VTK_RENDERING_BACKEND} #only for comp.data.disp.attr.
  TEST_DEPENDS
    vtkFiltersTexture
    vtkInteractionStyle
    vtkIOGeometry
    vtkIOPLY
    vtkIOXML
    vtkRenderingAnnotation
    vtkTestingCore
    vtkTestingRendering
  EXCLUDE_FROM_ALL
  )
