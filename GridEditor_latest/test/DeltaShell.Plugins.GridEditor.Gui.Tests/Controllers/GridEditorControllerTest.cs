using System.Collections.Generic;
using DelftTools.Utils.Collections.Extensions;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Native;
using DeltaShell.Plugins.GridEditor.Gui.Controllers;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using GeoAPI.Extensions.Coverages;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Coverages;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Extensions.Grids;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using Is = Rhino.Mocks.Constraints.Is;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.Controllers
{
    [TestFixture]
    public class GridEditorControllerTest
    {
        [Test]
        public void ResetUnstructuredGridsRefreshesRibbonUnstructuredGrids()
        {
            var unstructuredGrids = new List<UnstructuredGrid>
            {
                new UnstructuredGrid(),
                new UnstructuredGrid()
            };

            var mocks = new MockRepository();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var interactor = mocks.StrictMock<IMapInteractor>();

            interactor.Expect(i => i.GetUnstructuredGrids()).Return(unstructuredGrids);

            ribbonState.Expect(r => r.UnstructuredGrids = null).IgnoreArguments().WhenCalled((m) =>
            {
                Assert.AreEqual(1,m.Arguments.Length);
                var grids = (IList<UnstructuredGrid>) m.Arguments[0];

                Assert.AreEqual(unstructuredGrids[0], grids[0]);
                Assert.AreEqual(unstructuredGrids[1], grids[1]);

            });

            mocks.ReplayAll();

            var controller = new GridEditorController
            {
                RibbonState = ribbonState,
                MapInteractor = interactor
            };

            controller.ResetUnstructuredGrids();

            mocks.VerifyAll();
        }

        [Test]
        public void SetEditingCreatesADisposableMeshGeometry()
        {
            var grid = new UnstructuredGrid();
            
            var mocks = new MockRepository();

            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var interactor = mocks.StrictMock<IMapInteractor>();

            interactor.Expect(i => i.StartInteraction(null)).IgnoreArguments();
            interactor.Expect(i => i.StopInteraction()).IgnoreArguments();
            interactor.Expect(i => i.Renderer).Return(null);
            interactor.Expect(i => i.GetSplineGeometry = null).IgnoreArguments();
            interactor.Expect(i => i.InsertVertex = null).IgnoreArguments();
            interactor.Expect(i => i.InsertEdge = null).IgnoreArguments();
            interactor.Expect(i => i.GetVertexIndex = null).IgnoreArguments();
            interactor.Expect(i => i.DeleteVertex = null).IgnoreArguments();
            interactor.Expect(i => i.MergeVertices = null).IgnoreArguments();
            interactor.Expect(i => i.DeleteEdges = null).IgnoreArguments();
            interactor.Expect(i => i.MoveVertex = null).IgnoreArguments();

            ribbonState.Expect(r => r.SelectedGrid).Return(grid).Repeat.Any();
            ribbonState.Expect(r => r.RefreshState()).Repeat.Twice();

            mocks.ReplayAll();

            var controller = new GridEditorController
            {
                RibbonState = ribbonState,
                MapInteractor = interactor
            };


            Assert.IsNull(controller.GridEditorState.MeshGeometry, "No mesh by default");
            Assert.IsFalse(controller.IsEditing, "Default is not editing");

            // set editing on
            controller.IsEditing = true;

            Assert.NotNull(controller.GridEditorState.MeshGeometry, "Editing creates a mesh");
            
            // set editing off
            controller.IsEditing = false;

            Assert.IsNull(controller.GridEditorState.MeshGeometry, "Mesh is removed after leaving edit mode");

            mocks.VerifyAll();
        }

        [Test]
        public void SetNodesSetsGridDrawingMapToolShowPoints()
        {
            var mocks = new MockRepository();

            var interactor = mocks.StrictMock<IMapInteractor>();
            var renderer = mocks.StrictMock<IGridEditorStateRenderer>();
            var meshRenderer = mocks.StrictMock<IDisposableMeshGeometryRenderer>();

            interactor.Expect(i => i.Renderer).Return(renderer).Repeat.Any();
            renderer.Expect(r => r.MeshGeometryRenderer).Return(meshRenderer).Repeat.Any();

            meshRenderer.Expect(r => r.ShowGridVertices = true);
            meshRenderer.Expect(r => r.ShowGridVertices = false);

            mocks.ReplayAll();

            var controller = new GridEditorController
            {
                MapInteractor = interactor
            };

            controller.ShowVertices = true;

            // set editing off
            controller.ShowVertices = false;
            
            mocks.VerifyAll();
        }

        [Test]
        public void SetVertexSizeSetsGridDrawingMapToolPointSize()
        {
            var mocks = new MockRepository();
            var interactor = mocks.StrictMock<IMapInteractor>();
            var renderer = mocks.StrictMock<IGridEditorStateRenderer>();
            var meshRenderer = mocks.StrictMock<IDisposableMeshGeometryRenderer>();

            interactor.Expect(i => i.Renderer).Return(renderer).Repeat.Any();
            renderer.Expect(r => r.MeshGeometryRenderer).Return(meshRenderer).Repeat.Any();

            meshRenderer.Expect(r => r.PointSize = 10);
            
            mocks.ReplayAll();

            var controller = new GridEditorController
            {
                MapInteractor = interactor
            };

            controller.VertexSize = 10;

            mocks.VerifyAll();
        }

        [Test]
        public void DeleteSelectedVertices_ShouldNotTriggerApiDeleteVertex_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            
            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = ()=> api,
                GridEditorState = 
                {
                    SelectedVertices = new []{1,2}
                },
            };

            // Act
            gridEditorController.DeleteSelectedVertices();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void DeleteSelectedVertices_ShouldDisposeMeshGeometry_WhenNoVerticesAreSelected()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();

            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);
            api.Expect(a => a.CreateGridState()).Return(1);
            api.Expect(a => a.GetGridState(1)).Return(null);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Repeat.Twice().Return(true);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api,
                RibbonState = ribbonState,
                GridEditorState =
                {
                    SelectedVertices = null
                }
            };

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.DeleteSelectedVertices();

            // Assert
            Assert.IsNull(gridEditorController.GridEditorState.SelectedVertices, "Selected vertices should be removed after deleting");
            mocks.VerifyAll();
        }

        [Test]
        public void DeleteSelectedVertices_ShouldTriggerApiDeleteVertex_ForSelectedVertices()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();

            api.Expect(a => a.CreateGridState()).Return(1);
            api.Expect(a => a.GetGridState(1)).Return(null);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.DeleteVertex(1,0)).Return(true);
            api.Expect(a => a.DeleteVertex(1, 1)).Return(true);

            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api,
                GridEditorState = {SelectedVertices = new []{ 0, 1 } }
            };
            
            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.DeleteSelectedVertices();

            // Assert
            Assert.IsNull(gridEditorController.GridEditorState.SelectedVertices, "Selected vertices should be removed after deleting");

            mocks.VerifyAll();
        }

        [Test]
        public void HasSelectedVertices_ShouldWorkOnGridEditorControllerData()
        {
            //Arrange
            var gridEditorController = new GridEditorController
            {
                GridEditorState = 
                {
                    SelectedVertices = new int[] { 0, 1 }
                },
            };

            // Assert
            Assert.IsTrue(gridEditorController.HasSelectedVertices);

            // Act
            gridEditorController.GridEditorState.SelectedVertices = null;

            // Assert
            Assert.IsFalse(gridEditorController.HasSelectedVertices);
        }

        [Test]
        public void FlipEdges_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };
            
            // Act
            gridEditorController.FlipEdges();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void FlipEdges_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.FlipEdges(gridGeomId, false, false, 
                        ProjectToLandBoundaryOptions.NetBoundaryToLandBoundary)).Return(true);

            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);

            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api
            };

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.FlipEdges();

            // Assert
            Assert.AreEqual(disposableMeshGeometry, gridEditorController.GridEditorState.MeshGeometry);
            mocks.VerifyAll();
        }

        [Test]
        public void InsertEdges_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            var coordinate = new Coordinate(10, 10);
            gridEditorController.InsertVertex(coordinate);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void InsertEdges_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            int newIndex = 0;
            var verticesCoordinates = new[] { new Coordinate(10, 10), new Coordinate(10, 10) };
            var verticesIndexes = new[] { 0, 1 };
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.InsertVertex(gridGeomId, verticesCoordinates[0].X, verticesCoordinates[0].Y, verticesCoordinates[0].Z, ref newIndex)).Repeat.Twice().Return(true).IgnoreArguments();
            api.Expect(a => a.InsertEdge(gridGeomId, verticesIndexes[0], verticesIndexes[1], ref newIndex)).Repeat.Once().Return(true).IgnoreArguments();

            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);

            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);


            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api
            };
            gridEditorController.GridEditorState.SelectedVertices = verticesIndexes;

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.InsertVertex(verticesCoordinates[0]);
            gridEditorController.InsertVertex(verticesCoordinates[1]);
            gridEditorController.InsertEdge(0, 1);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MergeVertices_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.MergeVertices();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MergeVertices_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);

            api.Expect(a => a.MergeVertices(gridGeomId, null))
                .IgnoreArguments()
                .WhenCalled(m =>
                {
                    Assert.AreEqual(gridGeomId, m.Arguments[0]);
                    Assert.IsInstanceOf<DisposableGeometryList>(m.Arguments[1]);
                })
                .Return(true);

            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.MergeVertices();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeGrid_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.MakeGrid();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeGrid_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();

            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();
            var disposableGeometryList = new DisposableGeometryList();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.MakeGrid(gridGeomId, gridEditorController.MakeGridParameters, ref disposableGeometryList)).IgnoreArguments().Return(true);

            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.MakeGrid();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeGridFromSplines_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.MakeGridFromSplines();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeGridFromSplines_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            api.Expect(a => a.MakeGridFromSplines(Arg<int>.Is.Equal(gridGeomId), 
                    ref Arg<DisposableGeometryList>.Ref(Is.Anything(), null).Dummy, 
                    ref Arg<CurvilinearParameters>.Ref(Is.Anything(), null).Dummy))
                .IgnoreArguments()
                .WhenCalled(m =>
                {
                    Assert.AreEqual(gridGeomId, m.Arguments[0]);
                })
                .Return(true);

            mocks.ReplayAll();

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            gridEditorController.GridEditorState.Splines.AddRange(new[]
            {
                new Spline(),
                new Spline(),
                new Spline(),
                new Spline()
            });

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;
            gridEditorController.IsEditing = true;

            // Act
            gridEditorController.MakeGridFromSplines();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void Orthogonalize_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            api.Expect(a => a.Orthogonalize(gridGeomId, 
                    false, 
                    false,
                    0,
                    null,
                    null,
                    null))
                .IgnoreArguments()
                .WhenCalled(m =>
                {
                    Assert.AreEqual(gridGeomId, m.Arguments[0]);
                    Assert.NotNull(m.Arguments[4]);
                    Assert.IsInstanceOf<DisposableGeometryList>(m.Arguments[5]);
                    Assert.IsInstanceOf<DisposableGeometryList>(m.Arguments[6]);
                    Assert.AreEqual(gridEditorController.GridEditorState.OrthogonalizationParameters, m.Arguments[4]);
                })
                .Return(true);

            mocks.ReplayAll();

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;
            gridEditorController.IsEditing = true;

            // Act
            gridEditorController.Orthogonalize();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void Orthogonalize_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.Orthogonalize();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeGridFromSplinesOrthogonal_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.MakeGridFromSplinesOrthogonal();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeGridFromSplinesOrthogonal_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            api.Expect(a => a.MakeOrthogonalGridFromSplines( Arg<int>.Is.Equal(gridGeomId),
                    ref Arg<DisposableGeometryList>.Ref( Is.Anything(), null).Dummy,
                    ref Arg<CurvilinearParameters>.Ref(Is.Anything(), null).Dummy,
                    ref Arg<SplinesToCurvilinearParameters>.Ref(Is.Anything(), null).Dummy))
                .IgnoreArguments()
                .WhenCalled(m =>
                {
                    Assert.AreEqual(gridGeomId, m.Arguments[0]);
                })
                .Return(true);

            mocks.ReplayAll();

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            gridEditorController.GridEditorState.Splines.AddRange(new[]
            {
                new Spline(),
                new Spline(),
                new Spline(),
                new Spline()
            });

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;
            gridEditorController.IsEditing = true;

            // Act
            gridEditorController.MakeGridFromSplinesOrthogonal();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeTriangularGridInPolygon_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.MakeTriangularGridInPolygon();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeTriangularGridInPolygon_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            api.Expect(a => a.MakeTriangularGridInPolygon(Arg<int>.Is.Equal(gridGeomId),
                    ref Arg<DisposableGeometryList>.Ref(Is.Anything(), null).Dummy))
                .IgnoreArguments()
                .WhenCalled(m =>
                {
                    Assert.AreEqual(gridGeomId, m.Arguments[0]);
                })
                .Return(true);

            mocks.ReplayAll();

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;
            gridEditorController.IsEditing = true;

            gridEditorController.GridEditorState.Polygons.AddRange(new[]
            {
                new Feature(),
                new Feature(),
                new Feature()
            });

            // Act
            gridEditorController.MakeTriangularGridInPolygon();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GetMeshBoundaryPolygon_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            gridEditorController.GetMeshBoundaryPolygon();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GetMeshBoundaryPolygon_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var grid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();


            var gridGeomId = 1;
            var numberOfPolygonVertices = -1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Repeat.Any().Return(disposableMeshGeometry);
            api.Expect(a => a.CountMeshBoundaryPolygonVertices(gridGeomId, ref numberOfPolygonVertices)).WhenCalled((m) =>
                {
                    m.Arguments[1] = 10;
                }).Repeat.Once().Return(true);

            api.Expect(a => a.GetMeshBoundaryPolygon(Arg<int>.Is.Equal(gridGeomId),
                ref Arg<DisposableGeometryList>.Ref(Is.Anything(), null).Dummy)).Repeat.Once().Return(true);
            ribbonState.Expect(r => r.SelectedGrid).Return(grid).Repeat.Any();


            mocks.ReplayAll();

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;
            gridEditorController.IsEditing = true;

            // Act
            gridEditorController.GetMeshBoundaryPolygon();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void RefineSelectedPolygons_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var gridEditorController = new GridEditorController();
            var mocks = new MockRepository();
            var grid = new UnstructuredGrid();
            var coordinates = new[]
            {
                new Coordinate(0, 0, 0),
                new Coordinate(10, 0, 0),
                new Coordinate(5, 10, 0),
                new Coordinate(0, 0, 0)
            };
            var polygon = new Polygon(new LinearRing(coordinates));
            var polygonSelection = new PolygonSelection
            {
                SelectedIndices = new[] {0, 0},
                Feature = new Feature2DPolygon {Geometry = polygon}
            };

            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var interactor = mocks.StrictMock<IMapInteractor>();
            var renderer = mocks.StrictMock <IGridEditorStateRenderer>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();
            var DisposableGeometryList = new List<IGeometry> {polygonSelection.Feature.Geometry}.DisposableGeometryListFromGeometries();
            
            int gridGeomId = 1;
            int startIndex = 0;
            int endIndex = 0;
            double distance = 10.0;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Repeat.Any().Return(disposableMeshGeometry);
            api.Expect(a => a.CountVerticesRefinededPolygon(
                Arg<int>.Is.Equal(gridGeomId),
                ref Arg<DisposableGeometryList>.Ref(Is.Anything(), DisposableGeometryList).Dummy,
                Arg<int>.Is.Equal(startIndex),
                Arg<int>.Is.Equal(endIndex),
                Arg<double>.Is.Equal(distance), 
                ref Arg<int>.Ref(Is.Anything(), 0).Dummy)).WhenCalled((m) =>
            {
                m.Arguments[5] = 10;
            }).Repeat.Once().Return(true);

            api.Expect(a => a.GetRefinededPolygon(Arg<int>.Is.Equal(gridGeomId),
                ref Arg<DisposableGeometryList>.Ref(Is.Anything(), DisposableGeometryList).Dummy,
                Arg<int>.Is.Equal(startIndex),
                Arg<int>.Is.Equal(endIndex),
                Arg<double>.Is.Equal(distance),
                ref Arg<DisposableGeometryList>.Ref(Is.Anything(), DisposableGeometryList).Dummy)).Repeat.Once().Return(true);

            ribbonState.Expect(r => r.SelectedGrid).Return(grid).Repeat.Once();
            interactor.Expect(i => i.PolygonSelections).Return(new List<IPolygonSelection>{polygonSelection}).Repeat.Once();
            interactor.Expect(i => i.Renderer).Return(renderer).Repeat.Twice();
            renderer.Expect(r => r.MeshGeometryRenderer.MeshChanged()).Repeat.Once();
            renderer.Expect(r => r.Refresh()).Repeat.Once();
            
            mocks.ReplayAll();

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            gridEditorController.RibbonState = ribbonState;
            gridEditorController.CreateApi = () => api;
            gridEditorController.IsEditing = true;
            gridEditorController.MapInteractor = interactor;

            // Act
            gridEditorController.RefineSelectedPolygons();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void RefineGridBasedOnSamples_ShouldCallApi_WhenInEditMode()
        {
            // Arrange
            var mocks = new MockRepository();
            var grid = new UnstructuredGrid();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var interactor = mocks.StrictMock<IMapInteractor>();
            var renderer = mocks.StrictMock<IGridEditorStateRenderer>();
            var meshGeometryRenderer = mocks.StrictMock<IDisposableMeshGeometryRenderer>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();

            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();
            var interpolationParameters = new InterpolationParameters();
            var samplesRefineParameters = new SamplesRefineParameters();

            int gridGeomId = 1;
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.GetGridState(gridGeomId)).Repeat.Any().Return(disposableMeshGeometry);

            interactor.Expect(i => i.Renderer).Return(renderer).Repeat.Any();
            interactor.Expect(i => i.StartInteraction(null)).IgnoreArguments();
            interactor.Expect(i => i.GetSplineGeometry = null).IgnoreArguments();
            interactor.Expect(i => i.InsertVertex = null).IgnoreArguments();
            interactor.Expect(i => i.InsertEdge = null).IgnoreArguments();
            interactor.Expect(i => i.GetVertexIndex = null).IgnoreArguments();
            interactor.Expect(i => i.DeleteVertex = null).IgnoreArguments();
            interactor.Expect(i => i.MergeVertices = null).IgnoreArguments();
            interactor.Expect(i => i.DeleteEdges = null).IgnoreArguments();
            interactor.Expect(i => i.MoveVertex = null).IgnoreArguments();
            renderer.Expect(r => r.MeshGeometryRenderer).Return(meshGeometryRenderer).Repeat.Any();
            renderer.Expect(r => r.Refresh()).Repeat.Once();

            meshGeometryRenderer.Expect(m => m.MeshChanged());
            meshGeometryRenderer.Expect(m => m.GetSelectedVertices = null).IgnoreArguments();

            ribbonState.Expect(r => r.RefreshState()).Repeat.Once();
            ribbonState.Expect(r => r.SelectedGrid).Return(grid).Repeat.Any();

            api.Expect(a => a.RefineGridBasedOnSamples(
                Arg<int>.Is.Equal(gridGeomId),
                ref Arg<DisposableGeometryList>.Ref(Is.Anything(), null).Dummy,
                 Arg<InterpolationParameters>.Is.Equal(interpolationParameters),
                Arg<SamplesRefineParameters>.Is.Equal(samplesRefineParameters))).IgnoreArguments().Repeat.Any().Return(true);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api,
                GridEditorState =
                {
                    MeshGeometry = new DisposableMeshGeometry()
                },
                MapInteractor = interactor
            };

            gridEditorController.IsEditing = true;

            var pointValueList = new List<IPointValue>
            {
                new PointValue
                {
                    X = 10.0,
                    Y = 10.0,
                    Value = 20.0
                },
                new PointValue
                {
                    X = 20.0,
                    Y = 20.0,
                    Value = 30.0
                }
            };

            gridEditorController.GridEditorState.SamplePoints.PointValues = pointValueList;

            // Act
            gridEditorController.RefineGridBasedOnSamples();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void MakeOffsettedPolygon_ShouldCallApiWithSelectedPolygons_WhenInEditMode()
        {
            // Arrange
            var polygon1 = new Polygon(new LinearRing(new[] { new Coordinate(0, 0), new Coordinate(10, 10), new Coordinate(20, 0), new Coordinate(0, 0) }));
            var polygon2 = new Polygon(new LinearRing(new[] { new Coordinate(0, 0), new Coordinate(10, 10), new Coordinate(20, 0), new Coordinate(0, 0) }));

            var polygonFeature1 = new Feature { Geometry = polygon1 };
            var polygonFeature2 = new Feature { Geometry = polygon2 };

            var mocks = new MockRepository();
            var grid = new UnstructuredGrid();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var interactor = mocks.StrictMock<IMapInteractor>();
            
            interactor.Expect(i => i.PolygonSelections).Return(new IPolygonSelection[]
                {
                    new PolygonSelection {Feature = polygonFeature1, SelectedIndices = new[] {1}}
                })
                .Repeat.Any();

            interactor.Expect(i => i.StartInteraction(null)).IgnoreArguments();
            interactor.Expect(i => i.GetSplineGeometry = null).IgnoreArguments();
            interactor.Expect(i => i.InsertVertex = null).IgnoreArguments();
            interactor.Expect(i => i.InsertEdge = null).IgnoreArguments();
            interactor.Expect(i => i.GetVertexIndex = null).IgnoreArguments();
            interactor.Expect(i => i.DeleteVertex = null).IgnoreArguments();
            interactor.Expect(i => i.DeleteEdges = null).IgnoreArguments();
            interactor.Expect(i => i.MergeVertices = null).IgnoreArguments();
            interactor.Expect(i => i.MoveVertex = null).IgnoreArguments();
            interactor.Expect(i => i.Renderer).Return(null).Repeat.Any();

            ribbonState.Expect(r => r.SelectedGrid).Return(grid).Repeat.Any();
            ribbonState.Expect(r => r.RefreshState()).Repeat.Once();

            int gridGeomId = 1;

            // from edit mode
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);

            // expected api calls
            DisposableGeometryList disposableGeometryListIn = null;
            var numberOfPolygonVertices = 10;

            api.Expect(a => a.CountVerticesOffsettedPolygon(gridGeomId, ref disposableGeometryListIn, true, 0, ref numberOfPolygonVertices)).IgnoreArguments().Return(true).WhenCalled(m => m.Arguments[4] = 10);
            DisposableGeometryList disposableGeometryListOut = new DisposableGeometryList();
            api.Expect(a => a.GetOffsettedPolygon(gridGeomId, ref disposableGeometryListIn, true, 0, ref disposableGeometryListOut)).IgnoreArguments().Return(true);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api,
                MapInteractor = interactor,
                RibbonState = ribbonState
            };

            gridEditorController.IsEditing = true;

            gridEditorController.GridEditorState.Polygons.Add(polygonFeature1);
            gridEditorController.GridEditorState.Polygons.Add(polygonFeature2);

            // Act
            gridEditorController.AddOffsettedPolygon();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenGridEditorController_SettingPolygonRefinementDistance_UpdatesGridEditorState()
        {
            //Arrange
            var gridEditorController = new GridEditorController();

            Assert.AreEqual(10, gridEditorController.GridEditorState.PolygonRefinementDistance);

            // Act
            gridEditorController.PolygonRefinementDistance = 100;

            // Assert
            Assert.AreEqual(100, gridEditorController.GridEditorState.PolygonRefinementDistance);
        }

        [Test]
        public void GivenGridEditorController_SettingPolygonOffsetDistance_UpdatesGridEditorState()
        {
            //Arrange
            var gridEditorController = new GridEditorController();

            Assert.AreEqual(10, gridEditorController.GridEditorState.PolygonOffsetDistance);

            // Act
            gridEditorController.PolygonOffsetDistance = 100;

            // Assert
            Assert.AreEqual(100, gridEditorController.GridEditorState.PolygonOffsetDistance);
        }

        [Test]
        public void GivenGridEditorController_SettingDrawFilledInPolygons_UpdatesMapInteractorRendererPolygonRendererDrawFilledInPolygons()
        {
            //Arrange
            var mocks = new MockRepository();
            var interactor = mocks.StrictMock<IMapInteractor>();
            var renderer = mocks.StrictMock<IGridEditorStateRenderer>();
            var polygonRenderer = mocks.StrictMock<IGridEditorPolygonRenderer>();

            interactor.Expect(i => i.Renderer).Return(renderer).Repeat.Any();
            renderer.Expect(r => r.PolygonRenderer).Return(polygonRenderer);
            renderer.Expect(r => r.Refresh());
            
            polygonRenderer.Expect(r => r.DrawFilledInPolygons = false);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                MapInteractor = interactor
            };

            // Act
            gridEditorController.DrawFilledInPolygons = false;

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenGridEditorController_SettingDrawSelectedVertices_SetsMeshGeometryRendererDrawSelectedVertices()
        {
            var mocks = new MockRepository();
            var interactor = mocks.StrictMock<IMapInteractor>();
            var renderer = mocks.StrictMock<IGridEditorStateRenderer>();
            var meshRenderer = mocks.StrictMock<IDisposableMeshGeometryRenderer>();

            interactor.Expect(i => i.Renderer).Return(renderer).Repeat.Any();
            renderer.Expect(r => r.MeshGeometryRenderer).Return(meshRenderer).Repeat.Any();
            renderer.Expect(r => r.Refresh());

            meshRenderer.Expect(r => r.DrawSelectedVertices = true);
            
            mocks.ReplayAll();

            var controller = new GridEditorController
            {
                MapInteractor = interactor
            };

            controller.DrawSelectedVertices = true;

            mocks.VerifyAll();
        }

        [Test] public void GivenGridEditorController_EnableMapToolByType_ShouldCallMapInteractorActivateToolAndRefreshRibbon()
        {
            //Arrange
            var mocks = new MockRepository();
            var interactor = mocks.StrictMock<IMapInteractor>();
            var ribbon = mocks.StrictMock<IGridRibbonState>();

            interactor.Expect(i => i.SelectedMapToolType).Return(MapToolType.Polygon).Repeat.Any();
            interactor.Expect(i => i.ActivateTool(MapToolType.InsertVertices));

            ribbon.Expect(r => r.RefreshState());

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                MapInteractor = interactor,
                RibbonState = ribbon
            };

            // Act
            gridEditorController.EnableMapToolByType(MapToolType.InsertVertices);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GetVertexIndex_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            var coordinate = new Coordinate(10, 10);
            double searchRadius = 20.0;
            gridEditorController.GetVertexIndex(coordinate, searchRadius);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GetVertexIndex_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableGeometryList = new DisposableGeometryList();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            int index = 0;
            var verticesCoordinates = new[] { new Coordinate(0, 0), new Coordinate(10, 10) };
            var verticesIndexes = new[] { 0, 1 };
            double searchRadius = 10.0;
            
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.InsertVertex(gridGeomId, verticesCoordinates[0].X, verticesCoordinates[0].Y, verticesCoordinates[0].Z, ref index)).Repeat.Twice().Return(true).IgnoreArguments();
            api.Expect(a => a.InsertEdge(gridGeomId, verticesIndexes[0], verticesIndexes[1], ref index)).Repeat.Once().Return(true).IgnoreArguments();
            api.Expect(a => a.GetVertexIndex(gridGeomId, ref disposableGeometryList, searchRadius, ref index)).Repeat.Once().Return(true).IgnoreArguments();

            api.Expect(a => a.GetGridState(gridGeomId)).Repeat.Twice().Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api
            };

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.InsertVertex(verticesCoordinates[0]);
            gridEditorController.InsertVertex(verticesCoordinates[1]);
            gridEditorController.InsertEdge(0, 1);
            gridEditorController.GetVertexIndex(verticesCoordinates[0], searchRadius);

            // Assert
            mocks.VerifyAll();
        }


        [Test]
        public void DeleteVertex_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            int index = 0;
            gridEditorController.DeleteVertex(index);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void DeleteVertex_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();

            var gridGeomId = 1;
            int index = 0;
            var verticesCoordinate = new Coordinate(0, 0);
            
            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.InsertVertex(gridGeomId, verticesCoordinate.X, verticesCoordinate.Y, verticesCoordinate.Z, ref index)).Repeat.Once().Return(true).IgnoreArguments();
            api.Expect(a => a.DeleteVertex(gridGeomId,index)).Repeat.Once().Return(true).IgnoreArguments();

            api.Expect(a => a.GetGridState(gridGeomId)).Repeat.Once().Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api
            };

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.InsertVertex(verticesCoordinate);
            gridEditorController.DeleteVertex(0);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void DeleteEdges_ShouldNotCallApi_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                CreateApi = () => api
            };

            // Act
            var coordinate = new Coordinate(10, 10);
            double searchRadius = 20.0;
            gridEditorController.DeleteEdges(coordinate, searchRadius);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void DeleteEdges_ShouldCallApi_WhenInEditMode()
        {
            //Arrange
            var unstructuredGrid = new UnstructuredGrid();

            var mocks = new MockRepository();
            var api = mocks.StrictMock<IGridgeomStatefulApi>();
            var ribbonState = mocks.StrictMock<IGridRibbonState>();
            var disposableMeshGeometry = mocks.StrictMock<DisposableMeshGeometry>();
            var disposableGeometryList = new DisposableGeometryList();

            var gridGeomId = 1;
            int index = 0;
            var verticesCoordinate = new Coordinate(0, 0);
            double searchRadius = 10.0;

            api.Expect(a => a.CreateGridState()).Return(gridGeomId);
            api.Expect(a => a.SetGridState(0, null, false)).IgnoreArguments().Return(true);
            api.Expect(a => a.InsertVertex(gridGeomId, verticesCoordinate.X, verticesCoordinate.Y, verticesCoordinate.Z, ref index)).Repeat.Twice().Return(true).IgnoreArguments();
            api.Expect(a => a.DeleteEdge(gridGeomId, ref disposableGeometryList, index)).Repeat.Once().Return(true).IgnoreArguments();

            api.Expect(a => a.GetGridState(gridGeomId)).Repeat.Once().Return(disposableMeshGeometry);
            ribbonState.Expect(r => r.SelectedGrid).Return(unstructuredGrid);

            mocks.ReplayAll();

            var gridEditorController = new GridEditorController
            {
                RibbonState = ribbonState,
                CreateApi = () => api
            };

            Assert.IsNull(gridEditorController.GridEditorState.MeshGeometry);

            // Act
            gridEditorController.IsEditing = true;
            gridEditorController.InsertVertex(verticesCoordinate);
            gridEditorController.InsertVertex(verticesCoordinate);
            gridEditorController.DeleteEdges(verticesCoordinate, searchRadius);

            // Assert
            mocks.VerifyAll();
        }
    }
}
