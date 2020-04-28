using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Threading.Tasks;
using DelftTools.Utils.Collections;
using DelftTools.Utils.Collections.Extensions;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Native;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Extensions.Coverages;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using log4net;
using NetTopologySuite.Extensions.Coverages;
using NetTopologySuite.Extensions.Grids;
using NetTopologySuite.Geometries;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Controllers
{
    /// <summary>
    /// Class controlling the grid editing functionality
    /// </summary>
    internal class GridEditorController : IGridEditorController
    {
        private bool isEditing;
        private IGridgeomStatefulApi api;
        private int gridGeomId;
        private UnstructuredGrid selectedGrid;
        private bool invertVerticesSelection = false;
        private static readonly ILog log = LogManager.GetLogger(typeof(GridEditorController));
        private bool cancelRequested = false;

        public GridEditorController()
        {
            CreateApi = () => new GridgeomStatefulApiRemote();
        }
        
        /// <inheritdoc/>
        public bool IsEditing
        {
            get { return isEditing; }
            set
            {
                selectedGrid = RibbonState.SelectedGrid;
                if (selectedGrid == null) return;

                isEditing = value;
                
                if (isEditing)
                {
                    EnterEditMode();
                }
                else
                {
                    ExitEditMode();
                }
            }
        }

        /// <inheritdoc/>
        public float VertexSize
        {
            get { return MapInteractor?.Renderer?.MeshGeometryRenderer.PointSize ?? 4f; }
            set
            {
                if (MapInteractor?.Renderer == null) return;
                MapInteractor.Renderer.MeshGeometryRenderer.PointSize = value;
            }
        } 

        /// <inheritdoc/>
        public bool DrawSelectedVertices
        {
            get { return MapInteractor?.Renderer?.MeshGeometryRenderer.DrawSelectedVertices ?? false; }
            set
            {
                if (MapInteractor?.Renderer == null) return;
                MapInteractor.Renderer.MeshGeometryRenderer.DrawSelectedVertices = value;
                MapInteractor.Renderer.Refresh();
            }
        }

        /// <inheritdoc/>
        public bool DrawFilledInPolygons
        {
            get { return MapInteractor?.Renderer?.PolygonRenderer.DrawFilledInPolygons ?? false; }
            set
            {
                if (MapInteractor?.Renderer == null) return;
                MapInteractor.Renderer.PolygonRenderer.DrawFilledInPolygons = value;
                MapInteractor.Renderer.Refresh();
            }
        }

        /// <inheritdoc/>
        public IMapInteractor MapInteractor { get; set; }

        /// <inheritdoc/>
        public Func<bool> CommitState { get; set; }

        /// <inheritdoc/>
        public bool ShowVertices
        {
            get { return MapInteractor?.Renderer?.MeshGeometryRenderer.ShowGridVertices ?? false; }
            set
            {
                if (MapInteractor?.Renderer == null) return;
                MapInteractor.Renderer.MeshGeometryRenderer.ShowGridVertices = value;
            }
        }

        /// <inheritdoc/>
        public bool HasSelectedVertices
        {
            get
            {
                var selectedPointsLength = GridEditorState?.SelectedVertices?.Length;
                return selectedPointsLength.HasValue && selectedPointsLength > 0;
            }
        }

        /// <inheritdoc/>
        public bool HasSamples
        {
            get { return (GridEditorState?.SamplePoints?.PointValues?.Count ?? 0) > 0; }
        }

        /// <inheritdoc/>
        public GenerateGridWorkFlowType SelectedGenerateGridWorkFlowType { get; set; }

        /// <summary>
        /// Grid editor ribbon state (see:<see cref="IGridRibbonState"/>)
        /// </summary>
        public IGridRibbonState RibbonState { get; set; }

        /// <inheritdoc />
        public OrthogonalizationParameters OrthogonalizationParameters
        {
            get { return GridEditorState?.OrthogonalizationParameters; }
        }

        /// <inheritdoc />
        public MakeGridParameters MakeGridParameters
        {
            get { return GridEditorState?.MakeGridParameters; }
        }

        /// <inheritdoc />
        public CurvilinearParameters CurvilinearParameters 
        {
            get { return GridEditorState?.CurvilinearParameters; }
        }
        
        /// <inheritdoc />
        public InterpolationParameters InterpolationParameters
        {
            get { return GridEditorState?.InterpolationParameters; }
        }

        /// <inheritdoc />
        public SamplesRefineParameters SamplesRefineParameters
        {
            get { return GridEditorState?.SamplesRefineParameters; }
        }

        /// <inheritdoc />
        public SplinesToCurvilinearParameters SplinesToCurvilinearParameters
        {
            get { return GridEditorState?.SplinesToCurvilinearParameters; } 

        }

        /// <inheritdoc />
        public double PolygonRefinementDistance
        {
            get
            {
                return GridEditorState.PolygonRefinementDistance;
            }
            set
            {
                GridEditorState.PolygonRefinementDistance = value;
            }
        }

        /// <inheritdoc />
        public double PolygonOffsetDistance
        {
            get
            {
                return GridEditorState.PolygonOffsetDistance;
            }
            set
            {
                GridEditorState.PolygonOffsetDistance = value;
            }
        }

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public bool IsTriangulationRequired { get; set; }

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public bool IsAccountingForLandBoundariesRequired { get; set; }

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public ProjectToLandBoundaryOptions ProjectToLandBoundaryOption { get; set; } =
            ProjectToLandBoundaryOptions.NetBoundaryToLandBoundary;

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public DeleteMeshOptions DeleteMeshOption { get; set; } = DeleteMeshOptions.AllVerticesInside;

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public bool InvertVerticesSelection
        {
            get { return invertVerticesSelection; }
            set
            {
                invertVerticesSelection = value;
                GridEditorState.SelectedVertices = null;
                MapInteractor?.Renderer?.MeshGeometryRenderer?.MeshChanged();
            }
        }

        /// <summary>
        /// State of the grid editor
        /// </summary>
        public GridEditorState GridEditorState { get; } = new GridEditorState();

        /// <inheritdoc />
        public MapToolType SelectedMapToolType
        {
            get
            {
                return MapInteractor?.SelectedMapToolType ?? MapToolType.None;
            }
        }

        internal Func<IGridgeomStatefulApi> CreateApi { get; set; }

        /// <inheritdoc />
        public void ResetUnstructuredGrids()
        {
            RibbonState.UnstructuredGrids = MapInteractor?.GetUnstructuredGrids()?.ToList();
        }

        /// <inheritdoc />
        public void EnableMapToolByType(MapToolType mapToolType)
        {
            var selectedMapToolType = SelectedMapToolType != mapToolType 
                ? mapToolType 
                : MapToolType.None;

            MapInteractor.ActivateTool(selectedMapToolType);
            RibbonState?.RefreshState();
        }

        /// <inheritdoc/>
        public void DeleteSelectedVertices()
        {
            if (GridEditorState == null || !isEditing || api == null) return;

            var nodes = GridEditorState.SelectedVertices;
            if (nodes == null)
            {
                GridEditorState.MeshGeometry?.Dispose();
                GridEditorState.MeshGeometry = new DisposableMeshGeometry(new UnstructuredGrid());
                DoWithApiAndRefreshState(api => api.SetGridState(gridGeomId, GridEditorState.MeshGeometry,
                    GridEditorState.MeshCoordinateSystem?.IsGeographic ?? false), true);
            }
            else
            {
                bool successful = true;
                for (int i = 0; i < nodes.Length; i++)
                {
                    successful = DoWithApi(api => api.DeleteVertex(gridGeomId, nodes[i]));
                    if (!successful) break;
                }

                if (successful)
                {
                    RefreshGridStateMesh(true);
                }
            }

            
        }

        /// <inheritdoc/>
        public void DeleteMeshWithOptions()
        {

            if (!isEditing || api == null) return;
            var disposableGeometryList = GetSelectedPolygons()?.ToDisposableGeometryList();

            int deleteOption = (int) DeleteMeshOption;
            DoWithApiAndRefreshState(api => api.DeleteMeshWithOptions(gridGeomId, ref disposableGeometryList, ref deleteOption), true);
        }

        /// <inheritdoc/>
        public void DeleteSamplesInPolygons()
        {

            if (!isEditing || api == null) return;

            var disposableGeometryList = GetSelectedPolygons()?.ToDisposableGeometryList();
            var samples = GridEditorState?.SamplePoints?.PointValues ?? new List<IPointValue>();
            var inputPoints = samples.ToDisposableGeometryList();
            var selectedSamples = samples.ToDisposableGeometryList();

            DoWithApiAndRefreshState(api => api.PointsInPolygon(gridGeomId, ref disposableGeometryList, ref inputPoints, ref selectedSamples), false);

            PointCloud selectedCloud = new PointCloud();
            for (int i = 0; i < samples.Count; i++)
            {
                if (selectedSamples.ZCoordinates[i] < 1.0)
                {
                    selectedCloud.PointValues.Add(samples[i]);
                }
            }

            GridEditorState.SamplePoints.PointValues = selectedCloud.PointValues;
        }

        /// <inheritdoc/>
        public void FlipEdges()
        {
            if (!isEditing || api == null)
            {
                return;
            }
            DoWithApiAndRefreshState(api => api.FlipEdges(gridGeomId, IsTriangulationRequired, IsAccountingForLandBoundariesRequired, ProjectToLandBoundaryOption));

        }

        /// <inheritdoc/>
        public bool MergeVertices()
        {
            if (!isEditing || api == null) return false;
            using (var disposableGeometryList = GetSelectedPolygons()?.ToDisposableGeometryList())
            {
                DoWithApiAndRefreshState(api => api.MergeVertices(gridGeomId, disposableGeometryList), true);
            }

            return true;
        }

        /// <inheritdoc/>
        public async Task Orthogonalize()
        {
            if (!isEditing || api == null || GridEditorState == null) return;

            using (var polygons = GetSelectedPolygons()?.ToDisposableGeometryList())
            using (var landBoundaries = GridEditorState?.LandBoundaries?.ToDisposableGeometryList())
            {

                var stopwatch = new Stopwatch();

                int orthogonalizationOuterIterations = OrthogonalizationParameters.OuterIterations;
                int orthogonalizationBoundaryIterations = OrthogonalizationParameters.BoundaryIterations;
                int orthogonalizationInnerIterations = OrthogonalizationParameters.InnerIterations;
                bool state = api.OrthogonalizationInitialize(gridGeomId,
                     IsTriangulationRequired,
                     IsAccountingForLandBoundariesRequired,
                     ProjectToLandBoundaryOption,
                     OrthogonalizationParameters,
                     polygons,
                     landBoundaries);

                var totalProgressIterations = orthogonalizationOuterIterations * orthogonalizationBoundaryIterations;
                RibbonState.IsTaskRunning = true;
                RibbonState?.UpdateProgress(0, $"Starting");
                stopwatch.Start();

                for (int outerIter = 0; outerIter < orthogonalizationOuterIterations; outerIter++)
                {
                    api.OrthogonalizationPrepareOuterIteration(gridGeomId);

                    int elapsedBoundaryIter = 0;
                    for (int boundaryIter = 0; boundaryIter < orthogonalizationBoundaryIterations; boundaryIter++)
                    {
                        await Task.Run(() =>
                        {
                            for (int innerIter = 0; innerIter < orthogonalizationInnerIterations; innerIter++)
                            {
                                state = api.OrthogonalizationInnerIteration(gridGeomId);
                            }
                        });

                        double currentIteration = (outerIter * orthogonalizationBoundaryIterations) + boundaryIter;
                        RibbonState?.UpdateProgress((float)(currentIteration / totalProgressIterations * 100.0), $"Iteration {currentIteration} of {totalProgressIterations}");

                        // show intermediate results
                        elapsedBoundaryIter++;
                        if (stopwatch.Elapsed > new TimeSpan(0, 0, 0, 0, 500) && elapsedBoundaryIter > 5)
                        {
                            RefreshGridStateMesh();
                            stopwatch.Restart();
                            elapsedBoundaryIter = 0;
                        }

                        //RefreshGridStateMesh();

                        if (cancelRequested)
                        {
                            cancelRequested = false;
                            break;
                        }
                    }

                    //update mu
                    state = api.OrthogonalizationFinalizeOuterIteration(gridGeomId);

                }

                api.OrthogonalizationDelete(gridGeomId);

                RibbonState.IsTaskRunning = false;
                RibbonState?.UpdateProgress(100, $"Finished");

                stopwatch.Stop();
                RefreshGridStateMesh();

            }

        }

        public void RequestCancelation()
        {
            cancelRequested = true;
        }


        ///// <inheritdoc/>
        //public void MakeGridFromSplinesOrthogonal()
        //{

        //    if (!isEditing || api == null || GridEditorState.Splines.Count < 2) return;

        //    var disposableGeometryList = GridEditorState.Splines
        //        .Select(s => s.UserGeometry)
        //        .Where(g => g != null)
        //        .ToList()
        //        .DisposableGeometryListFromGeometries();

        //    var curvilinearParameters = GridEditorState.CurvilinearParameters;

        //    var splinesToCurvilinearParameters = GridEditorState.SplinesToCurvilinearParameters;

        //    DoWithApiAndRefreshState(api => api.MakeOrthogonalGridFromSplines(gridGeomId, ref disposableGeometryList, ref curvilinearParameters, ref splinesToCurvilinearParameters), true);

        //}


        /// <inheritdoc/>
        public async Task MakeGridFromSplinesOrthogonal()
        {

            if (!isEditing || api == null || GridEditorState.Splines.Count < 2) return;

            var disposableGeometryList = GridEditorState.Splines
                .Select(s => s.UserGeometry)
                .Where(g => g != null)
                .ToList()
                .DisposableGeometryListFromGeometries();

            var curvilinearParameters = GridEditorState.CurvilinearParameters;

            var splinesToCurvilinearParameters = GridEditorState.SplinesToCurvilinearParameters;

            var stopwatch = new Stopwatch();

            bool state = api.MakeOrthogonalGridFromSplinesInitialize(gridGeomId, ref disposableGeometryList, ref curvilinearParameters,
                ref splinesToCurvilinearParameters);

            int totalNumberOfIterations = curvilinearParameters.MRefinement;
            int elapsedBoundaryIter = 0;
            stopwatch.Start();
            for (int layer = 1; layer < totalNumberOfIterations; layer++)
            {
                await Task.Run(() => { state = api.MakeOrthogonalGridFromSplinesIteration(gridGeomId, layer); });

                RibbonState?.UpdateProgress((float)(layer / totalNumberOfIterations * 100.0), $"Iteration {layer} of {totalNumberOfIterations}");


                //if (stopwatch.Elapsed > new TimeSpan(0, 0, 0, 0, 500) && elapsedBoundaryIter > 1)
                //{
                    api.MakeOrthogonalGridFromSplinesRefreshMesh(gridGeomId);
                    RefreshGridStateMesh();
                    stopwatch.Restart();
                    elapsedBoundaryIter = 0;
                //}

                if (cancelRequested)
                {
                    cancelRequested = false;
                    break;
                }
            }

            api.MakeOrthogonalGridFromSplinesRefreshMesh(gridGeomId);

            RibbonState.IsTaskRunning = false;
            RibbonState?.UpdateProgress(100, $"Finished");

            stopwatch.Stop();
            RefreshGridStateMesh();

        }

        /// <inheritdoc/>
        public void MakeGrid()
        {
            if (!isEditing || api == null) return;

            var disposableGeometryList = GetSelectedPolygons()?.ToDisposableGeometryList();

            DoWithApiAndRefreshState(api=>api.MakeGrid(gridGeomId, GridEditorState.MakeGridParameters, ref disposableGeometryList), true);

            MapInteractor?.ZoomToGridExtent();
        }

        /// <inheritdoc/>
        public void MakeGridFromSplines()
        {
            if (!isEditing || api == null || GridEditorState.Splines.Count < 4) return;

            var disposableGeometryList = GridEditorState.Splines
                .Select(s => s.UserGeometry)
                .Where(g => g != null)
                .ToList()
                .DisposableGeometryListFromGeometries();

            var curvilinearParameters = GridEditorState.CurvilinearParameters;

            DoWithApiAndRefreshState(api=>api.MakeGridFromSplines(gridGeomId, ref disposableGeometryList, ref curvilinearParameters));
        }

        /// <inheritdoc/>
        public void MakeTriangularGridInPolygon()
        {
            var selectedPolygons = GetSelectedPolygons();

            if (!isEditing || api == null || selectedPolygons == null || selectedPolygons.Count < 1) return;

            var disposableGeometryList = selectedPolygons.ToDisposableGeometryList();
            DoWithApiAndRefreshState(api=>api.MakeTriangularGridInPolygon(gridGeomId, ref disposableGeometryList));
        }

        /// <inheritdoc/>
        public void MakeTriangularGridFromSamples()
        {
            if (!isEditing || api == null || GridEditorState.MeshGeometry == null) return;

            var samples = GridEditorState?.SamplePoints?.PointValues ?? new List<IPointValue>();
            var disposableGeometryList = samples.ToDisposableGeometryList();

            DoWithApiAndRefreshState(api=>api.MakeTriangularGridFromSamples(gridGeomId, ref disposableGeometryList));
        }

        /// <inheritdoc/>
        public void RefineSelectedPolygons()
        {
            var polygonSelections = MapInteractor.PolygonSelections.ToList();

            if (!isEditing || api == null || polygonSelections.Count < 1) return;

            foreach (var polygonSelection in polygonSelections)
            {
                var disposableGeometryList = new List<IGeometry> {polygonSelection.Feature.Geometry}.DisposableGeometryListFromGeometries();

                var distance = GridEditorState.PolygonRefinementDistance;
                var numberOfPolygonVertices = 0;
                int startIndex = 0;
                int endIndex = 0;
                if (polygonSelection.SelectedIndices.Length >= 2)
                {
                    startIndex = polygonSelection.SelectedIndices[0];
                    endIndex = polygonSelection.SelectedIndices[1];
                }

                bool successful = DoWithApi(api=>api.CountVerticesRefinededPolygon(gridGeomId, ref disposableGeometryList,
                     startIndex,  endIndex,  distance,  ref numberOfPolygonVertices));

                if (!successful) return;

                var geometryListOut = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(numberOfPolygonVertices);

                successful = DoWithApi(api => api.GetRefinededPolygon(gridGeomId, ref disposableGeometryList,
                     startIndex,  endIndex,  distance, ref geometryListOut));

                if (!successful) return;

                var feature = geometryListOut.ToFeatureList().FirstOrDefault();
                if (feature == null) continue;

                polygonSelection.Feature.Geometry = new Polygon(new LinearRing(feature.Geometry.Coordinates));
            }

            RefreshGridStateMesh();
        }

        /// <inheritdoc/>
        public void ShowOrthogonality()
        {
            if (!isEditing || api == null) return;

            var disposableMeshGeometryRenderer = MapInteractor?.Renderer?.MeshGeometryRenderer;
            if (disposableMeshGeometryRenderer != null)
            {

                disposableMeshGeometryRenderer.GetEdgeValues = () =>
                {
                    int numberOfEdges = GridEditorState.MeshGeometry.numberOfEdges;
                    var disposableGeometryListInOut = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(numberOfEdges);
                    DoWithApi(api => api.GetOrthogonality(gridGeomId, ref disposableGeometryListInOut));
                    return disposableGeometryListInOut;
                };
                disposableMeshGeometryRenderer.DrawEdgeValues = true;
            }
            
            RefreshGridStateMesh();
        }

        /// <inheritdoc/>
        public void ShowSmoothness()
        {
            if (!isEditing || api == null) return;

            var disposableMeshGeometryRenderer = MapInteractor?.Renderer?.MeshGeometryRenderer;
            if (disposableMeshGeometryRenderer != null)
            {

                disposableMeshGeometryRenderer.GetEdgeValues = () =>
                {
                    int numberOfEdges = GridEditorState.MeshGeometry.numberOfEdges;
                    var disposableGeometryListInOut = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(numberOfEdges);
                    DoWithApi(api => api.GetSmoothness(gridGeomId, ref disposableGeometryListInOut));
                    return disposableGeometryListInOut;
                };
                disposableMeshGeometryRenderer.DrawEdgeValues = true;
            }

            RefreshGridStateMesh();
        }

        /// <inheritdoc/>
        public void GenerateGridForSelectedWorkflow()
        {
            switch (SelectedGenerateGridWorkFlowType)
            {
                case GenerateGridWorkFlowType.CurveLinearTransfinite:
                    MakeGridFromSplines();
                    break;
                case GenerateGridWorkFlowType.CurveLinearOrthogonal:
                    MakeGridFromSplinesOrthogonal();
                    break;
                case GenerateGridWorkFlowType.Regular:
                    MakeGrid();
                    break;
                case GenerateGridWorkFlowType.TriangularWithPolygon:
                    MakeTriangularGridInPolygon();
                    break;
                case GenerateGridWorkFlowType.RegularWithPolygon:
                    MakeGrid();
                    break;
                case GenerateGridWorkFlowType.Samples:
                    MakeTriangularGridFromSamples();
                    break;
                case GenerateGridWorkFlowType.RefinementSamples:
                    RefineGridBasedOnSamples();
                    break;
                default:
                    throw new ArgumentOutOfRangeException($"Invalid value of {nameof(SelectedGenerateGridWorkFlowType)}");
            }
        }

        public void ClipSelectedPolygons()
        {
            throw new NotImplementedException();
        }

        public void MergeSelectedPolygons()
        {
            throw new NotImplementedException();
        }

        public void SubtractSelectedPolygons()
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public void ImportExport(ImportExportAction action, ImportExportType importExportType, string filePath)
        {
            if (filePath == null) return;

            GridEditorState.MeshCoordinateSystem = selectedGrid.CoordinateSystem;

            GridEditorState.ImportExport(action, importExportType, filePath);

            if (action == ImportExportAction.Import)
            {
                if (importExportType == ImportExportType.Grid)
                {
                    // reset state
                    DoWithApi(api=>api.SetGridState(gridGeomId, GridEditorState.MeshGeometry, GridEditorState.MeshCoordinateSystem?.IsGeographic ?? false));
                    MapInteractor?.ZoomToGridExtent();
                }

                if (importExportType == ImportExportType.Spline)
                {
                    // set spline geometry
                    GridEditorState.Splines.ForEach(s => s.Geometry = new LineString(GetSplineGeometry(s.UserGeometry.Coordinates)));
                    MapInteractor?.ZoomToSplinesExtent();
                }

                if (importExportType == ImportExportType.Samples)
                {
                    MapInteractor?.ZoomToSamplesExtent();
                }

                if (importExportType == ImportExportType.LandBoundary)
                {
                    MapInteractor?.ZoomToLandBoundariesExtent();
                }
            }

            RefreshGridStateMesh(true);
        }

        /// <inheritdoc/>
        public void GetMeshBoundaryPolygon()
        {
            if (!isEditing || api == null || GridEditorState.MeshGeometry == null) return;

            int numberOfPolygonVertices = -1;
            bool successful = DoWithApi(api=>api.CountMeshBoundaryPolygonVertices(gridGeomId, ref numberOfPolygonVertices));

            if (!successful) return;

            var disposableGeometryList = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(numberOfPolygonVertices);

            successful = DoWithApi(api=>api.GetMeshBoundaryPolygon(gridGeomId, ref disposableGeometryList));

            if (!successful) return;

            var featureList = disposableGeometryList.ToFeatureList();

            GridEditorState.Polygons.AddRange(featureList);

            GridEditorState.SelectedVertices = null;

            RefreshDrawing();
        }

        /// <inheritdoc/>
        public void RefineGridBasedOnSamples()
        {
            if (!isEditing || api == null || GridEditorState.MeshGeometry == null) return;

            var samples = GridEditorState?.SamplePoints?.PointValues ?? new List<IPointValue>();
            var disposableGeometryList = samples.ToDisposableGeometryList();

            DoWithApiAndRefreshState(api=>api.RefineGridBasedOnSamples(gridGeomId, ref disposableGeometryList, InterpolationParameters, SamplesRefineParameters));
        }

        /// <inheritdoc/>
        public void RefineGridBasedOnPolygon()
        {
            var selectedPolygons = GetSelectedPolygons();

            if (!isEditing || api == null || selectedPolygons == null || GridEditorState.MeshGeometry == null) return;

            var disposableGeometryListIn = GetSelectedPolygons()?.ToDisposableGeometryList();

            DoWithApiAndRefreshState(api=>api.RefineGridBasedOnPolygon(gridGeomId, ref disposableGeometryListIn, InterpolationParameters));
        }

        /// <inheritdoc/>
        public void AddOffsettedPolygon()
        {
            var selectedPolygons = GetSelectedPolygons();

            if (!isEditing || api == null || selectedPolygons == null || selectedPolygons.Count < 1) return;

            var disposableGeometryListIn = GetSelectedPolygons()?.ToDisposableGeometryList();

            int numberOfPolygonVertices = -1;
            bool innerPolygon = false;
            var distance = GridEditorState.PolygonOffsetDistance;
            bool successful = DoWithApi(api=>api.CountVerticesOffsettedPolygon(gridGeomId, ref disposableGeometryListIn, innerPolygon, distance, ref numberOfPolygonVertices));

            if (!successful) return;

            var disposableGeometryListOut = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(numberOfPolygonVertices);

            successful = DoWithApi(api=>api.GetOffsettedPolygon(gridGeomId, ref disposableGeometryListIn, innerPolygon, distance, ref disposableGeometryListOut));

            if (!successful) return;

            var featureList = disposableGeometryListOut.ToFeatureList();

            GridEditorState.Polygons.AddRange(featureList);

            GridEditorState.SelectedVertices = null;

            RefreshDrawing();
        }
        
        private void ExitEditMode()
        {
            if (CommitState == null || CommitState())
            {
                GridEditorState.MeshGeometry = api.GetGridStateWithCells(gridGeomId);
                CommitMeshToSelectedGrid(GridEditorState.MeshGeometry);
            }

            // reset GridEditorState and map-tools
            GridEditorState.Reset();

            // dispose api
            bool successful = DoWithApi(api => api.RemoveGridState(gridGeomId));

            if (!successful) return;

            api.Dispose();

            // stop map interaction
            MapInteractor?.StopInteraction();

            // refresh ribbon
            RibbonState.RefreshState();
        }

        private void EnterEditMode()
        {
            // prepare GridEditorState
            GridEditorState.MeshCoordinateSystem = selectedGrid?.CoordinateSystem;
            GridEditorState.MeshGeometry = new DisposableMeshGeometry(selectedGrid);

            // setup api
            InitializeApi();

            // update map-tools
            if (MapInteractor == null) return;

            // start map interaction
            MapInteractor.StartInteraction(GridEditorState);
            MapInteractor.GetSplineGeometry = GetSplineGeometry;
            MapInteractor.InsertVertex = InsertVertex;
            MapInteractor.InsertEdge = InsertEdge;
            MapInteractor.GetVertexIndex = GetVertexIndex;
            MapInteractor.DeleteVertex = DeleteVertex;
            MapInteractor.MergeVertices = MergeVertices;
            MapInteractor.DeleteEdges = DeleteEdges;
            MapInteractor.MoveVertex = MoveVertex;

            var disposableMeshGeometryRenderer = MapInteractor?.Renderer?.MeshGeometryRenderer;
            if (disposableMeshGeometryRenderer != null)
            {
                disposableMeshGeometryRenderer.GetSelectedVertices = () =>
                {
                    var selectedVertices = GridEditorState.SelectedVertices;
                    if (selectedVertices != null)
                    {
                        return selectedVertices;
                    }
                    
                    var disposableGeometryListIn = GridEditorState.Polygons.ToDisposableGeometryList();
                    int inside = InvertVerticesSelection? 0 : 1 ;
                    return (GridEditorState.SelectedVertices = api.GetSelectedVertices(gridGeomId, ref disposableGeometryListIn, inside));
                };
            }

            RibbonState.RefreshState();
        }

        private void InitializeApi()
        {
            api = CreateApi();
            gridGeomId = api.CreateGridState();

            api.SetGridState(gridGeomId, GridEditorState.MeshGeometry, GridEditorState.MeshCoordinateSystem?.IsGeographic ?? false);
        }

        public int InsertVertex(Coordinate coordinate)
        {
            int vertexIndex = -1;
            if (!isEditing || api == null) return vertexIndex;
            DoWithApi(api => api.InsertVertex(gridGeomId, coordinate.X, coordinate.Y, coordinate.Z, ref vertexIndex));
            return vertexIndex;
        }

        public int InsertEdge(int startVertexIndex, int endVertexIndex)
        {
            int edgeIndex = -1;
            if (!isEditing || api == null) return edgeIndex;
            DoWithApiAndRefreshState(api => api.InsertEdge(gridGeomId, startVertexIndex, endVertexIndex, ref edgeIndex));
            return edgeIndex;
        }

        public int GetVertexIndex(Coordinate coordinate, double searchRadius)
        {
            int vertexIndex = -1;
            if (!isEditing || api == null) return vertexIndex;

            var geometryListIn = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(1);
            geometryListIn.XCoordinates[0] = coordinate.X;
            geometryListIn.YCoordinates[0] = coordinate.Y;

            DoWithApiAndRefreshState(api => api.GetVertexIndex(gridGeomId, ref geometryListIn, searchRadius, ref vertexIndex));
            return vertexIndex;
        }

        public bool DeleteVertex(int vertexIndex)
        {
            if (!isEditing || api == null) return false;
            DoWithApiAndRefreshState(api => api.DeleteVertex(gridGeomId, vertexIndex));
            return true;
        }

        public bool DeleteEdges(Coordinate coordinate, double searchRadius)
        {
            if (!isEditing || api == null|| coordinate == null) return false;

            var geometryListIn = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(1);
            geometryListIn.XCoordinates[0] = coordinate.X;
            geometryListIn.YCoordinates[0] = coordinate.Y;

            DoWithApiAndRefreshState(api => api.DeleteEdge(gridGeomId, ref geometryListIn, searchRadius));
            return true;
        }

        public bool MoveVertex(Coordinate coordinate, int vertexIndex)
        {
            if (!isEditing || api == null || coordinate == null || vertexIndex < 0) return false;

            var geometryListIn = NativeStructConversionExtensions.CreateEmptyDisposableGeometryList(1);
            geometryListIn.XCoordinates[0] = coordinate.X;
            geometryListIn.YCoordinates[0] = coordinate.Y;

            DoWithApiAndRefreshState(api => api.MoveVertex(gridGeomId, ref geometryListIn, vertexIndex));
            return true;
        }

        private void CommitMeshToSelectedGrid(DisposableMeshGeometry meshGeometry)
        {
            selectedGrid.ResetState(meshGeometry.CreateVertices(), meshGeometry.CreateEdges(), meshGeometry.CreateCells());
        }

        private Coordinate[] GetSplineGeometry(Coordinate[] coordinates)
        {
            var emptyCoordinatesArray = Enumerable.Range(0, (GridEditorState.PointsBetweenSplineCornerPoints + 1) * coordinates.Length)
                .Select(i => new Coordinate(0, 0))
                .ToArray();

            var geometryListOut = emptyCoordinatesArray.ToDisposableGeometryList();
            using (var geometryListIn = coordinates.ToDisposableGeometryList())
            {
                var result = DoWithApi(api=>api.GetSplines(geometryListIn, ref geometryListOut, GridEditorState.PointsBetweenSplineCornerPoints));
                return result
                    ? Enumerable.Range(0, geometryListOut.NumberOfCoordinates)
                        .Select(i => new Coordinate(geometryListOut.XCoordinates[i], geometryListOut.YCoordinates[i], geometryListOut.ZCoordinates[i]))
                        .ToArray()
                    : coordinates;
            }
        }

        private void RefreshGridStateMesh(bool clearSelectedVertices = false)
        {
            if (clearSelectedVertices)
            {
                GridEditorState.SelectedVertices = null;
            }

            GridEditorState.MeshGeometry = api.GetGridState(gridGeomId);

            MapInteractor?.Renderer?.MeshGeometryRenderer.MeshChanged();

            RefreshDrawing();
        }

        private void RefreshDrawing()
        {
            MapInteractor?.Renderer?.Refresh();
        }   

        private IList<IFeature> GetSelectedPolygons()
        {
            var polygonSelections = MapInteractor?.PolygonSelections.ToList();
            
            return polygonSelections != null && polygonSelections.Any()
                ? polygonSelections.Select(s => s.Feature).ToList() 
                : GridEditorState.Polygons;
        }


        private bool DoWithApi(Func<IGridgeomStatefulApi, bool> apiFunction)
        {
            try
            {
                return apiFunction(api);
            }
            catch (InvalidOperationException)
            {
                log.Error("Api error, please try again.");

                api.Dispose();
                api = null;

                InitializeApi();

                return false;
            }
        }

        private void DoWithApiAndRefreshState(Func<IGridgeomStatefulApi, bool> apiFunction, bool clearSelectedVertices = false)
        {
            if (DoWithApi(apiFunction))
            {
                RefreshGridStateMesh(clearSelectedVertices);
            }
        }

    }
}