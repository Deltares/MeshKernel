using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using DelftTools.Utils;
using DelftTools.Utils.Collections;
using DelftTools.Utils.Collections.Extensions;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using DeltaShell.Plugins.GridEditor.Gui.MapTools;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Extensions.Grids;
using NetTopologySuite.Geometries;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.Layers;
using SharpMap.UI.Forms;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.Controllers
{
    internal class MapInteractor : IMapInteractor
    {
        private IMapControl mapControl;
        private Func<Coordinate[], Coordinate[]> getSplineGeometry;

        /// <summary>
        /// Function to obtain the <see cref="IMapControl"/> of the current view
        /// </summary>
        public Func<IMapControl> GetCurrentMapControl { get; set; }

        /// <inheritdoc />
        public IGridEditorStateRenderer Renderer { get; private set; }

        /// <summary>
        /// Map tool for creating polygons
        /// </summary>
        public IGridEditorMapTool PolygonMapTool { get; private set; }

        /// <summary>
        /// Map tool for editing splines
        /// </summary>
        public IGridEditorMapTool SplineMapTool { get; private set; }

        /// <summary>
        /// Map tool for inserting vertices
        /// </summary>
        public InsertVerticesMapTool InsertVerticesMapTool { get; private set; }

        /// <summary>
        /// Map tool for deleting edges
        /// </summary>
        public DeleteEdgesMapTool DeleteEdgesMapTool { get; private set; }

        /// <summary>
        /// Map tool for deleting a vertex
        /// </summary>
        public DeleteVerticesMapTool DeleteVerticesMapTool { get; private set; }

        /// <summary>
        /// Map tool for moving vertices
        /// </summary>
        public MoveVerticesMapTool MoveVerticesMapTool { get; private set; }

        /// <inheritdoc /> 
        public Func<Coordinate[], Coordinate[]> GetSplineGeometry
        {
            get { return getSplineGeometry;}
            set
            {
                getSplineGeometry = value;
                Renderer.SetGetSplineGeometryFunction(value);
            }
        }

        /// <inheritdoc />
        public Func<Coordinate, int> InsertVertex
        {
            get { return InsertVerticesMapTool.InsertVertex; }
            set
            {
                InsertVerticesMapTool.InsertVertex = value;
            }
        }

        /// <inheritdoc />
        public Func<int, int, int> InsertEdge
        {
            get { return InsertVerticesMapTool.InsertEdge; }
            set
            {
                InsertVerticesMapTool.InsertEdge = value;
            }
        }

        /// <inheritdoc />
        public Func<Coordinate, double, int> GetVertexIndex
        {
            get
            {
                return InsertVerticesMapTool.GetVertexIndex;
            }
            set
            {
                InsertVerticesMapTool.GetVertexIndex = value;
                DeleteVerticesMapTool.GetVertexIndex = value;
                MoveVerticesMapTool.GetVertexIndex = value;
            }
        }

        /// <inheritdoc />
        public Func<int, bool> DeleteVertex
        {
            get { return DeleteVerticesMapTool.DeleteVertex; }
            set
            {
                DeleteVerticesMapTool.DeleteVertex = value;
            }
        }

        /// <inheritdoc />
        public Func<bool> MergeVertices
        {
            get { return InsertVerticesMapTool.MergeVertices; }
            set
            {
                InsertVerticesMapTool.MergeVertices = value;
            }
        }

        /// <inheritdoc />
        public Func<Coordinate, double, bool> DeleteEdges
        {
            get { return DeleteEdgesMapTool.DeleteEdges; }
            set
            {
                DeleteEdgesMapTool.DeleteEdges = value;
            }
        }

        /// <inheritdoc />
        public Func<Coordinate, int, bool> MoveVertex
        {
            get { return MoveVerticesMapTool.MoveVertex; }
            set
            {
                MoveVerticesMapTool.MoveVertex = value;
            }
        }

        /// <inheritdoc />
        public IGridEditorMapTool LandBoundariesMapTool { get; private set; }

        /// <inheritdoc />
        public MapToolType SelectedMapToolType { get; private set; } = MapToolType.None;

        /// <inheritdoc />
        public void StartInteraction(GridEditorState gridEditorState)
        {
            mapControl = GetCurrentMapControl?.Invoke();

            Renderer = new GridEditorStateGroupRenderer { GridEditorState = gridEditorState };

            PolygonMapTool = new LineStringMapTool
            {
                LayerFilter = l => l?.CustomRenderers?.Contains((IFeatureRenderer) Renderer.PolygonRenderer) ?? false,
                GridEditorState = gridEditorState,
                LineEditorStrategies = new ILineStringMapToolStrategy[]
                {
                    new SelectLineStringMapToolStrategy(),
                    new AddNewPolygonLineStringMapToolStrategy
                    {
                        CommitGeometry = (g, s) =>
                        {
                            var updatedMultiPolygon = gridEditorState?.GetUpdatedMultiPolygon(g as IPolygon, PolygonEditAction.Add);
                            if (updatedMultiPolygon == null) return;

                            s.Polygons.Clear();
                            s.Polygons.AddRange(updatedMultiPolygon.Geometries.Select(ge => new Feature { Geometry = ge }));
                        }
                    },
                    new DraggingLineStringMapToolStrategy(), 
                    new EditingLineStringMapToolStrategy()
                }
            };

            SplineMapTool = new LineStringMapTool
            {
                LayerFilter = l => l?.Name == "Splines",
                GridEditorState = gridEditorState,
                LineEditorStrategies = new ILineStringMapToolStrategy[]
                {
                    new SelectLineStringMapToolStrategy(),
                    new AddNewLineStringMapToolStrategy
                    {
                        TransformPreviewCoordinates = (c) => c.Length < 2 ? c : GetSplineGeometry?.Invoke(c),
                        CommitGeometry = (g, s) =>
                        {
                            s.Splines.Add(new Spline
                            {
                                Geometry = new LineString(GetSplineGeometry?.Invoke(g.Coordinates)),
                                UserGeometry = g,
                                PointsBetweenSplineCornerPoints = gridEditorState.PointsBetweenSplineCornerPoints 
                            });
                        }
                    },
                    new DraggingLineStringMapToolStrategy(),
                    new EditingLineStringMapToolStrategy()
                }
            };

            LandBoundariesMapTool = new LineStringMapTool
            {
                LayerFilter = l => l?.Name == "Land boundaries",
                GridEditorState = gridEditorState,
                LineEditorStrategies = new ILineStringMapToolStrategy[]
                {
                    new SelectLineStringMapToolStrategy(),
                    new AddNewLineStringMapToolStrategy
                    {
                        CommitGeometry = (g, s) => { s.LandBoundaries.Add(new Feature {Geometry = g}); }
                    },
                    new DraggingLineStringMapToolStrategy(),
                    new EditingLineStringMapToolStrategy()
                }
            };
            
            InsertVerticesMapTool = new InsertVerticesMapTool { GridEditorState = gridEditorState };

            DeleteVerticesMapTool = new DeleteVerticesMapTool { GridEditorState = gridEditorState };

            DeleteEdgesMapTool = new DeleteEdgesMapTool { GridEditorState = gridEditorState };

            MoveVerticesMapTool = new MoveVerticesMapTool { GridEditorState = gridEditorState };

            mapControl?.Tools?.AddRange(new[]
            {
                (IMapTool) PolygonMapTool,
                (IMapTool) SplineMapTool,
                (IMapTool) LandBoundariesMapTool,
                InsertVerticesMapTool,
                DeleteVerticesMapTool,
                DeleteEdgesMapTool,
                MoveVerticesMapTool
            });

            mapControl?.Map?.Layers.Insert(0, (ILayer) Renderer);
            mapControl?.Map?.BringToFront((ILayer) Renderer);

            if (mapControl is MapControl control)
            {
                control.ToolActivated += OnToolActivated;
            }
        }

        /// <inheritdoc />
        public void StopInteraction()
        {
            // Clear and remove map tools
            var mapTools = new[] {PolygonMapTool, SplineMapTool, LandBoundariesMapTool};

            mapTools.ForEach(t => t.GridEditorState = null);
            mapTools.OfType<IDisposable>().ForEach(d => d.Dispose());
            mapTools.ForEach(t => mapControl?.Tools?.Remove((IMapTool) t));

            PolygonMapTool = null;
            SplineMapTool = null;
            LandBoundariesMapTool = null;
            InsertVerticesMapTool = null;

            // Clear and remove grid editor layer
            Renderer.GridEditorState = null;
            mapControl?.Map.Layers.Remove((ILayer) Renderer);
            Renderer = null;

            SelectedMapToolType = MapToolType.None;

            if (mapControl is MapControl control)
            {
                control.ToolActivated -= OnToolActivated;
            }

            mapControl = null;
        }

        /// <inheritdoc />
        public void ActivateTool(MapToolType mapToolType)
        {
            switch (mapToolType)
            {
                case MapToolType.Polygon:
                    mapControl?.ActivateTool((IMapTool) PolygonMapTool);
                    break;
                case MapToolType.Spline:
                    mapControl?.ActivateTool((IMapTool) SplineMapTool);
                    break;
                case MapToolType.LandBoundaries:
                    mapControl?.ActivateTool((IMapTool)LandBoundariesMapTool);
                    break;
                case MapToolType.InsertVertices:
                    mapControl?.ActivateTool(InsertVerticesMapTool);
                    break;
                case MapToolType.DeleteVertices:
                    mapControl?.ActivateTool(DeleteVerticesMapTool);
                    break;
                case MapToolType.DeleteEdges:
                    mapControl?.ActivateTool(DeleteEdgesMapTool);
                    break;
                case MapToolType.MoveVertices:
                    mapControl?.ActivateTool(MoveVerticesMapTool);
                    break;
                case MapToolType.None:
                    mapControl?.ActivateTool(mapControl.SelectTool);
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(mapToolType), mapToolType, null);
            }

            SelectedMapToolType = mapToolType;

            (mapControl as Control)?.Focus();

            mapControl?.Refresh();
        }

        /// <inheritdoc />
        public void ZoomToGridExtent()
        {
            GetCurrentMapControl?.Invoke()?.Map.ZoomToFit(Renderer?.MeshGeometryRenderer?.Envelope);
        }

        public void ZoomToSamplesExtent()
        {
            var pointValues = Renderer?.GridEditorState?.SamplePoints.PointValues;

            if (pointValues == null) return;

            double minX = pointValues.Select(s => s.X).Min();
            double maxX = pointValues.Select(s => s.X).Max();
            double minY = pointValues.Select(s => s.Y).Min();
            double maxY = pointValues.Select(s => s.Y).Max();

            var envelope = new Envelope(minX, maxX, minY, maxY);

            GetCurrentMapControl?.Invoke()?.Map.ZoomToFit(envelope);
        }

        public void ZoomToSplinesExtent()
        {
            var splines = Renderer?.GridEditorState?.Splines;

            var envelope = GetEnvelopeExtent(splines);

            if (envelope == null)
                return;

            GetCurrentMapControl?.Invoke()?.Map.ZoomToFit(envelope);
        }

        public void ZoomToLandBoundariesExtent()
        {
            var landBoundaries = Renderer?.GridEditorState?.LandBoundaries;

            var envelope = GetEnvelopeExtent(landBoundaries);

            if (envelope == null)
                return;

            GetCurrentMapControl?.Invoke()?.Map.ZoomToFit(envelope);
        }

        public IEnumerable<IPolygonSelection> PolygonSelections
        {
            get
            {
                var featureInteractors = mapControl?.SelectTool?.SelectedFeatureInteractors;
                if (featureInteractors == null) yield break;

                foreach (var featureInteractor in featureInteractors)
                {
                    yield return new PolygonSelection
                    {
                        Feature = featureInteractor.SourceFeature,
                        SelectedIndices = featureInteractor.Trackers
                            .Where(t => t.Selected)
                            .Select(t => t.Index)
                            .ToArray()
                    };
                }
            }
        }

        /// <inheritdoc />
        public IEnumerable<UnstructuredGrid> GetUnstructuredGrids()
        {
            return GetCurrentMapControl?.Invoke()?.Map?.GetAllVisibleLayers(false)?.OfType<UnstructuredGridLayer>()
                .Select(l => l.Grid).Where(g => g.IsEditable);
        }

        private void OnToolActivated(object sender, EventArgs<IMapTool> e)
        {
            SelectedMapToolType = MapToolType.None;
        }

        private Envelope GetEnvelopeExtent<T>(IList<T> features) where T : IFeature
        {

            if (features == null || features.Count == 0) 
                return null;

            double minX = features[0].Geometry.Envelope.Coordinates.Min(c => c.X);
            double maxX = features[0].Geometry.Envelope.Coordinates.Max(c => c.X);
            double minY = features[0].Geometry.Envelope.Coordinates.Min(c => c.Y);
            double maxY = features[0].Geometry.Envelope.Coordinates.Max(c => c.Y);

            foreach (var s in features)
            {

                double minSplineX = s.Geometry.Envelope.Coordinates.Min(c => c.X);
                double maxSplineX = s.Geometry.Envelope.Coordinates.Max(c => c.X);

                double minSplineY = s.Geometry.Envelope.Coordinates.Min(c => c.Y);
                double maxSplineY = s.Geometry.Envelope.Coordinates.Max(c => c.Y);

                minX = minX < minSplineX ? minX : minSplineX;
                maxX = maxX > maxSplineX ? maxX : maxSplineX;

                minY = minY < minSplineY ? minY : minSplineY;
                maxY = maxY > maxSplineY ? maxY : maxSplineY;
            }

            return new Envelope(minX, maxX, minY, maxY);

        }
    }
}