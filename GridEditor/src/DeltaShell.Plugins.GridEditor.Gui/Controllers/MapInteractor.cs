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

        /// <inheritdoc />
        public IChangePolygonMapTool ChangePolygonMapTool { get; private set; }

        /// <summary>
        /// Map tool for inserting edges
        /// </summary>
        public InsertEdgesMapTool InsertEdgesMapTool { get; private set; }

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
            get { return InsertEdgesMapTool.InsertVertex; }
            set
            {
                InsertEdgesMapTool.InsertVertex = value;
            }
        }

        /// <inheritdoc />
        public Func<int, int, int> InsertEdge
        {
            get { return InsertEdgesMapTool.InsertEdge; }
            set
            {
                InsertEdgesMapTool.InsertEdge = value;
            }
        }

        /// <inheritdoc />
        public Func<double, double, double, double, int> GetVertexIndex
        {
            get { return InsertEdgesMapTool.GetVertexIndex; }
            set
            {
                InsertEdgesMapTool.GetVertexIndex = value;
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
            
            InsertEdgesMapTool = new InsertEdgesMapTool { GridEditorState = gridEditorState };

            ChangePolygonMapTool = new ChangePolygonMapTool
            {
                MoveTool =
                {
                    LayerFilter = l => l?.CustomRenderers?.Contains((IFeatureRenderer) Renderer.PolygonRenderer) ?? false,
                }
            };

            mapControl?.Tools?.AddRange(new[]
            {
                (IMapTool) PolygonMapTool,
                (IMapTool) SplineMapTool,
                (IMapTool) LandBoundariesMapTool,
                ((ChangePolygonMapTool) ChangePolygonMapTool).MoveTool,
                InsertEdgesMapTool
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

            mapControl?.Tools?.Remove(((ChangePolygonMapTool) ChangePolygonMapTool).MoveTool);

            PolygonMapTool = null;
            SplineMapTool = null;
            LandBoundariesMapTool = null;
            ChangePolygonMapTool = null;
            InsertEdgesMapTool = null;

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
                case MapToolType.ChangePolygon:
                    mapControl?.ActivateTool(((ChangePolygonMapTool) ChangePolygonMapTool).MoveTool);
                    break;
                case MapToolType.Spline:
                    mapControl?.ActivateTool((IMapTool) SplineMapTool);
                    break;
                case MapToolType.LandBoundaries:
                    mapControl?.ActivateTool((IMapTool)LandBoundariesMapTool);
                    break;
                case MapToolType.InsertEdges:
                    mapControl?.ActivateTool(InsertEdgesMapTool);
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

        /// <inheritdoc />
        public IEnumerable<UnstructuredGrid> GetUnstructuredGrids()
        {
            return GetCurrentMapControl?.Invoke()?.Map?.GetAllVisibleLayers(false)?.OfType<UnstructuredGridLayer>()
                .Select(l => l.Grid);
        }

        private void OnToolActivated(object sender, EventArgs<IMapTool> e)
        {
            SelectedMapToolType = MapToolType.None;
        }
    }
}