using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors;
using GeoAPI.Extensions.CoordinateSystems;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Coverages;
using NetTopologySuite.Extensions.Features;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.Data.Providers;
using SharpMap.Editors.Interactors;
using SharpMap.Layers;
using SharpMap.Styles;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers
{
    public sealed class GridEditorStateGroupRenderer : GroupLayer, IGridEditorStateRenderer
    {
        private GridEditorState gridEditorState;
        private readonly DisposableMeshGeometryRenderer disposableMeshGeometryRenderer;
        private readonly VectorLayer polygonsLayer;
        private readonly VectorLayer splinesLayer;
        private readonly VectorLayer landBoundaryLayer;
        private readonly PointCloudLayer pointCloudLayer;

        public GridEditorStateGroupRenderer() : base("Grid editor")
        {
            disposableMeshGeometryRenderer = new DisposableMeshGeometryRenderer {Name = "Mesh"};

            PolygonRenderer = new GridEditorPolygonRenderer();
            polygonsLayer = new VectorLayer
            {
                Name = "Polygons",
                ReadOnly = false,
                DataSource = new FeatureCollection(),
                CustomRenderers = new List<IFeatureRenderer>(new []{(IFeatureRenderer) PolygonRenderer}),
                FeatureEditor = new GenericFeatureEditor<PolygonEditorFeatureInteractor>()
            };

            splinesLayer = new VectorLayer
            {
                Name = "Splines",
                DataSource = new FeatureCollection{FeatureType = typeof(Spline)},
                CustomRenderers = new List<IFeatureRenderer>(new IFeatureRenderer[] { new GridEditorSplineRenderer() }),
                FeatureEditor = new GenericFeatureEditor<SplineEditorFeatureInteractor>()
            };

            landBoundaryLayer = new VectorLayer
            {
                Name = "Land boundaries",
                Style = new VectorStyle
                {
                    Line = new Pen(Color.Purple, 2)
                    {
                        DashStyle = DashStyle.Dash
                    },
                    Outline = new Pen(Color.Transparent,0)
                },
                DataSource = new FeatureCollection {FeatureType = typeof(Feature)},
                FeatureEditor = new Feature2DEditor(null)
            };

            var pointCloudFeatureProvider = new PointCloudFeatureProvider();
            
            pointCloudLayer = (PointCloudLayer) SharpMapLayerFactory.CreateLayer(pointCloudFeatureProvider);
            pointCloudLayer.Name = "Samples";
            pointCloudLayer.AutoUpdateThemeOnDataSourceChanged = true;

            Layers.AddRange(new ILayer[]
            {
                disposableMeshGeometryRenderer,
                polygonsLayer,
                splinesLayer,
                landBoundaryLayer,
                pointCloudLayer
            });

            layersReadOnly = true;
        }

        /// <inheritdoc />
        public GridEditorState GridEditorState
        {
            get { return gridEditorState; }
            set
            {
                gridEditorState = value;

                disposableMeshGeometryRenderer.GetMeshGeometry = () => gridEditorState?.MeshGeometry;
                disposableMeshGeometryRenderer.GetCoordinateSystem = () => gridEditorState?.MeshCoordinateSystem;

                SyncLayerDataSource(polygonsLayer, gridEditorState?.MeshCoordinateSystem, (IList) gridEditorState?.Polygons ?? new List<Feature>());
                SyncLayerDataSource(splinesLayer, gridEditorState?.MeshCoordinateSystem, (IList)gridEditorState?.Splines ?? new List<Spline>());
                SyncLayerDataSource(landBoundaryLayer, gridEditorState?.MeshCoordinateSystem, (IList)gridEditorState?.LandBoundaries ?? new List<Feature>());

                var pointCloudFeatureProvider = pointCloudLayer.DataSource as PointCloudFeatureProvider;
                if (pointCloudFeatureProvider != null)
                {
                    pointCloudFeatureProvider.CoordinateSystem = gridEditorState?.MeshCoordinateSystem;
                    pointCloudFeatureProvider.PointCloud = gridEditorState?.SamplePoints ?? new PointCloud();
                }

                Refresh();
            }
        }

        private static void SyncLayerDataSource(ILayer vectorLayer, ICoordinateSystem coordinateSystem, IList features)
        {
            var featureCollection = vectorLayer?.DataSource as FeatureCollection;
            if (featureCollection == null) return;

            featureCollection.CoordinateSystem = coordinateSystem;
            featureCollection.Features = features;
        }

        /// <inheritdoc />
        public IGridEditorPolygonRenderer PolygonRenderer { get; }

        /// <inheritdoc />
        public IDisposableMeshGeometryRenderer MeshGeometryRenderer
        {
            get { return disposableMeshGeometryRenderer; }
        }

        /// <inheritdoc />
        public void Refresh()
        {
            RenderRequired = true;
        }

        public void SetGetSplineGeometryFunction(Func<Coordinate[], Coordinate[]> getSplineGeometryFunction)
        {
            var splineslayerFeatureEditor = (GenericFeatureEditor<SplineEditorFeatureInteractor>)splinesLayer?.FeatureEditor;
            if (splineslayerFeatureEditor == null) return;

            splineslayerFeatureEditor.AfterCreate = i => i.GetSplineCoordinates = getSplineGeometryFunction;
        }

        public override void Dispose()
        {
            base.Dispose();
            GridEditorState = null;
        }
    }
}