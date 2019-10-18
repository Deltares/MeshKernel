using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.UI.Forms;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool
{
    internal class LineStringMapTool : MapTool, IDisposable, IGridEditorMapTool
    {
        private readonly bool doRedraw;
        private readonly PolygonDrawingStyle style = new PolygonDrawingStyle();
        private readonly LineStringMapToolStateMachine stateMachine = new LineStringMapToolStateMachine();

        private ICollection<ILineStringMapToolStrategy> lineEditorStrategies;

        private Coordinate previousMouseMoveCoordinate;
        private IFeatureInteractor featureInteractor;
        private ILineStringMapToolStrategy currentStrategy;
        private ISelector selector;

        public LineStringMapTool(bool doRedraw = true)
        {
            this.doRedraw = doRedraw;
            GetSelector = () => new SelectToolSelectorFacade(MapControl.SelectTool);
        }

        /// <inheritdoc/>
        public GridEditorState GridEditorState { get; set; }

        /// <inheritdoc/>
        public override IMapControl MapControl
        {
            get { return base.MapControl; }
            set
            {
                if (base.MapControl != null && base.MapControl is Control previousControl)
                {
                    previousControl.PreviewKeyDown -= CurrentControlOnPreviewKeyDown;
                }

                base.MapControl = value;

                if (base.MapControl != null)
                {
                    selector = GetSelector?.Invoke();

                    if (base.MapControl is Control currentControl)
                    {
                        currentControl.PreviewKeyDown += CurrentControlOnPreviewKeyDown;
                    }
                }
                else
                {
                    selector = null;
                }
            }
        }

        internal Func<ISelector> GetSelector { get; set; }

        /// <summary>
        /// Strategies for handling commands of this tool
        /// </summary>
        public ICollection<ILineStringMapToolStrategy> LineEditorStrategies
        {
            get { return lineEditorStrategies; }
            set
            {
                lineEditorStrategies = value;
                currentStrategy = lineEditorStrategies?.FirstOrDefault();
            }
        }

        /// <inheritdoc/>
        public override void OnMouseDown(Coordinate worldPosition, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;
            var layer = GetLayer();
            if (layer == null) return;

            var feature = worldPosition.FindNearestFeature(layer);
            featureInteractor = selector?.SelectedFeatureInteractors?.FirstOrDefault(i => Equals(i.SourceFeature, feature));

            //var trackerFeature = featureInteractor?.GetTrackerAtCoordinate(worldPosition);
            var trackerFeature = worldPosition.SnappedTrackerAtCoordinate(featureInteractor);

            currentStrategy = GetNextStrategyForCommand(GetMouseDownState(feature, trackerFeature));

            currentStrategy?.HandleMouseDown(worldPosition, featureInteractor, feature, selector, trackerFeature);

            featureInteractor = selector?.SelectedFeatureInteractors?.FirstOrDefault(i => Equals(i.SourceFeature, feature));

            if (featureInteractor != null)
            {
                currentStrategy = GetNextStrategyForCommand(LineStringMapToolStateCommand.AfterMouseDownFeatureSelected);
            }

            ReDraw();
        }

        /// <inheritdoc/>
        public override void OnMouseUp(Coordinate worldPosition, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            currentStrategy?.HandleMouseUp(featureInteractor);

            currentStrategy = GetNextStrategyForCommand(LineStringMapToolStateCommand.AfterMouseUp);

            ReDraw();
        }

        /// <inheritdoc/>
        public override void OnMouseMove(Coordinate worldPosition, MouseEventArgs e)
        {
            if (previousMouseMoveCoordinate != null && worldPosition.Equals2D(previousMouseMoveCoordinate))
            {
                return;
            }

            currentStrategy?.HandleMouseMove(previousMouseMoveCoordinate, worldPosition, featureInteractor);

            previousMouseMoveCoordinate = worldPosition;

            ReDraw();
        }

        /// <inheritdoc/>
        public override void OnMouseDoubleClick(object sender, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            currentStrategy?.HandleMouseDoubleClick(previousMouseMoveCoordinate, GetLayer(), GridEditorState);

            currentStrategy = GetNextStrategyForCommand(LineStringMapToolStateCommand.AfterMouseDoubleClick);

            MapControl.Refresh();
        }

        /// <inheritdoc/>
        public override void OnKeyUp(KeyEventArgs e)
        {
            currentStrategy?.HandleKeyUp(featureInteractor);
        }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public override void OnDraw(Graphics graphics)
        {
            if (currentStrategy?.NewGeometryCoordinates != null)
            {
                PolygonDrawingHelper.DrawLinesFromCoordinates(graphics, currentStrategy.NewGeometryCoordinates, currentStrategy.UseLastSegmentStylePen, Map, style);
            }

            if (currentStrategy?.SnapResult != null)
            {
                var snappingImageCoordinate = Map.WorldToImage(currentStrategy.SnapResult.Location);
                const int halfPointSize = PolygonDrawingStyle.PointSize / 2;

                graphics.FillEllipse(style.SnappingPointBackgroundBrush, snappingImageCoordinate.X - halfPointSize,
                    snappingImageCoordinate.Y - halfPointSize, PolygonDrawingStyle.PointSize,
                    PolygonDrawingStyle.PointSize);
            }
        }

        /// <inheritdoc/>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            style?.Dispose();
        }

        private void ReDraw()
        {
            if (!doRedraw) return;

            StartDrawing();
            DoDrawing(true);
            StopDrawing();
        }

        private ILayer GetLayer()
        {
            return Layers?.FirstOrDefault();
        }

        private LineStringMapToolStateCommand GetMouseDownState(IFeature feature, TrackerFeature trackerFeature)
        {
            if (feature == null)
            {
                return selector.HasSelection
                    ? LineStringMapToolStateCommand.MouseDownNotOnFeatureHasSelection
                    : LineStringMapToolStateCommand.MouseDownNotOnFeatureNoSelection;
            }

            if (trackerFeature != null)
            {
                return LineStringMapToolStateCommand.MouseDownOnFeatureOnTracker;
            }

            return featureInteractor == null
                ? LineStringMapToolStateCommand.MouseDownOnOtherFeature
                : LineStringMapToolStateCommand.MouseDownOnSelectedFeature;
        }

        private ILineStringMapToolStrategy GetNextStrategyForCommand(LineStringMapToolStateCommand command)
        {
            var nextState = stateMachine.MoveNext(command);
            return lineEditorStrategies.FirstOrDefault(s => s.LineStringMapToolState == nextState) ?? currentStrategy;
        }

        private void CurrentControlOnPreviewKeyDown(object sender, PreviewKeyDownEventArgs e)
        {
            currentStrategy?.HandleKeyDown(IsCtrlPressed, e.KeyCode, featureInteractor);
        }
    }
}