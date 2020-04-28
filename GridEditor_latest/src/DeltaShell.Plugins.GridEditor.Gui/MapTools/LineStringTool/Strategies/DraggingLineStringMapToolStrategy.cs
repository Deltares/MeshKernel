using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal class DraggingLineStringMapToolStrategy : ILineStringMapToolStrategy
    {
        private Coordinate previousMouseDownCoordinate;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public LineStringMapToolState LineStringMapToolState { get; } = LineStringMapToolState.Dragging;

        /// <inheritdoc/>
        public Coordinate[] NewGeometryCoordinates { get; private set; }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public SnapResult SnapResult { get; } = null;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public bool UseLastSegmentStylePen { get; } = false;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleKeyUp(IFeatureInteractor featureInteractor) { /*No action*/  }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleKeyDown(bool isCtrlPressed, Keys key, IFeatureInteractor featureInteractor) { /*No action*/  }

        /// <inheritdoc/>
        public void HandleMouseDown(Coordinate currentMouseCoordinate,
            IFeatureInteractor featureInteractor, IFeature feature, ISelector selector,
            TrackerFeature trackerFeature)
        {
            previousMouseDownCoordinate = currentMouseCoordinate;

            featureInteractor?.Start();
            featureInteractor?.SetTrackerSelection(trackerFeature, true);
        }

        /// <inheritdoc/>
        public void HandleMouseUp(IFeatureInteractor featureInteractor)
        {
            if (featureInteractor==null) return;
            featureInteractor.Stop();
            featureInteractor.Layer.RenderRequired = true;
        }

        /// <inheritdoc/>
        public void HandleMouseMove(Coordinate previousMouseMoveCoordinate,
            Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor)
        {
            if (featureInteractor == null) return;
            var previousCoordinate = previousMouseMoveCoordinate ?? previousMouseDownCoordinate;

            var deltaX = currentMouseCoordinate.X - previousCoordinate.X;
            var deltaY = currentMouseCoordinate.Y - previousCoordinate.Y;

            foreach (var trackerFeature in featureInteractor.Trackers.Where(t => t.Selected))
            {
                featureInteractor.MoveTracker(trackerFeature, deltaX, deltaY);
            }

            NewGeometryCoordinates = featureInteractor?.TargetFeature?.Geometry?.Coordinates;
        }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleMouseDoubleClick(Coordinate currentMouseCoordinate, ILayer layer,
            GridEditorState gridEditorState)
        {
            // No action required
        }

        /// <inheritdoc/>
        public void ClearState()
        {
            previousMouseDownCoordinate = null;
            NewGeometryCoordinates = null;
        }
    }
}