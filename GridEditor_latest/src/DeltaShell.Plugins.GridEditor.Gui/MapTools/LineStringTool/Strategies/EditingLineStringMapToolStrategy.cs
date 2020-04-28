using System.Diagnostics.CodeAnalysis;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.Editors.Snapping;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal class EditingLineStringMapToolStrategy : ILineStringMapToolStrategy
    {
        private bool layerReadOnlySet;
        private bool previousLayerReadOnlyState;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public LineStringMapToolState LineStringMapToolState { get; } = LineStringMapToolState.Editing;

        /// <inheritdoc/>
        public Coordinate[] NewGeometryCoordinates { get; private set; }

        /// <inheritdoc/>
        public SnapResult SnapResult { get; private set; }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public bool UseLastSegmentStylePen { get; } = false;

        /// <inheritdoc/>
        public void HandleKeyUp(IFeatureInteractor featureInteractor)
        {
            if (!layerReadOnlySet) return;

            layerReadOnlySet = false;
            featureInteractor.Layer.ReadOnly = previousLayerReadOnlyState;
        }

        /// <inheritdoc/>
        public void HandleKeyDown(bool isCtrlPressed, Keys key, IFeatureInteractor featureInteractor)
        {
            if (key != Keys.Delete || SnapResult == null || SnapResult.SnapIndexPrevious != SnapResult.SnapIndexNext) return;

            layerReadOnlySet = true;
            previousLayerReadOnlyState = featureInteractor.Layer.ReadOnly;
            featureInteractor.Layer.ReadOnly = true; // prevent deletion of feature

            featureInteractor.Start();
            featureInteractor.RemoveTracker(featureInteractor.Trackers[SnapResult.SnapIndexPrevious]);
            featureInteractor.Stop();
        }

        /// <inheritdoc/>
        public void HandleMouseDown(Coordinate currentMouseCoordinate,
            IFeatureInteractor featureInteractor, IFeature feature, ISelector selector,
            TrackerFeature trackerFeature)
        {
            if (featureInteractor == null) return;
            if (SnapResult != null)
            {
                featureInteractor.Start();
                featureInteractor.InsertTracker(SnapResult.Location, SnapResult);
                featureInteractor.Stop();
            }

            featureInteractor.Layer.RenderRequired = true;
            NewGeometryCoordinates = featureInteractor?.SourceFeature?.Geometry?.Coordinates;
        }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleMouseUp(IFeatureInteractor featureInteractor) { /*No action*/  }

        /// <inheritdoc/>
        public void HandleMouseMove(Coordinate previousMouseMoveCoordinate,
            Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor)
        {
            if (featureInteractor == null) return;

            NewGeometryCoordinates = featureInteractor?.SourceFeature?.Geometry?.Coordinates;

            var maxDistance = 16 * featureInteractor.Layer.Map.PixelWidth;
            var snapResultTrackers = currentMouseCoordinate.GetTrackerSnapResult(maxDistance, featureInteractor.SourceFeature, featureInteractor.Layer);
            SnapResult = snapResultTrackers ?? currentMouseCoordinate.GetSnapResult(SnapRole.FreeAtObject, featureInteractor.SourceFeature, maxDistance, featureInteractor.Layer);
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
            NewGeometryCoordinates = null;
            SnapResult = null;
        }
    }
}