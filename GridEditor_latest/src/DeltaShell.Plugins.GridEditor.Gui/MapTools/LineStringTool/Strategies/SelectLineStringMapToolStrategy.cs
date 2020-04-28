using System.Diagnostics.CodeAnalysis;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal class SelectLineStringMapToolStrategy : ILineStringMapToolStrategy
    {
        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public LineStringMapToolState LineStringMapToolState { get; } = LineStringMapToolState.Selecting;

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
        public void HandleMouseDown(Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor,
            IFeature feature, ISelector selector, TrackerFeature trackerFeature)
        {
            if (selector.HasSelection)
            {
                selector.ClearSelection();
            }

            if (featureInteractor == null)
            {
                selector.Select(feature);
                NewGeometryCoordinates = feature?.Geometry?.Coordinates;
            }
        }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleMouseUp(IFeatureInteractor featureInteractor) { /*No action*/  }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleMouseMove(Coordinate previousMouseMoveCoordinate,
            Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor) { /*No action*/  }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleMouseDoubleClick(Coordinate currentMouseCoordinate, ILayer layer,
            GridEditorState gridEditorState){ /*No action*/  }

        /// <inheritdoc/>
        public void ClearState()
        {
            NewGeometryCoordinates = null;
        }
    }
}