using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Windows.Forms;
using DelftTools.Utils.Collections;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.CoordinateSystems.Transformations;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal class AddNewLineStringMapToolStrategy : ILineStringMapToolStrategy
    {
        private readonly IList<Coordinate> coordinates = new List<Coordinate>();
        private Coordinate previousMouseDownCoordinate;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public LineStringMapToolState LineStringMapToolState { get; } = LineStringMapToolState.AddingNew;

        /// <inheritdoc/>
        public Coordinate[] NewGeometryCoordinates { get; private set; }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage] // not used in this strategy
        public SnapResult SnapResult { get; } = null;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public bool UseLastSegmentStylePen { get; } = true;

        /// <summary>
        /// Action to commit the generated geometry to the <see cref="GridEditorState"/>
        /// </summary>
        public Action<IGeometry, GridEditorState> CommitGeometry { get; set; }

        /// <summary>
        /// Function for transforming coordinates
        /// </summary>
        public Func<Coordinate[], Coordinate[]> TransformPreviewCoordinates { get; set; } = c => c;

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleKeyUp(IFeatureInteractor featureInteractor) { /*No action*/ }

        /// <inheritdoc/>
        public void HandleKeyDown(bool isCtrlPressed, Keys key, IFeatureInteractor featureInteractor)
        {
            if (isCtrlPressed && key == Keys.Z && coordinates.Count > 0)
            {
                coordinates.RemoveAt(coordinates.Count - 1);
                NewGeometryCoordinates = TransformPreviewCoordinates(coordinates.ToArray());
            }

            if (isCtrlPressed && key == Keys.U)
            {
                coordinates.Clear();
            }
        }

        /// <inheritdoc/>
        public void HandleMouseDown(Coordinate currentMouseCoordinate,
            IFeatureInteractor featureInteractor, IFeature feature, ISelector selector,
            TrackerFeature trackerFeature)
        {
            if (Equals(currentMouseCoordinate, previousMouseDownCoordinate)) return;

            previousMouseDownCoordinate = currentMouseCoordinate;

            coordinates.Add(currentMouseCoordinate);
            NewGeometryCoordinates = TransformPreviewCoordinates(coordinates.ToArray());
        }

        /// <inheritdoc/>
        [ExcludeFromCodeCoverage]
        public void HandleMouseUp(IFeatureInteractor featureInteractor) { /*No action*/  }

        /// <inheritdoc/>
        public void HandleMouseMove(Coordinate previousMouseMoveCoordinate,
            Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor)
        {
            NewGeometryCoordinates = TransformPreviewCoordinates(coordinates.Plus(currentMouseCoordinate).ToArray());
        }

        /// <inheritdoc/>
        public void HandleMouseDoubleClick(Coordinate currentMouseCoordinate, ILayer layer, GridEditorState gridEditorState)
        {
            var coordinatesToUse = coordinates.Plus(currentMouseCoordinate).ToList().RemoveDuplicateCoordinates().ToArray();

            if (layer?.CoordinateTransformation != null)
            {
                coordinatesToUse = GeometryTransform.TransformLineString(new LineString(coordinatesToUse), layer.CoordinateTransformation.MathTransform.Inverse()).Coordinates;
            }

            CommitCoordinates(coordinatesToUse, gridEditorState);

            coordinates.Clear();
            NewGeometryCoordinates = null;
        }

        public void ClearState()
        {
            coordinates.Clear();
            NewGeometryCoordinates = null;
        }

        protected virtual void CommitCoordinates(Coordinate[] lineCoordinates, GridEditorState gridEditorState)
        {
            if (lineCoordinates.Length < 2) return;

            CommitGeometry?.Invoke(new LineString(lineCoordinates), gridEditorState);
        }
    }
}