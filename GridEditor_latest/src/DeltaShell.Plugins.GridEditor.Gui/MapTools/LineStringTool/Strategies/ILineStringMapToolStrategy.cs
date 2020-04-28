using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal interface ILineStringMapToolStrategy
    {
        /// <summary>
        /// <see cref="LineStringMapToolState"/> to use this strategy
        /// </summary>
        LineStringMapToolState LineStringMapToolState { get; }

        /// <summary>
        /// New geometry coordinates to render
        /// </summary>
        Coordinate[] NewGeometryCoordinates { get; }

        /// <summary>
        /// Snap result to render
        /// </summary>
        SnapResult SnapResult { get; }

        /// <summary>
        /// Render <see cref="NewGeometryCoordinates"/> with the last segment in a different style
        /// </summary>
        bool UseLastSegmentStylePen { get; }

        /// <summary>
        /// Handler for key up event
        /// </summary>
        /// <param name="featureInteractor"></param>
        void HandleKeyUp(IFeatureInteractor featureInteractor);

        /// <summary>
        /// Handler for a key combination
        /// </summary>
        /// <param name="isCtrlPressed">True if control is pressed</param>
        /// <param name="key"><see cref="Keys"/> that was pressed</param>
        /// <param name="featureInteractor"><see cref="IFeatureInteractor"/> for the currently edited feature</param>
        void HandleKeyDown(bool isCtrlPressed, Keys key, IFeatureInteractor featureInteractor);

        /// <summary>
        /// Handler for mouse down
        /// </summary>
        /// <param name="currentMouseCoordinate">The current mouse down coordinate</param>
        /// <param name="featureInteractor">Current feature interactor for the feature that was at the <see cref="currentMouseCoordinate"/></param>
        /// <param name="feature">Current feature that was at the <see cref="currentMouseCoordinate"/></param>
        /// <param name="selector">Select tool</param>
        /// <param name="trackerFeature">Tracker feature that was at <see cref="currentMouseCoordinate"/></param>
        void HandleMouseDown(Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor, IFeature feature, ISelector selector, TrackerFeature trackerFeature);

        /// <summary>
        /// Handler for mouse up
        /// </summary>
        /// <param name="featureInteractor"><see cref="IFeatureInteractor"/> for the currently edited feature</param>
        void HandleMouseUp(IFeatureInteractor featureInteractor);

        /// <summary>
        /// Handler for mouse move
        /// </summary>
        /// <param name="previousMouseMoveCoordinate">Previous mouse move coordinate</param>
        /// <param name="currentMouseCoordinate">Current mouse move coordinate</param>
        /// <param name="featureInteractor"><see cref="IFeatureInteractor"/> for the currently edited feature</param>
        void HandleMouseMove(Coordinate previousMouseMoveCoordinate, Coordinate currentMouseCoordinate, IFeatureInteractor featureInteractor);

        /// <summary>
        /// Handler for mouse double click
        /// </summary>
        /// <param name="currentMouseCoordinate"></param>
        /// <param name="layer"></param>
        /// <param name="gridEditorState"></param>
        void HandleMouseDoubleClick(Coordinate currentMouseCoordinate, ILayer layer, GridEditorState gridEditorState);

        /// <inheritdoc/>
        /// <summary>
        /// Clears map tool data
        /// </summary>
        void ClearState();
    }
}