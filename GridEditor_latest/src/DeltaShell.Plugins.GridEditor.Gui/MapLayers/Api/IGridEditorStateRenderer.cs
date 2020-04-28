using System;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Geometries;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api
{
    /// <summary>
    /// Renderer of a <see cref="GridEditorState"/>
    /// </summary>
    internal interface IGridEditorStateRenderer
    {
        /// <summary>
        /// <see cref="GridEditorState"/> to render
        /// </summary>
        GridEditorState GridEditorState { get; set; }

        /// <summary>
        /// Sub renderer for the <see cref="Data.GridEditorState.MeshGeometry"/>
        /// </summary>
        IDisposableMeshGeometryRenderer MeshGeometryRenderer { get; }

        /// <summary>
        /// Sub Renderer for rendering the <see cref="Data.GridEditorState.Polygons"/>
        /// </summary>
        IGridEditorPolygonRenderer PolygonRenderer { get; }

        /// <summary>
        /// Schedule a refresh action
        /// </summary>
        void Refresh();

        /// <summary>
        /// Spline function required to render splines when moving user geometry points
        /// </summary>
        void SetGetSplineGeometryFunction(Func<Coordinate[], Coordinate[]> getSplineGeometryFunction);
    }
}