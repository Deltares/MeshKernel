using System;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using GeoAPI.Geometries;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api
{
    public interface IDisposableMeshGeometryRenderer
    {
        /// <summary>
        /// Draw grid vertices
        /// </summary>
        bool ShowGridVertices { get; set; }

        /// <summary>
        /// Draw size of points
        /// </summary>
        float PointSize { get; set; }
        
        /// <summary>
        /// Draw selected vertices
        /// </summary>
        bool DrawSelectedVertices { get; set; }

        /// <summary>
        /// Draw edge values
        /// </summary>
        bool DrawEdgeValues { get; set; }

        /// <summary>
        /// Gets the envelope of the mesh
        /// </summary>
        Envelope Envelope { get; }

        /// <summary>
        /// Function for getting the selected vertices
        /// </summary>
        Func<int[]> GetSelectedVertices { get; set; }

        /// <summary>
        /// Function for getting edge values
        /// </summary>
        Func<DisposableGeometryList> GetEdgeValues { get; set; }

        /// <summary>
        /// Signal renderer that the mesh has changed
        /// </summary>
        void MeshChanged();
    }
}