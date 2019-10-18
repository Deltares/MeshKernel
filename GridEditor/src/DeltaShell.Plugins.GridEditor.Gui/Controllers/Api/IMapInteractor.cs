using System;
using System.Collections.Generic;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Grids;

namespace DeltaShell.Plugins.GridEditor.Gui.Controllers.Api
{
    internal interface IMapInteractor 
    {
        /// <summary>
        /// Renderer of data
        /// </summary>
        IGridEditorStateRenderer Renderer { get; }

        /// <summary>
        /// Currently selected <see cref="MapToolType"/>
        /// </summary>
        MapToolType SelectedMapToolType { get; }

        /// <summary>
        /// Map tool for changing a polygon 
        /// </summary>
        IChangePolygonMapTool ChangePolygonMapTool { get; }

        /// <summary>
        /// Function to calculate the spline points
        /// </summary>
        Func<Coordinate[], Coordinate[]> GetSplineGeometry { get; set; }

        /// <summary>
        /// Function to insert a vertex in the current grid
        /// </summary>
        Func<Coordinate, int> InsertVertex { get; set; }

        /// <summary>
        /// Function to insert an edge in the current grid
        /// </summary>
        Func<int, int, int> InsertEdge { get; set; }

        /// <summary>
        /// Function to get the index of the closest vertex
        /// </summary>
        Func<double, double, double, double, int> GetVertexIndex { get; set; }

        /// <summary>
        /// Gets all <see cref="UnstructuredGrid"/>s from the map
        /// </summary>
        /// <returns></returns>
        IEnumerable<UnstructuredGrid> GetUnstructuredGrids();

        /// <summary>
        /// Starts the map interaction
        /// </summary>
        /// <param name="gridEditorState"></param>
        void StartInteraction(GridEditorState gridEditorState);

        /// <summary>
        /// Stops the map interaction
        /// </summary>
        void StopInteraction();

        /// <summary>
        /// Activates a certain map tool (<see cref="MapToolType"/>)
        /// </summary>
        /// <param name="mapToolType"></param>
        void ActivateTool(MapToolType mapToolType);

        /// <summary>
        /// Zooms the map to the extent of the grid
        /// </summary>
        void ZoomToGridExtent();
    }
}