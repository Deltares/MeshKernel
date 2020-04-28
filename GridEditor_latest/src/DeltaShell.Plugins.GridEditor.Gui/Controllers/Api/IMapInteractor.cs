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
        Func<Coordinate, double, int> GetVertexIndex { get; set; }

        /// <summary>
        /// Function to delete a vertex at a given index
        /// </summary>
        Func<int, bool> DeleteVertex { get; set; }

        /// <summary>
        /// Function to merge the mesh vertices 
        /// </summary>
        Func<bool> MergeVertices { get; set; }

        /// <summary>
        /// Function to delete an edge within a search radius
        /// </summary>
        Func<Coordinate, double, bool> DeleteEdges { get; set; }

        /// <summary>
        /// Function to move a selected vertex to a new position
        /// </summary>
        Func<Coordinate, int, bool> MoveVertex { get; set; }

        /// <summary>
        /// List of currently selected polygons 
        /// </summary>
        IEnumerable<IPolygonSelection> PolygonSelections { get; }

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

        /// <summary>
        /// Zooms the map to the samples extent
        /// </summary>
        void ZoomToSamplesExtent();

        /// <summary>
        /// Zooms the map to spline extent
        /// </summary>
        void ZoomToSplinesExtent();

        /// <summary>
        /// Zooms the map to land boundaries extent
        /// </summary>
        void ZoomToLandBoundariesExtent();
    }
}