using System;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    public interface IGridgeomStatefulApi : IDisposable
    {
        /// <summary>
        /// Create a new grid state and return the generated gridStateId/>
        /// </summary>
        /// <returns>Generated gridStateId</returns>
        int CreateGridState();

        /// <summary>
        /// Deallocate grid state (collections of mesh arrays with auxiliary variables)
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <returns>If the operation succeeded</returns>
        bool RemoveGridState(int gridStateId);

        /// <summary>
        /// Synchronize provided grid (<param name="disposableMeshGeometry"/>) with the grid state with <param name="gridStateId"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableMeshGeometry">Grid state as <see cref="DisposableMeshGeometry"/> object</param>
        /// <param name="isGeographic">Coordinate reference system type (cartesian or sferic)</param>
        /// <returns>If the operation succeeded</returns>
        bool SetGridState(int gridStateId, DisposableMeshGeometry disposableMeshGeometry, bool isGeographic);

        /// <summary>
        /// Gets the grid state as a <see cref="MeshGeometry"/> structure excluding the cell information
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <returns><see cref="DisposableMeshGeometry"/> with the grid state</returns>
        DisposableMeshGeometry GetGridState(int gridStateId);

        /// <summary>
        /// Gets the grid state as a <see cref="MeshGeometry"/> structure including the cell information
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <returns><see cref="DisposableMeshGeometry"/> with the grid state</returns>
        DisposableMeshGeometry GetGridStateWithCells(int gridStateId);

        /// <summary>
        /// Deletes a node with specified <param name="vertexIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="vertexIndex">The index of the node to delete</param>
        /// <returns>If the operation succeeded</returns>
        bool DeleteVertex(int gridStateId, int vertexIndex);

        /// <summary>
        /// Flips the links
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account or not for existing land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <returns>If the operation succeeded</returns>
        bool FlipEdges(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired, ProjectToLandBoundaryOptions projectToLandBoundaryOption);

        /// <summary>
        /// Merges vertex <param name="startVertexIndex"/> to <param name="endVertexIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="startVertexIndex">The index of the first vertex to merge</param>
        /// <param name="endVertexIndex">The index of the second vertex to merge</param>
        /// <returns>If the operation succeeded</returns>
        bool MergeTwoVertices(int gridStateId, int startVertexIndex, int endVertexIndex);

        /// <summary>
        /// Merges vertices within a distance of 0.001 m, effectively removing small edges 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The polygon where to perform the operation</param>
        /// <returns>If the operation succeeded</returns>
        bool MergeVertices(int gridStateId, DisposableGeometryList disposableGeometryList);

        /// <summary>
        /// Orthogonalization
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells (if activated squares becomes triangles) </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account for land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <param name="orthogonalizationParameters">The structure containing the user defined orthogonalization parameters</param>
        /// <param name="geometryListNativePolygon">The polygon where to perform the orthogonalization</param>
        /// <param name="geometryListNativeLandBoundaries">The land boundaries to account for in the orthogonalization process</param>
        /// <returns>Error code</returns>
        bool Orthogonalize(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired,
            ProjectToLandBoundaryOptions projectToLandBoundaryOption, OrthogonalizationParameters orthogonalizationParameters, 
            DisposableGeometryList geometryListNativePolygon, DisposableGeometryList geometryListNativeLandBoundaries);

        /// <summary>
        /// Make a new grid
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="makeGridParameters">The structure containing the make grid parameters </param>
        /// <param name="disposableGeometryListIn"></param>
        /// <returns>Error code</returns>
        bool MakeGrid(int gridStateId, MakeGridParameters makeGridParameters, ref DisposableGeometryList disposableGeometryListIn);

        /// <summary>
        /// Get spline intermediate points 
        /// </summary>
        /// <param name="disposableGeometryListIn">The input corner vertices of the splines</param>
        /// <param name="disposableGeometryListOut">The output spline </param>
        /// <param name="numberOfPointsBetweenVertices">The number of spline vertices between the corners points</param>
        /// <returns>Error code</returns>
        bool GetSplines(DisposableGeometryList disposableGeometryListIn, ref DisposableGeometryList disposableGeometryListOut, int numberOfPointsBetweenVertices);

        /// <summary>
        /// Make curvilinear grid from splines 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The corner vertices of the splines</param>
        /// <param name="curvilinearParameters">The parameters for the generation of the curvilinear grid</param>
        /// <returns>Error code</returns>
        bool MakeGridFromSplines(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, ref CurvilinearParameters curvilinearParameters);

        /// <summary>
        /// Make curvilinear grid from splines with advancing front.
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The corner vertices of the splines</param>
        /// <param name="curvilinearParameters">The parameters for the generation of the curvilinear grid</param>
        /// <param name="splinesToCurvilinearParameters">The parameters of the advancing front algorithm</param>
        /// <returns></returns>
        bool MakeOrthogonalGridFromSplines(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,  ref CurvilinearParameters curvilinearParameters, ref SplinesToCurvilinearParameters splinesToCurvilinearParameters);

        bool MakeOrthogonalGridFromSplinesInitialize(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, ref CurvilinearParameters curvilinearParameters, ref SplinesToCurvilinearParameters splinesToCurvilinearParameters);
        bool MakeOrthogonalGridFromSplinesIteration(int gridStateId,  int layer);
        bool MakeOrthogonalGridFromSplinesRefreshMesh(int gridStateId);
        bool MakeOrthogonalGridFromSplinesDelete(int gridStateId);

        /// <summary>
        /// Make a triangular grid in a polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The polygon where to triangulate</param>
        /// <returns></returns>
        bool MakeTriangularGridInPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryList);


        /// <summary>
        /// Make a triangular grid from samples
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The samples where to triangulate</param>
        /// <returns></returns>
        bool MakeTriangularGridFromSamples(int gridStateId, ref DisposableGeometryList disposableGeometryList);

        /// <summary>
        /// Counts the number of polygon vertices contained in the mesh boundary polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="numberOfPolygonVertices">The number of polygon points</param>
        /// <returns></returns>
        bool CountMeshBoundaryPolygonVertices(int gridStateId, ref int numberOfPolygonVertices);

        /// <summary>
        /// Retrives the mesh boundary polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryList">The output network boundary polygon</param>
        /// <returns></returns>
        bool GetMeshBoundaryPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryList);

        /// <summary>
        /// Get the number of vertices of the offsetted polygon 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The input polygon</param>
        /// <param name="distance">The offset distance</param>
        /// <param name="innerPolygon">Compute inner polygon or not</param>
        /// <param name="numberOfPolygonVertices">The number of vertices of the offsetted polygon</param>
        /// <returns></returns>
        bool CountVerticesOffsettedPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, bool innerPolygon, double distance, ref int numberOfPolygonVertices);

        /// <summary>
        /// Get the offsetted polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The input polygon</param>
        /// <param name="distance">The offset distance</param>
        /// <param name="innerPolygon">Compute inner polygon or not</param>
        /// <param name="disposableGeometryListOut">The offsetted polygon</param>
        /// <returns></returns>
        bool GetOffsettedPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, bool innerPolygon, double distance, ref DisposableGeometryList disposableGeometryListOut);

        /// <summary>
        /// Count the number of polygon vertices after refinement
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The input polygon</param>
        /// <param name="firstIndex">The index of the first vertex</param>
        /// <param name="secondIndex">The index of the second vertex</param>
        /// <param name="distance">The refinement distance</param>
        /// <param name="numberOfPolygonVertices">The number of vertices after refinement </param>
        /// <returns></returns>
        bool CountVerticesRefinededPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, int firstIndex,  int secondIndex,  double distance,  ref int numberOfPolygonVertices);

        /// <summary>
        /// Gets the refined polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The input polygon</param>
        /// <param name="firstIndex">The index of the first vertex</param>
        /// <param name="secondIndex">The index of the second vertex</param>
        /// <param name="distance">The refinement distance</param>
        /// <param name="disposableGeometryListOut">The refined polygon</param>
        /// <returns></returns>
        bool GetRefinededPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, int firstIndex, int secondIndex, double distance, ref DisposableGeometryList disposableGeometryListOut);


        /// <summary>
        /// Refines a grid based on samples
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The input samples</param>
        /// <param name="interpolationParameters">The settings for the interpolation algorithm</param>
        /// <param name="samplesRefineParameters">The settings for the interpolation related to samples</param>
        /// <returns></returns>
        bool RefineGridBasedOnSamples(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, InterpolationParameters interpolationParameters, SamplesRefineParameters samplesRefineParameters);


        /// <summary>
        /// Refines a grid based on polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn">The closed polygon where to perform the refinement</param>
        /// <param name="interpolationParameters">The settings for the interpolation algorithm</param>
        /// <returns></returns>
        bool RefineGridBasedOnPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, InterpolationParameters interpolationParameters);

        /// <summary>
        /// Returns the vertices indexes inside selected polygons
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListIn"></param>
        /// <param name="inside"> Select inside (0) or outside (1) polygon</param>
        /// <returns></returns>
        int[] GetSelectedVertices(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, int inside);
        
        /// <summary>
        /// Get the edges orthogonality
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListOut">In ZCoordinates field the orthogonality values are stored</param>
        /// <returns></returns>
        bool GetOrthogonality(int gridStateId, ref DisposableGeometryList disposableGeometryListOut);

        /// <summary>
        /// Get the edges smoothness
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryListOut">In ZCoordinates field the smoothness values are stored</param>
        /// <returns></returns>
        bool GetSmoothness(int gridStateId, ref DisposableGeometryList disposableGeometryListOut);

        /// <summary>
        /// Inserts a new vertex
        /// </summary>
        /// <param name="gridGeomId">Id of the grid state</param>
        /// <param name="xCoordinate">x coordinate of the vertex</param>
        /// <param name="yCoordinate">y coordinate of the vertex</param>
        /// <param name="zCoordinate">z coordinate of the vertex</param>
        /// <param name="vertexIndex">the index of the new vertex</param>
        /// <returns></returns>
        bool InsertVertex(int gridGeomId, double xCoordinate, double yCoordinate, double zCoordinate, ref int vertexIndex);

        /// <summary>
        /// Insert a new edge
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="startVertexIndex">The index of the first vertex to connect</param>
        /// <param name="endVertexIndex">The index of the second vertex to connect</param>
        /// <param name="edgeIndex">The index of the new edge</param>
        /// <returns>If the operation succeeded</returns>
        bool InsertEdge(int gridStateId, int startVertexIndex, int endVertexIndex, ref int edgeIndex);

        /// <summary>
        /// Get the index of the closest existing vertex
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">Vertex coordinates</param>
        /// <param name="searchRadius">the radius where to search for the vertex</param>
        /// <param name="vertexIndex">the index of the closest vertex</param>
        /// <returns></returns>
        bool GetVertexIndex(int gridStateId, ref DisposableGeometryList geometryListIn, double searchRadius, ref int vertexIndex);


        /// <summary>
        /// Deletes the closest mesh edge within the search radius from the input point
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input point coordinates</param>
        /// <param name="searchRadius">The search radius</param>
        /// <returns> true if the edge has been deleted, false if not (the edge is outside the search radius) </returns>
        bool DeleteEdge(int gridStateId, ref DisposableGeometryList geometryListIn, double searchRadius);

        /// <summary>
        /// Deletes a mesh in a polygon using several options
        /// </summary>
        /// <param name="gridGeomId">Id of the grid state</param>
        /// <param name="disposableGeometryListOut">The polygon where to perform the operation</param>
        /// <param name="deletionOption">The deletion option (to be detailed)</param>
        /// <returns>If the operation succeeded</returns>
        bool DeleteMeshWithOptions(int gridGeomId, ref DisposableGeometryList disposableGeometryListOut, ref int deletionOption);

        /// <summary>
        /// Function to move a selected vertex to a new position
        /// </summary>
        /// <param name="gridGeomId">Id of the grid state</param>
        /// <param name="disposableGeometryLisIn">The new coordinate</param>
        /// <param name="vertexIndex">The vertex index (to be detailed)</param>
        /// <returns>If the operation succeeded</returns>
        bool MoveVertex(int gridGeomId, ref DisposableGeometryList disposableGeometryLisIn, int vertexIndex);

        /// <summary>
        /// Selects points in polygons
        /// </summary>
        /// <param name="gridGeomId">Id of the grid state</param>
        /// <param name="inputPolygon">The polygon(s) used for selection</param>
        /// <param name="inputPoints">The points to select</param>
        /// <param name="selectedPoints">The selected points in the zCoordinates field (0.0 not selected, 1.0 selected)</param>
        /// <returns>If the operation succeeded</returns>
        bool PointsInPolygon(int gridGeomId, ref DisposableGeometryList inputPolygon, ref DisposableGeometryList inputPoints, ref DisposableGeometryList selectedPoints);


        bool OrthogonalizationInitialize(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired,
            ProjectToLandBoundaryOptions projectToLandBoundaryOption,
            OrthogonalizationParameters orthogonalizationParameters, DisposableGeometryList geometryListNativePolygon,
            DisposableGeometryList geometryListNativeLandBoundaries);

        bool OrthogonalizationPrepareOuterIteration(int gridGeomId);

        bool OrthogonalizationInnerIteration(int gridGeomId);

        bool OrthogonalizationFinalizeOuterIteration(int gridGeomId);

        bool OrthogonalizationDelete(int gridGeomId);
    }
}