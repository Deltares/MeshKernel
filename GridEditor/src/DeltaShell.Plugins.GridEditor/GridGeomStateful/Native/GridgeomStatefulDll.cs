using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Runtime.InteropServices;
using DelftTools.Utils.Interop;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    // never hit by code coverage because tests use remoting and this only contains dll imports
    [ExcludeFromCodeCoverage]
    internal static class GridgeomStatefulDll
    {

        public struct TestStructure
        {
            public IntPtr array;
        }

        private const string GridGeomStatefulDllName = "gridgeomStateful_dll.dll";

        static GridgeomStatefulDll()
        {
            var dir = Path.GetDirectoryName(typeof(GridgeomStatefulDll).Assembly.Location);
            NativeLibrary.LoadNativeDll(GridGeomStatefulDllName, Path.Combine(dir, "Lib"));
        }

        #region State management

        /// <summary>
        /// Create a new grid state and return the generated <param name="gridStateId"/>
        /// </summary>
        /// <param name="gridStateId">Identifier for the created grid state</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_new_grid", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int CreateGridState([In, Out] ref int gridStateId);

        /// <summary>
        /// Deallocate grid state (collections of mesh arrays with auxiliary variables)
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_deallocate_state", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int RemoveGridState([In] ref int gridStateId);

        /// <summary>
        /// Gets the grid state as a <see cref="MeshGeometry"/> structure
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="meshGeometryDimensions">Grid dimensions</param>
        /// <param name="meshGeometry">Grid data</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_get_mesh", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetMeshState([In] ref int gridStateId, [In, Out] ref MeshGeometryDimensions meshGeometryDimensions, [In, Out] ref MeshGeometry meshGeometry);

        /// <summary>
        /// Synchronize provided grid (<param name="meshGeometryDimensions"/> and <param name="meshGeometry"/>) with the grid state with <param name="gridStateId"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="meshGeometryDimensions">Grid dimensions</param>
        /// <param name="meshGeometry">Grid data</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_set_state", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int SetState([In] ref int gridStateId, [In] ref MeshGeometryDimensions meshGeometryDimensions, [In] ref MeshGeometry meshGeometry, [In] ref bool IsGeographic);

        #endregion

        #region Node operations

        /// <summary>
        /// Deletes a node with specified <param name="nodeIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="nodeIndex">The nodeIndex to delete</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_delete_node", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int DeleteNode([In] ref int gridStateId, [In] ref int nodeIndex);

        #endregion

        #region Links operations  

        /// <summary>
        /// Flip the edges
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells (if activated squares becomes triangles) </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account for land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_flip_links", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int FlipEdges([In] ref int gridStateId, [In] ref int isTriangulationRequired, [In] ref int isAccountingForLandBoundariesRequired, [In] ref int projectToLandBoundaryOption);

        /// <summary>
        /// Insert a new edge connecting <param name="startVertexIndex"/> and <param name="endVertexIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="startVertexIndex">The index of the first vertex to connect</param>
        /// <param name="endVertexIndex">The index of the second vertex to connect</param>
        /// <param name="newEdgeIndex">The index of the new edge</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_insert_edge", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int InsertEdge([In] ref int gridStateId, [In] ref int startVertexIndex, [In] ref int endVertexIndex, [In, Out] ref int edgeIndex);
        #endregion

        #region Vertices operations  
        /// <summary>
        /// Merges vertex <param name="startVertexIndex"/> to <param name="endVertexIndex"/>
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="startVertexIndex">The index of the first vertex to merge</param>
        /// <param name="endVertexIndex">The index of the second vertex to merge</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_merge_two_nodes", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MergeTwoVertices([In] ref int gridStateId, [In] ref int startVertexIndex, [In] ref int endVertexIndex);

        /// <summary>
        /// Merges vertices within a distance of 0.001 m, effectively removing small edges 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The polygon where to perform the operation</param>
        /// <returns>If the operation succeeded</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_merge_nodes", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MergeVertices([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative);

        /// <summary>
        /// Inserts a new vertex
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The polygon where to perform the operation</param>
        /// <param name="vertexIndex">The index of the new vertex</param>
        /// <returns>If the operation succeeded</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_insert_node", CallingConvention = CallingConvention.Cdecl)]
        public static extern int InsertVertex([In] ref int gridStateId, [In] ref double xCoordinate, [In] ref double yCoordinate, [In] ref double zCoordinate, [In, Out] ref int vertexIndex);


        /// <summary>
        /// Deletes a mesh in a polygon using several options
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="disposableGeometryList">The polygon where to perform the operation</param>
        /// <param name="deletionOption">The deletion option (to be detailed)</param>
        /// <returns>If the operation succeeded</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_delete_mesh", CallingConvention = CallingConvention.Cdecl)]
        public static extern int DeleteMeshWithOptions([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative, [In] ref int deletionOption);


        #endregion

        #region Orthogonalization
        /// <summary>
        /// Orthogonalization
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="isTriangulationRequired">The option to triangulate also non triangular cells (if activated squares becomes triangles) </param>
        /// <param name="isAccountingForLandBoundariesRequired">The option to account for land boundaries</param>
        /// <param name="projectToLandBoundaryOption">The option to determine how to snap to land boundaries</param>
        /// <param name="orthogonalizationParametersNative">The structure containing the orthogonalization parameters</param>
        /// <param name="geometryListNativePolygon">The polygon where to perform the orthogonalization</param>
        /// <param name="geometryListNativeLandBoundaries">The land boundaries to account for in the orthogonalization process</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_orthogonalize", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int Orthogonalization([In] ref int gridStateId, [In] ref int isTriangulationRequired, [In] ref int isAccountingForLandBoundariesRequired, [In] ref int projectToLandBoundaryOption, [In] ref OrthogonalizationParametersNative orthogonalizationParametersNative, [In] ref GeometryListNative geometryListNativePolygon, [In] ref GeometryListNative geometryListNativeLandBoundaries);

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_orthogonalize_initialize", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int OrthogonalizationInitialize([In] ref int gridStateId);

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_orthogonalize_prepare_outer_iteration", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int OrthogonalizationPrepareOuterIteration([In] ref int gridStateId);

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_orthogonalize_inner_iteration", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int OrthogonalizationInnerIteration([In] ref int gridStateId);

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_orthogonalize_finalize_outer_iteration", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int OrthogonalizationFinalizeOuterIteration([In] ref int gridStateId);

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_orthogonalize_delete", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int OrthogonalizationDelete([In] ref int gridStateId);

        #endregion

        #region Make grid
        /// <summary>
        /// Make a new grid
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="makeGridParameters">The structure containing the make grid parameters </param>
        /// <param name="geometryListNative">The polygon to account for</param>
        /// <returns>Error code</returns>

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_make_net", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MakeGrid([In] ref int gridStateId, [In] ref MakeGridParametersNative makeGridParameters, [In] ref GeometryListNative geometryListNative);

        /// <summary>
        /// Make a triangular grid in a polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListNative">The polygon where to triangulate</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_mesh_from_polygon", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MakeTriangularGridFromPolygon([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative);


        /// <summary>
        /// Make a triangular grid from samples
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListNative">The samples where to triangulate</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_mesh_from_samples", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MakeTriangularGridFromSamples([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative);

        #endregion

        #region Curvilinear grids

        /// <summary>
        /// Get spline intermediate points 
        /// </summary>
        /// <param name="disposableGeometryListIn">The input corner vertices of the splines</param>
        /// <param name="disposableGeometryListOut">The output spline </param>
        /// <param name="numberOfPointsBetweenVertices">The number of spline vertices between the corners points</param>
        /// <returns>Error code</returns>

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_get_splines", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetSplines([In] ref GeometryListNative geometryListNativeIn, [In, Out] ref GeometryListNative geometryListNativeOut, [In] ref int numberOfPointsBetweenVertices);

        /// <summary>
        /// Get spline intermediate points 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListNativeIn">The input splines corners</param>
        /// <param name="curvilinearParametersNative">The input parameters to generate the curvilinear grid</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_curvilinear_mesh_from_splines", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MakeGridFromSplines([In] ref int gridStateId, [In] ref GeometryListNative geometryListNativeIn, [In] ref CurvilinearParametersNative curvilinearParametersNative);

        /// <summary>
        /// Make curvilinear grid from splines with an advancing front.
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListNative">The input splines corners</param>
        /// <param name="curvilinearParametersNative">The input parameters to generate the curvilinear grid</param> 
        /// <param name="splinesToCurvilinearParametersNative">The parameters of the advancing front algorithm</param>
        /// <returns>Error code</returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_curvilinear_mesh_from_splines_ortho", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int MakeOrthogonalGridFromSplines([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative, [In] ref CurvilinearParametersNative curvilinearParametersNative, [In] ref SplinesToCurvilinearParametersNative splinesToCurvilinearParametersNative);
        #endregion

        #region Mesh operations
        /// <summary>
        /// Counts the number of polygon vertices contained in the mesh boundary polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="numberOfPolygonVertices">The number of polygon points</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_copy_mesh_boundaries_to_polygon_count_edges", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int CountMeshBoundaryPolygonVertices([In] ref int gridStateId, [In, Out] ref int numberOfPolygonVertices);


        /// <summary>
        /// Retrives the network boundary polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryList">The output network boundary polygon</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_copy_mesh_boundaries_to_polygon", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetMeshBoundaryPolygon([In] ref int gridStateId, [In, Out] ref GeometryListNative geometryListNative);

        /// <summary>
        /// Refine a grid based on the samples contained in the geometry list
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListNative">The sample set</param>
        /// <param name="interpolationParametersNative">The interpolation parameters</param>
        /// <param name="sampleRefineParametersNative">The interpolation settings related to the samples</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_refine_mesh_based_on_samples", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int RefineGridBasedOnSamples([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative, [In] ref InterpolationParametersNative interpolationParametersNative, [In] ref SampleRefineParametersNative sampleRefineParametersNative);

        /// <summary>
        /// Refine a grid based on polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListNative">The closed polygon where to perform the refinement</param>
        /// <param name="interpolationParametersNative">The interpolation parameters</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_refine_mesh_based_on_polygon", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int RefineGridBasedOnPolygon([In] ref int gridStateId, [In] ref GeometryListNative geometryListNative, [In] ref InterpolationParametersNative interpolationParametersNative);

        #endregion

        #region Polygon operations

        /// <summary>
        /// Get the number of vertices of the offsetted polygon 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input polygon</param>
        /// <param name="xCoordOffsetPoint">The x coordinate of the offset point</param>
        /// <param name="yCoordOffsetPoint">The y coordinate of the offset point</param>
        /// <param name="innerPolygon">Compute inner polygon or not</param>
        /// <param name="numberOfPolygonVertices">The number of vertices of the offsetted polygon</param>
        /// <returns></returns>

        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_offsetted_polygon_count", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int CountVerticesOffsettedPolygon([In] ref int gridStateId, [In] ref GeometryListNative geometryListIn, [In] ref bool innerPolygon, [In] ref double distance, [In, Out] ref int numberOfPolygonVertices);

        /// <summary>
        /// Get the offsetted polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input polygon</param>
        /// <param name="xCoordOffsetPoint">The x coordinate of the offset point</param>
        /// <param name="yCoordOffsetPoint">The y coordinate of the offset point</param>
        /// <param name="innerPolygon">Compute inner polygon or not</param>
        /// <param name="geometryListOut">The offsetted polygon</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_offsetted_polygon", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetOffsettedPolygon([In] ref int gridStateId, [In] ref GeometryListNative geometryListIn, [In] ref bool innerPolygon, [In] ref double distance, [In, Out] ref GeometryListNative geometryListOut);

        /// <summary>
        /// Count the number of vertices after polygon refinment
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input polygon</param>
        /// <param name="firstIndex">The index of the first vertex</param>
        /// <param name="secondIndex">The index of the second vertex</param>
        /// <param name="distance">The refinement distance</param>
        /// <param name="numberOfPolygonVertices">The number of vertices after refinement </param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_refine_polygon_count", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int CountVerticesRefinededPolygon([In] ref int gridStateId, [In] ref GeometryListNative geometryListIn, [In] ref int firstIndex, [In] ref int secondIndex, [In] ref double distance, [In, Out] ref int numberOfPolygonVertices);

        /// <summary>
        /// Gets the refined polygon
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input polygons</param>
        /// <param name="firstIndex">The index of the first vertex</param>
        /// <param name="secondIndex">The index of the second vertex</param>
        /// <param name="distance">The refinement distance</param>
        /// <param name="geometryListOut"></param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_refine_polygon", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetRefinededPolygon([In] ref int gridStateId, [In] ref GeometryListNative geometryListIn, [In] ref int firstIndex, [In] ref int secondIndex, [In] ref double distance, [In, Out] ref GeometryListNative geometryListOut);

        #endregion

        /// <summary>
        /// Counts the number of vertices inside selected polygons
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input polygons</param>
        /// <param name="numberOfMeshVertices">The number of selected vertices</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_count_vertices_in_polygons", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int CountVerticesInPolygon(ref int gridStateId, ref GeometryListNative geometryListIn, ref int inside, ref int numberOfMeshVertices);

        /// <summary>
        /// Gets the selected vertexes indexes
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The input polygons</param>
        /// <param name="numberOfMeshVertices">The number of selected vertices</param>
        /// <param name="selectedVerticesPtr">The selected vertices indexes</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_vertices_in_polygons", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetSelectedVerticesInPolygon(ref int gridStateId, ref GeometryListNative geometryListIn, ref int inside, ref int numberOfMeshVertices, ref IntPtr selectedVerticesPtr);

        /// <summary>
        /// Gets the orthogonality
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The orthogonality values of each edge</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_get_orthogonality", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetOrthogonality(ref int gridStateId, ref GeometryListNative geometryListIn);

        /// <summary>
        /// Gets the smoothness 
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="geometryListIn">The smoothness values of each edge</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_get_smoothness", CallingConvention = CallingConvention.Cdecl)]
        internal static extern int GetSmoothness(ref int gridStateId, ref GeometryListNative geometryListIn);

        /// <summary>
        /// Get the index of the closest vertex
        /// </summary>
        /// <param name="gridStateId">Id of the grid state</param>
        /// <param name="xCoordinate">x coordinate of the vertex</param>
        /// <param name="yCoordinate">y coordinate of the vertex</param>
        /// <param name="zCoordinate">z coordinate of the vertex</param>
        /// <param name="searchRadius">the radius where to search for the vertex</param>
        /// <param name="vertexIndex">the index of the closest vertex</param>
        /// <returns></returns>
        [DllImport(GridGeomStatefulDllName, EntryPoint = "ggeo_get_vertex_index", CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetVertexIndex(ref int gridStateId, ref double xCoordinate, ref double yCoordinate, ref double zCoordinate, ref double searchRadius, ref int vertexIndex);
    }
}