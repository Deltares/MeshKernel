using System;
using DelftTools.Utils.Remoting;
using ProtoBufRemote;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    public class GridgeomStatefulApiRemote : IGridgeomStatefulApi
    {
        private readonly IGridgeomStatefulApi gridgeomStatefulApi;

        static GridgeomStatefulApiRemote()
        {
            RemotingTypeConverters.RegisterTypeConverter(new CoordinateConverter());
        }

        public GridgeomStatefulApiRemote()
        {
            gridgeomStatefulApi = RemoteInstanceContainer.CreateInstance<IGridgeomStatefulApi, GridgeomStatefulApi>();
        }

        /// <inheritdoc />
        public int CreateGridState()
        {
            return gridgeomStatefulApi.CreateGridState();
        }

        /// <inheritdoc />
        public bool RemoveGridState(int gridStateId)
        {
            return gridgeomStatefulApi.RemoveGridState(gridStateId);
        }

        /// <inheritdoc />
        public bool SetGridState(int gridStateId, DisposableMeshGeometry disposableMeshGeometry, bool isGeographic)
        {
            return gridgeomStatefulApi.SetGridState(gridStateId, disposableMeshGeometry, isGeographic);
        }

        /// <inheritdoc />
        public DisposableMeshGeometry GetGridState(int gridStateId)
        {
            return gridgeomStatefulApi.GetGridState(gridStateId);
        }

        /// <inheritdoc />
        public bool DeleteVertex(int gridStateId, int vertexIndex)
        {
            return gridgeomStatefulApi.DeleteVertex(gridStateId, vertexIndex);
        }

        /// <inheritdoc />
        public bool FlipEdges(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired, ProjectToLandBoundaryOptions projectToLandBoundaryOption)
        {
            return gridgeomStatefulApi.FlipEdges(gridStateId, isTriangulationRequired, isAccountingForLandBoundariesRequired, projectToLandBoundaryOption);
        }

        /// <inheritdoc />
        public bool MergeTwoVertices(int gridStateId, int startVertexIndex, int endVertexIndex)
        {
            return gridgeomStatefulApi.MergeTwoVertices(gridStateId, startVertexIndex, endVertexIndex);
        }

        /// <inheritdoc />
        public bool MergeVertices(int gridStateId, DisposableGeometryList disposableGeometryList)
        {
            return gridgeomStatefulApi.MergeVertices(gridStateId, disposableGeometryList);
        }

        /// <inheritdoc />
        public bool Orthogonalize(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired,
            ProjectToLandBoundaryOptions projectToLandBoundaryOption,
            OrthogonalizationParameters orthogonalizationParameters, DisposableGeometryList geometryListNativePolygon, DisposableGeometryList geometryListNativeLandBoundaries)
        {
            return gridgeomStatefulApi.Orthogonalize(gridStateId, isTriangulationRequired, isAccountingForLandBoundariesRequired,
                projectToLandBoundaryOption, orthogonalizationParameters, geometryListNativePolygon, geometryListNativeLandBoundaries);
        }

        /// <inheritdoc />
        public bool MakeGrid(int gridStateId, MakeGridParameters makeGridParameters, ref DisposableGeometryList disposableGeometryListIn)
        {
            return gridgeomStatefulApi.MakeGrid(gridStateId, makeGridParameters, ref disposableGeometryListIn);
        }

        /// <inheritdoc />
        public bool GetSplines(DisposableGeometryList disposableGeometryListIn, ref DisposableGeometryList disposableGeometryListOut, int numberOfPointsBetweenVertices)
        {
            return gridgeomStatefulApi.GetSplines(disposableGeometryListIn, ref disposableGeometryListOut, numberOfPointsBetweenVertices);
        }

        /// <inheritdoc />
        public bool MakeGridFromSplines(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
            ref CurvilinearParameters curvilinearParameters)
        {
            return gridgeomStatefulApi.MakeGridFromSplines(gridStateId, ref disposableGeometryListIn, ref curvilinearParameters);
        }

        /// <inheritdoc />
        public bool MakeOrthogonalGridFromSplines(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
            ref CurvilinearParameters curvilinearParameters, ref SplinesToCurvilinearParameters splinesToCurvilinearParameters)
        {
            return gridgeomStatefulApi.MakeOrthogonalGridFromSplines(gridStateId, ref disposableGeometryListIn, ref curvilinearParameters, ref splinesToCurvilinearParameters);
        }

        /// <inheritdoc />
        public bool MakeTriangularGridInPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryList)
        {
            return gridgeomStatefulApi.MakeTriangularGridInPolygon(gridStateId, ref disposableGeometryList);
        }

        public bool MakeTriangularGridFromSamples(int gridStateId, ref DisposableGeometryList disposableGeometryList)
        {
            return gridgeomStatefulApi.MakeTriangularGridFromSamples(gridStateId, ref disposableGeometryList);
        }

        /// <inheritdoc />
        public bool CountMeshBoundaryPolygonVertices(int gridStateId, ref int numberOfPolygonVertices)
        {
            return gridgeomStatefulApi.CountMeshBoundaryPolygonVertices(gridStateId, ref numberOfPolygonVertices);
        }

        /// <inheritdoc />
        public bool GetMeshBoundaryPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryList)
        {
            return gridgeomStatefulApi.GetMeshBoundaryPolygon(gridStateId, ref disposableGeometryList);
        }

        /// <inheritdoc />
        public bool CountVerticesOffsettedPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
            bool innerPolygon, double distance, ref int numberOfPolygonVertices)
        {
            return gridgeomStatefulApi.CountVerticesOffsettedPolygon(gridStateId, ref disposableGeometryListIn,
                innerPolygon, distance, ref numberOfPolygonVertices);
        }

        /// <inheritdoc />
        public bool GetOffsettedPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
            bool innerPolygon, double distance, ref DisposableGeometryList disposableGeometryListOut)
        {
            return gridgeomStatefulApi.GetOffsettedPolygon(gridStateId, ref disposableGeometryListIn,
                innerPolygon, distance, ref disposableGeometryListOut);
        }

        /// <inheritdoc />
        public bool CountVerticesRefinededPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,  int firstIndex,
             int secondIndex,  double distance,  ref int numberOfPolygonVertices)
        {
            return gridgeomStatefulApi.CountVerticesRefinededPolygon( gridStateId, ref disposableGeometryListIn,  firstIndex,
             secondIndex,  distance,  ref numberOfPolygonVertices);
        }

        /// <inheritdoc />
        public bool GetRefinededPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,  int firstIndex,
             int secondIndex,  double distance, ref DisposableGeometryList disposableGeometryListOut)
        {
            return gridgeomStatefulApi.GetRefinededPolygon(gridStateId, ref disposableGeometryListIn,  firstIndex,
                 secondIndex,  distance, ref disposableGeometryListOut);
        }

        /// <inheritdoc />
        public bool RefineGridBasedOnSamples(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, InterpolationParameters interpolationParameters,  SamplesRefineParameters samplesRefineParameters)
        {
            return gridgeomStatefulApi.RefineGridBasedOnSamples(gridStateId, ref disposableGeometryListIn, interpolationParameters, samplesRefineParameters);
        }

        /// <inheritdoc />
        public bool RefineGridBasedOnPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
            InterpolationParameters interpolationParameters)
        {
            return gridgeomStatefulApi.RefineGridBasedOnPolygon(gridStateId, ref disposableGeometryListIn, interpolationParameters);
        }

        /// <inheritdoc />
        public int[] GetSelectedVertices(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, int inside)
        {
            return gridgeomStatefulApi.GetSelectedVertices(gridStateId, ref disposableGeometryListIn, inside);
        }

        /// <inheritdoc />
        public bool GetOrthogonality(int gridStateId, ref DisposableGeometryList disposableGeometryListOut)
        {
            return gridgeomStatefulApi.GetOrthogonality(gridStateId, ref disposableGeometryListOut);
        }

        public bool GetSmoothness(int gridStateId, ref DisposableGeometryList disposableGeometryListOut)
        {
            return gridgeomStatefulApi.GetSmoothness(gridStateId, ref disposableGeometryListOut);
        }

        /// <inheritdoc />
        public bool InsertVertex(int gridGeomId, double xCoordinate, double yCoordinate, double zCoordinate, ref int vertexIndex)
        {
            return gridgeomStatefulApi.InsertVertex(gridGeomId, xCoordinate, yCoordinate, zCoordinate, ref vertexIndex);
        }

        /// <inheritdoc />
        public bool DeleteMeshWithOptions(int gridGeomId, ref DisposableGeometryList disposableGeometryListOut, ref int deletionOption)
        {
            return gridgeomStatefulApi.DeleteMeshWithOptions(gridGeomId, ref disposableGeometryListOut, ref deletionOption);
        }

        public bool OrthogonalizationInitialize(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired,
            ProjectToLandBoundaryOptions projectToLandBoundaryOption,
            OrthogonalizationParameters orthogonalizationParameters, DisposableGeometryList geometryListNativePolygon, DisposableGeometryList geometryListNativeLandBoundaries)
        {
            return gridgeomStatefulApi.OrthogonalizationInitialize(gridStateId, isTriangulationRequired, isAccountingForLandBoundariesRequired,
                projectToLandBoundaryOption, orthogonalizationParameters, geometryListNativePolygon, geometryListNativeLandBoundaries);
        }


        public bool OrthogonalizationPrepareOuterIteration(int gridGeomId)
        {
            return gridgeomStatefulApi.OrthogonalizationPrepareOuterIteration(gridGeomId);
        }

        public bool OrthogonalizationInnerIteration(int gridGeomId)
        {
            return gridgeomStatefulApi.OrthogonalizationInnerIteration(gridGeomId);
        }

        public bool OrthogonalizationFinalizeOuterIteration(int gridGeomId)
        {
            return gridgeomStatefulApi.OrthogonalizationFinalizeOuterIteration(gridGeomId);
        }

        public bool OrthogonalizationDelete(int gridGeomId)
        {
            return gridgeomStatefulApi.OrthogonalizationDelete(gridGeomId);
        }

        /// <inheritdoc />
        public bool InsertEdge(int gridStateId,  int startVertexIndex, int endVertexIndex, ref int edgeIndex)
        {
            return gridgeomStatefulApi.InsertEdge(gridStateId, startVertexIndex, endVertexIndex, ref edgeIndex);
        }

        /// <inheritdoc />
        public bool GetVertexIndex(int gridGeomId, double xCoordinate, double yCoordinate, double zCoordinate,
            double searchRadius, ref int vertexIndex)
        {
            return gridgeomStatefulApi.GetVertexIndex(gridGeomId, xCoordinate, yCoordinate, zCoordinate,
                searchRadius, ref vertexIndex);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (!disposing) return;

            gridgeomStatefulApi?.Dispose();

            if (RemoteInstanceContainer.NumInstances == 0) return;

            RemoteInstanceContainer.RemoveInstance(gridgeomStatefulApi);
        }
    }
}