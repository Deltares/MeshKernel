using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Runtime.InteropServices;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Native;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [ExcludeFromCodeCoverage] 
    // Excluded because it is tested through the GridgeomStatefulApiRemote
    // DotCover on the buildserver does not work correctly with remoting
    public sealed class GridgeomStatefulApi : IGridgeomStatefulApi
    {
        /// <inheritdoc />
        public int CreateGridState()
        {
            var gridStateId = 0;
            GridgeomStatefulDll.CreateGridState(ref gridStateId);
            return gridStateId;
        }

        /// <inheritdoc />
        public bool RemoveGridState(int gridStateId)
        {
            return GridgeomStatefulDll.RemoveGridState(ref gridStateId) == 0;
        }

        /// <inheritdoc />
        public bool SetGridState(int gridStateId, DisposableMeshGeometry disposableMeshGeometry, bool isGeographic)
        {
            var meshGeometryDimensions = disposableMeshGeometry.CreateMeshDimensions();
            var meshGeometry = disposableMeshGeometry.CreateMeshGeometry();

            var result = GridgeomStatefulDll.SetState(ref gridStateId, ref meshGeometryDimensions, ref meshGeometry, ref isGeographic);

            if (disposableMeshGeometry.IsMemoryPinned)
            {
                disposableMeshGeometry.UnPinMemory();
            }

            return result == 0;
        }
        
        /// <inheritdoc />
        public DisposableMeshGeometry GetGridState(int gridStateId)
        {
            var newMeshGeometry = new MeshGeometry();
            var newMeshDimensions = new MeshGeometryDimensions();

            GridgeomStatefulDll.GetMeshState(ref gridStateId, ref newMeshDimensions, ref newMeshGeometry);

            return new DisposableMeshGeometry
            {
                xNodes = newMeshGeometry.nodex.CreateValueArray<double>(newMeshDimensions.numnode),
                yNodes = newMeshGeometry.nodey.CreateValueArray<double>(newMeshDimensions.numnode),
                zNodes = newMeshGeometry.nodez.CreateValueArray<double>(newMeshDimensions.numnode),
                edgeNodes = newMeshGeometry.edge_nodes.CreateValueArray<int>(newMeshDimensions.numedge * 2).ToArray(),
                numberOfEdges = newMeshDimensions.numedge,
                numberOfFaces = newMeshDimensions.numface,
                numberOfNodes = newMeshDimensions.numnode,
                maxNumberOfFaceNodes = newMeshDimensions.maxnumfacenodes
            };
        }

        /// <inheritdoc />
        public bool DeleteVertex(int gridStateId, int vertexIndex)
        {
            return GridgeomStatefulDll.DeleteNode(ref gridStateId, ref vertexIndex) == 0;
        }

        public bool FlipEdges(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired, ProjectToLandBoundaryOptions projectToLandBoundaryOption)
        {
            int isTriangulationRequiredInt = isTriangulationRequired ? 1 : 0;
            int isAccountingForLandBoundariesRequiredInt = isAccountingForLandBoundariesRequired ? 1 : 0;
            int projectToLandBoundaryOptionInt = (int)projectToLandBoundaryOption;
            return GridgeomStatefulDll.FlipEdges(ref gridStateId, ref isTriangulationRequiredInt, ref isAccountingForLandBoundariesRequiredInt, ref projectToLandBoundaryOptionInt) == 0;
        }

        /// <inheritdoc />
        public bool InsertEdge(int gridStateId, int startVertexIndex, int endVertexIndex, ref int edgeIndex)
        {
            return GridgeomStatefulDll.InsertEdge(ref gridStateId, ref startVertexIndex, ref endVertexIndex, ref edgeIndex) == 0;
        }

        /// <inheritdoc />
        public bool MergeTwoVertices(int gridStateId, int startVertexIndex, int endVertexIndex)
        {
            return GridgeomStatefulDll.MergeTwoVertices(ref gridStateId, ref startVertexIndex, ref endVertexIndex) == 0;
        }

        /// <inheritdoc />
        public bool MergeVertices(int gridStateId, DisposableGeometryList disposableGeometryList)
        {
            var geometryList = disposableGeometryList.CreateGeometryListNative();
            return GridgeomStatefulDll.MergeVertices(ref gridStateId, ref geometryList) == 0;
        }
        
        public bool Orthogonalize(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired,
            ProjectToLandBoundaryOptions projectToLandBoundaryOption,
            OrthogonalizationParameters orthogonalizationParameters, DisposableGeometryList geometryListNativePolygon,
            DisposableGeometryList geometryListNativeLandBoundaries)
        {
            int isTriangulationRequiredInt = isTriangulationRequired ? 1 : 0;
            int isAccountingForLandBoundariesRequiredInt = isAccountingForLandBoundariesRequired ? 1 : 0;
            int projectToLandBoundaryOptionInt = (int)projectToLandBoundaryOption;

            var nativeOrthogonalizationParameters = orthogonalizationParameters.ToOrthogonalizationParametersNative();

            var geometryListPolygon = geometryListNativePolygon?.CreateGeometryListNative() ??
                               new GeometryListNative { numberOfCoordinates = 0 };
            var geometryListLandBoundaries = geometryListNativeLandBoundaries?.CreateGeometryListNative() ??
                                      new GeometryListNative { numberOfCoordinates = 0 };

            return GridgeomStatefulDll.Orthogonalization(ref gridStateId, ref isTriangulationRequiredInt,
                           ref isAccountingForLandBoundariesRequiredInt, ref projectToLandBoundaryOptionInt,
                           ref nativeOrthogonalizationParameters, ref geometryListPolygon,ref geometryListLandBoundaries) == 0;
        }

        public bool MakeGrid(int gridStateId, MakeGridParameters makeGridParameters, ref DisposableGeometryList disposableGeometryListIn)
        {
            var makeGridParametersNative =
                makeGridParameters.ToMakeGridParametersNative();
            var geometryListNative = disposableGeometryListIn.CreateGeometryListNative();

            return GridgeomStatefulDll.MakeGrid(ref gridStateId, ref makeGridParametersNative, ref geometryListNative) ==0;
        }

        public bool GetSplines(DisposableGeometryList disposableGeometryListIn, ref DisposableGeometryList disposableGeometryListOut, int numberOfPointsBetweenVertices)
        {
            var geometryListIn = disposableGeometryListIn.CreateGeometryListNative();
            var geometryListOut = disposableGeometryListOut.CreateGeometryListNative();
            if (GridgeomStatefulDll.GetSplines(ref geometryListIn, ref geometryListOut,
                    ref numberOfPointsBetweenVertices) != 0)
            {
                return false;
            }
            disposableGeometryListOut.NumberOfCoordinates = geometryListOut.numberOfCoordinates;
            return true;
        }

        public bool MakeGridFromSplines(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
            ref CurvilinearParameters curvilinearParameters)
        {
            var geometryListIn = disposableGeometryListIn.CreateGeometryListNative();
            var curvilinearParametersNative = curvilinearParameters.ToCurvilinearParametersNative();
            return GridgeomStatefulDll.MakeGridFromSplines(ref gridStateId, ref geometryListIn,
                       ref curvilinearParametersNative) == 0;
        }

        public bool MakeOrthogonalGridFromSplines(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,
             ref CurvilinearParameters curvilinearParameters, ref SplinesToCurvilinearParameters splinesToCurvilinearParameters)
        {
            var geometryListNative = disposableGeometryListIn.CreateGeometryListNative();
            var curvilinearParametersNative = curvilinearParameters.ToCurvilinearParametersNative();
            var splinesToCurvilinearParametersNative = splinesToCurvilinearParameters.ToSplinesToCurvilinearParametersNative();
            return GridgeomStatefulDll.MakeOrthogonalGridFromSplines(ref gridStateId, ref geometryListNative,
                        ref curvilinearParametersNative, ref splinesToCurvilinearParametersNative) == 0;
        }

        public bool MakeTriangularGridInPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryList)
        {
            var geometryListNative = disposableGeometryList.CreateGeometryListNative();
            return GridgeomStatefulDll.MakeTriangularGridFromPolygon(ref gridStateId, ref geometryListNative) == 0;
        }

        public bool MakeTriangularGridFromSamples(int gridStateId, ref DisposableGeometryList disposableGeometryList)
        {
            var geometryListNative = disposableGeometryList.CreateGeometryListNative();
            return GridgeomStatefulDll.MakeTriangularGridFromSamples(ref gridStateId, ref geometryListNative) == 0;
        }

        public bool CountMeshBoundaryPolygonVertices(int gridStateId, ref int numberOfPolygonVertices)
        {
            return GridgeomStatefulDll.CountMeshBoundaryPolygonVertices(ref gridStateId, ref numberOfPolygonVertices) == 0;
        }

        public bool GetMeshBoundaryPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryList)
        {
            var geometryListNative = disposableGeometryList.CreateGeometryListNative();
            return GridgeomStatefulDll.GetMeshBoundaryPolygon(ref gridStateId, ref geometryListNative) == 0;
        }

        public bool CountVerticesOffsettedPolygon( int gridStateId,  ref DisposableGeometryList disposableGeometryListIn, bool innerPolygon, double distance,ref int numberOfPolygonVertices)
        {
            var geometryListNativeIn = disposableGeometryListIn.CreateGeometryListNative();
            return GridgeomStatefulDll.CountVerticesOffsettedPolygon(ref gridStateId, ref geometryListNativeIn, ref innerPolygon, ref distance, ref numberOfPolygonVertices) == 0;
        }

        public bool GetOffsettedPolygon( int gridStateId,  ref DisposableGeometryList disposableGeometryListIn, bool innerPolygon, double distance, ref DisposableGeometryList disposableGeometryListOut)
        {
            var geometryListNativeIn = disposableGeometryListIn.CreateGeometryListNative();
            var geometryListNativeOut = disposableGeometryListOut.CreateGeometryListNative();
            return GridgeomStatefulDll.GetOffsettedPolygon(ref gridStateId, ref geometryListNativeIn, ref innerPolygon, ref distance, ref geometryListNativeOut) == 0;
        }

        public bool CountVerticesRefinededPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,  int firstIndex,  int secondIndex,  double distance,  ref int numberOfPolygonVertices)
        {
            var geometryListInNative = disposableGeometryListIn.CreateGeometryListNative();
            return GridgeomStatefulDll.CountVerticesRefinededPolygon(ref gridStateId, ref geometryListInNative, ref firstIndex, ref secondIndex, ref distance, ref numberOfPolygonVertices) == 0;
        }

        public bool GetRefinededPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn,  int firstIndex,
             int secondIndex,  double distance, ref DisposableGeometryList disposableGeometryListOut)
        {
            var geometryListNativeIn = disposableGeometryListIn.CreateGeometryListNative();
            var geometryListNativeOut = disposableGeometryListOut.CreateGeometryListNative();
            return GridgeomStatefulDll.GetRefinededPolygon(ref gridStateId, ref geometryListNativeIn, ref firstIndex, ref secondIndex, ref distance, ref geometryListNativeOut) == 0;
        }

        public bool RefineGridBasedOnSamples(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, InterpolationParameters interpolationParameters,  SamplesRefineParameters samplesRefineParameters)
        {
            var disposableGeometryListInNative = disposableGeometryListIn.CreateGeometryListNative();
            var interpolationParametersNative = interpolationParameters.ToInterpolationParametersNative();
            var samplesRefineParametersNative = samplesRefineParameters.ToSampleRefineParametersNative();
            return GridgeomStatefulDll.RefineGridBasedOnSamples(ref gridStateId, ref disposableGeometryListInNative, ref interpolationParametersNative, ref samplesRefineParametersNative) == 0;
        }

        public bool RefineGridBasedOnPolygon(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, InterpolationParameters interpolationParameters)
        {
            var disposableGeometryListInNative = disposableGeometryListIn.CreateGeometryListNative();
            var interpolationParametersNative = interpolationParameters.ToInterpolationParametersNative();
            return GridgeomStatefulDll.RefineGridBasedOnPolygon(ref gridStateId, ref disposableGeometryListInNative, ref interpolationParametersNative) == 0;
        }

        public int[] GetSelectedVertices(int gridStateId, ref DisposableGeometryList disposableGeometryListIn, int inside)
        {
            var geometryListNativeIn = disposableGeometryListIn.CreateGeometryListNative();
            int numberOfMeshVertices = -1;
            GridgeomStatefulDll.CountVerticesInPolygon(ref gridStateId, ref geometryListNativeIn, ref inside, ref numberOfMeshVertices);
            IntPtr selectedVerticesPtr = Marshal.AllocCoTaskMem(Marshal.SizeOf(typeof(int)) * numberOfMeshVertices);
            GridgeomStatefulDll.GetSelectedVerticesInPolygon(ref gridStateId, ref geometryListNativeIn, ref inside, ref numberOfMeshVertices, ref selectedVerticesPtr);
            int[] selectedVertices = new int[numberOfMeshVertices];
            Marshal.Copy(selectedVerticesPtr, selectedVertices, 0, numberOfMeshVertices);
            return selectedVertices;
        }

        public bool GetOrthogonality(int gridStateId, ref DisposableGeometryList disposableGeometryListOut)
        {
            var geometryListNativeInOut = disposableGeometryListOut.CreateGeometryListNative();
            return GridgeomStatefulDll.GetOrthogonality(ref gridStateId, ref geometryListNativeInOut) == 0;
        }

        public bool GetSmoothness(int gridStateId, ref DisposableGeometryList disposableGeometryListOut)
        {
            var geometryListNativeInOut = disposableGeometryListOut.CreateGeometryListNative();
            return GridgeomStatefulDll.GetSmoothness(ref gridStateId, ref geometryListNativeInOut) == 0;
        }

        public bool InsertVertex(int gridGeomId, double xCoordinate, double yCoordinate, double zCoordinate, ref int vertexIndex)
        {
            return GridgeomStatefulDll.InsertVertex(ref gridGeomId, ref xCoordinate, ref yCoordinate, ref zCoordinate, ref vertexIndex) == 0;
        }

        public bool GetVertexIndex(int gridGeomId, double xCoordinate, double yCoordinate, double zCoordinate,
            double searchRadius, ref int vertexIndex)
        {
            return GridgeomStatefulDll.GetVertexIndex(ref gridGeomId, ref xCoordinate, ref yCoordinate, ref zCoordinate, ref searchRadius, ref vertexIndex) == 0;
        }

        public bool DeleteMeshWithOptions(int gridGeomId, ref DisposableGeometryList disposableGeometryListOut, ref int deletionOption)
        {
            var geometryListNativeIn = disposableGeometryListOut.CreateGeometryListNative();
            return GridgeomStatefulDll.DeleteMeshWithOptions(ref gridGeomId, ref geometryListNativeIn, ref  deletionOption) == 0;
        }

        public bool OrthogonalizationInitialize(int gridStateId, bool isTriangulationRequired, bool isAccountingForLandBoundariesRequired,
            ProjectToLandBoundaryOptions projectToLandBoundaryOption,
            OrthogonalizationParameters orthogonalizationParameters, DisposableGeometryList geometryListNativePolygon,
            DisposableGeometryList geometryListNativeLandBoundaries)
        {
            int isTriangulationRequiredInt = isTriangulationRequired ? 1 : 0;
            int isAccountingForLandBoundariesRequiredInt = isAccountingForLandBoundariesRequired ? 1 : 0;
            int projectToLandBoundaryOptionInt = (int)projectToLandBoundaryOption;

            var nativeOrthogonalizationParameters = orthogonalizationParameters.ToOrthogonalizationParametersNative();

            var geometryListPolygon = geometryListNativePolygon?.CreateGeometryListNative() ??
                                      new GeometryListNative { numberOfCoordinates = 0 };
            var geometryListLandBoundaries = geometryListNativeLandBoundaries?.CreateGeometryListNative() ??
                                             new GeometryListNative { numberOfCoordinates = 0 };

            return GridgeomStatefulDll.OrthogonalizationInitialize(ref gridStateId, ref isTriangulationRequiredInt,
                       ref isAccountingForLandBoundariesRequiredInt, ref projectToLandBoundaryOptionInt,
                       ref nativeOrthogonalizationParameters, ref geometryListPolygon, ref geometryListLandBoundaries) == 0;
        }

        public bool OrthogonalizationPrepareOuterIteration(int gridGeomId)
        {
            return GridgeomStatefulDll.OrthogonalizationPrepareOuterIteration(ref gridGeomId) == 0;
        }

        public bool OrthogonalizationInnerIteration(int gridGeomId)
        {
            return GridgeomStatefulDll.OrthogonalizationInnerIteration(ref gridGeomId) == 0;
        }

        public bool OrthogonalizationFinalizeOuterIteration(int gridGeomId)
        {
            return GridgeomStatefulDll.OrthogonalizationFinalizeOuterIteration(ref gridGeomId) == 0;
        }

        public bool OrthogonalizationDelete(int gridGeomId)
        {
            return GridgeomStatefulDll.OrthogonalizationDelete(ref gridGeomId) == 0;
        }

        public void Dispose()
        {
            // Do nothing because no remoting is used
        }
    }
}