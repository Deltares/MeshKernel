using System;
using System.Threading.Tasks;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using DeltaShell.Plugins.GridEditor.Helpers;

namespace DeltaShell.Plugins.GridEditor.Gui.Controllers.Api
{
    internal interface IGridEditorController
    {
        /// <summary>
        /// Is in edit state (mode)
        /// </summary>
        bool IsEditing { get; set; }

        /// <summary>
        /// Sets the size of the points for rendering
        /// </summary>
        float VertexSize { get; set; }

        /// <summary>
        /// Enable vertices in rendering
        /// </summary>
        bool ShowVertices { get; set; }

        /// <summary>
        /// Subselection of the orthogonalization parameters that can be set
        /// </summary>
        OrthogonalizationParameters OrthogonalizationParameters { get; }

        /// <summary>
        /// The parameters for the make grid algorithm
        /// </summary>
        MakeGridParameters MakeGridParameters { get; }

        /// <summary>
        /// The parameters for curvilinear grids
        /// </summary>
        CurvilinearParameters CurvilinearParameters { get; }

        /// <summary>
        /// The parameters of the advancing front algorithm
        /// </summary>
        SplinesToCurvilinearParameters SplinesToCurvilinearParameters { get; }

        /// <summary>
        /// The parameters for interpolating samples
        /// </summary>
        InterpolationParameters InterpolationParameters { get; }

        /// <summary>
        /// The parameters for refining the grid based on samples
        /// </summary>
        SamplesRefineParameters SamplesRefineParameters { get; }

        /// <summary>
        /// The distance used when refining the polygon
        /// </summary>
        double PolygonRefinementDistance { get; set; }

        /// <summary>
        /// The distance used for offsetting a polygon
        /// </summary>
        double PolygonOffsetDistance { get; set; }

        /// <summary>
        /// Triangulate all cell option for flip edges and orthogonalization
        /// </summary>
        bool  IsTriangulationRequired { get; set; }

        /// <summary>
        /// Account for land boundaries in flip edges and orthogonalization
        /// </summary>
        bool IsAccountingForLandBoundariesRequired { get; set; }

        /// <summary>
        /// Project to land boundary option used in flip edges and orthogonalization
        /// </summary>
        ProjectToLandBoundaryOptions ProjectToLandBoundaryOption { get; set; }

        /// <summary>
        /// Option for mesh deletion
        /// </summary>
        DeleteMeshOptions DeleteMeshOption { get; set; }

        /// <summary>
        /// Currently selected <see cref="MapToolType"/>
        /// </summary>
        MapToolType SelectedMapToolType { get; }

        /// <summary>
        /// Draws the <see cref="GridEditorState.SelectedVertices"/>
        /// </summary>
        bool DrawSelectedVertices { get; set; }

        /// <summary>
        /// Inverts the selection of the vertices <see cref="GridEditorState.SelectedVertices"/>
        /// </summary>
        bool InvertVerticesSelection { get; set; }

        /// <summary>
        /// Draw the <see cref="GridEditorState.Polygons"/> filled in
        /// </summary>
        bool DrawFilledInPolygons { get; set; }

        /// <summary>
        /// <see cref="IMapInteractor"/> for doing map related actions
        /// </summary>
        IMapInteractor MapInteractor { get; set; }

        /// <summary>
        /// Function to check if the changes to the grid should be committed
        /// </summary>
        Func<bool> CommitState { get; set; }

        /// <summary>
        /// Has sample values declared
        /// </summary>
        bool HasSamples { get; }

        /// <summary>
        /// Currently selected workflow for generating grids
        /// </summary>
        GenerateGridWorkFlowType SelectedGenerateGridWorkFlowType { get; set; }

        /// <summary>
        /// Resets the unstructured grids for ribbon state
        /// </summary>
        void ResetUnstructuredGrids();

        /// <summary>
        /// Enables the map tool of the specified type <see cref="MapToolType"/>
        /// </summary>
        /// <param name="mapToolType"></param>
        void EnableMapToolByType(MapToolType mapToolType);

        /// <summary>
        /// Deletes the currently selected vertices
        /// </summary>
        void DeleteSelectedVertices();

        /// <summary>
        /// Deletes mesh with options
        /// </summary>
        void DeleteMeshWithOptions();

        /// <summary>
        /// Deletes samples inside a polygon
        /// </summary>
        void DeleteSamplesInPolygons();

        /// <summary>
        /// Flips the edges 
        /// </summary>
        void FlipEdges();

        /// <summary>
        /// Merge all vertices in a polygon or all vertices of a grid
        /// </summary>
        bool MergeVertices();

        /// <summary>
        /// Orthogonalize unstructured grid
        /// </summary>
        Task Orthogonalize();

        /// <summary>
        /// Makes a new grid
        /// </summary>
        void MakeGrid();

        /// <summary>
        /// Makes a new grid from splines
        /// </summary>
        void MakeGridFromSplines();

        /// <summary>
        /// Make curvilinear grid from splines with advancing front.
        /// </summary>
        Task MakeGridFromSplinesOrthogonal();

        /// <summary>
        /// Makes a triangular grid inside a polygon
        /// </summary>
        void MakeTriangularGridInPolygon();

        /// <summary>
        /// Makes a triangular grid from samples
        /// </summary>
        void MakeTriangularGridFromSamples();

        /// <summary>
        /// Retrieves the mesh boundary polygon
        /// </summary>
        void GetMeshBoundaryPolygon();

        /// <summary>
        /// Refines a mesh based on samples
        /// </summary>
        void RefineGridBasedOnSamples();

        /// <summary>
        /// Refines a mesh based on samples
        /// </summary>
        void RefineGridBasedOnPolygon();

        /// <summary>
        /// Given a polygon, create another one with an offset
        /// </summary>
        void AddOffsettedPolygon();

        /// <summary>
        /// Imports or exports (depending on <paramref name="action"/>) the defined type (<paramref name="importExportType"/>) to
        /// <paramref name="filePath"/>
        /// </summary>
        /// <param name="action">Import or export action</param>
        /// <param name="importExportType">Type of data to import/export</param>
        /// <param name="filePath">Path to import/export towards</param>
        void ImportExport(ImportExportAction action, ImportExportType importExportType, string filePath);

        /// <summary>
        /// Refines the vertices of the selected polygons
        /// </summary>
        void RefineSelectedPolygons();

        /// <summary>
        /// Show orthogonality 
        /// </summary>
        void ShowOrthogonality();

        /// <summary>
        /// Show smoothness 
        /// </summary>
        void ShowSmoothness();

        /// <summary>
        /// Generates a grid for the <see cref="SelectedGenerateGridWorkFlowType"/>
        /// </summary>
        void GenerateGridForSelectedWorkflow();

        /// <summary>
        /// Clip first selected polygon with second polygon
        /// </summary>
        void ClipSelectedPolygons();

        /// <summary>
        /// Merge first selected polygon with second polygon
        /// </summary>
        void MergeSelectedPolygons();

        /// <summary>
        /// Subtract first selected polygon with second polygon
        /// </summary>
        void SubtractSelectedPolygons();

        void RequestCancelation();
    }
}