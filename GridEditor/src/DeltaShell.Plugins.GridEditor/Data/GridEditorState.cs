using System.Collections.Generic;
using DelftTools.Utils.Collections.Generic;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using GeoAPI.Extensions.CoordinateSystems;
using GeoAPI.Extensions.Feature;
using NetTopologySuite.Extensions.Coverages;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.Data
{
    public class GridEditorState
    {
        private readonly EventedList<IFeature> polygons;

        public GridEditorState()
        {
            polygons = new EventedList<IFeature>();
            polygons.CollectionChanged += (sender, args) => { SelectedVertices = null; };
            polygons.PropertyChanged += (sender, args) =>
            {
                if (args.PropertyName == nameof(IFeature.Geometry))
                {
                    SelectedVertices = null;
                }
            };
        }

        /// <summary>
        /// A cloud of sample points
        /// </summary>
        public PointCloud SamplePoints { get; } = new PointCloud();

        /// <summary>
        /// List of land boundaries
        /// </summary>
        public IList<IFeature> LandBoundaries { get; } = new List<IFeature>();
        
        /// <summary>
        /// List of splines
        /// </summary>
        public IList<Spline> Splines { get; } = new List<Spline>();

        /// <summary>
        /// Selection polygon(s)
        /// </summary>
        public IList<IFeature> Polygons
        {
            get { return polygons; }
        }

        /// <summary>
        /// Selected points of the grid
        /// </summary>
        public int[] SelectedVertices { get; set; }

        /// <summary>
        /// <see cref="DisposableMeshGeometry"/> representing the current grid
        /// </summary>
        public DisposableMeshGeometry MeshGeometry { get; set; }

        /// <summary>
        /// Coordinate system of the mesh
        /// </summary>
        public ICoordinateSystem MeshCoordinateSystem { get; set; }

        /// <summary>
        /// Sub-selection of the orthogonalization parameters that can be set
        /// </summary>
        public OrthogonalizationParameters OrthogonalizationParameters { get; } = OrthogonalizationParameters.CreateDefault();

        /// <summary>
        /// The parameters for the make grid algorithm
        /// </summary>
        public MakeGridParameters MakeGridParameters { get; } =  MakeGridParameters.CreateDefault();

        /// <summary>
        /// The parameters for curvilinear grids
        /// </summary>
        public CurvilinearParameters CurvilinearParameters { get; } =  CurvilinearParameters.CreateDefault();
        
        /// <summary>
        /// The parameters for interpolating samples
        /// </summary>
        public InterpolationParameters InterpolationParameters { get; } = InterpolationParameters.CreateDefault();

        /// <summary>
        /// The parameters for refining the grid based on samples
        /// </summary>
        public SamplesRefineParameters SamplesRefineParameters { get; } = SamplesRefineParameters.CreateDefault();

        /// <summary>
        /// The distance used when refining the polygon
        /// </summary>
        public double PolygonRefinementDistance { get; set; } = 10.0;

        /// <summary>
        /// The distance for offsetting a polygon
        /// </summary>
        public double PolygonOffsetDistance { get; set; } = 10.0;

        /// <summary>
        /// The number of points between splines corner points
        /// </summary>
        public int PointsBetweenSplineCornerPoints { get; } = 20;

        /// <summary>
        /// The parameters of the advancing front algorithm
        /// </summary>
        public SplinesToCurvilinearParameters SplinesToCurvilinearParameters { get; } = SplinesToCurvilinearParameters.CreateDefault();

        /// <summary>
        /// Clears all set data
        /// </summary>
        public void Reset()
        {
            MeshGeometry?.Dispose();
            MeshGeometry = null;

            Polygons.Clear();
            LandBoundaries.Clear();
            Splines.Clear();
            SelectedVertices = null;

            SamplePoints.PointValues.Clear();
        }
    }
}