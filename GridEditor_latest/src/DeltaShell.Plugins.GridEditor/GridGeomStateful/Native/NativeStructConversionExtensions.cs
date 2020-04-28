using System;
using System.Collections.Generic;
using System.Linq;
using DelftTools.Utils.Collections;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using GeoAPI.Extensions.Coverages;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Native
{
    public static class NativeStructConversionExtensions
    {
        private const double GeometrySeparator = -999.0;
        private const double InnerOuterSeparator = -998.0;

        internal static OrthogonalizationParametersNative ToOrthogonalizationParametersNative(this OrthogonalizationParameters orthogonalizationParameters)
        {
            return new OrthogonalizationParametersNative
            {
                OuterIterations = orthogonalizationParameters.OuterIterations,
                BoundaryIterations = orthogonalizationParameters.BoundaryIterations,
                InnerIterations = orthogonalizationParameters.InnerIterations,
                OrthogonalizationToSmoothingFactor = 0.975,
                AtpfB = 1,
                Circumormasscenter = 1,
                Smoothorarea = 1.0,
                AdaptMethod = 1,
                AdaptBeta = 0,
                AdaptNiterU = 0,
                AdaptNiterG = 4,
                OrthoPure = 0.5
            };
        }

        internal static MakeGridParametersNative ToMakeGridParametersNative(this MakeGridParameters makeGridParameters)
        {
            return new MakeGridParametersNative
            {
                GridType = (int)makeGridParameters.GridType,
                NumberOfColumns = makeGridParameters.NumberOfColumns,
                NumberOfRows = makeGridParameters.NumberOfRows,
                GridAngle = makeGridParameters.GridAngle,
                GridBlockSize = makeGridParameters.GridBlockSize,
                OriginXCoordinate = makeGridParameters.OriginXCoordinate,
                OriginYCoordinate = makeGridParameters.OriginYCoordinate,
                OriginZCoordinate = makeGridParameters.OriginZCoordinate,
                XGridBlockSize = makeGridParameters.XGridBlockSize,
                YGridBlockSize = makeGridParameters.YGridBlockSize
            };
        }

        internal static CurvilinearParametersNative ToCurvilinearParametersNative(this CurvilinearParameters curvilinearParameters)
        {
            return new CurvilinearParametersNative
            {
                MRefinement = curvilinearParameters.MRefinement,
                NRefinement = curvilinearParameters.NRefinement,
                SmoothingIterations = curvilinearParameters.SmoothingIterations,
                SmoothingParameter = curvilinearParameters.SmoothingParameter,
                AttractionParameter = curvilinearParameters.AttractionParameter
            };
        }

        internal static SplinesToCurvilinearParametersNative ToSplinesToCurvilinearParametersNative(this SplinesToCurvilinearParameters splinesToCurvilinearParameters)
        {
            return new SplinesToCurvilinearParametersNative
            {
                AspectRatio = splinesToCurvilinearParameters.AspectRatio,
                AspectRatioGrowFactor = splinesToCurvilinearParameters.AspectRatioGrowFactor,
                AverageWidth = splinesToCurvilinearParameters.AverageWidth,
                CurvatureAdapetedGridSpacing = splinesToCurvilinearParameters.CurvatureAdaptedGridSpacing,
                GrowGridOutside = splinesToCurvilinearParameters.GrowGridOutside,
                MaximumNumberOfGridCellsInTheUniformPart = splinesToCurvilinearParameters.MaximumNumberOfGridCellsInTheUniformPart,
                GridsOnTopOfEachOtherTolerance = splinesToCurvilinearParameters.GridsOnTopOfEachOtherTolerance,
                MinimumCosineOfCrossingAngles = splinesToCurvilinearParameters.MinimumCosineOfCrossingAngles,
                CheckFrontCollisions = splinesToCurvilinearParameters.CheckFrontCollisions,
                UniformGridSize = splinesToCurvilinearParameters.UniformGridSize,
                RemoveSkinnyTriangles = splinesToCurvilinearParameters.RemoveSkinnyTriangles
            };
        }

        internal static InterpolationParametersNative ToInterpolationParametersNative(this InterpolationParameters interpolationParameters)
        {
            return new InterpolationParametersNative
            {
                InterpolationType = interpolationParameters.InterpolationType,
                DisplayInterpolationProcess = interpolationParameters.DisplayInterpolationProcess,
                MaxNumberOfRefinementIterations = interpolationParameters.MaxNumberOfRefinementIterations,
                AveragingMethod = interpolationParameters.AveragingMethod,
                MinimumNumberOfPoints = interpolationParameters.MinimumNumberOfPoints,
                RelativeSearchRadius = interpolationParameters.RelativeSearchRadius,
                InterpolateTo = interpolationParameters.InterpolateTo
            };
        }

        internal static SampleRefineParametersNative ToSampleRefineParametersNative(this SamplesRefineParameters samplesRefineParameters)
        {
            return new SampleRefineParametersNative
            {
                SampleVectorDimension = samplesRefineParameters.SampleVectorDimension,
                MinimumCellSize = samplesRefineParameters.MinimumCellSize,
                DirectionalRefinement = samplesRefineParameters.DirectionalRefinement,
                RefinementType = samplesRefineParameters.RefinementType,
                ConnectHangingNodes = samplesRefineParameters.ConnectHangingNodes,
                MaximumTimeStepInCourantGrid = samplesRefineParameters.MaximumTimeStepInCourantGrid,
                AccountForSamplesOutside = samplesRefineParameters.AccountForSamplesOutside
            };
        }

        /// <summary>
        /// Converts a multi-polygon into a <see cref="DisposableGeometryList"/>
        /// </summary>
        /// <param name="features">Multi-polygon to convert</param>
        /// <returns><see cref="DisposableGeometryList"/> with the multi-polygon information</returns>
        public static DisposableGeometryList ToDisposableGeometryList(this IList<IFeature> features)
        {
            return DisposableGeometryListFromGeometries(features?.Select(f => f.Geometry).ToList());
        }

        /// <summary>
        /// Converts a coordinate array into a <see cref="DisposableGeometryList"/>
        /// </summary>
        /// <param name="coordinates">Coordinate array to convert</param>
        /// <returns><see cref="DisposableGeometryList"/> with the coordinate array information</returns>
        public static DisposableGeometryList ToDisposableGeometryList(this Coordinate[] coordinates)
        {
            return DisposableGeometryListFromGeometries(new IGeometry[] { new LineString(coordinates) });
        }

        public static DisposableGeometryList ToDisposableGeometryList(this IList<IPointValue> pointValues)
        {
            return DisposableGeometryListFromGeometries(new IGeometry[]
            {
                new LineString(pointValues.Select(v => new Coordinate(v.X, v.Y, v.Value)).ToArray())
            });
        }

        /// <summary>
        /// Converts a list of geometries into a <see cref="DisposableGeometryList"/>
        /// </summary>
        /// <param name="geometries"></param>
        /// <returns><see cref="DisposableGeometryList"/> with the coordinate array information</returns>
        public static DisposableGeometryList DisposableGeometryListFromGeometries(this IList<IGeometry> geometries)
        {
            var xCoordinates = new List<double>();
            var yCoordinates = new List<double>();
            var zCoordinates = new List<double>();

            if (geometries == null)
            {
                return new DisposableGeometryList
                {
                    GeometrySeparator = GeometrySeparator,
                    InnerOuterSeparator = InnerOuterSeparator,
                    NumberOfCoordinates = xCoordinates.Count,
                    XCoordinates = xCoordinates.ToArray(),
                    YCoordinates = yCoordinates.ToArray(),
                    ZCoordinates = zCoordinates.ToArray()
                };
            }

            var coordinateArrays = new[] { xCoordinates, yCoordinates, zCoordinates };

            for (int i = 0; i < geometries.Count; i++)
            {
                var geometry = geometries[i];
                if (geometry == null) continue;

                if (i != 0)
                {
                    coordinateArrays.ForEach(a => a.Add(GeometrySeparator));
                }

                if (geometry is IPolygon polygon)
                {
                    AddCoordinatesToArrays(polygon.ExteriorRing.Coordinates, xCoordinates, yCoordinates, zCoordinates);

                    for (int j = 0; j < polygon.InteriorRings.Length; j++)
                    {
                        coordinateArrays.ForEach(a => a.Add(InnerOuterSeparator));

                        var interiorRing = polygon.InteriorRings[j];
                        AddCoordinatesToArrays(interiorRing.Coordinates, xCoordinates, yCoordinates, zCoordinates);
                    }
                }
                else
                {
                    AddCoordinatesToArrays(geometry.Coordinates, xCoordinates, yCoordinates, zCoordinates);
                }
            }

            return new DisposableGeometryList
            {
                GeometrySeparator = GeometrySeparator,
                InnerOuterSeparator = InnerOuterSeparator,
                NumberOfCoordinates = xCoordinates.Count,
                XCoordinates = xCoordinates.ToArray(),
                YCoordinates = yCoordinates.ToArray(),
                ZCoordinates = zCoordinates.ToArray()
            };
        }

        public static DisposableGeometryList CreateEmptyDisposableGeometryList(int length)
        {
            return new DisposableGeometryList
            {
                XCoordinates = new double[length],
                YCoordinates = new double[length],
                ZCoordinates = new double[length],
                GeometrySeparator = GeometrySeparator,
                InnerOuterSeparator = InnerOuterSeparator,
                NumberOfCoordinates = length
            };
        }

        public static ICollection<IFeature> ToFeatureList(this DisposableGeometryList disposableGeometryList)
        {
            var features = new List<IFeature>();

            if (disposableGeometryList == null)
                return features;

            var innerOuterSeparator = -998.0;
            var geometrySeparator = disposableGeometryList.GeometrySeparator;

            var coordinatesOuter = new List<Coordinate>();
            var coordinatesInner = new List<List<Coordinate>>();
            var currentCoordinateList = coordinatesOuter;

            for (int i = 0; i < disposableGeometryList.NumberOfCoordinates; i++)
            {
                var xCoordinate = disposableGeometryList.XCoordinates[i];
                var yCoordinate = disposableGeometryList.YCoordinates[i];

                if (Math.Abs(xCoordinate - innerOuterSeparator) < 1e-10)
                {
                    // add coordinate list for inner rings
                    currentCoordinateList = new List<Coordinate>();
                    coordinatesInner.Add(currentCoordinateList);
                    continue;
                }

                if (Math.Abs(xCoordinate - geometrySeparator) < 1e-10 ||
                    i == disposableGeometryList.NumberOfCoordinates - 1)
                {
                    // add first coordinate to close the linear 
                    coordinatesOuter.Add(new Coordinate(coordinatesOuter[0]));
                    coordinatesInner.ForEach(c => c.Add(c[0]));

                    // create polygon
                    var shell = new LinearRing(coordinatesOuter.ToArray());
                    var linearRings = coordinatesInner
                        .Select(c => (ILinearRing) new LinearRing(c.ToArray()))
                        .ToArray();

                    var polygon = new Polygon(shell, linearRings);
                    features.Add(new Feature { Geometry = polygon });
                    coordinatesOuter.Clear();
                    coordinatesInner.Clear();
                    currentCoordinateList = coordinatesOuter;
                    continue;
                }


                currentCoordinateList.Add(new Coordinate(xCoordinate, yCoordinate));
            }

            return features;
        }

        private static void AddCoordinatesToArrays(Coordinate[] interiorRingCoordinates, 
            ICollection<double> xCoordinates,
            ICollection<double> yCoordinates,
            ICollection<double> zCoordinates)
        {
            for (var index = 0; index < interiorRingCoordinates.Length; index++)
            {
                var coordinate = interiorRingCoordinates[index];
                xCoordinates.Add(coordinate.X);
                yCoordinates.Add(coordinate.Y);
                zCoordinates.Add(double.IsNaN(coordinate.Z) ? 0.0 : coordinate.Z);
            }
        }
    }
}