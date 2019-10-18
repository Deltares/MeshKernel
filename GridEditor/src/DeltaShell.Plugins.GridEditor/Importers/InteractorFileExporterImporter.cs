using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using DelftTools.Utils;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using log4net;
using static System.String;

namespace DeltaShell.Plugins.GridEditor.Importers
{
    public static class InteractorFileExporterImporter
    {
        private static readonly ILog log = LogManager.GetLogger(typeof(InteractorFileExporterImporter));

        public static void Export<T>(ref IEnumerable<T> features, string path) where T : IFeature
        {
            var directory = Path.GetDirectoryName(path);
            if (!IsNullOrEmpty(directory) && !Directory.Exists(directory))
            {
                Directory.CreateDirectory(directory);
            }

            using (var writer = new StreamWriter(path))
            using (CultureUtils.SwitchToInvariantCulture())
            {
                var i = 1;
                const int numColumns = 2; // X, Y
                foreach (var feature in features)
                {
                    //here we should use the feature name
                    writer.WriteLine(nameof(feature) + "_" + i++);
                    writer.WriteLine($"    {feature.Geometry.Coordinates.Length}    {numColumns}");

                    foreach (var coordinate in feature.Geometry.Coordinates)
                    {
                        writer.WriteLine($"{coordinate.X,24}{coordinate.Y,24}");
                    }
                }
            }
        }

        public static IEnumerable<T> Import<T>(string path, Func<Coordinate[],T> createFeature) where T: IFeature
        {
            if (!File.Exists(path))
            {
                log.Error($"Can not import {path}, because it does not exist.");
                return Enumerable.Empty<T>();
            }

            var features = new List<T>();

            using (var reader = new StreamReader(path))
            using (CultureUtils.SwitchToInvariantCulture())
            {
                int lineCount = 0;
                
                // Read the comments first
                var line = "";
                var splineName = "";
                var numPoints = -1;
                var numColumns = -1;
                var points = new List<Coordinate>();
                var numPointsToRead = 0;

                // Start reading each spline
                while (line != null)
                {
                    line = reader.ReadLine()?.Trim();
                    lineCount++;

                    if (IsNullOrEmpty(line) ||
                        line[0] == '#' ||
                        line[0] == '!' ||
                        line[0] == '*')
                    {
                        continue;
                    }

                    if (IsNullOrEmpty(splineName))
                    {
                        // start of the spline, skip the name
                        splineName = line;
                        continue;
                    }

                    var lineFields = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);

                    if (numPoints == -1 || numColumns == -1)
                    {
                        if (lineFields.Length != 2)
                        {
                            throw new FormatException($"Invalid numpoints/numcolums {lineFields} in file {path}");
                        }

                        if (!(int.TryParse(lineFields[0], out numPoints) &&
                              int.TryParse(lineFields[1], out numColumns) && numColumns == 2))
                        {
                            throw new FormatException($"Invalid format of numPoints or numColumns : {lineFields} in file {path}");
                        }

                        numPointsToRead = numPoints;
                        continue;
                    }

                    numPointsToRead--;
                    
                    if (lineFields.Length != 2)
                    {
                        throw new FormatException($"Invalid point row on line {lineCount} in file {path}");
                    }

                    if (!(double.TryParse(lineFields[0], out var x) &&
                          double.TryParse(lineFields[1], out var y)))
                    {
                        throw new FormatException($"Invalid format of x or y coordinate : {lineFields} in file {path}");
                    }

                    points.Add(new Coordinate(x, y));

                    if (numPointsToRead != 0)
                    {
                        continue;
                    }

                    features.Add(createFeature(points.ToArray()));

                    // reset temp variables
                    splineName = "";
                    numPoints = -1;
                    numColumns = -1;
                    points.Clear();
                }
            }

            return features;
        }
    }
}
