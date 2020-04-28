using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using DelftTools.Utils;
using DelftTools.Utils.IO;
using GeoAPI.Extensions.Coverages;
using log4net;
using NetTopologySuite.Extensions.Coverages;

namespace DeltaShell.Plugins.GridEditor.Importers
{
    public static class AscReaderFast
    {
        private static readonly ILog log = LogManager.GetLogger(typeof(AscReaderFast));

        /// <summary>
        /// Reads an asc file and converts it to a pointValue list
        /// </summary>
        /// <param name="path">Path to the asc file</param>
        /// <param name="skipNoData">If true, removes points that have no data value</param>
        /// <returns>List of asc file values as <see cref="IPointValue"/></returns>
        public static IList<IPointValue> Read(string path, bool skipNoData = false)
        {
            if (!File.Exists(path))
            {
                log.Error($"Can not find {path}");
                return new List<IPointValue>();
            }

            if (FileUtils.IsFileLocked(path))
            {
                log.Error($"Can open {path} because it is locked");
                return new List<IPointValue>();
            }

            try
            {
                using (var reader = new StreamReader(new FileStream(path, FileMode.Open)))
                {
                    var cultureInfo = CultureInfo.InvariantCulture;
                    var ncols = GetIntValueFromParameter(reader, "ncols", cultureInfo);
                    var nrows = GetIntValueFromParameter(reader, "nrows", cultureInfo);
                    var xllcorner = GetDoubleValueFromParameter(reader, "xllcorner", cultureInfo);
                    var yllcorner = GetDoubleValueFromParameter(reader, "yllcorner", cultureInfo);
                    var cellsize = GetDoubleValueFromParameter(reader, "cellsize", cultureInfo);
                    var NODATA_value = GetDoubleValueFromParameter(reader, "NODATA_value", cultureInfo);

                    var points = new ChunkedList<IPointValue>(ncols * nrows);

                    for (int curYIndex = nrows - 1; curYIndex >= 0; curYIndex--)
                    {
                        var readLine = reader.ReadLine();
                        if (string.IsNullOrEmpty(readLine)) continue;

                        var values = readLine
                            .Split(new[] {" "}, StringSplitOptions.RemoveEmptyEntries)
                            .Select(v => Convert.ToDouble(v, cultureInfo))
                            .ToArray();

                        for (int curXIndex = 0; curXIndex < ncols; curXIndex++)
                        {
                            var value = values[curXIndex];
                            double x = xllcorner + cellsize * curXIndex;
                            double y = yllcorner + cellsize * curYIndex;

                            if (skipNoData && value == NODATA_value) continue;

                            var pointValue = new PointValue
                            {
                                X = x,
                                Y = y,
                                Value = value
                            };

                            points.Add(pointValue);
                        }
                    }

                    return points;
                }
            }
            catch (FormatException)
            {
                log.Error($"Error reading asc file {path}");
            }
            catch (UnauthorizedAccessException)
            {
                log.Error($"Could not access the file {path}");
            }

            return new List<IPointValue>();
        }

        private static double GetDoubleValueFromParameter(TextReader reader, string parameterName, IFormatProvider cultureInfo)
        {
            var readLine = reader.ReadLine();

            return readLine != null 
                ? Convert.ToDouble(readLine.Remove(0, parameterName.Length), cultureInfo) 
                : default(double);
        }

        private static int GetIntValueFromParameter(TextReader reader, string parameterName, IFormatProvider cultureInfo)
        {
            var readLine = reader.ReadLine();

            return readLine != null
                ? Convert.ToInt32(readLine.Remove(0, parameterName.Length), cultureInfo)
                : default(int);
        }
    }
}