using System;
using GeoAPI.Geometries;
using ProtoBufRemote;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    public class CoordinateConverter : ITypeToProtoConverter
    {
        public object ToProtoObject(object original)
        {
            var coordinate = original as Coordinate;
            return coordinate == null
                ? null
                : new[] { coordinate.X, coordinate.Y, coordinate.Z };
        }

        public object FromProtoObject(object protoObject)
        {
            if (!(protoObject is double[] array)) return null;
            return new Coordinate(array[0], array[1], array[2]);
        }

        public Type GetProtoType()
        {
            return typeof(double[]);
        }

        public Type GetSourceType()
        {
            return typeof(Coordinate);
        }
    }
}