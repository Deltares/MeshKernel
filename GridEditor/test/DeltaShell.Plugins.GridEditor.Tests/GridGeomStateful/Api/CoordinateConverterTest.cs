using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using GeoAPI.Geometries;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.GridGeomStateful.Api
{
    [TestFixture]
    public class CoordinateConverterTest
    {
        [Test]
        public void ConvertCoordinateToDoubleArray()
        {
            var x = 10;
            var y = 100;
            var z = 2;

            var coordinate = new Coordinate(x, y, z);
            var converter = new CoordinateConverter();

            var convertedValue = converter.ToProtoObject(coordinate);

            Assert.AreEqual(typeof(double[]), converter.GetProtoType());

            Assert.AreEqual(new []{x,y,z},convertedValue);

            convertedValue = converter.ToProtoObject(null);
            Assert.IsNull(convertedValue);
        }

        [Test]
        public void ConvertDoubleArrayToCoordinate()
        {
            var x = 10;
            var y = 100;
            var z = 2;

            var doubles = new double[] { x, y, z };
            var converter = new CoordinateConverter();

            var convertedValue = converter.FromProtoObject(doubles);
            Assert.AreEqual(typeof(Coordinate), converter.GetSourceType());

            Assert.AreEqual(new Coordinate(x, y, z), convertedValue);

            convertedValue = converter.FromProtoObject(null);
            Assert.IsNull(convertedValue);
        }
    }
}