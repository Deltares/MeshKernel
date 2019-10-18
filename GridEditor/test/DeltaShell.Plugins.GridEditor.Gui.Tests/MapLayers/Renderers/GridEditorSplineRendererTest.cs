using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Runtime.InteropServices;
using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using GeoAPI.CoordinateSystems.Transformations;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.Api.Layers;
using Point = NetTopologySuite.Geometries.Point;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapLayers.Renderers
{
    [TestFixture]
    public class GridEditorSplineRendererTest
    {
        [Test]
        public void GridEditorSplineRenderer_ShouldUseCoordinateTransformation_WhenReturningGetRenderedFeatureGeometry()
        {
            //Arrange
            var mocks = new MockRepository();
            var layer = mocks.StrictMock<ILayer>();
            var coordinateTransformation = mocks.StrictMock<ICoordinateTransformation>();
            var mathTransform = mocks.StrictMock<IMathTransform>();

            layer.Expect(l => l.CoordinateTransformation).Return(coordinateTransformation).Repeat.Any();
            coordinateTransformation.Expect(t => t.MathTransform).Return(mathTransform).Repeat.Any();
            mathTransform.Expect(t => t.Transform(new[] {10.0, 20.0})).Return(new []{100.0, 100.0});

            mocks.ReplayAll();

            var renderer = new GridEditorSplineRenderer();
            var feature = new Feature {Geometry = new Point(10, 20)};

            // Act

            renderer.GetRenderedFeatureGeometry(feature, layer);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GridEditorSplineRenderer_ShouldReferGetFeaturesToLayer_WhenGetFeaturesIsCalled()
        {
            //Arrange
            var mocks = new MockRepository();
            var layer = mocks.StrictMock<ILayer>();

            layer.Expect(l => l.GetFeatures(null, false))
                .IgnoreArguments()
                .Return(new IFeature[]{}).Repeat.Any();
            
            mocks.ReplayAll();

            var renderer = new GridEditorSplineRenderer();

            // Act
            renderer.GetFeatures(new LineString(new Coordinate[]{ }), layer);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GridEditorSplineRenderer_ShouldRenderFeatureCorrectly_WhenRenderIsCalled()
        {
            //Arrange
            var mocks = new MockRepository();
            var layer = mocks.StrictMock<ILayer>();
            var map = mocks.StrictMock<IMap>();

            map.Expect(m => m.WorldLeft).Return(-20).Repeat.Any();
            map.Expect(m => m.WorldTop).Return(60).Repeat.Any();
            map.Expect(m => m.PixelWidth).Return(0.5).Repeat.Any();
            map.Expect(m => m.PixelHeight).Return(0.5).Repeat.Any();

            layer.Expect(l => l.CoordinateTransformation).Return(null).Repeat.Any();
            layer.Expect(l => l.Map).Return(map).Repeat.Any();

            mocks.ReplayAll();

            var renderer = new GridEditorSplineRenderer();
            var spline = new Spline
            {
                Geometry = new LineString(new[]
                {
                    new Coordinate(10, 10),
                    new Coordinate(13, 13),
                    new Coordinate(17, 17),
                    new Coordinate(20, 20)
                }),
                UserGeometry = new LineString(new []{new Coordinate(10,10), new Coordinate(20,20)})
            };

            // Act
            var referenceImagePath = TestHelper.GetTestFilePath("GridEditorSplineRendererReferenceImage.png");
            using (var bitmap = new Bitmap(150,150))
            using (var g = Graphics.FromImage(bitmap))
            {
                renderer.Render(spline, g, layer);

                // Assert
                Assert.IsTrue(GuiTestHelper.CompareBitmapWithReferenceFile(bitmap, referenceImagePath));
            }
        }
    }
}