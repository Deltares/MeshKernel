using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using NUnit.Framework;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapLayers.Renderers
{
    [TestFixture]
    public class DisposableMeshGeometryRendererTest
    {
        [Test]
        [TestCase(nameof(DisposableMeshGeometryRenderer.ShowGridVertices), false, true)]
        [TestCase(nameof(DisposableMeshGeometryRenderer.DrawEdgeValues), false, true)]
        [TestCase(nameof(DisposableMeshGeometryRenderer.PointSize), 5, 15)]
        [TestCase(nameof(DisposableMeshGeometryRenderer.DrawSelectedVertices), false, true)]
        public void GivenDisposableMeshGeometryRenderer_SettingDrawingProperties_SetsRenderRequired(string propertyName, object valueBefore, object valueToSet)
        {
            //Arrange
            var render = new DisposableMeshGeometryRenderer{RenderRequired = false};

            Assert.IsFalse(render.RenderRequired);

            // Act
            var propertyInfo = render.GetType().GetProperty(propertyName);

            Assert.NotNull(propertyInfo);
            Assert.AreEqual(valueBefore,propertyInfo.GetValue(render));

            propertyInfo?.SetValue(render, valueToSet);

            // Assert
            Assert.IsTrue(render.RenderRequired);
            Assert.AreEqual(valueToSet, propertyInfo.GetValue(render));
        }

        [Test]
        public void GivenDisposableMeshGeometryRenderer_GettingEnvelope_UsesGridExtend()
        {
            //Arrange
            var render = new DisposableMeshGeometryRenderer();

            // Act
            Assert.IsTrue(render.Envelope.IsNull, "Render envelope should be null without mesh");

            render.GetMeshGeometry = () => new DisposableMeshGeometry
            {
                xNodes = new double[] { 10, 10, 20 },
                yNodes = new double[] { 0, 10, 20 },
                numberOfNodes = 3
            };

            // Assert
            var renderEnvelope = render.Envelope;

            Assert.IsFalse(renderEnvelope.IsNull, "Render envelope should be set with mesh");
            Assert.AreEqual(10, renderEnvelope.Width);
            Assert.AreEqual(20, renderEnvelope.Height);
        }
    }
}