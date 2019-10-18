using DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapLayers.FeatureEditors
{
    [TestFixture]
    public class GenericFeatureEditorTest
    {
        [Test]
        public void GenericFeatureEditor_ShouldCreateProvidedInteractorType_WhenCallingCreateInteractor()
        {
            //Arrange
            var afterCreateCalled = false;
            var genericFeatureEditor = new GenericFeatureEditor<PolygonEditorFeatureInteractor>{AfterCreate = i =>
            {
                afterCreateCalled = true;
            }};
            
            // Act
            var interactor = genericFeatureEditor.CreateInteractor(null, null);

            // Assert
            Assert.IsTrue(afterCreateCalled);
            Assert.IsInstanceOf<PolygonEditorFeatureInteractor>(interactor);
        }
    }
}