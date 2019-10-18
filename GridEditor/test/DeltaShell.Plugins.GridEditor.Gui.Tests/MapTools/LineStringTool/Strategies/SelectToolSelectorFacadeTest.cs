using System;
using System.Collections.Generic;
using System.Linq;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using GeoAPI.Extensions.Feature;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.Api.Editors;
using SharpMap.Layers;
using SharpMap.UI.Forms;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool.Strategies
{
    [TestFixture]
    public class SelectToolSelectorFacadeTest
    {
        [Test]
        public void GivenSelectToolSelectorFacade_HasSelection_ShouldCheckForSelectedFeatureInteractors()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();

            mocks.ReplayAll();
            
            var selectTool = new SelectTool();
            var facade = new SelectToolSelectorFacade(selectTool);

            // Act & Assert
            Assert.IsFalse(facade.HasSelection);

            selectTool.SelectedFeatureInteractors.Add(featureInteractor);

            Assert.IsTrue(facade.HasSelection);
            mocks.VerifyAll();
        }

        [Test]
        public void GivenSelectToolSelectorFacade_SelectedFeatureInteractors_ShouldReturnSelectoolSelectedFeatureInteractors()
        {
            //Arrange
            var selectTool = new SelectTool();
            var facade = new SelectToolSelectorFacade(selectTool);

            // Act & Assert
            Assert.AreEqual(selectTool.SelectedFeatureInteractors,facade.SelectedFeatureInteractors);
        }

        [Test]
        public void GivenSelectToolSelectorFacade_ClearSelection_ShouldClearSelectoolSelectedFeatureInteractors()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();

            mocks.ReplayAll();

            var selectTool = new SelectTool();
            var facade = new SelectToolSelectorFacade(selectTool);

            // Act & Assert
            selectTool.SelectedFeatureInteractors.Add(featureInteractor);

            Assert.AreEqual(1, facade.SelectedFeatureInteractors.Count());

            facade.ClearSelection();

            Assert.AreEqual(0, selectTool.SelectedFeatureInteractors.Count);

            mocks.VerifyAll();
        }

        [Test]
        public void GivenSelectToolSelectorFacade_Select_ShouldAddFeatureInteractorToSelectool()
        {
            //Arrange
            var vectorLayer = new VectorLayer();
            var mocks = new MockRepository();

            var featureInteractor = mocks.StrictMultiMock<IFeatureInteractor>(typeof(IDisposable));
            var featureEditor = mocks.StrictMock<IFeatureEditor>();
            var feature = mocks.StrictMock<IFeature>();
            var map = mocks.StrictMock<IMap>();
            var mapControl = mocks.StrictMock<IMapControl>();

            vectorLayer.FeatureEditor = featureEditor;

            featureInteractor.Expect(i => i.Trackers).Return(new List<TrackerFeature>()).Repeat.AtLeastOnce();
            featureInteractor.Expect(i => i.SourceFeature).Return(feature).Repeat.AtLeastOnce();

            featureEditor.Expect(e => e.CreateInteractor(vectorLayer, feature)).Return(featureInteractor).Repeat.AtLeastOnce();

            mapControl.Expect(c => c.Map).Return(map).Repeat.AtLeastOnce();
            mapControl.Expect(c => c.SelectedFeatures = null).IgnoreArguments().Repeat.AtLeastOnce();

            map.Expect(m => m.IsDisposing).Return(false).Repeat.AtLeastOnce();
            map.Expect(m => m.CoordinateSystem).Return(null);
            map.Expect(m => m.GetLayerByFeature(feature)).Return(vectorLayer);

            mocks.ReplayAll();

            var selectTool = new SelectTool{MapControl = mapControl};
            var facade = new SelectToolSelectorFacade(selectTool);

            // Act & Assert
            Assert.AreEqual(0, facade.SelectedFeatureInteractors.Count());

            facade.Select(feature);

            Assert.AreEqual(1, selectTool.SelectedFeatureInteractors.Count);

            mocks.VerifyAll();
        }
    }
}