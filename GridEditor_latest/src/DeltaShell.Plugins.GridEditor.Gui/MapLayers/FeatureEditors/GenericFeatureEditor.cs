using System;
using GeoAPI.Extensions.Feature;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.Editors;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors
{
    internal class GenericFeatureEditor<T> : FeatureEditor where T : IFeatureInteractor
    {
        public override IFeatureInteractor CreateInteractor(ILayer layer, IFeature feature)
        {
            var featureInteractor = (T) Activator.CreateInstance(typeof(T), layer, feature, null, null);
            AfterCreate?.Invoke(featureInteractor);
            return featureInteractor;
        }

        internal Action<T> AfterCreate { get; set; }
    }
}