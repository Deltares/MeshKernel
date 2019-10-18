namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool
{
    internal enum LineStringMapToolStateCommand
    {
        MouseDownNotOnFeatureNoSelection,
        MouseDownNotOnFeatureHasSelection,
        MouseDownOnOtherFeature,
        MouseDownOnSelectedFeature,
        MouseDownOnFeatureOnTracker,

        AfterMouseDownFeatureSelected,

        AfterMouseDoubleClick,
        AfterMouseUp,
    }
}