using System.ComponentModel;

namespace DeltaShell.Plugins.GridEditor.Helpers
{
    public enum PolygonEditAction
    {
        [Description("Create new polygon and remove current selection")]
        Replace,

        [Description("Add polygon to current selection")]
        Add,

        [Description("Merge polygon with current selection")]
        Merge,

        [Description("Subtract polygon from current selection")]
        Subtract,

        [Description("Add polygon and subtract current selection")]
        AddAndClip
    }
}