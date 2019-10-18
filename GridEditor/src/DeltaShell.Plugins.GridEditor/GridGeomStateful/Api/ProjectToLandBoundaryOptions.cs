using System.ComponentModel;
using DelftTools.Utils;

namespace DeltaShell.Plugins.GridEditor.GridGeomStateful.Api
{
    [TypeConverter(typeof(EnumDescriptionAttributeTypeConverter))]
    public enum ProjectToLandBoundaryOptions
    {
        [Description("Do not project to land boundary")]
        No = 0,

        [Description("To original net-boundary")]
        ToOriginalNetBoundary = 1,

        [Description("Net-boundary to landBoundary")]
        NetBoundaryToLandBoundary = 2,

        [Description("Net-boundary to land boundary and inner to land")]
        NetBoundaryToLandBoundaryAndInnerToLand = 3,

        [Description("Whole net")]
        WholeNet = 4,

        [Description("Ok")]
        Ok = 5
    }
}