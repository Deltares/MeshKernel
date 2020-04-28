using System.ComponentModel;
using DelftTools.Utils;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon
{
    [TypeConverter(typeof(EnumDescriptionAttributeTypeConverter))]
    public enum GenerateGridWorkFlowType
    {
        [Description("Curvelinear grid (transfinite interp.)")]
        CurveLinearTransfinite,
        [Description("Curvelinear grid (orthogonal)")]
        CurveLinearOrthogonal,
        [Description("Rectangular grid")]
        Regular,
        [Description("Triangular grid within polygon")]
        TriangularWithPolygon,
        [Description("Rectangular grid within polygon")]
        RegularWithPolygon,
        [Description("Grid from samples")]
        Samples,
        [Description("Grid refinement based on samples")]
        RefinementSamples
    }
}