using System.ComponentModel;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon
{
    public enum GenerateGridWorkFlowType
    {
        [Description("Curvelinear (transfinite interp.)")]
        CurveLinearTransfinite,
        [Description("Curvelinear (orthogonal)")]
        CurveLinearOrthogonal,
        [Description("Regular")]
        Regular,
        [Description("Triangular with polygon")]
        TriangularWithPolygon,
        [Description("Regular with polygon")]
        RegularWithPolygon,
        [Description("Samples based")]
        Samples,
        [Description("RefinementSamples based on samples")]
        RefinementSamples
    }
}