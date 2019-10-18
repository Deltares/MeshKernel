using System;
using System.Drawing;
using System.Drawing.Drawing2D;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers
{
    /// <summary>
    /// Polygon drawing style. Collection of <see cref="Pen"/> and <see cref="Brush"/> objects
    /// used in drawing polygons
    /// </summary>
    public sealed class PolygonDrawingStyle : IDisposable
    {
        private const int LineThickness = 2;

        /// <summary>
        /// Size of the vertices
        /// </summary>
        public const int PointSize = 6;

        /// <summary>
        /// Default color for selected line
        /// </summary>
        private static Color SelectedLineColor = Color.CornflowerBlue;

        /// <summary>
        /// Default line color for selected point
        /// </summary>
        private static Color SelectedPointLineColor = Color.CornflowerBlue;

        /// <summary>
        /// Default background color for selected point
        /// </summary>
        private static Color SelectedPointBackgroundColor = Color.White;

        /// <summary>
        /// Polygon fill brush
        /// </summary>
        public Brush PolygonBrush { get; } = new SolidBrush(Color.FromArgb(125, Color.LightSteelBlue));

        /// <summary>
        /// Polygon border pen
        /// </summary>
        public Pen PolygonBorderPen { get; } = new Pen(Color.Black, LineThickness);

        /// <summary>
        /// Point background brush during editing
        /// </summary>
        public Brush EditingPointBackgroundBrush { get; } = new SolidBrush(SelectedPointBackgroundColor);

        /// <summary>
        /// Point border pen during editing
        /// </summary>
        public Pen EditingPointBorderPen { get; } = new Pen(SelectedPointLineColor, LineThickness);

        /// <summary>
        /// Snapped point background
        /// </summary>
        public Brush SnappingPointBackgroundBrush { get; } = new SolidBrush(SelectedPointLineColor);

        /// <summary>
        /// Line pen during editing
        /// </summary>
        public Pen EditingLinePen { get; } = new Pen(SelectedLineColor, LineThickness);

        /// <summary>
        /// Last segment pen during polygon creation
        /// </summary>
        public Pen CreatingPolygonLastSegmentLinePen { get; } = new Pen(SelectedLineColor, LineThickness) {DashStyle = DashStyle.Dash};

        /// <summary>
        /// First point of polygon border pen
        /// </summary>
        public Pen FirstPolygonPointPen { get; } = new Pen(Color.Black, LineThickness);

        /// <summary>
        /// <inheritdoc cref="IDisposable"/>
        /// </summary>
        public void Dispose()
        {
            PolygonBrush?.Dispose();
            PolygonBorderPen?.Dispose();
            EditingPointBackgroundBrush?.Dispose();
            EditingPointBorderPen?.Dispose();
            SnappingPointBackgroundBrush?.Dispose();
            EditingLinePen?.Dispose();
            CreatingPolygonLastSegmentLinePen?.Dispose();
            FirstPolygonPointPen?.Dispose();
        }
    }   
}