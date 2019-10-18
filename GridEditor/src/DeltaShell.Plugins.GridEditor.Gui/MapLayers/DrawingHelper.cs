using System.Drawing;
using System.Drawing.Drawing2D;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers
{
    internal static class DrawingHelper
    {
        internal static Bitmap GenerateEllipseBitmap(Color penColor, Color brushColor, int width, int height)
        {
            using (var pen = new Pen(penColor))
            using (var brush = new SolidBrush(brushColor))
            {
                return GenerateEllipseBitmap(pen, brush, width, height);
            }
        }

        internal static Bitmap GenerateEllipseBitmap(Pen pen, Brush brush, int width, int height)
        {
            var penWidth = (int)pen.Width;
            width += 1;
            height += 1;
            var bm = new Bitmap(width + penWidth, height + penWidth);

            using (var graphics = Graphics.FromImage(bm))
            {
                graphics.SmoothingMode = SmoothingMode.HighQuality;

                var halfPenWidth = penWidth / 2.0f;
                var rect = new RectangleF(halfPenWidth, halfPenWidth, width -1, height -1);

                graphics.FillEllipse(brush, rect);
                graphics.DrawEllipse(pen, rect);
            }

            return bm;
        }
    }
}